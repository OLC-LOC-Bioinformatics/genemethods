#!/usr/bin/env python
from olctools.accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject, run_subprocess, \
    SetupLogging
from genemethods.assemblypipeline.legacy_vtyper import epcr_primers, Filer
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser
from threading import Thread
from csv import DictReader
from queue import Queue
import multiprocessing
import logging
import shutil
import os


def fasta_primers(primerfile, forward_dict, reverse_dict):
    """
    Read in the FASTA-formatted primer file, and populate forward and reverse primer dictionaries
    :param primerfile: Name and absolute path to the FASTA-formatted primer file
    :param forward_dict: Dictionary of primer name: forward primer sequence
    :param reverse_dict: Dictionary of primer name: reverse primer sequence
    :return: Populated forward_dict and reverse_dict
    """
    logging.info('Populating primer dictionaries')
    for record in SeqIO.parse(primerfile, 'fasta'):
        # Set the primer name and sequence
        primer_name = record.id
        primer = str(record.seq)
        # Split the base name of the target from the direction
        # e.g. vtx1a-F1 is split in vtx1a and F1
        basename, direction = primer_name.split('-')
        # Populate the dictionaries of forward and reverse primers based on the direction determined above
        if direction.startswith('F'):
            # Add the current primer sequence to the dictionary
            forward_dict[basename] = primer
        else:
            reverse_dict[basename] = primer
    return forward_dict, reverse_dict


def epcr_primer_file(formattedprimers, forward_dict, reverse_dict, min_amplicon_size, max_amplicon_size):
    """
    Create the ipcress-compatible primer file from the dictionaries of primer combinations
    :param formattedprimers: Name and absolute path of the formatted primer file
    :param forward_dict: Dictionary of primer name: forward primer sequence
    :param reverse_dict: Dictionary of primer name: reverse primer sequence
    :param min_amplicon_size: Integer of the minimum amplicon size allowed
    :param max_amplicon_size: Integer of the maximum amplicon size allowed
    """
    logging.info('Creating ipcress-compatible primer file')
    with open(formattedprimers, 'w') as formatted:
        # Iterate through all the targets
        for primer_name in sorted(forward_dict):
            # Create the string to write to the ePCR-compatible primer file
            # e.g. vtx1a	CCTTTCCAGGTACAACAGCGGTT	GGAAACTCATCAGATGCCATTCTGG   0   1500
            output_string = '{pn}\t{fp}\t{rp}\t{min}\t{max}\n' \
                .format(pn=primer_name,
                        fp=forward_dict[primer_name],
                        rp=reverse_dict[primer_name],
                        min=min_amplicon_size,
                        max=max_amplicon_size)
            # Write the string to file
            formatted.write(output_string)


def ipcress_threads(metadata, analysistype, formattedprimers, mismatches, ipcress_queue, threads):
    """
    Run ipcress in a multi-threaded fashion
    :param metadata: List of metadata objects for each sample
    :param analysistype: String of the current analysis type
    :param formattedprimers: Name and absolute path of the formatted primer file
    :param mismatches: Integer of the current maximum number of mismatches allowed for the analyses
    :param ipcress_queue: Queue object
    :param threads: Integer of the number of threads to create
    :return: metadata: Updated list of metadata objects
    """
    # Create a thread targeting the ipcress function supplying ipcress queue and the analysis type
    for _ in range(threads):
        thread = Thread(target=ipcress, args=(ipcress_queue, analysistype))
        thread.setDaemon(True)
        thread.start()
    # Create the threads for the ePCR analysis
    logging.info('Running ipcress analyses')
    for sample in metadata:
        setattr(sample, analysistype, GenObject())
        # Get the primers ready
        sample[analysistype].primers = formattedprimers
        sample[analysistype].report = '{of}.txt' \
            .format(of=os.path.join(sample.general.outputdirectory, sample.name))
        sample.commands.ipcress = 'ipcress -i {pf} -P -m {mismatches} -s {fasta} > {report}' \
            .format(pf=formattedprimers,
                    mismatches=mismatches,
                    fasta=sample.general.bestassemblyfile,
                    report=sample[analysistype].report)
        if sample.general.bestassemblyfile != 'NA':
            ipcress_queue.put(sample)
    # Join the threads
    ipcress_queue.join()
    return metadata


def ipcress(ipcress_queue, analysistype):
    """
    Run the ipcress command
    :param ipcress_queue: Queue object
    :param analysistype: String of the current analysis type
    """
    while True:
        # sample = ipcress_queue.get()
        sample = ipcress_queue.get()
        # Make the output path
        make_path(sample.general.outputdirectory)
        # Run the commands if the ePCR output file doesn't exist
        if not os.path.isfile(sample[analysistype].report):
            run_subprocess(sample.commands.ipcress)
        # Signal that the thread is complete
        ipcress_queue.task_done()


def ipcress_parse(metadata, analysistype):
    """
    Parse the ipcress outputs to extract relevant results
    :param metadata: List of metadata objects
    :param analysistype: String of the current analysis type
    :return: metadata: Updated list of metadata objects
    """
    logging.info('Parsing ipcress outputs')
    for sample in metadata:
        # Initialise the results GenObject
        sample[analysistype].results = GenObject()
        if sample.general.bestassemblyfile != 'NA':
            # Open the output file
            with open(sample[analysistype].report, 'r') as ipcress_report:
                # Initialise the experiment name
                experiment = str()
                # The output will be multiple blocks starting with 'Ipcress result' e.g.
                '''
                Ipcress result
                --------------
                 Experiment: ES10
                    Primers: A B
                    Target: CP022407.1:filter(unmasked) Escherichia coli O121:H19 strain 16-9255 chromosome, complete genome
                    Matches: 21/21 19/19
                    Product: 103 bp (range 0-1500)
                Result type: forward
            
                ...GGAAAGAATACTGGACCAGTC............................. # forward
                   |||||||||||||||||||||-->
                5'-GGAAAGAATACTGGACCAGTC-3' 3'-GGTCATGGACACTTAGTCC-5' # primers
                                            <--|||||||||||||||||||
                ...............................GGTCATGGACACTTAGTCC... # revcomp
                --
                ipcress: CP022407.1:filter(unmasked) ES10 103 A 5388922 0 B 5389006 0 forward
                >ES10_product_1 seq CP022407.1:filter(unmasked) start 5388922 length 103
                GGAAAGAATACTGGACCAGTCGCTGGAATCTGCAACCGTTACTGCAAAGTGCTCAGTTGACAGGAATGAC
                TGTCACAATCAAATCCAGTACCTGTGAATCAGG
                '''
                for line in ipcress_report:
                    # Allow for multiple experiments with the same name by appending an iterator to the end e.g. ES10_0
                    iterator = 0
                    # Find the 'Experiment' line e.g. Experiment: ES10
                    if 'Experiment:' in line:
                        # Extract the experiment name e.g. ES10 and append the iterator
                        experiment = line.split()[-1].rstrip() + f'_{iterator}'
                        # Determine if the current experiment name has been initialised in the GenObject
                        try:
                            # While the experiment attribute exists, iterate through the experiment names
                            while sample[analysistype].results[experiment]:
                                # Increase the iterator
                                iterator += 1
                                # Recreate the experiment name
                                experiment = line.split()[-1].rstrip() + f'_{iterator}'
                            # Create the experiment name GenObject
                            sample[analysistype].results[experiment] = GenObject()
                        # Otherwise, the GenObject has not been initialise yet
                        except AttributeError:
                            sample[analysistype].results[experiment] = GenObject()
                    if 'Target:' in line:
                        contig = line.split('Target:')[-1].split()[0].replace(':filter(unmasked)', '')
                        try:
                            sample[analysistype].results[experiment][contig].contig = contig
                        except AttributeError:
                            sample[analysistype].results[experiment][contig] = GenObject()
                            sample[analysistype].results[experiment][contig].contig = contig
                    # Determine if the primers annealed to the forward or the reverse complement strand
                    # e.g. Result type: forward
                    if 'Result type:' in line:
                        # Set the boolean appropriately
                        if 'revcomp' in line:
                            direction = False
                        else:
                            direction = True
                    # Lines with the query sequence start with '...'
                    if line.startswith('...'):
                        # The top query sequence only has three periods
                        # e.g. ...GGAAAGAATACTGGACCAGTC............................. # forward
                        if not line.startswith('....'):
                            # If forward, populate the forward_query attribute with the cleaned line
                            if direction:
                                sample[analysistype].results[experiment][contig].forward_query = line.replace('.', '') \
                                    .replace(' # forward\n', '')
                            # Reverse
                            else:
                                sample[analysistype].results[experiment][contig].reverse_query = \
                                    line.replace('.', '').replace(' # forward\n', '')
                        # The bottom query sequence has many periods
                        # e.g. ...............................GGTCATGGACACTTAGTCC... # revcomp
                        else:
                            # Reverse
                            if direction:
                                sample[analysistype].results[experiment][contig].reverse_query = \
                                    revcomp(line.replace('.', '').replace(' # revcomp\n', ''))
                            # Forward
                            else:
                                sample[analysistype].results[experiment][contig].forward_query = \
                                    revcomp(seq_string=line.replace('.', '').replace(' # revcomp\n', ''))
                    # The line with the primers starts with 5'
                    # e.g. 5'-GGAAAGAATACTGGACCAGTC-3' 3'-GGTCATGGACACTTAGTCC-5' # primers
                    if line.startswith('5\''):
                        # Forward reverse
                        if direction:
                            # Split the cleaned line to create the forward_ref and reverse_ref attributes
                            sample[analysistype].results[experiment][contig].forward_ref, \
                                sample[analysistype].results[experiment][contig].reverse_ref = \
                                line.replace('5', '').replace('\'', '').replace('3', '').replace('-', '') \
                                    .replace(' # primers\n', '').split()
                            sample[analysistype].results[experiment][contig].reverse_ref = \
                                revcomp(seq_string=sample[analysistype].results[experiment][contig].reverse_ref)
                        # Reverse forward
                        else:
                            reverse_ref, forward_ref = \
                                line.replace('5', '').replace('\'', '').replace('3', '').replace('-', '') \
                                    .replace(' # primers\n', '').split()
                            sample[analysistype].results[experiment][contig].reverse_ref = reverse_ref
                            sample[analysistype].results[experiment][contig].forward_ref = \
                                revcomp(seq_string=forward_ref)
                    # Lines starting with 'ipcress:' contain results in an easily parsable format
                    # e.g. ipcress: CP022407.1:filter(unmasked) ES10 103 A 5388922 0 B 5389006 0 forward
                    # The eleven fields are contig (CP022407.1:filter(unmasked)), primer set (ES10), amplicon length (103),
                    # forward primer (A) - affected by reverse complement, position of forward primer (5388922), mismatches
                    # in forward primer (0), reverse primer (B), position of reverse primer (5389006), mismatches in reverse
                    # primer, analysis direction (forward) - will be 'revcomp' for experiments in reverse complement
                    if line.startswith('ipcress:'):
                        # Forward
                        if direction:
                            ipcress_str, \
                                sample[analysistype].results[experiment][contig].contig, \
                                sample[analysistype].results[experiment][contig].primer_set, \
                                sample[analysistype].results[experiment][contig].amplicon_length, \
                                sample[analysistype].results[experiment][contig].forward_primer, \
                                sample[analysistype].results[experiment][contig].forward_pos, \
                                sample[analysistype].results[experiment][contig].forward_mismatch, \
                                sample[analysistype].results[experiment][contig].reverse_primer, \
                                sample[analysistype].results[experiment][contig].reverse_pos, \
                                sample[analysistype].results[experiment][contig].reverse_mismatch, \
                                sample[analysistype].results[experiment][contig].direction = line.rstrip().split()
                            # Clean up the name of the contig e.g. CP022407.1:filter(unmasked) becomes CP022407.1:
                            sample[analysistype].results[experiment][contig].contig = \
                                sample[analysistype].results[experiment][contig].contig.replace(':filter(unmasked)', '')
                        # Reverse complement
                        else:
                            ipcress_str, \
                                sample[analysistype].results[experiment][contig].contig, \
                                sample[analysistype].results[experiment][contig].primer_set, \
                                sample[analysistype].results[experiment][contig].amplicon_length, \
                                sample[analysistype].results[experiment][contig].reverse_primer, \
                                sample[analysistype].results[experiment][contig].reverse_pos, \
                                sample[analysistype].results[experiment][contig].reverse_mismatch, \
                                sample[analysistype].results[experiment][contig].forward_primer, \
                                sample[analysistype].results[experiment][contig].forward_pos, \
                                sample[analysistype].results[experiment][contig].forward_mismatch, \
                                sample[analysistype].results[experiment][contig].direction = line.rstrip().split()
                            # Clean up the name of the contig e.g. CP022407.1:filter(unmasked) becomes CP022407.1:
                            sample[analysistype].results[experiment][contig].contig = \
                                sample[analysistype].results[experiment][contig].contig.replace(':filter(unmasked)', '')
                    # Find the extracted amplicon header
                    # e.g. >ES10_product_1 seq CP022407.1:filter(unmasked) start 5388922 length 103
                    if line.startswith('>'):
                        # Set the header attribute with the with cleaned the header string
                        sample[analysistype].results[experiment][contig].header = \
                            line.rstrip().lstrip('>').replace(':filter(unmasked)', '')
                        # Initialise a string to store the amplicon sequence
                        sample[analysistype].results[experiment][contig].sequence = str()
                        for subline in ipcress_report:
                            # There's an empty line after the sequence - break when it is encountered
                            if subline == '\n' or subline.startswith('--') or not subline:
                                break
                            # Add the sequence data to the growing string
                            sample[analysistype].results[experiment][contig].sequence += subline.rstrip()
    return metadata


def revcomp(seq_string, blast_primers=False):
    """
    Use SeqIO to revert primers in a reverse complement experiment to the forward orientation
    :param seq_string: String of primer sequence
    :param blast_primers: Boolean of whether the sequence to manipulate is from BLAST outputs
    :return: rev_comp: String of the calculated forward primer sequence
    """
    # Create a Seq object from the sequence string
    seq_record = Seq(seq_string)
    if not blast_primers:
        # Calculate the forward sequence by first finding the reverse complement, then the complement
        # e.g. seq_string: TTTAATGGTTACAGTCAT -> ATGACTGTAACCATTAAA -> TACTGACATTGGTAATTT
        rev_comp = str(seq_record.reverse_complement().complement())
    # BLAST outputs only need to be reverse complemented
    # e.g. seq_string: ACATCAAAAGTATATCGTTC -> GAACGATATACTTTTGATGT
    else:
        rev_comp = str(seq_record.reverse_complement())
    return rev_comp


def ipcress_mismatches(metadata, analysistype, iupac):
    """
    Find the location(s) within the primer of any mismatches
    :param metadata: List of metadata objects for all samples
    :param analysistype: String of the current analysis type
    :param iupac: Dictionary of all IUPAC codes
    :return: metadata: Updated list of metadata objects
    """
    for sample in metadata:
        # Iterate through all the ipcress experiments parsed above
        for experiment in sample[analysistype].results.datastore:
            for contig in sample[analysistype].results[experiment].datastore:
                # Calculate the total number of mismatches
                sample[analysistype].results[experiment][contig].total_mismatch = \
                    int(sample[analysistype].results[experiment][contig].forward_mismatch) + \
                    int(sample[analysistype].results[experiment][contig].reverse_mismatch)
                # Determine the range of the amplicon - make sure that the direction is accounted for
                if sample[analysistype].results[experiment][contig].direction == 'forward':
                    sample[analysistype].results[experiment][contig].amplicon_range = \
                        range(int(sample[analysistype].results[experiment][contig].forward_pos),
                              int(sample[analysistype].results[experiment][contig].reverse_pos))
                else:
                    sample[analysistype].results[experiment][contig].amplicon_range = \
                        range(int(sample[analysistype].results[experiment][contig].reverse_pos),
                              int(sample[analysistype].results[experiment][contig].forward_pos))
                # Determine if the forward primer has any mismatches
                if int(sample[analysistype].results[experiment][contig].forward_mismatch) > 0:
                    # Calculate details of the mismatches
                    sample[analysistype].results[experiment][contig].forward_mismatch_details = \
                        determine_mismatch_locations(
                            ref_primer=sample[analysistype].results[experiment][contig].forward_ref,
                            query_primer=sample[analysistype].results[experiment][contig].forward_query,
                            iupac=iupac)
                # If no mismatches, initialise the attributes as an empty string
                else:
                    sample[analysistype].results[experiment][contig].forward_mismatch_details = str()
                # Mismatches in reverse primer
                if int(sample[analysistype].results[experiment][contig].reverse_mismatch) > 0:
                    sample[analysistype].results[experiment][contig].reverse_mismatch_details = \
                        determine_mismatch_locations(
                            ref_primer=sample[analysistype].results[experiment][contig].reverse_ref,
                            query_primer=sample[analysistype].results[experiment][contig].reverse_query,
                            iupac=iupac)
                else:
                    sample[analysistype].results[experiment][contig].reverse_mismatch_details = str()
    return metadata


def determine_mismatch_locations(ref_primer, query_primer, iupac):
    """
    Find the position of mismatches within the primer
    :param ref_primer: String of the primer sequence
    :param query_primer: String of the query sequence
    :param iupac: Dictionary of IUPAC codes
    :return: mismatch_string: Formatted string of mismatches
    """
    # Initialise a string to hold mismatch details
    mismatch_string = str()
    # Iterate through every base of the primer sequence
    for pos, ref_base in enumerate(ref_primer.upper()):
        # Use the iterator to extract the corresponding base of the query sequence
        query_base = query_primer[pos]
        # Determine if the ref base and the query base do not match. If the primer has an IUPAC degenerate base e.g.
        # 'B': ['C', 'G', 'T'], check to see if the query base is coded by the degenerate base
        if query_base != ref_base and query_base not in iupac[ref_base]:
            # If the mismatch string already exists, add a semicolon delimiter
            if mismatch_string:
                mismatch_string += ';'
            # Add the current position (+1 due to 0-based indexing) plus the details of the mismatch
            # e.g. 12T>G indicates that at position 12 of the query sequence, a G is present instead of the G in the
            # primer sequence
            mismatch_string += '{pos}{ref}>{query}'.format(pos=pos + 1,
                                                           ref=ref_base,
                                                           query=query_base)
    # Return the formatted mismatch string
    return mismatch_string


def best_hit(metadata, analysistype, range_buffer=0):
    """
    Filter and return the best hit for a region
    :param metadata: List of metadata objects for all samples
    :param analysistype: String of current analysis type
    :param range_buffer: Integer of number of bases to add to the range when finding overlaps - default is 0 (ranges
    must actually overlap), increasing the buffer will allow for the overlaps to increase
    :return: Updated list of metadata objects
    """
    for sample in metadata:
        # Initialise a dictionary to store amplicon range: additional details
        overlaps = dict()
        # Initialise the .best_hits attribute as required
        if not GenObject.isattr(sample[analysistype], 'best_hits'):
            sample[analysistype].best_hits = dict()
        for experiment in sample[analysistype].results.datastore:
            for contig in sample[analysistype].results[experiment].datastore:
                # Ensure that the amplicon_range attribute exists - happens when the method is used for singleton BLAST
                # analyses
                try:
                    amplicon_range = sample[analysistype].results[experiment][contig].amplicon_range
                    # Ensure that the dictionary has been created
                    if overlaps:
                        # Create a boolean to track if the current primer set overlaps with any other in the dictionary
                        overlap = False
                        # Iterate through all the previous ranges: dictionaries in the overlaps dictionary
                        try:
                            for ranges, overlap_dict in overlaps[contig].items():
                                # Increase the range by the provided buffer value
                                buffered_min = amplicon_range[0] - range_buffer if amplicon_range[0] - \
                                                                                   range_buffer >= 0 else 0
                                buffered_max = amplicon_range[-1] + range_buffer
                                # Check if the current range overlaps with a previously stored range
                                if range(max(buffered_min, ranges[0]),
                                         min(buffered_max, ranges[-1]) + 1):
                                    # Set the overlap boolean to True
                                    overlap = True
                                    # If the current primer set has fewer mismatches than the previous best primer set
                                    # for this range, replace the details in the overlap dictionary with the details
                                    # from the current set
                                    if sample[analysistype].results[experiment][contig].total_mismatch < \
                                            overlap_dict['mismatches']:
                                        overlap_dict['mismatches'] = sample[analysistype].results[experiment][contig].\
                                            total_mismatch
                                        overlap_dict['experiment'] = experiment
                            # If the current range does not overlap with any previous primer sets, add this range and
                            # details to the dictionary
                            if not overlap:
                                overlaps[contig][amplicon_range] = {
                                    'mismatches': sample[analysistype].results[experiment][contig].total_mismatch,
                                    'experiment': experiment
                                }
                        except KeyError:
                            overlaps[contig] = dict()
                            # Set the amplicon range as the key for an additional dictionary containing total
                            # mismatches and experiment name
                            overlaps[contig][amplicon_range] = {
                                'mismatches': sample[analysistype].results[experiment][contig].total_mismatch,
                                'experiment': experiment
                            }
                    # If this is the first sample, populate the dictionary
                    else:
                        overlaps[contig] = dict()
                        # Set the amplicon range as the key for an additional dictionary containing total mismatches and
                        # experiment name
                        overlaps[contig][amplicon_range] = {
                            'mismatches': sample[analysistype].results[experiment][contig].total_mismatch,
                            'experiment': experiment
                        }
                    # Create the best_hits attribute using the overlaps dictionary
                    sample[analysistype].best_hits = overlaps
                except AttributeError:
                    pass
    return metadata


def empty_results(metadata, analysistype):
    """
    Determine if any samples did not have any hits
    :param metadata: List of metadata objects for all samples
    :param analysistype: String of current analysis type
    :return: empty: Boolean of whether samples did not have hits
    """
    # Initialise the boolean to be False
    empty = False
    # Iterate through all the samples
    for sample in metadata:
        # Check to see if the results attribute is empty
        if not sample[analysistype].results.datastore:
            empty = True
    return empty


def create_fasta_format(forward_dict, reverse_dict, fastaprimerfile):
    """
    Create a FASTA-formatted primer file to allow for BLAST analyses
    :param forward_dict: Dictionary of primer name: forward primer sequence
    :param reverse_dict: Dictionary of primer name: reverse primer sequence
    :param fastaprimerfile: Name and absolute path of FASTA-formatted primer file to create
    """
    logging.info('Creating FASTA-formatted primer file')
    records = list()
    with open(fastaprimerfile, 'w') as formatted:
        # Iterate through all the targets
        for primer_name in sorted(forward_dict):
            # Create SeqRecords for the forward and reverse primers
            records.append(SeqRecord(seq=Seq(forward_dict[primer_name]),
                                     id='{pn}-{direction}'
                                     .format(pn=primer_name, direction='F'),
                                     name='',
                                     description=''))
            records.append(SeqRecord(seq=Seq(reverse_dict[primer_name]),
                                     id='{pn}-{direction}'
                                     .format(pn=primer_name, direction='R'),
                                     name='',
                                     description=''))
        SeqIO.write(records, formatted, 'fasta')


def make_blastdb(formattedprimers):
    """
    Create a BLAST database of the primer file
    """
    logging.info('Creating BLAST database file from FASTA-formatted primers')
    # Remove the path and the file extension for easier future globbing
    db = os.path.splitext(formattedprimers)[0]
    nhr = '{db}.nhr'.format(db=db)  # add nhr for searching
    if not os.path.isfile(str(nhr)):
        # Create the databases
        command = 'makeblastdb -in {primerfile} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {outfile}'\
            .format(primerfile=formattedprimers,
                    outfile=db)
        run_subprocess(command)


def run_blast(metadata, analysistype, formattedprimers, blastheader, threads):
    """
    Run BLASTn analyses of the query file against the primer database
    :param metadata: List of metadata objects for all samples
    :param analysistype: String of the current analysis type
    :param formattedprimers: String of name and absolute path to formatted primer file
    :param blastheader: List of all column headers used in the BLAST analyses
    :param threads: Integer of the number of threads for BLAST to use
    :return:
    """
    logging.info('Running BLAST analyses')
    for sample in metadata:
        sample[analysistype].blastresults = '{of}_blast_results.tsv' \
                .format(of=os.path.join(sample.general.outputdirectory, sample.name))
        # If a report was created, but no results entered - program crashed, or no sequences passed thresholds,
        # remove the report, and run the blast analyses again
        try:
            size = os.path.getsize(sample[analysistype].blastresults)
            if size == 0:
                os.remove(sample[analysistype].blastresults)
        except FileNotFoundError:
            pass
        # Check to see if the results attribute is empty
        if not sample[analysistype].results.datastore and not os.path.isfile(sample[analysistype].blastresults):
            db = os.path.splitext(formattedprimers)[0]
            # BLAST command line call. Note the high number of alignments.
            # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
            # alignments are reported. Also note the custom outfmt. Using very permissive BLAST settings (word_size: 4,
            # dust:'no', penalty: -1, and task: 'blastn-short' to allow BLAST to return results even with the short
            # primer input sequences)
            blastn = NcbiblastnCommandline(query=sample.general.bestassemblyfile,
                                           db=db,
                                           evalue=1,
                                           task='blastn-short',
                                           word_size=4,
                                           dust='no',
                                           penalty=-2,
                                           reward=3,
                                           gapopen=5,
                                           gapextend=5,
                                           num_alignments=1000000,
                                           num_threads=threads,
                                           outfmt='6 qseqid sseqid positive mismatch gaps evalue bitscore slen '
                                                  'length qstart qend qseq sstart send sseq',
                                           out=sample[analysistype].blastresults)
            # Save the blast command in the metadata
            sample[analysistype].blastcommand = str(blastn)
            if sample.general.bestassemblyfile != 'NA':
                # Run the blastn command
                blastn()
                # Add a header to the report
                with open(sample[analysistype].blastresults, 'r+') as f:
                    # Read in the information from the blastresults file
                    content = f.read()
                    # Go back to the start of the file
                    f.seek(0, 0)
                    # Write the formatted header (\n) followed by the content to the file
                    f.write('\t'.join(blastheader) + '\n' + content)
    return metadata


def parse_blast(metadata, analysistype, fieldnames, forward_dict, reverse_dict, mismatches, min_hits=2):
    """
    Parse the BLAST outputs, and populate metadata objects in a similar fashion to the ipcress-processed samples
    :param metadata: List of metadata objects of all samples
    :param analysistype: String of current analysis type
    :param fieldnames: List of BLAST headers used in the analysis
    :param forward_dict: Dictionary of primer name: forward primer sequence
    :param reverse_dict: Dictionary of primer name: reverse primer sequences
    :param mismatches: Integer of the maximum number of mismatches allowed per primer
    :param min_hits: Integer of the minimum number of primers (out of forward and reverse) that require hits. Useful
    if searching for singletons
    :return: List of updated metadata objects
    """
    logging.info('Parsing BLAST outputs')
    for sample in metadata:
        # Only process the samples that were analysed with BLAST
        if os.path.isfile(sample[analysistype].blastresults):
            # Open blast output csv file
            csvfile = open(sample[analysistype].blastresults)
            # Skip header
            csvfile.readline()
            # Open the sequence profile file as a dictionary
            blastdict = DictReader(csvfile, fieldnames=fieldnames, dialect='excel-tab')
            hit_dict = dict()
            # Go through each BLAST result
            for row in blastdict:
                # Extract the primer name and the direction (F/R) from the subject_id
                # e.g. Lin-R yields a primer_name of Lin, and a direction of R
                primer_name = row['subject_id'][:-2]
                direction = row['subject_id'][-1]
                positives = int(row['positives'])
                contig = row['query_id']
                # Initialise the primer name in the dictionary as required
                if primer_name not in hit_dict:
                    hit_dict[primer_name] = dict()
                if contig not in hit_dict[primer_name]:
                    hit_dict[primer_name][contig] = dict()
                # Populate the hit_dict for the forward primer
                if direction == 'F':
                    # Find hits that do not exceed the maximum number of mismatches; calculated by subtracting the
                    # number of positives from the total length of the primer sequence
                    if positives >= (len(forward_dict[primer_name]) - mismatches):
                        # Initialise the 'forward' key in the dictionary as required
                        if 'forward' not in hit_dict[primer_name]:
                            hit_dict[primer_name][contig]['forward'] = dict()
                        # If the dictionary is still empty, add the number of positives: blast details
                        if not hit_dict[primer_name][contig]['forward']:
                            hit_dict[primer_name][contig]['forward'][positives] = row
                        # Extract the number of positives from the previous best result in the dictionary
                        previous_best = list(hit_dict[primer_name][contig]['forward'].keys())[0]
                        # Check if the current number of positives is better than the previous best
                        if positives > previous_best:
                            # Add the current positives: blast details to the dictionary
                            hit_dict[primer_name][contig]['forward'][positives] = row
                            # Delete the previous best positives key from the dictionary
                            del(hit_dict[primer_name][contig]['forward'][previous_best])
                # Same as above, but with the reverse primer
                else:
                    if positives >= (len(reverse_dict[primer_name]) - mismatches):
                        if 'reverse' not in hit_dict[primer_name][contig]:
                            hit_dict[primer_name][contig]['reverse'] = dict()
                        if not hit_dict[primer_name][contig]['reverse']:
                            hit_dict[primer_name][contig]['reverse'][positives] = row
                        previous_best = list(hit_dict[primer_name][contig]['reverse'].keys())[0]
                        if positives > previous_best:
                            hit_dict[primer_name][contig]['reverse'][positives] = row
                            del (hit_dict[primer_name][contig]['reverse'][previous_best])
            # Iterate through hit_dict
            for experiment, contig_dict in hit_dict.items():
                for contig, detail_dict in contig_dict.items():
                    # Only continue if the detail_dict has a length of two (positive hits for both forward and reverse
                    # primers), unless searching for singletons
                    if len(detail_dict) >= min_hits:
                        # Initialise the experiment GenObject
                        if not GenObject.isattr(sample[analysistype].results, experiment):
                            sample[analysistype].results[experiment] = GenObject()
                        # Forward primer
                        if 'forward' in detail_dict:
                            for positives, forward_details in detail_dict['forward'].items():
                                if not GenObject.isattr(sample[analysistype].results[experiment], contig):
                                    sample[analysistype].results[experiment][contig] = GenObject()
                                # Check to see if the .contig attribute already exists
                                try:
                                    if sample[analysistype].results[experiment][contig].contig:
                                        # Append a semicolon plus 'forward_' and the extracted query_id
                                        sample[analysistype].results[experiment][contig].contig += \
                                            ';forward_' + forward_details['query_id']
                                    else:
                                        sample[analysistype].results[experiment][contig].contig = \
                                            'forward_' + forward_details['query_id']
                                # Otherwise create and populate the attribute with 'forward_' and the extracted query_id
                                except AttributeError:
                                    sample[analysistype].results[experiment][contig].contig = \
                                        'forward_' + forward_details['query_id']
                                sample[analysistype].results[experiment][contig].primer_set = experiment
                                # Set the amplicon length to zero
                                sample[analysistype].results[experiment][contig].amplicon_length = 0
                                # Forward primers are referred to as 'A' in ipcress outputs
                                sample[analysistype].results[experiment][contig].forward_primer = 'A'
                                # Create an attribute to store the range over which the primer had hits
                                # e.g. 4-20 indicates that the first three bases were not covered
                                sample[analysistype].results[experiment][contig].forward_pos_range = [
                                    min(int(forward_details['subject_start']),
                                        int(forward_details['subject_end'])),
                                    max(int(forward_details['subject_start']),
                                        int(forward_details['subject_end'])),
                                    ]
                                # Create an attribute to store the range over which the primer hit the query sequence
                                sample[analysistype].results[experiment][contig].forward_range = \
                                    range(
                                        min(int(forward_details['query_start']),
                                            int(forward_details['query_end'])),
                                        max(int(forward_details['query_start']),
                                            int(forward_details['query_end']))
                                      )
                                # Calculate the number of mismatches by subtracting the total number of positives from
                                # the total length of the primer sequence
                                sample[analysistype].results[experiment][contig].forward_mismatch = \
                                    len(forward_dict[experiment]) - int(forward_details['positives'])
                                # Determine if the hit was on the forward or reverse strand. Hits on the forward strand
                                # will have a subject_start before the subject_end (opposite for hits on the
                                # rev strand) e.g. start: 4, end 20 for hits on the forward strand
                                direction = 'forward' if int(forward_details['subject_start']) < int(
                                    forward_details['subject_end']) else 'revcomp'
                                sample[analysistype].results[experiment][contig].forward_direction = direction
                                # Either create or update the .direction attribute with the calculated direction
                                try:
                                    if sample[analysistype].results[experiment][contig].direction:
                                        sample[analysistype].results[experiment][contig].direction += \
                                            ';forward_' + direction
                                    else:
                                        sample[analysistype].results[experiment][contig].direction = \
                                            'forward_' + direction
                                except AttributeError:
                                    sample[analysistype].results[experiment][contig].direction = \
                                        'forward_' + direction
                                # The forward_ref attribute is the sequence of the forward primer extracted from
                                # forward_dict
                                sample[analysistype].results[experiment][contig].forward_ref = forward_dict[experiment]
                                # Extract the sequence of the query. Use the extracted 'query_sequence' if the hit is
                                # on the forward strand, otherwise, calculate the reverse complement
                                sample[analysistype].results[experiment][contig].forward_query = \
                                    forward_details['query_sequence'] if \
                                    direction == 'forward' else \
                                    revcomp(seq_string=forward_details['query_sequence'],
                                            blast_primers=True)

                        # Reverse primer - same as above
                        if 'reverse' in detail_dict:
                            for positives, reverse_details in detail_dict['reverse'].items():
                                if not GenObject.isattr(sample[analysistype].results[experiment], contig):
                                    sample[analysistype].results[experiment][contig] = GenObject()
                                try:
                                    if sample[analysistype].results[experiment][contig].contig:
                                        sample[analysistype].results[experiment][contig].contig += \
                                            ';reverse_' + reverse_details['query_id']
                                    else:
                                        sample[analysistype].results[experiment][contig].contig = \
                                            'reverse_' + reverse_details['query_id']
                                except AttributeError:
                                    sample[analysistype].results[experiment][contig].contig = \
                                        'reverse_' + reverse_details['query_id']
                                sample[analysistype].results[experiment][contig].primer_set = experiment
                                sample[analysistype].results[experiment][contig].reverse_primer = 'B'
                                sample[analysistype].results[experiment][contig].reverse_pos_range = [
                                    min(int(reverse_details['subject_start']),
                                        int(reverse_details['subject_end'])),
                                    max(int(reverse_details['subject_start']),
                                        int(reverse_details['subject_end'])),
                                ]
                                sample[analysistype].results[experiment][contig].reverse_range = \
                                    range(
                                        min(int(reverse_details['query_start']),
                                            int(reverse_details['query_end'])),
                                        max(int(reverse_details['query_start']),
                                            int(reverse_details['query_end']))
                                    )
                                sample[analysistype].results[experiment][contig].reverse_mismatch = \
                                    len(reverse_dict[experiment]) - int(reverse_details['positives'])
                                direction = 'forward' if int(reverse_details['subject_start']) < int(
                                    reverse_details['subject_end']) else 'revcomp'
                                sample[analysistype].results[experiment][contig].reverse_direction = direction
                                try:
                                    if sample[analysistype].results[experiment][contig].direction:
                                        sample[analysistype].results[experiment][contig].direction \
                                            += ';reverse_' + direction
                                    else:
                                        sample[analysistype].results[experiment][contig].direction \
                                            = 'reverse_' + direction
                                except AttributeError:
                                    sample[analysistype].results[experiment][contig].direction \
                                        = 'reverse_' + direction
                                sample[analysistype].results[experiment][contig].reverse_ref = reverse_dict[experiment]
                                sample[analysistype].results[experiment][contig].reverse_query = \
                                    reverse_details['query_sequence'] if \
                                    direction == 'forward' else \
                                    revcomp(seq_string=reverse_details['query_sequence'],
                                            blast_primers=True)
    return metadata


def blast_primer_mismatch_details(metadata, analysistype, iupac):
    """
    Determine mismatch details for BLAST analyses
    :param metadata: List of metadata objects for all samples
    :param analysistype: String of current analysis type
    :param iupac: Dictionary of IUPAC base: corresponding list of bases
    :return: List of updated metadata objects
    """
    for sample in metadata:
        # Only process samples that were investigated with BLAST
        if os.path.isfile(sample[analysistype].blastresults):
            for experiment in sample[analysistype].results.datastore:
                for contig in sample[analysistype].results[experiment].datastore:
                    # Only samples with the .contig attribute have BLAST results for both forward and reverse primers
                    try:

                        # Calculate the mismatch details of the forward primer
                        try:
                            sample[analysistype].results[experiment][contig].forward_mismatch_details, \
                                sample[analysistype].results[experiment][contig].forward_query = \
                                blast_mismatches(ref_primer=sample[analysistype].results[experiment][contig]
                                                 .forward_ref,
                                                 query_primer=sample[analysistype].results[experiment][contig]
                                                 .forward_query,
                                                 primer_pos_range=sample[analysistype].results[experiment][contig]
                                                 .forward_pos_range,
                                                 query_pos_range=sample[analysistype].results[experiment][contig]
                                                 .forward_range,
                                                 iupac=iupac,
                                                 sample=sample,
                                                 contig=contig,
                                                 direction=sample[analysistype].results[experiment][contig]
                                                 .forward_direction)
                        except AttributeError:
                            pass
                        # Reverse primer
                        try:
                            sample[analysistype].results[experiment][contig].reverse_mismatch_details, \
                                sample[analysistype].results[experiment][contig].reverse_query = \
                                blast_mismatches(ref_primer=sample[analysistype].results[experiment][contig]
                                                 .reverse_ref,
                                                 query_primer=sample[analysistype].results[experiment][contig]
                                                 .reverse_query,
                                                 primer_pos_range=sample[analysistype].results[experiment][contig]
                                                 .reverse_pos_range,
                                                 query_pos_range=sample[analysistype].results[experiment][contig]
                                                 .reverse_range,
                                                 iupac=iupac,
                                                 sample=sample,
                                                 contig=contig,
                                                 direction=sample[analysistype].results[experiment][contig]
                                                 .reverse_direction)
                        except AttributeError:
                            pass
                    except AttributeError:
                        pass
    return metadata


def blast_mismatches(ref_primer, query_primer, primer_pos_range, query_pos_range, iupac, sample, contig, direction):
    """
    Calculate mismatch details of BLAST hits
    :param ref_primer: String of primer sequence
    :param query_primer: String of query sequence hitting the primer sequence
    :param primer_pos_range: Range of bases in the primer covered by the hit to the query sequence
    :param query_pos_range: Range of bases in the query sequence
    :param iupac: Dictionary of IUPAC base: corresponding list of bases
    :param sample: Metadata object of current sample
    :param contig: Name of contig with the hit
    :param direction: Direction of the hit
    :return: mismatch_string: String of the mismatch details
    :return: query_primer: Query primer padded with '-' to be the same length as the ref_primer
    """
    front_padding = 0
    end_padding = 0
    # Determine if the BLAST hit starts at the beginning of the ref_primer
    if not primer_pos_range[0] == 1:
        # Pad the start of query_primer with the number of missing bases
        front_padding = primer_pos_range[0] - 1
    # Check if the query primer is shorter than the reference primer
    if len(query_primer) < len(ref_primer):
        # Pad the end of the primer with the number of missing bases
        end_padding = len(ref_primer) - len(query_primer) - front_padding
    query_primer = extract_sequence(sample=sample,
                                    contig=contig,
                                    query_pos_range=query_pos_range,
                                    front_padding=front_padding,
                                    end_padding=end_padding,
                                    direction=direction)
    # Calculate the mismatch details
    mismatch_string = determine_mismatch_locations(ref_primer=ref_primer,
                                                   query_primer=query_primer,
                                                   iupac=iupac)
    return mismatch_string, query_primer


def extract_sequence(sample, contig, query_pos_range, front_padding, end_padding, direction):
    """
    Extract the sequence from the query genome using the range of the hit
    :param sample: Metadata object for the current query
    :param contig: Name of contig with the hit
    :param query_pos_range: The range of the hit on the contig
    :param front_padding: Integer of the number of bases missing from the 5' end of the match
    :param end_padding: Integer of the number of bases missing from the 3' end of the match
    :param direction: Direction of the hit
    :return: query_sequence: Extracted query sequence
    """
    # Initialise the string of the query sequence
    query_sequence = str()
    # Read in the query genome using SeqIO
    for record in SeqIO.parse(sample.general.bestassemblyfile, 'fasta'):
        # Ensure that the contig is the one with the hit
        if record.id == contig:
            # Set the start and stop positions of the amplicon sequence depending on the direction
            if direction == 'revcomp':
                # Reverse hits are the beginning of the range minus the number of missing bases to pad the end (minus
                # one) due to zero-based indexing
                start = query_pos_range[0] - end_padding - 1
                # The end position is the the final position of the range plus the number of missing bases at the front
                end = query_pos_range[-1] + front_padding + 1
            else:
                start = query_pos_range[0] - front_padding - 1
                end = query_pos_range[-1] + end_padding + 1
            # Extract the query sequence
            query_sequence = str(record.seq)[start: end]
            # Use the reverse complement of the sequence if the direction is reverse
            if direction == 'revcomp':
                seq_o = Seq(query_sequence)
                query_sequence = seq_o.reverse_complement()
    # Ensure that the sequence is in uppercase
    query_sequence = query_sequence.upper()
    return query_sequence


def amplicon_write(metadata, analysistype, reportpath,):
    """
    Write the amplicons to file
    :param metadata: List of metadata objects
    :param analysistype: Sting of current analysis type
    :param reportpath: Name and absolute path in which the alleles are to be written
    """
    logging.info('Creating amplicon files')
    with open(os.path.join(reportpath, 'alleles.fasta'), 'w') as alleles:
        for sample in metadata:
            for experiment in sample[analysistype].results.datastore:
                for contig in sample[analysistype].results[experiment].datastore:
                    # Ensure that the sequence attribute exists
                    try:
                        # Create a SeqRecord from the sequence
                        seq_record = SeqRecord(seq=Seq(sample[analysistype].results[experiment][contig].sequence),
                                               id='{sn}_{header}'
                                               .format(sn=sample.name,
                                                       header=sample[analysistype].results[experiment][contig].header),
                                               name='',
                                               description='')
                        # Use SeqIO to write the SeqRecord to file
                        SeqIO.write(seq_record, alleles, 'fasta')
                    except AttributeError:
                        pass


def best_hits_singletons(metadata, analysistype):
    """
    Find the best hits by single primers, and remove duplicates
    :param metadata: List of metadata objects for all samples
    :param analysistype: String of current analysis type
    :return: Update list of metadata objects
    """
    for sample in metadata:
        # If the best_hits dictionary is empty, there were no amplicons produced by any of the primer sets for this
        # sample
        if not sample[analysistype].best_hits:
            # Initialise a dictionary to store the hit information
            overlaps = dict()
            for experiment in sample[analysistype].results.datastore:
                for contig in sample[analysistype].results[experiment].datastore:
                    # Process the forward and reverse primers to allow for the attribute names to not exist
                    if sample[analysistype].results[experiment][contig].direction == 'forward':
                        # As the total_mismatch attribute was initialised as an empty string, typecasting it to an
                        # int will cause a ValueError to be raised
                        try:
                            # Use the .forward_mismatch attribute to set the .total_match
                            sample[analysistype].results[experiment][contig].total_mismatch = \
                                int(sample[analysistype].results[experiment][contig].forward_mismatch)
                            # Use the .forward_range to set the amplicon_range
                            sample[analysistype].results[experiment][contig].amplicon_range = \
                                sample[analysistype].results[experiment][contig].forward_range
                            # Update the overlap dictionary, and the overlap boolean
                            overlaps = \
                                best_matcher(overlaps=overlaps,
                                             amplicon_range=sample[analysistype].results[experiment][
                                                      contig].amplicon_range,
                                             total_mismatch=sample[analysistype].results[
                                                      experiment][contig].total_mismatch,
                                             experiment=experiment,
                                             contig=contig)
                        except ValueError:
                            pass
                    # Reverse primer - same as above, but with the appropriately named attributes
                    else:
                        try:
                            sample[analysistype].results[experiment][contig].total_mismatch = \
                                int(sample[analysistype].results[experiment][contig].reverse_mismatch)
                            sample[analysistype].results[experiment][contig].amplicon_range = \
                                sample[analysistype].results[experiment][contig].reverse_range
                            overlaps = \
                                best_matcher(overlaps=overlaps,
                                             amplicon_range=sample[analysistype].results[experiment][
                                                      contig].amplicon_range,
                                             total_mismatch=sample[analysistype].results[
                                                      experiment][contig].total_mismatch,
                                             experiment=experiment,
                                             contig=contig)
                        except ValueError:
                            pass
                    # Update the best_hits dictionary if overlaps exists
                    if overlaps:
                        sample[analysistype].best_hits = overlaps
    return metadata


def best_matcher(overlaps, amplicon_range, total_mismatch, experiment, contig):
    """
    Determine the best primer hit to a contig
    :param overlaps: Dictionary of details on the best primer hit in a region
    :param amplicon_range: Range of the current primer hit
    :param total_mismatch: Integer of the total number of mismatches between the primer and the query
    :param experiment: String of the name of the primer set
    :param contig: String of the name of current contig
    :return:
    """
    # Create a boolean to store whether an overlap is encountered between the current primer hit range and the
    # previous best primer hit range
    overlap = False
    # Ensure that the dictionary has been created
    if overlaps:
        try:
            # Iterate through all the previous ranges: dictionaries in the overlaps dictionary
            for ranges, overlap_dict in overlaps[contig].items():
                # Check if the current range overlaps with a previously stored range
                if range(max(amplicon_range[0], ranges[0]),
                         min(amplicon_range[-1], ranges[-1]) + 1):
                    # Set the overlap boolean to True
                    overlap = True
                    # If the current primer set has fewer mismatches than the previous best primer set for this
                    # range, replace the details in the overlap dictionary with the details from the current set
                    if total_mismatch < overlap_dict['mismatches']:
                        overlap_dict['mismatches'] = total_mismatch
                        overlap_dict['experiment'] = experiment
            # If the current range does not overlap with any previous primer sets, add this range and details
            # to the dictionary
            if not overlap:
                overlaps[contig][amplicon_range] = {
                    'mismatches': total_mismatch,
                    'experiment': experiment
                }
        except KeyError:
            overlaps[contig] = dict()
            # Set the amplicon range as the key for an additional dictionary containing total mismatches and
            # experiment name
            overlaps[contig][amplicon_range] = {
                'mismatches': total_mismatch,
                'experiment': experiment
            }
    # If this is the first sample, populate the dictionary
    else:
        # Set the amplicon range as the key for an additional dictionary containing total mismatches and
        # experiment name
        overlaps[contig] = dict()
        # Set the amplicon range as the key for an additional dictionary containing total mismatches and
        # experiment name
        overlaps[contig][amplicon_range] = {
            'mismatches': total_mismatch,
            'experiment': experiment
        }
    # Create the best_hits attribute using the overlaps dictionary
    return overlaps


class VtyperIP(object):

    def vtyper(self):
        self.forward_dict, self.reverse_dict = epcr_primers(primerfile=self.primerfile,
                                                            forward_dict=self.forward_dict,
                                                            reverse_dict=self.reverse_dict,
                                                            formattedprimers=self.formattedprimers,
                                                            copyfile=False)
        if not os.path.isfile(self.formattedprimers):
            epcr_primer_file(forward_dict=self.forward_dict,
                             reverse_dict=self.reverse_dict,
                             formattedprimers=self.formattedprimers,
                             min_amplicon_size=0,
                             max_amplicon_size=1500)
        self.metadata = ipcress_threads(metadata=self.metadata,
                                        analysistype=self.analysistype,
                                        formattedprimers=self.formattedprimers,
                                        mismatches=self.mismatches,
                                        ipcress_queue=self.ipcress_queue,
                                        threads=self.threads)
        self.metadata = ipcress_parse(metadata=self.metadata,
                                      analysistype=self.analysistype)
        self.metadata = ipcress_mismatches(metadata=self.metadata,
                                           analysistype=self.analysistype,
                                           iupac=self.iupac)
        self.metadata = best_hit(metadata=self.metadata,
                                 analysistype=self.analysistype,
                                 range_buffer=self.range_buffer)
        empty = empty_results(metadata=self.metadata,
                              analysistype=self.analysistype)
        if empty:
            create_fasta_format(forward_dict=self.forward_dict,
                                reverse_dict=self.reverse_dict,
                                fastaprimerfile=self.fastaprimers)
        # Create a BLAST database from the primer file
        make_blastdb(formattedprimers=self.fastaprimers)
        self.metadata = run_blast(metadata=self.metadata,
                                  analysistype=self.analysistype,
                                  formattedprimers=self.fastaprimers,
                                  blastheader=self.fieldnames,
                                  threads=self.threads)
        self.metadata = parse_blast(metadata=self.metadata,
                                    analysistype=self.analysistype,
                                    fieldnames=self.fieldnames,
                                    forward_dict=self.forward_dict,
                                    reverse_dict=self.reverse_dict,
                                    mismatches=self.mismatches + 1)
        self.metadata = blast_primer_mismatch_details(metadata=self.metadata,
                                                      analysistype=self.analysistype,
                                                      iupac=self.iupac)
        self.metadata = best_hits_singletons(metadata=self.metadata,
                                             analysistype=self.analysistype)
        self.toxin_profile()
        self.vtyper_report()

    def toxin_profile(self):
        """
        Create a sample-specific list of all toxins
        """
        for sample in self.metadata:
            sample[self.analysistype].toxinprofile = 'ND'
            # Create a set to store all the unique results
            toxin_set = set()
            for contig, amplicon_dict in sample[self.analysistype].best_hits.items():
                # Iterate over all the best primer hits for a range
                for amplicon_range, hit_dict in amplicon_dict.items():
                    # Remove the primer numbering information, and add the cleaned toxin to the set
                    # e.g. vtx1a_0_0;vtx1c_0_0;vtx2d_1_2;vtx2d_1_3 becomes {vtx1a, vtx1c, vtx2d}
                    if not os.path.isfile(sample[self.analysistype].blastresults):
                        toxin_set.add(hit_dict['experiment'].split('_')[0])
                    # Add an asterisk if the hit was derived from BLAST-based analyses
                    else:
                        toxin_set.add(hit_dict['experiment'].split('_')[0] + '*')
                # Create a string of the entries in the sorted list of toxins joined with ";"
                sample[self.analysistype].toxinprofile = ";".join(sorted(list(toxin_set))) if toxin_set else 'ND'

    def vtyper_report(self):
        """
        Create a report of the ePCR-calculated toxin profiles
        """
        logging.info('Creating {at} report'.format(at=self.analysistype))
        with open(os.path.join(self.reportpath, '{at}.csv'.format(at=self.analysistype)), 'w') as report:
            data = 'Strain,ToxinProfile\n'
            for sample in self.metadata:
                data += '{sn},{tp}\n'.format(sn=sample.name,
                                             tp=sample[self.analysistype].toxinprofile)
            # Write the data to the report
            report.write(data)

    def __init__(self, metadataobject, analysistype, reportpath, mismatches=3):
        self.metadata = metadataobject
        self.analysistype = analysistype
        self.reportpath = reportpath
        self.mismatches = mismatches
        make_path(self.reportpath)
        # Extract the path of the current script from the full path + file name
        self.homepath = os.path.split(os.path.abspath(__file__))[0]
        self.primerfile = os.path.join(self.homepath, 'ssi_subtyping_primers.txt')
        self.formattedprimers = os.path.join(self.homepath, 'ipcress_format_subtyping_primers_degenerate.txt')
        self.fastaprimers = os.path.join(os.path.dirname(self.primerfile), 'formattedprimers.fa')
        self.forward_dict = dict()
        self.reverse_dict = dict()
        self.range_buffer = 150
        self.threads = multiprocessing.cpu_count() - 1
        self.ipcress_queue = Queue()
        # Dictionary of IUPAC codes
        self.iupac = {
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'S': ['C', 'G'],
            'W': ['A', 'T'],
            'K': ['G', 'T'],
            'M': ['A', 'C'],
            'B': ['C', 'G', 'T'],
            'D': ['A', 'G', 'T'],
            'H': ['A', 'C', 'T'],
            'V': ['A', 'C', 'G'],
            'N': ['A', 'C', 'G', 'T'],
            'A': ['A'],
            'C': ['C'],
            'G': ['G'],
            'T': ['T']
        }
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'subject_length', 'alignment_length',
                           'query_start', 'query_end', 'query_sequence',
                           'subject_start', 'subject_end', 'subject_sequence']


class CustomIP(object):

    def main(self):
        if self.primer_format == 'fasta':
            self.forward_dict, self.reverse_dict = fasta_primers(primerfile=self.primerfile,
                                                                 forward_dict=self.forward_dict,
                                                                 reverse_dict=self.reverse_dict)
            epcr_primer_file(forward_dict=self.forward_dict,
                             reverse_dict=self.reverse_dict,
                             formattedprimers=self.formattedprimers,
                             min_amplicon_size=self.min_amplicon_size,
                             max_amplicon_size=self.max_amplicon_size)
        else:
            self.forward_dict, self.reverse_dict = epcr_primers(primerfile=self.primerfile,
                                                                forward_dict=self.forward_dict,
                                                                reverse_dict=self.reverse_dict,
                                                                formattedprimers=self.formattedprimers)
        self.metadata = ipcress_threads(metadata=self.metadata,
                                        analysistype=self.analysistype,
                                        formattedprimers=self.formattedprimers,
                                        mismatches=self.mismatches,
                                        ipcress_queue=self.ipcress_queue,
                                        threads=self.threads)
        self.metadata = ipcress_parse(metadata=self.metadata,
                                      analysistype=self.analysistype)
        self.metadata = ipcress_mismatches(metadata=self.metadata,
                                           analysistype=self.analysistype,
                                           iupac=self.iupac)
        if self.contigbreaks:
            empty = empty_results(metadata=self.metadata,
                                  analysistype=self.analysistype)
            if empty and self.primer_format == 'epcr':
                create_fasta_format(forward_dict=self.forward_dict,
                                    reverse_dict=self.reverse_dict,
                                    fastaprimerfile=self.fastaprimers)
            # Create a BLAST database from the primer file
            make_blastdb(formattedprimers=self.fastaprimers)
            make_blastdb(formattedprimers=self.formattedprimers)
            self.metadata = run_blast(metadata=self.metadata,
                                      analysistype=self.analysistype,
                                      formattedprimers=self.fastaprimers,
                                      blastheader=self.fieldnames,
                                      threads=self.threads)
            self.metadata = parse_blast(metadata=self.metadata,
                                        analysistype=self.analysistype,
                                        fieldnames=self.fieldnames,
                                        forward_dict=self.forward_dict,
                                        reverse_dict=self.reverse_dict,
                                        mismatches=self.mismatches + 1)
            self.metadata = blast_primer_mismatch_details(metadata=self.metadata,
                                                          analysistype=self.analysistype,
                                                          iupac=self.iupac)

        self.ipcress_report()
        if self.export_amplicons:
            amplicon_write(metadata=self.metadata,
                           analysistype=self.analysistype,
                           reportpath=self.reportpath)

    def ipcress_report(self):
        logging.info('Creating summary report')
        with open(self.report, 'w') as report:
            data = 'Sample,Gene,Contig,GenomeLocation,AmpliconSize,Orientation,ForwardMismatches,' \
                   'ForwardMismatchDetails,ForwardLength,ReverseMismatches,ReverseMismatchDetails,ReverseLength,' \
                   'ForwardPrimer,ForwardQuery,ReversePrimer,ReverseQuery\n'
            for sample in self.metadata:
                results = False
                for experiment in sample[self.analysistype].results.datastore:
                    for contig in sample[self.analysistype].results[experiment].datastore:
                        results = True
                        location = str()
                        # Ensure that the .blastresults attribute has been created
                        sample[self.analysistype].blastresults = '{of}_blast_results.tsv' \
                            .format(of=os.path.join(sample.general.outputdirectory, sample.name))
                        if not os.path.isfile(sample[self.analysistype].blastresults):
                            if sample[self.analysistype].results[experiment][contig].direction == 'forward':
                                location = '{forward}-{reverse}'.format(
                                    forward=sample[self.analysistype].results[experiment][contig].forward_pos,
                                    reverse=str(int(sample[self.analysistype].results[experiment][contig].reverse_pos) +
                                                len(sample[self.analysistype].results[experiment][contig]
                                                    .reverse_query))
                                )
                            else:
                                location = '{reverse}-{forward}'.format(
                                    reverse=sample[self.analysistype].results[experiment][contig].reverse_pos,
                                    forward=str(int(sample[self.analysistype].results[experiment][contig].forward_pos) +
                                                len(sample[self.analysistype].results[experiment][contig]
                                                    .forward_query))
                                )
                        else:
                            if sample[self.analysistype].results[experiment][contig].direction == 'forward':
                                if location:
                                    location += ';'
                                location += 'forward_{forward}-{reverse}'.format(
                                    forward=sample[self.analysistype].results[experiment][contig].forward_range[0],
                                    reverse=sample[self.analysistype].results[experiment][contig].forward_range[-1]
                                )
                            else:
                                if location:
                                    location += ';'
                                location += 'reverse_{forward}-{reverse}'.format(
                                    forward=sample[self.analysistype].results[experiment][contig].reverse_range[0],
                                    reverse=sample[self.analysistype].results[experiment][contig].reverse_range[-1]
                                )
                        data += \
                            '{sn},{gene},{contig},{loc},{amplicon_size},{orientation},{f_mm},{f_md},{f_l},{r_mm},' \
                            '{r_md},{r_l},{fp},{fq},{rp},{rq}\n' \
                            .format(sn=sample.name,
                                    gene=sample[self.analysistype].results[experiment][contig].primer_set,
                                    contig=sample[self.analysistype].results[experiment][contig].contig.replace(
                                        '(unmasked)', ''),
                                    loc=location,
                                    amplicon_size=sample[self.analysistype].results[experiment][contig].amplicon_length,
                                    orientation=sample[self.analysistype].results[experiment][contig].direction,
                                    f_mm=sample[self.analysistype].results[experiment][contig].forward_mismatch,
                                    f_md=sample[self.analysistype].results[experiment][contig].forward_mismatch_details,
                                    f_l=len(sample[self.analysistype].results[experiment][contig].forward_ref),
                                    r_mm=sample[self.analysistype].results[experiment][contig].reverse_mismatch,
                                    r_md=sample[self.analysistype].results[experiment][contig].reverse_mismatch_details,
                                    r_l=len(sample[self.analysistype].results[experiment][contig].reverse_ref),
                                    fp=sample[self.analysistype].results[experiment][contig].forward_ref,
                                    fq=sample[self.analysistype].results[experiment][contig].forward_query,
                                    rp=sample[self.analysistype].results[experiment][contig].reverse_ref,
                                    rq=sample[self.analysistype].results[experiment][contig].reverse_query)
                    if not results:
                        # If there were no amplicons, add the sample name and nothing else
                        data += '{sn}\n'.format(sn=sample.name)
            report.write(data)

    def __init__(self, metadataobject, sequencepath, reportpath, primerfile, min_amplicon_size, max_amplicon_size,
                 primer_format, mismatches=2, export_amplicons=False, contigbreaks=False, range_buffer=0):
        self.metadata = metadataobject
        self.analysistype = 'custom_epcr'
        if sequencepath.startswith('~'):
            self.sequencepath = os.path.abspath(os.path.expanduser(os.path.join(sequencepath)))
        else:
            self.sequencepath = os.path.abspath(os.path.join(sequencepath))
        assert os.path.isdir(self.sequencepath), 'Cannot locate the supplied sequence path: {sp}. Please ensure that ' \
                                                 'the folder exists, and you typed in the path correctly'
        if reportpath.startswith('~'):
            self.reportpath = os.path.abspath(os.path.expanduser(os.path.join(reportpath)))
        else:
            self.reportpath = os.path.abspath(os.path.join(reportpath))
        self.report = os.path.join(self.reportpath, '{at}_report.csv'.format(at=self.analysistype))
        if primerfile.startswith('~'):
            self.primerfile = os.path.abspath(os.path.expanduser(os.path.join(primerfile)))
        else:
            self.primerfile = os.path.abspath(os.path.join(primerfile))
        if not self.primerfile:
            assert False, 'Please include the absolute path to the primer file'
        assert os.path.isfile(self.primerfile), 'Cannot locate the primer file: {pf}' \
            .format(pf=self.primerfile)
        self.min_amplicon_size = min_amplicon_size
        self.max_amplicon_size = max_amplicon_size
        self.primer_format = primer_format
        self.mismatches = mismatches
        self.formattedprimers = os.path.join(os.path.dirname(self.primerfile), 'epcr_formatted_primers',
                                             'formatted_primers.txt')
        self.fastaprimers = os.path.join(os.path.dirname(self.primerfile), 'formattedprimers.fa') if \
            self.primer_format != 'fasta' else self.primerfile
        make_path(os.path.dirname(self.formattedprimers))
        try:
            shutil.rmtree(self.reportpath)
        except FileNotFoundError:
            pass
        make_path(self.reportpath)
        self.export_amplicons = export_amplicons
        self.contigbreaks = contigbreaks
        self.range_buffer = range_buffer
        self.forward_dict = dict()
        self.reverse_dict = dict()
        self.threads = multiprocessing.cpu_count() - 1
        self.ipcress_queue = Queue()
        # Dictionary of IUPAC codes
        self.iupac = {
            'R': ['A', 'G'],
            'Y': ['C', 'T'],
            'S': ['C', 'G'],
            'W': ['A', 'T'],
            'K': ['G', 'T'],
            'M': ['A', 'C'],
            'B': ['C', 'G', 'T'],
            'D': ['A', 'G', 'T'],
            'H': ['A', 'C', 'T'],
            'V': ['A', 'C', 'G'],
            'N': ['A', 'C', 'G', 'T'],
            'A': ['A'],
            'C': ['C'],
            'G': ['G'],
            'T': ['T']
        }
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'subject_length', 'alignment_length',
                           'query_start', 'query_end', 'query_sequence',
                           'subject_start', 'subject_end', 'subject_sequence']


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Use ipcress in the exonerate package to perform '
                                        'in silico PCR on FASTA files')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path to folder containing FASTA files')
    parser.add_argument('-m', '--mismatches',
                        default=2,
                        help='Number of mismatches to allow for ipcress searches. Default is 2')
    parser.add_argument('-a', '--analysistype',
                        default='custom',
                        choices=['vtyper', 'custom'],
                        help='Either perform the standard vtyper analysis using the included primer files, or '
                             'supply your own FASTA-formatted (multi-)primer file with the following format: '
                             '>primer1-F\n'
                             'seq\n'
                             '>primer1-R\n'
                             'seq\n'
                             '>primer2-F\n'
                             'etc.')
    parser.add_argument('-pf', '--primerfile',
                        help='Absolute path and name of the primer file to test')
    parser.add_argument('-f', '--primer_format',
                        choices=['epcr', 'fasta'],
                        default='fasta',
                        type=str,
                        help='Format of the supplied primer file. Choices are "epcr" and "fasta". Default is '
                             'fasta. epcr format describes one white-space delimited PCR reaction per line. '
                             'The following fields be included: name of primer pair, sequence of forward '
                             'primer, sequence of reverse primer, min amplicon size, max amplicon size e.g.\n'
                             'vtx1a CCTTTCCAGGTACAACAGCGGTT GGAAACTCATCAGATGCCATTCTGG 0 1500\n'
                             'vtx1c CCTTTCCTGGTACAACTGCGGTT CAAGTGTTGTACGAAATCCCCTCTGA 0 1500\n'
                             '.......\n'
                             'fasta format must have every primer on a separate line AND -F/-R following the '
                             'name e.g.\n'
                             '>vtx1a-F\n'
                             'CCTTTCCAGGTACAACAGCGGTT\n'
                             '>vtx1a-R\n'
                             'GGAAACTCATCAGATGCCATTCTGG\n'
                             '>vtx1c-F\n'
                             '.......\n')
    parser.add_argument('-mas', '--max_amplicon_size',
                        default=1500,
                        type=int,
                        help='Maximum size of amplicons. Default is 1500')
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        help='Enable debug-level messages')
    parser.add_argument('-e', '--export_amplicons',
                        action='store_true',
                        help='Export the sequence of the calculated amplicons. Default is False')
    parser.add_argument('-rb', '--range_buffer',
                        type=int,
                        default=0,
                        help='Increase the buffer size between amplicons in a sequence. Useful if you have overlapping '
                             'primer sets, and you are looking for the best one. Default is 0 "custom" analyses, and '
                             '150 for "vtyper" (vtx1c and vtx1d sets don\'t quite overlap, but the script requires '
                             'an overlap for finding the best primer set)')
    parser.add_argument('-cb', '--contigbreaks',
                        action='store_true',
                        help='Use BLAST to find hits that may be on separate contigs')
    # Get the arguments into an object
    arguments = parser.parse_args()
    SetupLogging(debug=arguments.debug)
    arguments.reportpath = os.path.join(arguments.sequencepath, 'reports')
    arguments.runmetadata = MetadataObject()
    # Create metadata objects for the samples
    arguments.runmetadata.samples = Filer.filer(arguments)
    if arguments.analysistype == 'vtyper':
        epcr = VtyperIP(metadataobject=arguments.runmetadata.samples,
                        analysistype=arguments.analysistype,
                        reportpath=os.path.join(arguments.sequencepath, 'reports'))
        epcr.vtyper()
    else:
        arguments.runmetadata = MetadataObject()
        epcr = CustomIP(metadataobject=arguments.runmetadata,
                        sequencepath=arguments.sequencepath,
                        reportpath=os.path.join(arguments.sequencepath, 'reports'),
                        primerfile=arguments.primerfile,
                        min_amplicon_size=arguments.min_amplicon_size,
                        max_amplicon_size=arguments.max_amplicon_size,
                        primer_format=arguments.primer_format,
                        mismatches=arguments.mismatches,
                        export_amplicons=arguments.export_amplicons,
                        range_buffer=arguments.range_buffer,
                        contigbreaks=arguments.contigbreaks)
        epcr.main()


if __name__ == '__main__':
    # Run the script
    cli()
