#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject, run_subprocess, \
    SetupLogging
from Bio import SeqIO, Seq
from argparse import ArgumentParser
from itertools import product
from threading import Thread
from queue import Queue
import multiprocessing
from glob import glob
import logging
import shutil
import os
import re

__author__ = 'adamkoziol'


def fasta_primers(primerfile, forward_dict, reverse_dict):
    """
    Read in the primer file, and create a properly formatted output file that takes any degenerate bases
    into account
    :param primerfile: Name and path of FASTA-formatted primer file to parse
    :param forward_dict: Dictionary of name: sequence of all forward primers
    :param reverse_dict: Dictionary of name: sequence of all reverse primers
    :return: Populated forward_dict and reverse_dict
    """
    logging.info('Populating primer dictionaries')
    # Dictionary of degenerate IUPAC codes
    iupac = {
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
    primerlist = list()
    for record in SeqIO.parse(primerfile, 'fasta'):
        # from https://stackoverflow.com/a/27552377 - find any degenerate bases in the primer sequence, and
        # create all possibilities as a list
        try:
            primerlist = [''.join(base) for base in product(*[iupac[seq] for seq in str(record.seq).upper()])]
        except TypeError:
            logging.error("Invalid Primer Sequence: {seq}".format(seq=str(record.seq)))
            quit()
        # As the record.id is being updated in the loop below, set the name of the primer here so that will
        # be able to be recalled when setting the new record.ids
        primername = record.id
        # Iterate through all the possible primers created from any degenerate bases
        for primer in primerlist:
            # Split the base name of the target from the direction
            # e.g. vtx1a-F1 is split in vtx1a and F1
            basename, direction = primername.split('-')
            # Populate the dictionaries of forward and reverse primers based on the direction determined above
            if direction.startswith('F'):
                # Attempt to add the current primer sequence to the dictionary
                try:
                    forward_dict[basename].append(primer)
                # On a key error, initialise the list of primers
                except KeyError:
                    forward_dict[basename] = list()
                    forward_dict[basename].append(primer)
            else:
                try:
                    reverse_dict[basename].append(primer)
                except KeyError:
                    reverse_dict[basename] = list()
                    reverse_dict[basename].append(primer)
    return forward_dict, reverse_dict


def epcr_primer_file(formattedprimers, forward_dict, reverse_dict, minampliconsize=0, maxampliconsize=0):
    """
    Create the ePCR-compatible primer file from the dictionaries of primer combinations
    :param formattedprimers: Name and path of ePCR-formatted primer file to be created
    :param forward_dict: Dictionary of name: sequence of all forward primers
    :param reverse_dict: Dictionary of name: sequence of all reverse primers
    :param minampliconsize: Integer of the minimum amplicon size to search
    :param maxampliconsize: Integer of the maximum amplicon size to search
    """
    logging.info('Creating re-PCR-compatible primer file')
    with open(formattedprimers, 'w') as formatted:
        # Iterate through all the targets
        for basename in sorted(forward_dict):
            # Use enumerate to number the iterations for each forward and reverse primer in the lists
            for forward_index, forward_primer in enumerate(forward_dict[basename]):
                for reverse_index, reverse_primer in enumerate(reverse_dict[basename]):
                    # Set the name of the primer using the target name, and the indices of the primers
                    # e.g. vtx1a_0_0
                    primer_name = '{bn}_{fi}_{ri}'.format(bn=basename,
                                                          fi=forward_index,
                                                          ri=reverse_index)
                    # Create the string to write to the ePCR-compatible primer file
                    # e.g. vtx1a_0_0	CCTTTCCAGGTACAACAGCGGTT	GGAAACTCATCAGATGCCATTCTGG
                    output_string = '{pn}\t{fp}\t{rp}'.format(pn=primer_name,
                                                              fp=forward_primer,
                                                              rp=reverse_primer)
                    # If desired, add the minimum and maximum amplicon size information
                    if maxampliconsize:
                        output_string += '\t{min}\t{max}\n'.format(min=minampliconsize,
                                                                   max=maxampliconsize)
                    else:
                        output_string += '\n'
                    # Write the string to file
                    formatted.write(output_string)


def epcr_primers(primerfile, forward_dict, reverse_dict, formattedprimers, copyfile=True):
    """
    Read in the ePCR-formatted primer file, and populate forward and reverse primer dictionaries
    :param primerfile: Name and absolute path to the primer file to parse
    :param forward_dict: Dictionary of all forward primer sequences
    :param reverse_dict: Dictionary of all reverse primer sequences
    :param formattedprimers: Name and absolute path of ePCR-formatted primer file
    :param copyfile: Boolean of whether the primer file should be copied to the formatted primer file path
    :return: Updated forward_dict and reverse_dict
    """
    logging.info('Populating primer dictionaries')
    with open(primerfile, 'r') as primer_file:
        for line in primer_file:
            if len(line.split()) == 5:
                basename, forward_seq, reverse_seq, min_amplicon, max_amplicon = line.rstrip().split()
            elif len(line.split()) == 4:
                basename, forward_seq, reverse_seq, max_amplicon = line.rstrip().split()
            else:
                basename, forward_seq, reverse_seq = line.rstrip().split()
            # Add the current primer sequence to the dictionary
            forward_dict[basename] = forward_seq
            reverse_dict[basename] = reverse_seq
    # Copy the primer file to the location store in self.formattedprimers - makes code reuse easier
    if copyfile:
        try:
            shutil.copyfile(src=primerfile,
                            dst=formattedprimers)
        except shutil.SameFileError:
            pass
    # Return the updated dictionaries
    return forward_dict, reverse_dict


class Vtyper(object):

    def vtyper(self):
        self.primer_prep()
        self.epcr_threads(formattedprimers=self.formattedprimers)
        self.epcr_parse()
        self.epcr_report()

    def primer_prep(self):
        self.forward_dict, self.reverse_dict = fasta_primers(primerfile=self.primerfile,
                                                             forward_dict=self.forward_dict,
                                                             reverse_dict=self.reverse_dict)

        epcr_primer_file(formattedprimers=self.formattedprimers,
                         forward_dict=self.forward_dict,
                         reverse_dict=self.reverse_dict)

    def epcr_threads(self, formattedprimers, ampliconsize=1500):
        """
        Run ePCR in a multi-threaded fashion
        """
        # Create the threads for the ePCR analysis
        for _ in range(self.cpus):
            threads = Thread(target=self.epcr, args=())
            threads.setDaemon(True)
            threads.start()
        logging.info('Running ePCR analyses')
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                setattr(sample, self.analysistype, GenObject())
                # Get the primers ready
                sample[self.analysistype].primers = formattedprimers
                # Make the output path
                sample[self.analysistype].reportdir = os.path.join(sample.general.outputdirectory,
                                                                   self.analysistype)
                make_path(sample[self.analysistype].reportdir)
                outfile = os.path.join(sample[self.analysistype].reportdir, sample.name)
                # Set the hashing and mapping commands
                sample.commands.famap = '{famap} -b {outfile}.famap {fasta}'\
                    .format(famap=os.path.join(self.homepath, 'ePCR', 'famap'),
                            outfile=outfile,
                            fasta=sample.general.bestassemblyfile)
                sample.commands.fahash = '{fahash} -b {outfile}.hash {outfile}.famap'\
                    .format(fahash=os.path.join(self.homepath, 'ePCR', 'fahash'),
                            outfile=outfile)
                # re-PCR uses the subtyping primers list to search the contigs file using the following parameters
                # -S {hash file} (Perform STS lookup using hash-file), -r + (Enable/disable reverse STS lookup)
                # -m 1500 (Set variability for STS size for lookup), this very large, as I don't necessarily know
                # the size of the amplicon
                # -n 1 (Set max allowed mismatches per primer pair for lookup)
                # -g 0 (Set max allowed indels per primer pair for lookup),
                # -G (Print alignments in comments)
                # -o {output file}
                sample.commands.epcr = \
                    '{rePCR} -S {outfile}.hash -r + -d 1-{ampsize} -n {mismatches} -g 0 -G -q ' \
                    '-o {outfile}.txt {primers}'\
                    .format(rePCR=os.path.join(self.homepath, 'ePCR', 're-PCR'),
                            outfile=outfile,
                            ampsize=ampliconsize,
                            mismatches=self.mismatches,
                            primers=sample[self.analysistype].primers)
                sample[self.analysistype].resultsfile = '{of}.txt'.format(of=outfile)
                # Add the sample object and the output file to the queue
                self.epcrqueue.put((sample, outfile))
        # Join the threads
        self.epcrqueue.join()

    def epcr(self):
        while True:
            sample, linkfile = self.epcrqueue.get()
            # Run the commands if the ePCR output file doesn't exist
            if not os.path.isfile('{lf}.txt'.format(lf=linkfile)):
                run_subprocess(sample.commands.famap)
                run_subprocess(sample.commands.fahash)
                run_subprocess(sample.commands.epcr)
            # Clean up the temporary files
            try:
                os.remove('{lf}.famap'.format(lf=linkfile))
            except FileNotFoundError:
                pass
            try:
                os.remove('{lf}.hash'.format(lf=linkfile))
            except FileNotFoundError:
                pass
            # Signal that the thread is complete
            self.epcrqueue.task_done()

    def epcr_parse(self):
        """
        Parse the ePCR outputs
        """
        logging.info('Parsing ePCR outputs')
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Create a set to store all the unique results
                toxin_set = set()
                if os.path.isfile(sample[self.analysistype].resultsfile):
                    with open(sample[self.analysistype].resultsfile) as epcrresults:
                        for result in epcrresults:
                            # Only the lines without a # contain results
                            if "#" not in result:
                                # Split on \t
                                data = result.split('\t')
                                # The subtyping primer pair is the first entry on lines with results
                                vttype = data[0].split('_')[0]
                                # Add the verotoxin subtype to the set of detected subtypes
                                toxin_set.add(vttype)
                # Create a string of the entries in the sorted list of toxins joined with ";"
                sample[self.analysistype].toxinprofile = ";".join(sorted(list(toxin_set))) if toxin_set else 'ND'
            else:
                setattr(sample, self.analysistype, GenObject())
                sample[self.analysistype].toxinprofile = 'ND'

    def epcr_report(self):
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

    def __init__(self, inputobject, analysistype, mismatches=2):
        self.metadata = inputobject.runmetadata.samples
        self.analysistype = analysistype
        self.reportpath = inputobject.reportpath
        self.primer_format = inputobject.primer_format
        self.mismatches = mismatches
        make_path(self.reportpath)
        self.devnull = open(os.devnull, 'wb')
        self.cpus = multiprocessing.cpu_count() - 1
        self.epcrqueue = Queue()
        # Extract the path of the current script from the full path + file name
        self.homepath = os.path.split(os.path.abspath(__file__))[0]
        self.formattedprimers = os.path.join(self.homepath, 'ssi_subtyping_primers.txt')
        self.primerfile = os.path.join(self.homepath, 'primers.txt')
        self.forward_dict = dict()
        self.reverse_dict = dict()


class Custom(object):
    def main(self):
        if self.primer_format == 'fasta':
            # if not os.path.isfile(self.formattedprimers):
            logging.info('Extracting primer sequences')
            fasta_primers(primerfile=self.primerfile,
                          forward_dict=self.forward_dict,
                          reverse_dict=self.reverse_dict)
            logging.info('Creating ePCR-compatible primer file')
            epcr_primer_file(formattedprimers=self.formattedprimers,
                             forward_dict=self.forward_dict,
                             reverse_dict=self.reverse_dict)
        else:
            self.forward_dict, self.reverse_dict = epcr_primers(primerfile=self.primerfile,
                                                                forward_dict=self.forward_dict,
                                                                reverse_dict=self.reverse_dict,
                                                                formattedprimers=self.formattedprimers)
        self.vtyper_object.epcr_threads(formattedprimers=self.formattedprimers,
                                        ampliconsize=self.ampliconsize)
        logging.info('Parsing ePCR outputs'.format())
        self.parse_epcr()
        logging.info('Creating {at} report'.format(at=self.analysistype))
        self.create_epr_report()

    def parse_epcr(self):
        """
        Parse the ePCR output file. Populate dictionary of results. For alleles, find the best result based on the
        number of mismatches before populating dictionary
        """
        # Use the metadata object from the vtyper_object
        for sample in self.vtyper_object.metadata:
            # Initialise the dictionary
            sample[self.analysistype].result_dict = dict()
            # Read in the output file
            with open(sample[self.analysistype].resultsfile) as epcrresults:
                for result in epcrresults:
                    # Only the lines without a # contain results
                    if "#" not in result:
                        # Split on \t
                        # vtx2a_0_0 2014-SEQ-0121_127_length_1407_cov_50.7797_ID_10924  -   228   576 2   0 349/100-350
                        # primer_set: vtx2a_0_0, contig: 2014-SEQ-0121_127_length_1407_cov_50.7797_ID_10924, strand: -,
                        # start: 228, stop: 576, number of forward mismatches: 2, number of reverse mismatches: 2
                        # amplicon_combo: 349/100-350
                        primer_set, contig, strand, start, stop, total_mismatches, indels, amplicon_combo = \
                            result.rstrip().split('\t')
                        # Set the mismatches to be an int
                        total_mismatches = int(total_mismatches)
                        # Set the position of the amplicon on the contig. Ensure that the lower value is first
                        genome_pos = '{min}-{max}'.format(min=min([int(start), int(stop)]),
                                                          max=max([int(start), int(stop)]))
                        # Extract the gene name from the modified name used when creating the primer file: LMhlyA_0_0
                        # becomes LMhlyA
                        gene_re = re.search(r'([\w-]+)_(\d{1,3})_(\d{1,3})', primer_set)
                        gene = gene_re.groups()[0]
                        # Split the amplicon length from amplicon_combo: 349/100-350 -> 349
                        amplicon_length = amplicon_combo.split('/')[0]
                        # Populate the dictionary if the 'total_mismatches' key doesn't exist, or if the current number
                        # of mismatches is better than the previous 'best' number of mismatches
                        try:
                            if total_mismatches < sample[self.analysistype].result_dict[gene]['total_mismatches']:
                                self.populate_results_dict(sample=sample,
                                                           gene=gene,
                                                           total_mismatches=total_mismatches,
                                                           genome_pos=genome_pos,
                                                           amplicon_length=amplicon_length,
                                                           contig=contig,
                                                           primer_set=primer_set)
                        except KeyError:
                            self.populate_results_dict(sample=sample,
                                                       gene=gene,
                                                       total_mismatches=total_mismatches,
                                                       genome_pos=genome_pos,
                                                       amplicon_length=amplicon_length,
                                                       contig=contig,
                                                       primer_set=primer_set)

    def populate_results_dict(self, sample, gene, total_mismatches, genome_pos, amplicon_length, contig, primer_set):
        """
        Populate the results dictionary with the required key: value pairs
        :param sample: type MetadataObject: Current metadata sample to process
        :param gene: type STR: Gene of interest
        :param total_mismatches: type INT: Number of mismatches between primer pairs and subject sequence
        :param genome_pos: type STR: Positions of 5' and 3' ends of the amplicon
        :param amplicon_length: type INT: Total length of the amplicon
        :param contig: type STR: Contig name
        :param primer_set: type STR: Name of primer set from the ePCR-formatted file used in the analyses
        """
        sample[self.analysistype].result_dict[gene] = {
            'total_mismatches': total_mismatches,
            'genome_pos': genome_pos,
            'amplicon_length': amplicon_length,
            'contig': contig,
            'primer_set': primer_set
        }

    def create_epr_report(self):
        """
        Parse the results dictionaries, and create a final report
        """
        # Open the report as a .csv file
        with open(os.path.join(self.reportpath, 'ePCR_report.csv'), 'w') as report:
            # Initialise a string to store the header
            results = 'Sample,Gene,GenomeLocation,AmpliconSize,Contig,TotalMismatches,PrimerSet\n'
            for sample in self.vtyper_object.metadata:
                # Check to see if there are strain-specific results
                if sample[self.analysistype].result_dict:
                    for gene, result_dict in sample[self.analysistype].result_dict.items():
                        # Populate the string with the appropriate values extracted from the dictionary
                        results += '{sn},{gene},{genomelocation},{ampliconsize},{contig},{nm},{ps}\n'\
                            .format(sn=sample.name,
                                    gene=gene,
                                    genomelocation=result_dict['genome_pos'],
                                    ampliconsize=result_dict['amplicon_length'],
                                    contig=result_dict['contig'],
                                    nm=result_dict['total_mismatches'],
                                    ps=result_dict['primer_set'])
                        if self.export_amplicons:
                            self.ampliconfile(sample=sample,
                                              contig=result_dict['contig'],
                                              amplicon_range=result_dict['genome_pos'].split('-'),
                                              primer_set=result_dict['primer_set'])
                else:
                    results += '{sn}\n'.format(sn=sample.name)
            # Write the complete string to the report
            report.write(results)

    def ampliconfile(self, sample, contig, amplicon_range, primer_set):
        """
        Extracts amplicon sequence from contig file
        :param sample: sample metadata object
        :param contig: name of the contig hit by primers
        :param amplicon_range: type LIST: range of the amplicon within the contig
        :param primer_set: type STR: Name of primer set from the ePCR-formatted primer file used in the analyses
        """
        # Open the file
        sample[self.analysistype].ampliconfile = os.path.join(self.reportpath, '{sn}_amplicons.fa'
                                                              .format(sn=sample.name))
        with open(sample[self.analysistype].ampliconfile, 'a') as ampliconfile:
            # Load the records from the assembly into the dictionary
            for record in SeqIO.parse(sample.general.bestassemblyfile, 'fasta'):
                if record.id == contig:
                    try:
                        # Extract the start and end positions of the supplied range
                        start, end = amplicon_range
                        # Slice the gene sequence from the sequence record - remember to subtract one to
                        # allow for zero-based indexing
                        genesequence = str(record.seq)[int(start) - 1:int(end)]
                        # Set the record.id to be the sample name, the contig name,
                        # the range, and the primers
                        record.id = '{sn}_{contig}_{genomepos}_{primer_set}' \
                            .format(sn=sample.name,
                                    contig=contig,
                                    genomepos='_'.join(amplicon_range),
                                    primer_set=primer_set)
                        # Clear the record.description
                        record.description = ''
                        # Create a seq record from the sliced genome sequence
                        record.seq = Seq.Seq(genesequence)
                        # Write the amplicon to file
                        SeqIO.write(record, ampliconfile, 'fasta')
                    except AttributeError:
                        pass

    def __init__(self, inputobject, analysistype, primerfile, ampliconsize, primer_format, mismatches=2,
                 export_amplicons=False):
        self.runmetadata = inputobject.runmetadata
        self.analysistype = analysistype
        self.reportpath = inputobject.reportpath
        self.primerfile = primerfile
        if not self.primerfile:
            assert False, 'Please include the absolute path to the FASTA-formatted primer file'
        assert os.path.isfile(self.primerfile), 'Cannot locate the specified FASTA-formatted primer file: {pf}'\
            .format(pf=self.primerfile)
        self.ampliconsize = ampliconsize
        self.primer_format = primer_format
        self.mismatches = mismatches
        # self.primer_format = inputobject.primer_format
        self.formattedprimers = os.path.join(os.path.dirname(self.primerfile), 'epcr_formatted_primers',
                                             'formatted_primers.txt')
        make_path(os.path.dirname(self.formattedprimers))
        try:
            shutil.rmtree(self.reportpath)
        except FileNotFoundError:
            pass
        make_path(self.reportpath)
        self.devnull = open(os.devnull, 'wb')
        self.epcrqueue = Queue(maxsize=multiprocessing.cpu_count())
        self.export_amplicons = export_amplicons
        self.forward_dict = dict()
        self.reverse_dict = dict()
        # Create an object, so that the script can call methods from the Vtyper class
        self.vtyper_object = Vtyper(inputobject=self,
                                    analysistype=self.analysistype,
                                    mismatches=self.mismatches)


class Filer(object):

    @staticmethod
    def filer(args):
        """
        Create metadata objects with necessary attributes for each FASTA file found in the sequence path
        :param args: Argument parser object with necessary variables
        :return: samples: List of metadata objects
        """
        # List to store all the metadata objects
        samples = list()
        if args.sequencepath.startswith('~'):
            args.sequencepath = os.path.abspath(os.path.expanduser(os.path.join(args.sequencepath)))
        else:
            args.sequencepath = os.path.abspath(os.path.join(args.sequencepath))
        # Find all the sequence files in the path
        fastas = sorted(glob(os.path.join(args.sequencepath, '*.fa*')))
        for fasta in fastas:
            # Create a metadata object for each sample
            metadata = MetadataObject()
            # Populate the metadata object with the required attributes
            metadata.name = os.path.splitext(os.path.basename(fasta))[0]
            metadata.general = GenObject()
            metadata.commands = GenObject()
            metadata.general.bestassemblyfile = fasta
            metadata.general.outputdirectory = os.path.join(args.sequencepath, metadata.name)
            samples.append(metadata)
        return samples


if __name__ == '__main__':

    def argument_parser():

        # Parser for arguments
        parser = ArgumentParser(description='Perform verotoxin sub-typing on FASTA files')
        parser.add_argument('-s', '--sequencepath',
                            required=True,
                            help='Path to folder containing sequencing reads')
        parser.add_argument('-m', '--mismatches',
                            default=2,
                            help='Number of mismatches to allow for ePCR searches. Default is 2')
        parser.add_argument('-a', '--analysistype',
                            default='vtyper',
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
        parser.add_argument('-mas', '--maxampliconsize',
                            default=1500,
                            type=int,
                            help='Maximum size of amplicons. Default is 1500')
        parser.add_argument('-d', '--debug',
                            action='store_true',
                            help='Enable debug-level messages')
        parser.add_argument('-e', '--export_amplicons',
                            action='store_true',
                            help='Export the sequence of the calculated amplicons. Default is False')
        # Get the arguments into an object
        arguments = parser.parse_args()
        SetupLogging(debug=arguments.debug)
        arguments.reportpath = os.path.join(arguments.sequencepath, 'reports')
        arguments.runmetadata = MetadataObject()
        # Create metadata objects for the samples
        arguments.runmetadata.samples = Filer.filer(arguments)
        if arguments.analysistype == 'vtyper':
            # Perform vtx typing
            vtyper = Vtyper(inputobject=arguments,
                            analysistype='vtyper_legacy',
                            mismatches=arguments.mismatches)
            vtyper.vtyper()
        else:
            epcr = Custom(inputobject=arguments,
                          analysistype='custom_epcr',
                          primerfile=arguments.primerfile,
                          ampliconsize=arguments.maxampliconsize,
                          primer_format=arguments.primer_format,
                          mismatches=arguments.mismatches,
                          export_amplicons=arguments.export_amplicons)
            epcr.main()

    # Run the script
    argument_parser()
