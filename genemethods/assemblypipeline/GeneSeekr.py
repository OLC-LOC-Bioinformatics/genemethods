#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import combinetargets, GenObject, make_dict, make_path, \
    MetadataObject, printtime, run_subprocess, write_to_logfile
from genemethods.typing_classes.typing_classes import ResistanceNotes
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Application import ApplicationError
from Bio.pairwise2 import format_alignment
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser
from collections import defaultdict
from threading import Thread
from csv import DictReader
from queue import Queue
from glob import glob
import xlsxwriter
import threading
import time
import csv
import sys
import os
import re

__author__ = 'adamkoziol'


class GeneSeekr(object):

    def geneseekr(self):
        # Make blast databases (if necessary)
        printtime('Creating {} blast databases as required'.format(self.analysistype), self.start)
        self.makedbthreads()
        # Run the blast analyses
        printtime('Running {} blast analyses'.format(self.analysistype), self.start)
        self.blastnthreads()
        if self.unique:
            self.filterunique()
        if self.analysistype == 'resfinder':
            self.resfinderreporter()
        elif self.analysistype == 'virulence':
            self.virulencefinderreporter()
        else:
            self.reporter()
        # Remove the attributes from the object; they take up too much room on the .json report
        for sample in self.metadata:
            delattr(sample[self.analysistype], "targetnames")
            delattr(sample[self.analysistype], "targets")
        printtime('{} analyses complete'.format(self.analysistype), self.start)

    def filterunique(self):
        """
        Filters multiple BLAST hits in a common region of the genome. Leaves only the best hit
        """
        for sample in self.metadata:
            # Initialise variables
            sample[self.analysistype].blastresults = list()
            resultdict = dict()
            rowdict = dict()
            try:
                # Iterate through all the contigs, which had BLAST hits
                for contig in sample[self.analysistype].queryranges:
                    # Find all the locations in each contig that correspond to the BLAST hits
                    for location in sample[self.analysistype].queryranges[contig]:
                        # Extract the BLAST result dictionary for the contig
                        for row in sample[self.analysistype].results[contig]:
                            # Initialise variable to reduce the number of times row['value'] needs to be typed
                            contig = row['query_id']
                            high = row['high']
                            low = row['low']
                            percentidentity = row['percentidentity']
                            # Join the two ranges in the location list with a comma
                            locstr = ','.join([str(x) for x in location])
                            # Create a set of the location of all the base pairs between the low and high (-1) e.g.
                            # [6, 10] would give 6, 7, 8, 9, but NOT 10. This turns out to be useful, as there are
                            # genes located back-to-back in the genome e.g. strB and strA, with locations of 2557,3393
                            # and 3393,4196, respectively. By not including 3393 in the strB calculations, I don't
                            # have to worry about this single bp overlap
                            loc = set(range(low, high))
                            # Use a set intersection to determine whether the current result overlaps with location
                            # This will allow all the hits to be grouped together based on their location
                            if loc.intersection(set(range(location[0], location[1]))):
                                # Populate the grouped hits for each location
                                try:
                                    resultdict[contig][locstr].append(percentidentity)
                                    rowdict[contig][locstr].append(row)
                                # Initialise and populate the lists of the nested dictionary
                                except KeyError:
                                    try:
                                        resultdict[contig][locstr] = list()
                                        resultdict[contig][locstr].append(percentidentity)
                                        rowdict[contig][locstr] = list()
                                        rowdict[contig][locstr].append(row)
                                    # As this is a nested dictionary, it needs to be initialised here
                                    except KeyError:
                                        resultdict[contig] = dict()
                                        resultdict[contig][locstr] = list()
                                        resultdict[contig][locstr].append(percentidentity)
                                        rowdict[contig] = dict()
                                        rowdict[contig][locstr] = list()
                                        rowdict[contig][locstr].append(row)
            except KeyError:
                pass
            # Find the best hit for each location based on percent identity
            for contig in resultdict:
                # Do not allow the same gene to be added to the dictionary more than once
                genes = list()
                for location in resultdict[contig]:
                    # Initialise a variable to determine whether there is already a best hit found for the location
                    multiple = False
                    # Iterate through the BLAST results to find the best hit
                    for row in rowdict[contig][location]:
                        # Add the best hit to the .blastresults attribute of the object
                        if row['percentidentity'] == max(resultdict[contig][location]) and not multiple \
                                and row['subject_id'] not in genes:
                            sample[self.analysistype].blastresults.append(row)
                            genes.append(row['subject_id'])
                            multiple = True

    def makedbthreads(self):
        """
        Setup and create threads for class
        """
        # Find all the target folders in the analysis and add them to the targetfolders set
        for sample in self.metadata:
            if sample[self.analysistype].combinedtargets != 'NA':
                self.targetfolders.add(sample[self.analysistype].targetpath)
        # Create and start threads for each fasta file in the list
        for i in range(len(self.targetfolders)):
            # Send the threads to makeblastdb
            threads = Thread(target=self.makeblastdb, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        # Make blast databases for MLST files (if necessary)
        for targetdir in self.targetfolders:
            # List comprehension to remove any previously created database files from list
            self.targetfiles = glob(os.path.join(targetdir, '*.fasta'))
            try:
                _ = self.targetfiles[0]
            except IndexError:
                self.targetfiles = glob(os.path.join(targetdir, '*.fasta'))
            for targetfile in self.targetfiles:
                # Read the sequences from the target file to a dictionary
                self.records[targetfile] = SeqIO.to_dict(SeqIO.parse(targetfile, 'fasta'))
                # Add the fasta file to the queue
                self.dqueue.put(targetfile)
        self.dqueue.join()  # wait on the dqueue until everything has been processed

    def makeblastdb(self):
        """Makes blast database files from targets as necessary"""
        while True:  # while daemon
            fastapath = self.dqueue.get()  # grabs fastapath from dqueue
            # remove the path and the file extension for easier future globbing
            db = os.path.splitext(fastapath)[0]
            nhr = '{}.nhr'.format(db)  # add nhr for searching
            # fnull = open(os.devnull, 'w')  # define /dev/null
            if not os.path.isfile(str(nhr)):  # if check for already existing dbs
                # Create the databases
                # TODO use MakeBLASTdb class
                threadlock = threading.Lock()
                command = 'makeblastdb -in {} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {}'.format(fastapath, db)
                # subprocess.call(shlex.split('makeblastdb -in {} -parse_seqids -max_file_sz 2GB -dbtype nucl -out {}'
                #                            .format(fastapath, db)), stdout=fnull, stderr=fnull)
                out, err = run_subprocess(command)
                threadlock.acquire()
                write_to_logfile(command, command, self.logfile, None, None, None, None)
                write_to_logfile(out, err, self.logfile, None, None, None, None)
                threadlock.release()
            self.dqueue.task_done()  # signals to dqueue job is done

    def blastnthreads(self):
        """Setup and create  threads for blastn and xml path"""
        # Create the threads for the BLAST analysis
        for i in range(self.cpus):
            threads = Thread(target=self.runblast, args=())
            threads.setDaemon(True)
            threads.start()
        # Populate threads for each gene, genome combination
        for sample in self.metadata:
            make_path(sample[self.analysistype].reportdir)
            sample[self.analysistype].report = os.path.join(
                sample[self.analysistype].reportdir, '{}.csv'.format(sample.name))
            if sample[self.analysistype].combinedtargets != 'NA':
                # Add each fasta file combination to the threads
                self.blastqueue.put((sample.general.bestassemblyfile, sample[self.analysistype].combinedtargets,
                                     sample))
        # Join the threads
        self.blastqueue.join()

    def runblast(self):
        while True:  # while daemon
            (assembly, target, sample) = self.blastqueue.get()  # grabs fastapath from dqueue
            # Run the BioPython BLASTn module with the genome as query, fasta(target gene) as db.
            # Do not re-perform the BLAST search each time
            try:
                size = os.path.getsize(sample[self.analysistype].report)
                # If a report was created, but no results entered - program crashed, or no sequences passed thresholds,
                # remove the report, and run the blast analyses again
                if size == 0:
                    os.remove(sample[self.analysistype].report)
            except FileNotFoundError:
                pass
            # Split the extension from the file path
            db = os.path.splitext(target)[0]
            # BLAST command line call. Note the mildly restrictive evalue, and the high number of alignments.
            # Due to the fact that all the targets are combined into one database, this is to ensure that all potential
            # alignments are reported. Also note the custom outfmt: the doubled quotes are necessary to get it work
            blastn = NcbiblastnCommandline(query=assembly, db=db, evalue='1E-5', num_alignments=1000000,
                                           num_threads=12,
                                           outfmt="'6 qseqid sseqid positive mismatch gaps "
                                                  "evalue bitscore slen length qstart qend qseq sstart send sseq'",
                                           out=sample[self.analysistype].report)
            # Save the blast command in the metadata
            sample[self.analysistype].blastcommand = str(blastn)
            # Only run blast if the report doesn't exist
            if not os.path.isfile(sample[self.analysistype].report):
                try:
                    blastn()
                except ApplicationError:
                    self.blastqueue.task_done()
                    self.blastqueue.join()
                    try:
                        os.remove(sample[self.analysistype].report)
                    except (IOError, ApplicationError):
                        pass
                    raise
            # Parse the output depending on whether unique results are desired
            if self.unique:
                self.uniqueblastparser(sample[self.analysistype].report, sample)
            else:
                # Run the blast parsing module
                self.blastparser(sample[self.analysistype].report, sample)
            self.blastqueue.task_done()  # signals to dqueue job is done

    def blastparser(self, report, sample):
        """
        Parse the blast results, and store necessary data in dictionaries in sample object
        :param report: Name of the blast output report being parsed
        :param sample: sample object
        """
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(report), fieldnames=self.fieldnames, dialect='excel-tab')
        resultdict = dict()
        # Initialise a dictionary to store all the target sequences
        sample[self.analysistype].targetsequence = dict()
        # Go through each BLAST result
        for row in blastdict:
            # Calculate the percent identity and extract the bitscore from the row
            # Percent identity is the (length of the alignment - number of mismatches) / total subject length
            percentidentity = float('{:0.2f}'.format((float(row['positives']) - float(row['gaps'])) /
                                                     float(row['subject_length']) * 100))
            target = row['subject_id']
            # If the percent identity is greater than the cutoff
            if percentidentity >= self.cutoff:
                # Update the dictionary with the target and percent identity
                resultdict.update({target: percentidentity})
                # Determine if the orientation of the sequence is reversed compared to the reference
                if int(row['subject_end']) < int(row['subject_start']):
                    # Create a sequence object using Biopython
                    seq = Seq(row['query_sequence'])
                    # Calculate the reverse complement of the sequence
                    querysequence = str(seq.reverse_complement())
                # If the sequence is not reversed, use the sequence as it is in the output
                else:
                    querysequence = row['query_sequence']
                # Add the sequence in the correct orientation to the sample
                sample[self.analysistype].targetsequence[target] = querysequence
            # Add the percent identity to the object
            sample[self.analysistype].blastresults = resultdict
        # Populate missing results with 'NA' values
        if len(resultdict) == 0:
            sample[self.analysistype].blastresults = 'NA'

    def uniqueblastparser(self, report, sample):
        """
        Find the best hit at a location, and discard any other matches
        :param report: Name of the blast output report being parsed
        :param sample: sample object
        """
        # Encountering the following error: # _csv.Error: field larger than field limit (131072)
        # According to https://stackoverflow.com/a/15063941, increasing the field limit should fix the issue
        csv.field_size_limit(sys.maxsize)
        # Open the sequence profile file as a dictionary
        blastdict = DictReader(open(report), fieldnames=self.fieldnames, dialect='excel-tab')
        # Initialise a dictionary to store all the target sequences
        sample[self.analysistype].targetsequence = dict()
        sample[self.analysistype].queryranges = dict()
        sample[self.analysistype].querypercent = dict()
        sample[self.analysistype].queryscore = dict()
        sample[self.analysistype].results = dict()
        # Go through each BLAST result
        for row in blastdict:
            # Calculate the percent identity and extract the bitscore from the row
            # Percent identity is the (length of the alignment - number of mismatches) / total subject length
            percentidentity = float('{:0.2f}'.format((float(row['positives'])) /
                                                     float(row['subject_length']) * 100))
            target = row['subject_id']
            contig = row['query_id']
            high = max([int(row['query_start']), int(row['query_end'])])
            low = min([int(row['query_start']), int(row['query_end'])])
            score = row['bit_score']
            # Create new entries in the blast results dictionaries with the calculated variables
            row['percentidentity'] = percentidentity
            row['low'] = low
            row['high'] = high
            row['alignment_fraction'] = float('{:0.2f}'.format(float(float(row['alignment_length']) /
                                                                     float(row['subject_length']) * 100)))
            # If the percent identity is greater than the cutoff
            if percentidentity >= self.cutoff:
                try:
                    sample[self.analysistype].results[contig].append(row)
                    # Boolean to store whether the list needs to be updated
                    append = True
                    # Iterate through all the ranges in the list - if the new range is different than any of the ranges
                    # seen before, append it. Otherwise, update the previous ranges with the new, longer range as
                    # necessary e.g. [2494, 3296] will be updated to [2493, 3296] with [2493, 3293], and
                    # [2494, 3296] will become [[2493, 3296], [3296, 4132]] with [3296, 4132]
                    for spot in sample[self.analysistype].queryranges[contig]:
                        # Update the low value if the new low value is slightly lower than before
                        if 1 <= (spot[0] - low) <= 100:
                            # Update the low value
                            spot[0] = low
                            # It is not necessary to append
                            append = False
                        # Update the previous high value if the new high value is slightly higher than before
                        elif 1 <= (high - spot[1]) <= 100:
                            # Update the high value in the list
                            spot[1] = high
                            # It is not necessary to append
                            append = False
                        # Do not append if the new low is slightly larger than before
                        elif 1 <= (low - spot[0]) <= 100:
                            append = False
                        # Do not append if the new high is slightly smaller than before
                        elif 1 <= (spot[1] - high) <= 100:
                            append = False
                        # Do not append if the high and low are the same as the previously recorded values
                        elif low == spot[0] and high == spot[1]:
                            append = False
                    # If the result appears to be in a new location, add the data to the object
                    if append:
                        sample[self.analysistype].queryranges[contig].append([low, high])
                        sample[self.analysistype].querypercent[contig] = percentidentity
                        sample[self.analysistype].queryscore[contig] = score
                # Initialise and populate the dictionary for each contig
                except KeyError:
                    sample[self.analysistype].queryranges[contig] = list()
                    sample[self.analysistype].queryranges[contig].append([low, high])
                    sample[self.analysistype].querypercent[contig] = percentidentity
                    sample[self.analysistype].queryscore[contig] = score
                    sample[self.analysistype].results[contig] = list()
                    sample[self.analysistype].results[contig].append(row)
                    sample[self.analysistype].targetsequence[target] = dict()
                # Determine if the query sequence is in a different frame than the subject, and correct
                # by setting the query sequence to be the reverse complement
                if int(row['subject_end']) < int(row['subject_start']):
                    # Create a sequence object using Biopython
                    seq = Seq(row['query_sequence'])
                    # Calculate the reverse complement of the sequence
                    querysequence = str(seq.reverse_complement())
                # If the sequence is not reversed, use the sequence as it is in the output
                else:
                    querysequence = row['query_sequence']
                # Add the sequence in the correct orientation to the sample
                sample[self.analysistype].targetsequence[target] = querysequence

    def reporter(self):
        """
        Creates .xlsx reports using xlsxwriter
        """
        # Create a workbook to store the report. Using xlsxwriter rather than a simple csv format, as I want to be
        # able to have appropriately sized, multi-line cells
        workbook = xlsxwriter.Workbook(os.path.join(self.reportpath, '{}.xlsx'.format(self.analysistype)))
        # New worksheet to store the data
        worksheet = workbook.add_worksheet()
        # Add a bold format for header cells. Using a monotype font size 10
        bold = workbook.add_format({'bold': True, 'font_name': 'Courier New', 'font_size': 10})
        # Format for data cells. Monotype, size 10, top vertically justified
        courier = workbook.add_format({'font_name': 'Courier New', 'font_size': 10})
        courier.set_align('top')
        # Initialise the position within the worksheet to be (0,0)
        row = 0
        # A dictionary to store the column widths for every header
        columnwidth = dict()
        for sample in self.metadata:
            # Reset the column to zero
            col = 0
            # Initialise a list to store all the data for each strain
            data = list()
            # Initialise a list of all the headers with 'Strain'
            headers = ['Strain']
            if sample[self.analysistype].targetnames != 'NA':
                # Append the sample name to the data list only if the script could find targets
                data.append(sample.name)
                if sample[self.analysistype].blastresults != 'NA':
                    for target in sorted(sample[self.analysistype].targetnames):
                        # Add the name of the gene to the header
                        headers.append(target)
                        try:
                            # Append the percent identity to the data list
                            data.append(str(sample[self.analysistype].blastresults[target]))
                            # Only if the alignment option is selected, for inexact results, add alignments
                            if self.align and sample[self.analysistype].blastresults[target] != 100.00:
                                # Align the protein (and nucleotide) sequences to the reference
                                self.alignprotein(sample, target)
                                # Add the appropriate headers
                                headers.extend(['{}_aa_Identity'.format(target),
                                                '{}_aa_Alignment'.format(target),
                                                '{}_aa_SNP_location'.format(target),
                                                '{}_nt_Alignment'.format(target),
                                                '{}_nt_SNP_location'.format(target)
                                                ])
                                # Add the alignment, and the location of mismatches for both nucleotide and amino
                                # acid sequences
                                data.extend([sample[self.analysistype].aaidentity[target],
                                             sample[self.analysistype].aaalign[target],
                                             sample[self.analysistype].aaindex[target],
                                             sample[self.analysistype].ntalign[target],
                                             sample[self.analysistype].ntindex[target],
                                             ])
                        # If there are no blast results for the target, add a '-'
                        except (KeyError, TypeError):
                            data.append('-')
                        # If there are no blast results at all, add a '-'
                        else:
                            data.append('-')
            # Write the header to the spreadsheet
            for header in headers:
                worksheet.write(row, col, header, bold)
                # Set the column width based on the longest header
                try:
                    columnwidth[col] = len(header)if len(header) > columnwidth[col] else columnwidth[col]
                except KeyError:
                    columnwidth[col] = len(header)
                worksheet.set_column(col, col, columnwidth[col])
                col += 1
            # Increment the row and reset the column to zero in preparation of writing results
            row += 1
            col = 0
            # List of the number of lines for each result
            totallines = list()
            # Write out the data to the spreadsheet
            for results in data:
                worksheet.write(row, col, results, courier)
                try:
                    # Counting the length of multi-line strings yields columns that are far too wide, only count
                    # the length of the string up to the first line break
                    alignmentcorrect = len(results.split('\n')[0])
                    # Count the number of lines for the data
                    lines = results.count('\n') if results.count('\n') >= 1 else 1
                    # Add the number of lines to the list
                    totallines.append(lines)
                # If there are no newline characters, set the width to the length of the string
                except AttributeError:
                    alignmentcorrect = len(results)
                    lines = 1
                    # Add the number of lines to the list
                    totallines.append(lines)
                # Increase the width of the current column, if necessary
                try:
                    columnwidth[col] = alignmentcorrect if alignmentcorrect > columnwidth[col] else columnwidth[col]
                except KeyError:
                    columnwidth[col] = alignmentcorrect
                worksheet.set_column(col, col, columnwidth[col])
                col += 1
            # Set the width of the row to be the number of lines (number of newline characters) * 12
            if len(totallines) != 0:
                worksheet.set_row(row, max(totallines) * 12)
            else:
                worksheet.set_row(row, 1)
            # Increase the row counter for the next strain's data
            row += 1
        # Close the workbook
        workbook.close()

    def alignprotein(self, sample, target):
        """
        Create alignments of the sample nucleotide and amino acid sequences to the reference sequences
        """
        # Initialise dictionaries
        sample[self.analysistype].dnaseq = dict()
        sample[self.analysistype].protseq = dict()
        sample[self.analysistype].ntindex = dict()
        sample[self.analysistype].aaindex = dict()
        sample[self.analysistype].ntalign = dict()
        sample[self.analysistype].aaalign = dict()
        sample[self.analysistype].aaidentity = dict()
        # Remove any gaps incorporated into the sequence
        sample[self.analysistype].targetsequence[target] = \
            sample[self.analysistype].targetsequence[target].replace('-', '')
        # In order to properly translate the nucleotide sequence, BioPython requests that the sequence is a multiple of
        # three - not partial codons. Trim the sequence accordingly
        remainder = 0 - len(sample[self.analysistype].targetsequence[target]) % 3
        seq = sample[self.analysistype].targetsequence[target] if remainder == 0 \
            else sample[self.analysistype].targetsequence[target][:remainder]
        # Set the DNA and protein sequences of the target in the sample
        sample[self.analysistype].dnaseq[target] = Seq(seq)
        # Translate the nucleotide sequence
        sample[self.analysistype].protseq[target] = str(sample[self.analysistype].dnaseq[target].translate())
        for targetfile in self.targetfiles:
            # Trim the reference sequence to multiples of three
            refremainder = 0 - len(self.records[targetfile][target].seq) % 3
            refseq = str(self.records[targetfile][target].seq) if refremainder % 3 == 0 \
                else str(self.records[targetfile][target].seq)[:refremainder]
            # Translate the nucleotide sequence of the reference sequence
            refdna = Seq(refseq)
            refprot = str(refdna.translate())
            # Use pairwise2 to perform a local alignment with the following parameters:
            # x     No match parameters. Identical characters have score of 1, otherwise 0.
            # s     Same open (-1)  and extend (-.1) gap penalties for both sequences
            ntalignments = pairwise2.align.localxs(seq, refseq, -1, -.1)
            # Use format_alignment to create a formatted alignment that is subsequently split on newlines e.g.
            '''
            ACCGT
            | ||
            A-CG-
            Score=3
            '''
            ntformat = (str(format_alignment(*ntalignments[0])).split('\n'))
            # Align the nucleotide sequence of the reference (ntalignments[2]) to the sample (ntalignments[0]).
            # If the corresponding bases match, add a |, otherwise a space
            ntalignment = ''.join(map(lambda x: '|' if len(set(x)) == 1 else ' ',
                                      zip(ntformat[0], ntformat[2])))
            # Create the nucleotide alignment: the sample sequence, the (mis)matches, and the reference sequence
            sample[self.analysistype].ntalign[target] = self.interleaveblastresults(ntformat[0], ntformat[2])
            # Regex to determine location of mismatches in the sequences
            count = 0
            sample[self.analysistype].ntindex[target] = str()
            for snp in re.finditer(' ', ntalignment):
                # If there are many SNPs, then insert line breaks for every 10 SNPs
                if count <= 10:
                    sample[self.analysistype].ntindex[target] += str(snp.start()) + ';'
                else:
                    sample[self.analysistype].ntindex[target] += '\n' + str(snp.start()) + ';'
                    count = 0
                count += 1
            # Perform the same steps, except for the amino acid sequence
            aaalignments = pairwise2.align.localxs(sample[self.analysistype].protseq[target], refprot, -1, -.1)
            aaformat = (str(format_alignment(*aaalignments[0])).split('\n'))
            aaalignment = ''.join(map(lambda x: '|' if len(set(x)) == 1 else ' ',
                                      zip(aaformat[0], aaformat[2])))
            sample[self.analysistype].aaidentity[target] = '{:.2f}'\
                .format(float(aaalignment.count('|')) / float(len(aaalignment)) * 100)
            sample[self.analysistype].aaalign[target] = self.interleaveblastresults(aaformat[0], aaformat[2])
            count = 0
            sample[self.analysistype].aaindex[target] = str()
            for snp in re.finditer(' ', aaalignment):
                if count <= 10:
                    sample[self.analysistype].aaindex[target] += str(snp.start()) + ';'
                else:
                    sample[self.analysistype].aaindex[target] += '\n' + str(snp.start()) + ';'
                    count = 0
                count += 1

    def resfinderreporter(self):
        """
        Custom reports for ResFinder analyses. These reports link the gene(s) found to their resistance phenotypes
        """

        target_dir = str()
        for folder in self.targetfolders:
            target_dir = folder
        genedict, altgenedict = ResistanceNotes.notes(target_dir)
        # Create a workbook to store the report. Using xlsxwriter rather than a simple csv format, as I want to be
        # able to have appropriately sized, multi-line cells
        workbook = xlsxwriter.Workbook(os.path.join(self.reportpath, '{}.xlsx'.format(self.analysistype)))
        # New worksheet to store the data
        worksheet = workbook.add_worksheet()
        # Add a bold format for header cells. Using a monotype font size 10
        bold = workbook.add_format({'bold': True, 'font_name': 'Courier New', 'font_size': 8})
        # Format for data cells. Monotype, size 10, top vertically justified
        courier = workbook.add_format({'font_name': 'Courier New', 'font_size': 8})
        courier.set_align('top')
        # Initialise the position within the worksheet to be (0,0)
        row = 0
        col = 0
        # A dictionary to store the column widths for every header
        columnwidth = dict()
        extended = False
        headers = ['Strain', 'Gene', 'Allele', 'Resistance', 'PercentIdentity', 'PercentCovered', 'Contig', 'Location',
                   'nt_sequence']
        for sample in self.metadata:
            sample[self.analysistype].sampledata = list()
            # Process the sample only if the script could find targets
            if sample[self.analysistype].blastresults != 'NA' and sample[self.analysistype].blastresults:
                for result in sample[self.analysistype].blastresults:
                    # Set the name to avoid writing out the dictionary[key] multiple times
                    name = result['subject_id']
                    # Use the ResistanceNotes gene name extraction method to get the necessary variables
                    gname, genename, accession, allele = ResistanceNotes.gene_name(name)
                    # Initialise a list to store all the data for each strain
                    data = list()
                    # Determine the name of the gene to use in the report and the resistance using the resistance
                    # method
                    finalgene, resistance = ResistanceNotes.resistance(gname, genename, genedict, altgenedict)
                    # Append the necessary values to the data list
                    data.append(finalgene)
                    data.append(allele)
                    data.append(resistance)
                    percentid = result['percentidentity']
                    data.append(percentid)
                    data.append(result['alignment_fraction'])
                    data.append(result['query_id'])
                    data.append('...'.join([str(result['low']), str(result['high'])]))
                    try:
                        # Only if the alignment option is selected, for inexact results, add alignments
                        if self.align and percentid != 100.00:

                            # Align the protein (and nucleotide) sequences to the reference
                            self.alignprotein(sample, name)
                            if not extended:
                                # Add the appropriate headers
                                headers.extend(['aa_Identity',
                                                'aa_Alignment',
                                                'aa_SNP_location',
                                                'nt_Alignment',
                                                'nt_SNP_location'
                                                ])
                                extended = True
                            # Create a FASTA-formatted sequence output of the query sequence
                            record = SeqRecord(sample[self.analysistype].dnaseq[name],
                                               id='{}_{}'.format(sample.name, name),
                                               description='')

                            # Add the alignment, and the location of mismatches for both nucleotide and amino
                            # acid sequences
                            data.extend([record.format('fasta'),
                                         sample[self.analysistype].aaidentity[name],
                                         sample[self.analysistype].aaalign[name],
                                         sample[self.analysistype].aaindex[name],
                                         sample[self.analysistype].ntalign[name],
                                         sample[self.analysistype].ntindex[name]
                                         ])
                        else:
                            record = SeqRecord(Seq(result['query_sequence']),
                                               id='{}_{}'.format(sample.name, name),
                                               description='')
                            data.append(record.format('fasta'))
                            if self.align:
                                # Add '-'s for the empty results, as there are no alignments for exact matches
                                data.extend(['-', '-', '-', '-', '-'])
                    # If there are no blast results for the target, add a '-'
                    except (KeyError, TypeError):
                        data.append('-')
                    sample[self.analysistype].sampledata.append(data)
        if 'nt_sequence' not in headers:
            headers.append('nt_sequence')
            # Write the header to the spreadsheet
        for header in headers:
            worksheet.write(row, col, header, bold)
            # Set the column width based on the longest header
            try:
                columnwidth[col] = len(header) if len(header) > columnwidth[col] else columnwidth[
                    col]
            except KeyError:
                columnwidth[col] = len(header)
            worksheet.set_column(col, col, columnwidth[col])
            col += 1
            # Increment the row and reset the column to zero in preparation of writing results
        row += 1
        col = 0
        # Write out the data to the spreadsheet
        for sample in self.metadata:
            if not sample[self.analysistype].sampledata:
                worksheet.write(row, col, sample.name, courier)
                # Increment the row and reset the column to zero in preparation of writing results
                row += 1
                col = 0
                # Set the width of the row to be the number of lines (number of newline characters) * 12
                worksheet.set_row(row)
                worksheet.set_column(col, col, columnwidth[col])
            for data in sample[self.analysistype].sampledata:
                columnwidth[col] = len(sample.name) + 2
                worksheet.set_column(col, col, columnwidth[col])
                worksheet.write(row, col, sample.name, courier)
                col += 1
                # List of the number of lines for each result
                totallines = list()
                for results in data:
                    #
                    worksheet.write(row, col, results, courier)
                    try:
                        # Counting the length of multi-line strings yields columns that are far too wide, only count
                        # the length of the string up to the first line break
                        alignmentcorrect = len(str(results).split('\n')[1])
                        # Count the number of lines for the data
                        lines = results.count('\n') if results.count('\n') >= 1 else 1
                        # Add the number of lines to the list
                        totallines.append(lines)
                    except IndexError:
                        try:
                            # Counting the length of multi-line strings yields columns that are far too wide, only count
                            # the length of the string up to the first line break
                            alignmentcorrect = len(str(results).split('\n')[0])
                            # Count the number of lines for the data
                            lines = results.count('\n') if results.count('\n') >= 1 else 1
                            # Add the number of lines to the list
                            totallines.append(lines)
                        # If there are no newline characters, set the width to the length of the string
                        except AttributeError:
                            alignmentcorrect = len(str(results))
                            lines = 1
                            # Add the number of lines to the list
                            totallines.append(lines)
                    # Increase the width of the current column, if necessary
                    try:
                        columnwidth[col] = alignmentcorrect if alignmentcorrect > columnwidth[col] else \
                            columnwidth[col]
                    except KeyError:
                        columnwidth[col] = alignmentcorrect
                    worksheet.set_column(col, col, columnwidth[col])
                    col += 1
                # Set the width of the row to be the number of lines (number of newline characters) * 12
                worksheet.set_row(row, max(totallines) * 11)
                # Increase the row counter for the next strain's data
                row += 1
                col = 0
        # Close the workbook
        workbook.close()

    def virulencefinderreporter(self):
        with open(os.path.join(self.reportpath, 'virulence.csv'), 'w') as report:
            header = 'Strain,Gene,PercentIdentity,PercentCovered,Contig,Location,Sequence\n'
            data = ''
            for sample in self.metadata:
                if sample.general.bestassemblyfile != 'NA':
                    if sample[self.analysistype].blastresults:
                        data += '{},'.format(sample.name)
                        #
                        multiple = False
                        for result in sample[self.analysistype].blastresults:
                            if self.analysistype == 'virulence':
                                gene = result['subject_id'].split(':')[0]
                            else:
                                gene = result['subject_id']
                            if multiple:
                                data += ','
                            data += '{},{},{},{},{}..{},{}\n' \
                                .format(gene, result['percentidentity'], result['alignment_fraction'],
                                        result['query_id'], result['low'], result['high'], result['query_sequence'])
                            # data += '\n'
                            multiple = True
                    else:
                        data += '{}\n'.format(sample.name)
                else:
                    data += '{}\n'.format(sample.name)
            report.write(header)
            report.write(data)

    @staticmethod
    def interleaveblastresults(query, subject):
        """
        Creates an interleaved string that resembles BLAST sequence comparisons
        :param query: Query sequence
        :param subject: Subject sequence
        :return: Properly formatted BLAST-like sequence comparison
        """
        # Initialise strings to hold the matches, and the final BLAST-formatted string
        matchstring = ''
        blaststring = ''
        # Iterate through the query
        for i, bp in enumerate(query):
            # If the current base in the query is identical to the corresponding base in the reference, append a '|'
            # to the match string, otherwise, append a ' '
            if bp == subject[i]:
                matchstring += '|'
            else:
                matchstring += ' '
        # Set a variable to store the progress through the sequence
        prev = 0
        # Iterate through the query, from start to finish in steps of 60 bp
        for j in range(0, len(query), 60):
            # BLAST results string. The components are: current position (padded to four characters), 'OLC', query
            # sequence, \n, matches, \n, 'ref', subject sequence. Repeated until all the sequence data are present.
            """
            0000 OLC ATGAAGAAGATATTTGTAGCGGCTTTATTTGCTTTTGTTTCTGTTAATGCAATGGCAGCT
                     ||||||||||| ||| | |||| ||||||||| || ||||||||||||||||||||||||
                 ref ATGAAGAAGATGTTTATGGCGGTTTTATTTGCATTAGTTTCTGTTAATGCAATGGCAGCT
            0060 OLC GATTGTGCAAAAGGTAAAATTGAGTTCTCTAAGTATAATGAGAATGATACATTCACAGTA
                     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
                 ref GATTGTGCAAAAGGTAAAATTGAGTTCTCTAAGTATAATGAGAATGATACATTCACAGTA
            """
            blaststring += '{} OLC {}\n         {}\n     ref {}\n' \
                .format('{:04d}'.format(j), query[prev:j + 60], matchstring[prev:j + 60], subject[prev:j + 60])
            # Update the progress variable
            prev = j + 60
        # Return the properly formatted string
        return blaststring

    def __init__(self, inputobject):

        self.metadata = inputobject.runmetadata.samples
        self.cutoff = inputobject.cutoff
        self.start = inputobject.start
        self.analysistype = inputobject.analysistype
        self.reportpath = inputobject.reportdir
        self.targetfolders = set()
        self.targetfiles = list()
        self.records = dict()
        self.pipeline = inputobject.pipeline
        self.referencefilepath = inputobject.referencefilepath
        self.cpus = inputobject.threads
        self.align = inputobject.align
        self.logfile = inputobject.logfile
        # self.resfinder = inputobject.resfinder
        # self.virulencefinder = inputobject.virulencefinder
        # If CGE-based analyses are specified, set self.unique to True, otherwise, use the arguments
        if self.analysistype == 'resfinder':
            self.unique = True
        elif self.analysistype == 'virulence':
            self.unique = True
            self.cutoff = 80
        else:
            self.unique = inputobject.unique
        # Fields used for custom outfmt 6 BLAST output:
        self.fieldnames = ['query_id', 'subject_id', 'positives', 'mismatches', 'gaps',
                           'evalue', 'bit_score', 'subject_length', 'alignment_length',
                           'query_start', 'query_end', 'query_sequence',
                           'subject_start', 'subject_end', 'subject_sequence']
        self.plusdict = defaultdict(make_dict)
        self.dqueue = Queue(maxsize=self.cpus)
        self.blastqueue = Queue(maxsize=self.cpus)


def sequencenames(contigsfile):
    """
    Takes a multifasta file and returns a list of sequence names
    :param contigsfile: multifasta of all sequences
    :return: list of all sequence names
    """
    sequences = list()
    for record in SeqIO.parse(open(contigsfile, "rU", encoding="iso-8859-15"), "fasta"):
        sequences.append(record.id)
    return sequences


if __name__ == '__main__':

    class Parser(object):

        def strainer(self):
            # Get the sequences in the sequences folder into a list. Note that they must have a file extension that
            # begins with .fa
            self.strains = sorted(glob(os.path.join(self.sequencepath, '*.fa*'.format(self.sequencepath))))
            self.targets = sorted(glob(os.path.join(self.targetpath, '*.tfa')))
            try:
                self.combinedtargets = glob(os.path.join(self.targetpath, '*.fasta'))[0]
            except IndexError:
                combinetargets(self.targets, self.targetpath)
                self.combinedtargets = glob(os.path.join(self.targetpath, '*.fasta'))[0]
            # Populate the metadata object. This object will be populated to mirror the objects created in the
            # genome assembly pipeline. This way this script will be able to be used as a stand-alone, or as part
            # of a pipeline
            assert self.strains, 'Could not find any files with an extension starting with "fa" in {}. Please check' \
                                 'to ensure that your sequence path is correct'.format(self.sequencepath)
            assert self.targets, 'Could not find any files with an extension starting with "fa" in {}. Please check' \
                                 'to ensure that your target path is correct'.format(self.targetpath)
            for sample in self.strains:
                # Create the object
                metadata = MetadataObject()
                # Set the base file name of the sequence. Just remove the file extension
                filename = os.path.splitext(os.path.split(sample)[1])[0]
                # Set the .name attribute to be the file name
                metadata.name = filename
                # Create the .general attribute
                metadata.general = GenObject()
                # Create the .mlst attribute
                setattr(metadata, self.analysistype, GenObject())
                # Set the .general.bestassembly file to be the name and path of the sequence file
                metadata.general.bestassemblyfile = sample
                metadata[self.analysistype].targets = self.targets
                metadata[self.analysistype].combinedtargets = self.combinedtargets
                metadata[self.analysistype].targetpath = self.targetpath
                metadata[self.analysistype].targetnames = sequencenames(self.combinedtargets)
                metadata[self.analysistype].reportdir = self.reportpath
                # Append the metadata for each sample to the list of samples
                self.samples.append(metadata)

        def __init__(self):
            parser = ArgumentParser(description='Use to find markers for any bacterial genome')
            parser.add_argument('--version',
                                action='version',
                                version='%(prog)s v0.5')
            parser.add_argument('-s', '--sequencepath',
                                required=True,
                                help='Specify input fasta folder')
            parser.add_argument('-t', '--targetpath',
                                required=True,
                                help='Specify folder of targets')
            parser.add_argument('-r', '--reportpath',
                                required=True,
                                help='Specify output folder for csv')
            parser.add_argument('-c', '--cutoff',
                                type=int,
                                default=70, help='Threshold for maximum unique bacteria for a single antibiotic')
            parser.add_argument('-n', '--numthreads',
                                type=int,
                                default=24,
                                help='Specify number of threads')
            parser.add_argument('-a', '--align',
                                action='store_true',
                                help='Optionally output alignments of genes with less than 100% identity to reference '
                                     'genes. This alignment will use amino acid sequences for both query and reference')
            parser.add_argument('-u', '--unique',
                                action='store_true',
                                help='Do not report multiple hits at the same location in a contig. Instead, store the'
                                     'best hit, and ignore the rest')
            parser.add_argument('-R', '--resfinder',
                                action='store_true',
                                help='Perform ResFinder-like analyses ')
            parser.add_argument('-v', '--virulencefinder',
                                action='store_true',
                                help='Perform VirulenceFinder-like analyses')
            args = parser.parse_args()
            self.sequencepath = os.path.join(args.sequencepath)
            assert os.path.isdir(self.sequencepath), 'Cannot locate sequence path as specified: {}'\
                .format(self.sequencepath)
            self.targetpath = os.path.join(args.targetpath)
            assert os.path.isdir(self.targetpath), 'Cannot locate target path as specified: {}'\
                .format(self.targetpath)
            self.reportpath = os.path.join(args.reportpath)
            make_path(self.reportpath)
            assert os.path.isdir(self.reportpath), 'Cannot locate report path as specified: {}'\
                .format(self.reportpath)
            self.cutoff = args.cutoff
            self.threads = args.numthreads
            self.align = args.align
            self.unique = args.unique
            self.resfinder = args.resfinder
            self.virulencefinder = args.virulencefinder
            self.strains = list()
            self.targets = list()
            self.combinedtargets = str()
            self.samples = list()
            self.logfile = os.path.join(self.sequencepath, 'log.txt')
            if self.resfinder:
                self.analysistype = 'resfinder'
            elif self.virulencefinder:
                self.analysistype = 'virulence'
            elif self.resfinder and self.virulencefinder:
                printtime('Cannot perform ResFinder and VirulenceFinder simultaneously. Please choose only one '
                          'of the -R and -v flags', self.start)
            else:
                self.analysistype = 'geneseekr'
            self.start = time.time()
            self.strainer()

    class MetadataInit(object):
        def __init__(self):
            # Run the parser
            self.runmetadata = Parser()
            # Get the appropriate variables from the metadata file
            self.start = self.runmetadata.start
            self.analysistype = self.runmetadata.analysistype
            self.cutoff = self.runmetadata.cutoff
            self.threads = int(self.runmetadata.threads)
            self.reportdir = self.runmetadata.reportpath
            self.pipeline = False
            self.referencefilepath = str()
            self.align = self.runmetadata.align
            self.unique = self.runmetadata.unique
            self.logfile = self.runmetadata.logfile
            # Run the analyses
            geneseekr = GeneSeekr(self)
            geneseekr.geneseekr()

    # Run the class
    MetadataInit()


class PipelineInit(object):
    def strainer(self):
        for sample in self.runmetadata.samples:
            if sample.general.bestassemblyfile != 'NA':
                setattr(sample, self.analysistype, GenObject())
                if self.genusspecific:
                    try:
                        genus = sample.general.closestrefseqgenus
                    except AttributeError:
                        genus = sample.general.referencegenus
                    # Allow Shigella to use the same targets as Escherichia
                    genus = genus if genus != 'Shigella' else 'Escherichia'
                    targetpath = os.path.join(self.referencefilepath, self.analysistype, genus)
                else:
                    targetpath = os.path.join(self.referencefilepath, self.analysistype)
                targets = glob(os.path.join(targetpath, '*.tfa'))
                targetcheck = glob(os.path.join(targetpath, '*.tfa'))
                if targetcheck:
                    try:
                        combinedtargets = glob(os.path.join(targetpath, '*.fasta'))[0]
                    except IndexError:
                        combinetargets(targets, targetpath)
                        combinedtargets = glob(os.path.join(targetpath, '*.fasta'))[0]
                    sample[self.analysistype].targets = targets
                    sample[self.analysistype].combinedtargets = combinedtargets
                    sample[self.analysistype].targetpath = targetpath
                    sample[self.analysistype].targetnames = sequencenames(combinedtargets)
                    sample[self.analysistype].reportdir = os.path.join(sample.general.outputdirectory,
                                                                       self.analysistype)
                else:
                    # Set the metadata file appropriately
                    sample[self.analysistype].targets = 'NA'
                    sample[self.analysistype].combinedtargets = 'NA'
                    sample[self.analysistype].targetpath = 'NA'
                    sample[self.analysistype].targetnames = 'NA'
                    sample[self.analysistype].reportdir = 'NA'
                    sample[self.analysistype].blastresults = 'NA'
            else:
                # Set the metadata file appropriately
                setattr(sample, self.analysistype, GenObject())
                sample[self.analysistype].targets = 'NA'
                sample[self.analysistype].combinedtargets = 'NA'
                sample[self.analysistype].targetpath = 'NA'
                sample[self.analysistype].targetnames = 'NA'
                sample[self.analysistype].reportdir = 'NA'
                sample[self.analysistype].blastresults = 'NA'

    def __init__(self, inputobject, analysistype, genusspecific, cutoff, unique):
        self.runmetadata = inputobject.runmetadata
        self.analysistype = analysistype
        self.path = inputobject.path
        self.start = inputobject.starttime
        self.referencefilepath = inputobject.reffilepath
        self.threads = inputobject.cpus
        self.reportdir = inputobject.reportpath
        self.cutoff = cutoff
        self.logfile = inputobject.logfile
        self.pipeline = True
        self.genusspecific = genusspecific
        self.align = False
        self.unique = unique
        # Get the alleles and profile into the metadata
        self.strainer()
