#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import GenObject, logstr, make_path, MetadataObject, \
    run_subprocess, write_to_logfile
from Bio.Sequencing.Applications import SamtoolsIndexCommandline
from Bio.Application import ApplicationError
from Bio import SeqIO
from threading import Lock, Thread
from argparse import ArgumentParser
from click import progressbar
from io import StringIO
import multiprocessing
from glob import glob
from time import time
import logging
import shutil
import os
import re

__author__ = 'adamkoziol'


class AssemblyEvaluation(object):

    def main(self):
        """
        Run the methods in the correct order
        """
        self.quast()
        self.parse_quast_report()
        self.sort_bam()
        self.indexing()
        self.qualimapper()
        self.parse_qualimap_report()
        self.clean_quast()
        self.extract_insert_size()
        self.calculate_weighted_insert_size()
        self.pilon()
        self.filter()
        self.clear()

    def quast(self):
        """
        Run quast on the samples
        """
        logging.info('Running Quast on assemblies')
        with progressbar(self.metadata) as bar:
            for sample in bar:
                # Create and populate the quast GenObject
                sample.quast = GenObject()
                sample.quast.outputdir = os.path.join(sample.general.outputdirectory, 'quast')
                sample.quast.report = os.path.join(sample.quast.outputdir, 'report.tsv')
                if sample.general.bestassemblyfile != "NA":
                    make_path(sample.quast.outputdir)
                    # Allow for non-paired samples
                    if len(sample.general.trimmedcorrectedfastqfiles) == 2:
                        sample.quast.cmd = 'quast --pe1 {forward} --pe2 {reverse}'\
                            .format(forward=sample.general.trimmedcorrectedfastqfiles[0],
                                    reverse=sample.general.trimmedcorrectedfastqfiles[1])
                    else:
                        sample.quast.cmd = 'quast --single {single}'\
                            .format(single=sample.general.trimmedcorrectedfastqfiles[0])
                    # Both paired and unpaired samples share the rest of the system call
                    # --debug is specified, as certain temporary files are either used for downstream analyses
                    # (BAM file), or parsed (insert size estimation)
                    sample.quast.cmd += ' -t {threads} --k-mer-stats --circos --rna-finding ' \
                                        '--conserved-genes-finding -o {outputdir} --debug {assembly}'\
                        .format(threads=self.cpus,
                                outputdir=sample.quast.outputdir,
                                assembly=sample.general.assemblyfile)
                    # Run the quast system call if the final quast report doesn't already exist
                    if not os.path.isfile(sample.quast.report):
                        out, err = run_subprocess(sample.quast.cmd)
                        # Write the appropriate information to the logfile
                        write_to_logfile(out='{cmd}\n{out}'.format(cmd=sample.quast.cmd,
                                                                   out=out),
                                         err=err,
                                         logfile=self.logfile,
                                         samplelog=sample.general.logout,
                                         sampleerr=sample.general.logerr,
                                         analysislog=None,
                                         analysiserr=None)

    def parse_quast_report(self):
        """
        Parse the quast report, and populate the metadata object with the extracted key: value pairs
        """
        for sample in self.metadata:
            if os.path.isfile(sample.quast.report):
                # Read in the report
                with open(sample.quast.report, 'r') as quast_report:
                    for line in quast_report:
                        # Sanitise the tab-delimited key: value pairs
                        key, value = self.analyze(line)
                        # Use the sanitised pair to populate the metadata object
                        setattr(sample.quast, key, value)

    def sort_bam(self):
        """
        Use samtools sort to sort the BAM files created by quast
        """
        logging.info('Sorting BAM files')
        with progressbar(self.metadata) as bar:
            for sample in bar:
                sample.quast.sortedbam = os.path.join(sample.quast.outputdir, '{sn}_sorted.bam'
                                                      .format(sn=sample.name))
                if sample.general.bestassemblyfile != 'NA' and not os.path.isfile(sample.quast.sortedbam):
                    try:
                        sample.quast.unfilteredbam = glob(os.path.join(sample.quast.outputdir, 'reads_stats',
                                                                       'temp_output', '*_unfiltered.bam'))[0]
                        sample.quast.samtools_sort_cmd = 'samtools sort --output-fmt BAM --reference {reference} ' \
                                                         '-@ {threads} -o {output} {input}'\
                            .format(reference=sample.general.assemblyfile,
                                    threads=self.cpus,
                                    output=sample.quast.sortedbam,
                                    input=sample.quast.unfilteredbam)
                        if not os.path.isfile(sample.quast.sortedbam):
                            out, err = run_subprocess(command=sample.quast.samtools_sort_cmd)
                            write_to_logfile(out='{cmd}\n{out}'.format(cmd=sample.quast.samtools_sort_cmd,
                                                                       out=out),
                                             err=err,
                                             logfile=self.logfile,
                                             samplelog=sample.general.logout,
                                             sampleerr=sample.general.logerr,
                                             analysislog=None,
                                             analysiserr=None)
                    except IndexError:
                        sample.quast.sortedbam = 'ND'
                        sample.quast.samtools_sort_cmd = str()

    def indexing(self):
        """
        Use samtools index to index the sorted BAM files
        """
        logging.info('Indexing sorted bam files')
        for i in range(self.cpus):
            threads = Thread(target=self.index, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                bamindex = SamtoolsIndexCommandline(input=sample.quast.sortedbam)
                sample.quast.sortedbai = sample.quast.sortedbam + '.bai'
                sample.quast.bamindex = str(bamindex)
                self.indexqueue.put((sample, bamindex))
        self.indexqueue.join()

    def index(self):
        while True:
            try:
                sample, bamindex = self.indexqueue.get()
                # Only make the call if the .bai file doesn't already exist
                if not os.path.isfile(sample.quast.sortedbai):
                    # Use cStringIO streams to handle output
                    stdout, stderr = map(StringIO, bamindex(cwd=sample.quast.outputdir))
                    if stderr:
                        # Write the standard error to log
                        with open(os.path.join(sample.quast.outputdir,
                                               'indexing_samtools_bam_index.log'), 'a+') as log:
                            log.writelines(logstr(bamindex, stderr.getvalue(), stdout.getvalue()))
                    stderr.close()
            except ApplicationError:
                pass
            self.indexqueue.task_done()

    def clean_quast(self):
        """
        Remove all the unnecessary temporary files created by quast
        """
        for sample in self.metadata:
            # Use os.walk to find all the files in the quast directory
            for path, directories, files in os.walk(sample.quast.outputdir):
                for quast_file in files:
                    # Set the absolute path of the file
                    file_path = os.path.join(path, quast_file)
                    # Remove all large files that do not have .err or .html extensions
                    if os.path.getsize(file_path) > 100000 and '_sorted.bam' not in quast_file and '.err' \
                            not in quast_file and '.html' not in quast_file:
                        os.remove(file_path)

    def extract_insert_size(self):
        """
        Parse the bwa index log information to extract insert size estimations
        """
        logging.info('Calculating insert size')
        for sample in self.metadata:
            sample.quast.reads_stats_file = os.path.join(sample.quast.outputdir, 'reads_stats', 'reads_stats.err')
            if os.path.isfile(sample.quast.reads_stats_file):
                # Initialise attributes for the insert size estimation
                sample.quast.total_reads = 0
                sample.quast.insert_mean = list()
                sample.quast.insert_std = list()
                sample.quast.read_blocks = list()
                current_reads = 0
                # Open the report
                with open(sample.quast.reads_stats_file, 'r') as read_stats:
                    for line in read_stats:
                        # BWA estimates the insert size distribution per 256*1024 read pairs. Extract the number of
                        # reads present in the current block being processed
                        # e.g. # candidate unique pairs for (FF, FR, RF, RR): (46, 226102, 14, 28)
                        if '# candidate unique pairs for' in line:
                            # Using the example above, current_reads will be 226102
                            current_reads = int(line.rstrip().replace(',', '').replace('(', '').replace(')', '')
                                                .split()[-3])
                            # Add current_reads to the total number of reads
                            sample.quast.total_reads += current_reads
                        # Continue parsing to find the FR section of the current block
                        elif 'analyzing insert size distribution for orientation FR' in line:
                            for sub_line in read_stats:
                                # Extract the mean and standard deviation of the insert size for this block
                                #  [M::mem_pestat] mean and std.dev: (487.88, 246.14)
                                if '[M::mem_pestat] mean and std.dev:' in sub_line:
                                    split_line = sub_line.rstrip().replace(',', '').replace('(', '').replace(')', '')\
                                        .split()
                                    mean = float(split_line[-2])
                                    std = float(split_line[-1])
                                    sample.quast.insert_mean.append(mean)
                                    sample.quast.insert_std.append(std)
                                    sample.quast.read_blocks.append(current_reads)
                                    # Break out of this loop
                                    break
            else:
                sample.quast.insert_mean = 'ND'
                sample.quast.insert_std = 'ND'

    def calculate_weighted_insert_size(self):
        """
         Calculate the weighted mean and standard deviation of the insert size from the extracted bwa blocks
        """
        for sample in self.metadata:
            # Initialise attributes to store the calculated weighted mean and standard deviation of insert sizes
            sample.quast.mean_insert = float()
            sample.quast.std_insert = float()
            if sample.quast.insert_mean != 'ND':
                # Iterate through all the read blocks present in the sample (read_block is the number of reads present
                # in the current block)
                for i, read_block in enumerate(sample.quast.read_blocks):
                    # Calculate the weight of the current block by dividing it (current number of reads) by the total
                    # number of reads in the sample
                    weight = read_block / sample.quast.total_reads
                    # Multiply the mean for this block to obtain the weighted mean, and add it to the total mean
                    sample.quast.mean_insert += sample.quast.insert_mean[i] * weight
                    # Same calculation, but for standard deviation
                    sample.quast.std_insert += sample.quast.insert_std[i] * weight
            # Set the attributes to floats with two decimal places
            sample.quast.mean_insert = float('{:.2f}'.format(sample.quast.mean_insert))
            sample.quast.std_insert = float('{:.2f}'.format(sample.quast.std_insert))

    def qualimapper(self):
        """
        Create threads and commands for performing reference mapping for qualimap analyses
        """
        logging.info('Running qualimap on samples')
        for i in range(self.cpus):
            # Send the threads to the merge method. :args is empty as I'm using
            threads = Thread(target=self.qualimap, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Create and populate the qualimap attribute
            sample.qualimap = GenObject()
            sample.qualimap.outputdir = os.path.join(sample.general.outputdirectory, 'qualimap')
            make_path(sample.qualimap.outputdir)
            sample.qualimap.reportfile = os.path.join(sample.qualimap.outputdir, 'genome_results.txt')
            sample.qualimap.length = dict()
            sample.qualimap.bases = dict()
            sample.qualimap.coverage = dict()
            sample.qualimap.stddev = dict()
            if sample.general.bestassemblyfile != "NA":
                self.qualimap_queue.put(sample)
        self.qualimap_queue.join()

    def qualimap(self):
        """
        Run qualimap
        """
        while True:
            sample = self.qualimap_queue.get()
            if sample.general.bestassemblyfile != "NA":

                # Define the Qualimap call
                qualimapcall = 'qualimap bamqc -bam {sorted_bam} -outdir {outdir}' \
                    .format(sorted_bam=sample.quast.sortedbam,
                            outdir=sample.qualimap.outputdir)
                sample.commands.qualimap = qualimapcall
                # If the report file doesn't exist, run Qualimap, and print logs to the log file
                if not os.path.isfile(sample.qualimap.reportfile):
                    out, err = run_subprocess(sample.commands.qualimap)
                    write_to_logfile(out=sample.commands.qualimap,
                                     err=sample.commands.qualimap,
                                     logfile=self.logfile,
                                     samplelog=sample.general.logout,
                                     sampleerr=sample.general.logerr,
                                     analysislog=None,
                                     analysiserr=None)
                    write_to_logfile(out=out,
                                     err=err,
                                     logfile=self.logfile,
                                     samplelog=sample.general.logout,
                                     sampleerr=sample.general.logerr,
                                     analysislog=None,
                                     analysiserr=None)
            self.qualimap_queue.task_done()

    def parse_qualimap_report(self):
        """
        Parse the qualimap report
        """
        for sample in self.metadata:
            if sample.general.bestassemblyfile != "NA":
                # Initialise a dictionary to hold the Qualimap results
                qdict = dict()
                try:
                    with open(sample.qualimap.reportfile) as report:
                        # Read the report
                        for line in report:
                            # Sanitise the keys and values using self.analyze
                            key, value = self.qualimap_analyze(line)
                            # If the keys and values exist, enter them into the dictionary
                            if (key, value) != (None, None):
                                # Only keep two decimal places
                                if type(value) is float:
                                    value = float('{:.2f}'.format(value))
                                qdict[key] = value
                            if 'Coverage per contig' in line:
                                for contigline in report:
                                    try:
                                        _, name, length, bases, coverage, stddev = contigline.rstrip().split('\t')
                                        sample.qualimap.length.update({name: length})
                                        sample.qualimap.bases.update({name: bases})
                                        sample.qualimap.coverage.update({name: coverage})
                                        sample.qualimap.stddev.update({name: stddev})
                                    except ValueError:
                                        pass

                except (IOError, FileNotFoundError):
                    pass
                # If there are values in the dictionary
                if qdict:
                    # Make new category for Qualimap results and populate this category with the report data
                    for attribute in qdict:
                        # Remove the 'X' from the depth values e.g. 40.238X
                        setattr(sample.qualimap, attribute, qdict[attribute].rstrip('X'))

    def pilon(self):
        """
        Run pilon to fix any misassemblies in the contigs - will look for SNPs and indels
        """
        logging.info('Improving quality of assembly with pilon')
        for i in range(self.threads):
            # Send the threads to the merge method. :args is empty
            threads = Thread(target=self.pilonthreads, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            sample.pilon = GenObject()
            if sample.general.bestassemblyfile != 'NA':
                # Set the name of the unfiltered assembly output file
                sample.general.contigsfile = sample.general.assemblyfile
                sample.pilon.outdir = os.path.join(sample.quast.outputdir, 'pilon')
                make_path(sample.quast.outputdir)
                # Create the command line command
                sample.pilon.cmd = 'pilon --genome {raw_assembly} --bam {sorted_bam} --fix bases ' \
                                   '--threads {threads} --outdir {outdir} --changes --mindepth 0.25' \
                    .format(raw_assembly=sample.general.contigsfile,
                            sorted_bam=sample.quast.sortedbam,
                            threads=self.threads,
                            outdir=sample.pilon.outdir)
                self.pilonqueue.put(sample)
        self.pilonqueue.join()

    def pilonthreads(self):
        while True:
            sample = self.pilonqueue.get()
            # Only perform analyses if the output file doesn't already exist
            if not os.path.isfile(sample.general.contigsfile):
                command = sample.quast.cmd
                out, err = run_subprocess(command)
                write_to_logfile(out=command,
                                 err=command,
                                 logfile=self.logfile,
                                 samplelog=sample.general.logout,
                                 sampleerr=sample.general.logerr,
                                 analysislog=None,
                                 analysiserr=None)
                write_to_logfile(out=out,
                                 err=err,
                                 logfile=self.logfile,
                                 samplelog=sample.general.logout,
                                 sampleerr=sample.general.logerr,
                                 analysislog=None,
                                 analysiserr=None)
            self.pilonqueue.task_done()

    def filter(self):
        """
        Filter contigs based on depth and length
        """
        logging.info('Filtering contigs')
        for i in range(self.cpus):
            # Send the threads to the filter method
            threads = Thread(target=self.filterthreads, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Set the name of the unfiltered assembly output file
            if sample.general.bestassemblyfile != 'NA':
                sample.general.contigsfile = sample.general.assemblyfile
                self.filterqueue.put(sample)
        self.filterqueue.join()

    def filterthreads(self):
        while True:
            sample = self.filterqueue.get()
            # Only run on samples that have been assembled
            if os.path.isfile(sample.general.contigsfile) and not os.path.isfile(sample.general.filteredfile):
                # Create a list to store all the records of contigs that pass the minimum depth filtering
                passdepth = list()
                for record in SeqIO.parse(open(sample.general.contigsfile, "rU"), "fasta"):
                    # Extract the values for the mean coverage depth
                    coveragemean = float(sample.qualimap.sample.qualimap.coverage[record.id])
                    coveragestd = float(sample.qualimap.stddev[record.id])
                    # Remove the _pilon added to the contig name in order to allow the contig name to match the original
                    # name used as the key in the sample.qualimap.coverage dictionary
                    contig = record.id.split('_pilon')[0]
                    # Only include contigs with a depth greater or equal to the mean coverage minus 1.5 times the
                    # coverage standard deviation
                    if float(sample.qualimap.coverage[contig]) > (coveragemean - coveragestd * 1.5) \
                            and len(record.seq) > 500:
                        # Replace 'NODE' in the fasta header with the sample name
                        # >NODE_1_length_705814_cov_37.107_ID_4231
                        newid = re.sub("Contig", sample.name, record.id)
                        record.id = str(record.id).replace('Contig', sample.name)
                        record.id = newid
                        # Clear the name and description attributes of the record
                        record.name = ''
                        record.description = ''
                        # Add this record to our list
                        passdepth.append(record)
                # Only create the file if there are contigs that pass the depth filter
                if passdepth:
                    # Open the filtered assembly file
                    with open(sample.general.filteredfile, 'w') as formatted:
                        # Write the records in the list to the file
                        SeqIO.write(passdepth, formatted, 'fasta')
            # If the filtered file was successfully created, copy it to the BestAssemblies folder
            if os.path.isfile(sample.general.filteredfile):
                # Set the assemblies path
                sample.general.bestassembliespath = os.path.join(self.path, 'BestAssemblies')
                # Set the name of the file in the best assemblies folder
                bestassemblyfile = os.path.join(sample.general.bestassembliespath, '{sn}.fasta'.format(sn=sample.name))
                # Add the name and path of the best assembly file to the metadata
                sample.general.bestassemblyfile = bestassemblyfile
                # Copy the filtered file to the BestAssemblies folder
                if not os.path.isfile(bestassemblyfile):
                    shutil.copyfile(sample.general.filteredfile, bestassemblyfile)
            else:
                sample.general.bestassemblyfile = 'NA'
            self.filterqueue.task_done()

    def clear(self):
        """
        Clear out large attributes from the metadata objects
        """
        for sample in self.metadata:
            try:
                delattr(sample.qualimap, 'bases')
                delattr(sample.qualimap, 'coverage')
                delattr(sample.qualimap, 'length')
                delattr(sample.qualimap, 'stddev')
            except AttributeError:
                pass

    @staticmethod
    def qualimap_analyze(line):
        # Split on ' = '
        if ' = ' in line:
            key, value = line.split(' = ')
            # Replace occurrences of
            key = key.replace('number of ', "").replace("'", "").title().replace(" ", "")
            # Should we keep comma separation?
            value = value.replace(",", "").replace(" ", "").rstrip()
        # Otherwise set the keys and values to None
        else:
            key, value = None, None
        return key, value

    @staticmethod
    def analyze(line):
        key, value = line.rstrip().split('\t')
        key = key.replace(' (%)', '').replace(' ', '_').replace('#', 'num').replace('(>=', 'greater_than')\
            .replace(')', '').replace('.', '').replace('\'', '')
        return key, value

    def __init__(self, inputobject):
        from queue import Queue
        self.metadata = inputobject.runmetadata.samples
        self.start = inputobject.starttime
        self.cpus = inputobject.cpus
        self.threadlock = Lock()
        try:
            self.threads = int(self.cpus / len(self.metadata)) if self.cpus / len(self.metadata) > 1 else 1
        except TypeError:
            self.threads = self.cpus
        self.logfile = inputobject.logfile
        self.path = inputobject.path
        # Initialise queues
        self.qualimap_queue = Queue(maxsize=self.cpus)
        self.pilonqueue = Queue(maxsize=self.cpus)
        self.indexqueue = Queue(maxsize=self.cpus)
        self.filterqueue = Queue(maxsize=self.cpus)


if __name__ == '__main__':
    class Parser(object):

        def associate(self):
            # Get the sequences in the sequences folder into a list. Note that they must have a file extension that
            # begins with .fa
            self.strains = [fasta for fasta in sorted(glob(os.path.join(self.assemblypath, '*.fa*')))
                            if '.fastq' not in fasta]
            for strain in self.strains:
                # Extract the name of the strain from the path and file extension
                strainname = os.path.split(strain)[1].split('.')[0]
                # Find the corresponding fastq files for each strain
                fastq = sorted(glob(os.path.join(self.fastqpath, '{sn}*fastq*'.format(sn=strainname))))
                # Ensure that fastq files are present for each assembly
                assert fastq, 'Cannot find fastq files for strain {sn}'.format(sn=strainname)
                # Create the object
                metadata = MetadataObject()
                # Set the .name attribute to be the file name
                metadata.name = strainname
                # Create the .general attribute
                metadata.general = GenObject()
                # Set the path of the assembly file
                metadata.general.bestassembliespath = self.assemblypath
                # Populate the .fastqfiles category of :self.metadata
                metadata.general.trimmedfastqfiles = fastq
                # Create the output directory path
                metadata.general.outputdirectory = os.path.join(self.path, strainname)
                metadata.mapping = GenObject()
                # Append the metadata for each sample to the list of samples
                self.samples.append(metadata)

        def __init__(self):
            parser = ArgumentParser(description='Calculates coverage depth by mapping FASTQ reads against assemblies')
            parser.add_argument('-p', '--path',
                                default=os.getcwd(),
                                help='Specify the path of the folder that either contains the files of interest, or'
                                     'will be used to store the outputs')
            parser.add_argument('-a', '--assemblies',
                                help='Path to a folder of assemblies. If not provided, the script will look for .fa'
                                     'or .fasta files in the path')
            parser.add_argument('-f', '--fastq',
                                help='Path to a folder of fastq files. If not provided, the script will look for '
                                     'fastq or .fastq.gz files in the path')
            parser.add_argument('-t', '--threads',
                                help='Number of threads. Default is the number of cores in the system')
            # Get the arguments into an object
            args = parser.parse_args()
            # Define variables from the arguments - there may be a more streamlined way to do this
            # Add trailing slashes to the path variables to ensure consistent formatting (os.path.join)
            self.path = os.path.join(args.path, '')
            self.assemblypath = os.path.join(args.assemblies, '') if args.assemblies else self.path
            self.fastqpath = os.path.join(args.fastq, '') if args.fastq else self.path
            # Use the argument for the number of threads to use, or default to the number of cpus in the system
            self.cpus = args.threads if args.threads else multiprocessing.cpu_count()
            # Initialise variables
            self.strains = list()
            self.samples = list()
            self.logfile = os.path.join(self.path, 'logfile.txt')

            # Associate the assemblies and fastq files in a metadata object
            self.associate()

    class MetadataInit(object):
        def __init__(self, start):
            # Run the parser
            self.runmetadata = Parser()
            # Get the appropriate variables from the metadata file
            self.path = self.runmetadata.path
            self.assemblypath = self.runmetadata.assemblypath
            self.fastqpath = self.runmetadata.fastqpath
            self.starttime = start
            self.cpus = self.runmetadata.cpus
            self.logfile = self.runmetadata.logfile
            # Run the analyses
            AssemblyEvaluation(self)

    # Run the class
    starttime = time()
    MetadataInit(start=starttime)
