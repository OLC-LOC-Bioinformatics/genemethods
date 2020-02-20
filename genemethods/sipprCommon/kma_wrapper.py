#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject, run_subprocess, \
    write_to_logfile
from genemethods.geneseekr.parser import Parser
from click import progressbar
import logging
import csv
import os
__author__ = 'adamkoziol'


class KMA(object):

    def main(self):
        self.targets()
        self.index_targets()
        self.load_kma_db()
        self.run_kma_mem_mode()
        self.unload_kma_db()
        self.parse_kma_outputs()

    def targets(self):
        """
        Use the parser class from geneseekr to find all the appropriate target files
        """
        parse = Parser(args=self,
                       pipeline=self.pipeline)
        if not self.genus_specific:
            parse.target_find()
        parse.metadata_populate()

    def index_targets(self):
        """
        Index the combined target FASTA file for KMA analyses if necessary
        """
        logging.info('Indexing databases for KMA analyses')
        for sample in self.metadata:
            # Set the name of the database by stripping of the extension of the combined FASTA file
            sample[self.analysistype].db_no_ext = os.path.splitext(sample[self.analysistype].combinedtargets)[0]
            if not os.path.isfile(sample[self.analysistype].db_no_ext + '.length.b'):
                kma_index_cmd = 'kma index -i {combined} -o {kma_db} -k {kmer} -k_t {kmer} -k_i {kmer}'\
                    .format(combined=sample[self.analysistype].combinedtargets,
                            kma_db=sample[self.analysistype].db_no_ext,
                            kmer=self.kmer)
                out, err = run_subprocess(command=kma_index_cmd)
                write_to_logfile(out='{cmd}\n{out}'.format(cmd=kma_index_cmd,
                                                           out=out),
                                 err=err,
                                 logfile=self.logfile)

    def load_kma_db(self):
        """
        Load all KMA databases into memory
        """
        logging.info('Loading {at} databases into memory for KMA analyses'.format(at=self.analysistype))
        for sample in self.metadata:
            # For genus-specific databases, each db needs to be loaded into memory
            if sample[self.analysistype].db_no_ext not in self.loaded_dbs:
                # Add the database to the list of loaded databases
                self.loaded_dbs.append(sample[self.analysistype].db_no_ext)
                # Load the database into memory at the appropriate level:
                # DB piece: Flag
                # *.comp.b: 1
                # *.decon.comp.b: 2
                # *.length.b: 4
                # *.seq.b *.index.b: 8
                # *.name: 16
                kma_db_load_cmd = 'kma shm -t_db {kma_db} -shmLvl {level}'\
                    .format(kma_db=sample[self.analysistype].db_no_ext,
                            level=self.level)
                run_subprocess(command=kma_db_load_cmd)

    def run_kma_mem_mode(self):
        """
        Use KMA with memory mode enabled
        """
        logging.info('Running {at} analyses with KMA with memory mode enabled'.format(at=self.analysistype))
        with progressbar(self.metadata) as bar:
            for sample in bar:
                if sample.general.bestassemblyfile != 'NA':
                    sample[self.analysistype].outputdir = os.path.join(sample.general.outputdirectory,
                                                                       self.analysistype)
                    sample[self.analysistype].output_prefix_mem_mode = \
                        os.path.join(sample[self.analysistype].outputdir, '{sn}_{at}_mem_mode'
                                     .format(sn=sample.name,
                                             at=self.analysistype))
                    sample[self.analysistype].kma_report_mem_mode = \
                        sample[self.analysistype].output_prefix_mem_mode + '.res'
                    sample[self.analysistype].kma_fasta_mem_mode = \
                        sample[self.analysistype].output_prefix_mem_mode + '.fsa'
                    # Determine whether the input FASTQ files are paired
                    prefix = self.prefix(sample=sample)
                    # Create the system call to KMA
                    # -ck: Count kmers instead of pseudo alignment
                    # -ID {percid}: Minimum percent identity is self.cutoff
                    # -ConClave 2: ConClave 2 allows for several closely related templates to present, which limits the
                    # assumptions taken by ConClave 1 at the cost of false positives.
                    # -1t1: Skip HMM, Force each query sequence to match to only one template, strictly global
                    # -ex_mode: Searh kmers exhaustively
                    # -boot: Bootstrap the query sequences, by subsampling from them
                    # -mem_mode	Use kmers to choose best template, and save memory. *.index and *.seq are not loaded
                    # into memory, which enables one to map against larger databases. Templates are chosen using k-mer
                    # counting.
                    # -mp 5: Minimum phred score of 5
                    # -shm {level}: Force KMA to use shared memory (db loaded above) at the requested level.
                    # Default is 4
                    sample[self.analysistype].kma_cmd_mem_mode = \
                        prefix + ' -o {output} -t {cpus} -t_db {db} -ck -ID {percid} -ConClave 2 -1t1 -ex_mode -boot ' \
                                 '-mem_mode -mp 5 -shm {level}'.format(
                            output=sample[self.analysistype].output_prefix_mem_mode,
                            cpus=self.threads,
                            percid=self.cutoff,
                            db=sample[self.analysistype].db_no_ext,
                            level=self.level)
                    # Add the additional arguments to the KMA call
                    if self.kma_kwargs:
                        sample[self.analysistype].kma_cmd_mem_mode += self.kma_kwargs
                    if not os.path.isfile(sample[self.analysistype].kma_report_mem_mode):
                        out, err = run_subprocess(command=sample[self.analysistype].kma_cmd_mem_mode)
                        write_to_logfile(out='{cmd}\n{out}'.format(cmd=sample[self.analysistype].kma_cmd_mem_mode,
                                                                   out=out),
                                         err=err,
                                         logfile=self.logfile,
                                         samplelog=sample.general.logout,
                                         sampleerr=sample.general.logerr,
                                         analysislog=sample[self.analysistype].log,
                                         analysiserr=sample[self.analysistype].log)

    def run_kma(self):
        """
        Run KMA on the samples with memory mode disabled
        """
        logging.info('Running {at} analyses with KMA'.format(at=self.analysistype))
        with progressbar(self.metadata) as bar:
            for sample in bar:
                if sample.general.bestassemblyfile != 'NA':
                    sample[self.analysistype].outputdir = os.path.join(sample.general.outputdirectory,
                                                                       self.analysistype)
                    sample[self.analysistype].output_prefix = os.path.join(sample[self.analysistype].outputdir,
                                                                           '{sn}_{at}'
                                                                           .format(sn=sample.name,
                                                                                   at=self.analysistype))
                    sample[self.analysistype].kma_report = sample[self.analysistype].output_prefix + '.res'
                    # Determine whether the FASTQ files are paired, and retrieve the prefix of the KMA system call
                    prefix = self.prefix(sample=sample)
                    # Create the system call. Note that I removed -ConClave 2 from the system call, as it was causing
                    # segmentation faults  some of the time when I tried to process raw FASTQ files
                    sample[self.analysistype].kma_cmd = prefix + ' -o {output} -t {cpus} -t_db {db} -mp 5 ' \
                                                                 '-ID {percid} -1t1 -ex_mode -boot' \
                        .format(
                                output=sample[self.analysistype].output_prefix,
                                cpus=self.threads,
                                percid=self.cutoff,
                                db=sample[self.analysistype].db_no_ext)
                    # Add the additional arguments to the KMA call
                    if self.kma_kwargs:
                        sample[self.analysistype].kma_cmd += self.kma_kwargs
                    # Run the analyses if the report does not exist
                    if not os.path.isfile(sample[self.analysistype].kma_report):
                        out, err = run_subprocess(command=sample[self.analysistype].kma_cmd)
                        write_to_logfile(out='{cmd}\n{out}'.format(cmd=sample[self.analysistype].kma_cmd, out=out),
                                         err=err, logfile=self.logfile, samplelog=sample.general.logout,
                                         sampleerr=sample.general.logerr, analysislog=sample[self.analysistype].log,
                                         analysiserr=sample[self.analysistype].log)

    def prefix(self, sample):
        """
        Determine whether the input FASTQ files are paired, and create the prefix of the system call to KMA to reflect
        the paired/unpaired nature
        """
        make_path(sample[self.analysistype].outputdir)
        sample[self.analysistype].log = os.path.join(sample[self.analysistype].outputdir, 'log')
        if self.datatype == 'raw':
            if len(sample.general.trimmedcorrectedfastqfiles) == 2:
                prefix = 'kma -ipe {forward} {reverse}' \
                    .format(forward=sample.general.trimmedcorrectedfastqfiles[0],
                            reverse=sample.general.trimmedcorrectedfastqfiles[1])
            else:
                prefix = 'kma -i {forward}' \
                    .format(forward=sample.general.trimmedcorrectedfastqfiles[0])
        else:
            prefix = 'kma -i {forward}' \
                    .format(forward=sample.general.bestassemblyfile)
        return prefix

    def unload_kma_db(self):
        """
        Remove all loaded KMA databases from memory
        """
        logging.info('Removing {at} databases from memory'.format(at=self.analysistype))
        for db in self.loaded_dbs:
            # Add the -destroy option to remove the database from memory
            unload_kma_db_cmd = 'kma shm -t_db {kma_db} -shmLvl {level} -destroy'\
                .format(kma_db=db,
                        level=self.level)
            run_subprocess(command=unload_kma_db_cmd)

    def parse_kma_outputs(self):
        """
        Parse the KMA outputs, and populate the GenObject to mimic the outputs from BLAST/sippr-based analyses
        """
        for sample in self.metadata:
            sample[self.analysistype].kmaresults = dict()
            sample[self.analysistype].kmadepthresults = dict()
            if sample.general.bestassemblyfile != 'NA':
                # Try to run the parsing method on both the memory mode-enabled and disabled analyses
                try:
                    if os.path.isfile(sample[self.analysistype].kma_report_mem_mode):
                        self.populate_outputs(sample=sample,
                                              report=sample[self.analysistype].kma_report_mem_mode)
                except AttributeError:
                    sample[self.analysistype].kma_report_mem_mode = str()
                try:
                    if os.path.isfile(sample[self.analysistype].kma_report):
                        self.populate_outputs(sample=sample,
                                              report=sample[self.analysistype].kma_report)
                except AttributeError:
                    sample[self.analysistype].kma_report = str()

    def populate_outputs(self, sample, report):
        """
        Parse the supplied KMA report (memory mode either enabled or disabled )for the sample
        """
        # Load the file into a dictionary using the csv library
        csv_dict = csv.DictReader(open(report), delimiter='\t')
        if sample.name not in self.kma_outputs:
            self.kma_outputs[sample.name] = list()
        genes_present = list()
        # Iterate through all the hits
        for kma_dict in csv_dict:
            if not self.headers:
                self.headers = [key for key in kma_dict]
            # Extract the necessary variables from the dictionary
            gene_allele = kma_dict['#Template']
            depth = float(kma_dict['Depth'].lstrip().rstrip())
            percent_identity = float(kma_dict['Template_Identity'].lstrip().rstrip())
            if self.analysistype == 'cgmlst':
                split_name = gene_allele.split('_')
                del split_name[-1]
                gene = '_'.join(split_name)
            else:
                gene = str()
            # Ensure that the hit passes the supplied depth and identity thresholds
            if percent_identity >= self.cutoff and depth >= self.averagedepth and gene not in genes_present:
                self.kma_outputs[sample.name].append(kma_dict)
                # Initialise the dictionary as necessary
                if gene_allele not in sample[self.analysistype].kmaresults:
                    sample[self.analysistype].kmaresults[gene_allele] = percent_identity
                    sample[self.analysistype].kmadepthresults[gene_allele] = {percent_identity: depth}
                else:
                    if percent_identity >= sample[self.analysistype].kmaresults[gene_allele]:
                        sample[self.analysistype].kmaresults[gene_allele] = percent_identity
                        sample[self.analysistype].kmadepthresults[gene_allele] = {percent_identity: depth}
                # Add the gene to the list of genes only for cgMLST analyses
                if self.analysistype == 'cgmlst':
                    genes_present.append(gene)

    def kma_report(self):
        """
        Write the raw KMA results to a summary file
        """
        data = str()
        with open(os.path.join(self.reportpath, '{at}_kma_outputs.csv'.format(at=self.analysistype)), 'w') as report:
            for sample, kma_dict in self.kma_outputs.items():
                for result in kma_dict:
                    data += '{sn},'.format(sn=sample)
                    for header in self.headers:
                        data += '{value},'.format(value=result[header])
                    data += '\n'
            report.write('Sample,{header}\n'.format(header=','.join(header.replace('#', '') for header in
                                                                    self.headers)))
            report.write(data)

    def __init__(self, args, pipeline, analysistype, datatype='raw', cutoff=98, averagedepth=2, k=16, level=4,
                 kma_kwargs=None):
        logging.info('Running {at} analyses with KMA'.format(at=analysistype))
        self.metadata = args.runmetadata.samples
        self.path = args.path
        self.sequencepath = args.path
        self.reportpath = args.reportpath
        make_path(self.reportpath)
        self.analysistype = analysistype
        self.targetpath = args.reffilepath
        self.pipeline = pipeline
        self.start = args.starttime
        self.threads = args.cpus
        self.cpus = self.threads
        self.homepath = args.homepath
        self.datatype = datatype
        self.cutoff = cutoff
        self.averagedepth = averagedepth
        self.kmer = k
        self.level = level
        self.kma_kwargs = kma_kwargs
        self.taxonomy = {'Escherichia': 'coli', 'Listeria': 'monocytogenes', 'Salmonella': 'enterica'}
        if self.analysistype == 'mlst':
            self.targetpath = os.path.join(self.targetpath, 'MLST')
            self.genus_specific = True
        elif self.analysistype == 'rmlst':
            self.targetpath = os.path.join(self.targetpath, 'rMLST')
            self.genus_specific = False
        elif self.analysistype == 'cgmlst':
            self.targetpath = os.path.join(self.targetpath, 'cgMLST')
            self.genus_specific = True
        else:
            self.targetpath = os.path.join(self.targetpath, analysistype)
            self.genus_specific = False
        self.logfile = os.path.join(self.path, 'log')
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = self.metadata
        self.loaded_dbs = list()
        self.kma_outputs = dict()
        self.headers = list()
