#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject, run_subprocess, \
    write_to_logfile
from genemethods.MLSTsippr.mlst import GeneSippr as MLSTSippr
from genemethods.sipprCommon.kma_wrapper import KMA
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from glob import glob
import hashlib
import logging
import os

__author__ = 'adamkoziol'


class KMAMLST(KMA):

    def main(self):
        self.targets()
        self.index_targets()
        self.load_kma_db()
        self.run_kma_mem_mode()
        self.unload_kma_db()
        self.parse_kma_outputs()
        self.typing()

    def typing(self):
        # Create the typing object
        typing = MLST(args=self,
                      pipelinecommit='',
                      startingtime=self.start,
                      scriptpath='',
                      analysistype=self.analysistype,
                      cutoff=1.0,
                      pipeline=self.pipeline)
        # Perform typing, and create reports
        typing.reporter()

    def __init__(self, args, pipeline, analysistype='mlst', datatype='raw', cutoff=100, averagedepth=2, k=16, level=4,
                 kma_kwargs=None):
        logging.info('Running {at} analyses with KMA'.format(at=analysistype))
        self.metadata = args.runmetadata.samples
        self.path = args.path
        self.sequencepath = args.path
        self.reportpath = args.reportpath
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
        self.logfile = os.path.join(self.path, 'log')
        self.runmetadata = MetadataObject()
        self.runmetadata.samples = self.metadata
        self.loaded_dbs = list()
        self.kma_outputs = dict()
        self.headers = list()
        self.new_allele_dict = dict()


class MLST(MLSTSippr):

    def reporter(self):
        if 'rmlst' in self.analysistype.lower():
            analysistype = 'rmlst'
        elif 'cgmlst' in self.analysistype.lower():
            analysistype = 'cgmlst'
        else:
            analysistype = 'mlst'
        # Populate self.plusdict in order to reuse parsing code from an assembly-based method
        for sample in self.runmetadata.samples:
            self.plusdict[sample.name] = dict()
            self.matchdict[sample.name] = dict()

            sample[self.analysistype].alleles = sorted(list(set(allele.split('_')[0]
                                                                for allele in sample[self.analysistype]
                                                                .targetnames)))
            # In order to work with the Enterobase cgMLST scheme that has underscores in the gene names (e.g.
            # AEJV01_03887, check for multiple underscores in the allele name, and act appropriately
            if len(sample[self.analysistype].alleles) > 53:
                allele_set = set()
                for allele in sample[self.analysistype].targetnames:
                    allele_set.add(allele)
                sample[self.analysistype].alleles = sorted(list(allele_set))
            # Allele names attribute is apparently the same as the alleles attribute
            sample[self.analysistype].allelenames = sample[self.analysistype].alleles
            try:
                sample[self.analysistype].profile = glob(os.path.join(sample[self.analysistype].targetpath,
                                                                      '*.txt'))[0]
            except IndexError:
                sample[self.analysistype].profile = 'ND'
            if sample.general.bestassemblyfile != 'NA':
                for gene in sample[analysistype].targetnames:
                    self.plusdict[sample.name][gene] = dict()
                    for allele, depth_dict in sample[self.analysistype].kmadepthresults.items():
                        for percent_id, depth in depth_dict.items():
                            if gene in allele:
                                # Split the allele number from the gene name using the appropriate delimiter
                                if '_' in allele:
                                    splitter = '_'
                                elif '-' in allele:
                                    splitter = '-'
                                else:
                                    splitter = ''
                                self.matchdict[sample.name].update({gene: allele.split(splitter)[-1]})
                                try:
                                    self.plusdict[sample.name][gene][allele.split(splitter)[-1]][float(percent_id)] \
                                        = depth
                                except KeyError:
                                    self.plusdict[sample.name][gene][allele.split(splitter)[-1]] = dict()
                                    self.plusdict[sample.name][gene][allele.split(splitter)[-1]][float(percent_id)] \
                                        = depth
                    if gene not in self.matchdict[sample.name]:
                        self.matchdict[sample.name].update({gene: 'N'})
        self.profiler()
        self.sequencetyper()
        self.new_alleles()
        self.mlstreporter()

    def new_alleles(self):
        """
        Parse samples to determine if new alleles need to be written to a file. Update the samples as required.
        """
        for sample in self.runmetadata.samples:
            # Only process the FASTA files if new alleles have been identified
            if sample[self.analysistype].new_alleles:
                # Read in the KMA outputs with SeqIO
                record_dict = SeqIO.to_dict(SeqIO.parse(sample[self.analysistype].kma_fasta_mem_mode, 'fasta'))
                # Retrieve all the new alleles for the sample
                for gene_allele in sample[self.analysistype].new_alleles:
                    # Extract the gene sequence from the dictionary
                    gene_sequence = str(record_dict[gene_allele].seq).upper()
                    # Split the string of gene_allele on underscores
                    split_name = gene_allele.split('_')
                    # Remove the allele from the list (this allows genes with underscores in the name to not be lost)
                    del split_name[-1]
                    # Recreate the gene name by joining the remaining items in the list with underscores
                    gene = '_'.join(split_name)
                    # Ensure that there are no gaps or 'N's in the allele sequence
                    if set(gene_sequence) == {"A", "T", "C", "G"}:
                        # Create a hash of the allele sequence to be used as the allele identifier. This way, if the
                        # allele is found multiple times, it will always produce the same name
                        # Uses logic from https://bitbucket.org/genomicepidemiology/cgmlstfinder/src/master/cgMLST.py
                        hash_str = hashlib.md5(gene_sequence.encode('utf-8')).hexdigest()
                        # Write the new allele to file
                        self.write_new_alleles(sample=sample,
                                               hash_str=hash_str,
                                               gene=gene,
                                               gene_seq=gene_sequence)
                        # Update the GenObject with the new allele information
                        for seqtype in self.resultprofile[sample.name]:
                            self.resultprofile[sample.name][seqtype][
                                sample[self.analysistype].matchestosequencetype][gene] = {hash_str: 100.00}

    def write_new_alleles(self, sample, hash_str, gene, gene_seq):
        """
        Write the new alleles to file in both the reports folder, as well as the target path
        :param sample: Metadata object of the current sample
        :param hash_str: hashlib-formatted string of the allele sequence to be used as part of the identifier
        :param gene: Name of the gene with the novel allele
        :param gene_seq: The sequence of the novel allele
        """
        # Set the name of the FASTA output files in both the report and target paths
        report_alleles = os.path.join(self.reportpath, 'new_cgmlst_alleles_{genus}.fasta'
                                      .format(genus=sample.general.closestrefseqgenus))
        db_alleles = os.path.join(sample[self.analysistype].targetpath, 'novel_alleles.fna')
        # Create the SeqRecord of the novel allele. Use Seq to prep the allele sequence. Set the ID to be the gene name
        # and the hash string
        record = SeqRecord(Seq(gene_seq),
                           id='{gene}_{hash_str}'.format(gene=gene,
                                                         hash_str=hash_str),
                           description=str())
        # Create and write the FASTA file in the report folder
        if not os.path.isfile(report_alleles):
            with open(report_alleles, 'w') as report_allele:
                SeqIO.write(record, report_allele, 'fasta')
        # If the report output file already exists, ensure that the novel allele isn't already present. If not, write
        # the novel allele to file
        else:
            record_dict = SeqIO.to_dict(SeqIO.parse(report_alleles, 'fasta'))
            if record.id not in record_dict:
                with open(report_alleles, 'a+') as report_allele:
                    SeqIO.write(record, report_allele, 'fasta')
        # Same logic as with the report alleles: create and populate the FASTA file file if it doesn't exist, otherwise,
        # ensure that the allele isn't already present in the file before writing
        if not os.path.isfile(db_alleles):
            with open(db_alleles, 'w') as db_allele:
                SeqIO.write(record, db_allele, 'fasta')
        else:
            record_dict = SeqIO.to_dict(SeqIO.parse(db_alleles, 'fasta'))
            if record.id not in record_dict:
                with open(db_alleles, 'a+') as db_allele:
                    SeqIO.write(record, db_allele, 'fasta')
