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
        quit()

    def new_alleles(self):
        for sample in self.runmetadata.samples:
            if sample[self.analysistype].new_alleles:
                record_dict = SeqIO.to_dict(SeqIO.parse(sample[self.analysistype].kma_fasta_mem_mode, 'fasta'))
                for gene_allele in sample[self.analysistype].new_alleles:
                    gene_sequence = str(record_dict[gene_allele].seq).upper()
                    split_name = gene_allele.split('_')
                    del split_name[-1]
                    gene = '_'.join(split_name)
                    if set(gene_sequence) == {"A", "T", "C", "G"}:

                        hash_str = hashlib.md5(gene_sequence.encode('utf-8')).hexdigest()
                        print(sample.name, gene_allele, hash_str)
                        self.write_new_alleles(sample=sample,
                                               hash_str=hash_str,
                                               gene_seq=gene_sequence)
                        for seqtype in self.resultprofile[sample.name]:
                            # for allele, percent_id in self.resultprofile[sample.name][seqtype][sample[self.analysistype].matchestosequencetype][gene].items():
                            # print(seqtype, gene, allele, percent_id, sample[self.analysistype].targetpath)
                            # print(self.resultprofile[sample.name][seqtype][sample[self.analysistype].matchestosequencetype][gene])

                            # self.resultprofile[sample.name][seqtype] += 1
                            # sample[self.analysistype].matchestosequencetype += 1
                            self.resultprofile[sample.name][seqtype][
                                sample[self.analysistype].matchestosequencetype][gene] = {hash_str: 100.00}

                    else:
                        print('Illegal character in {sn} {split} {gene} {bad}'.format(sn=sample.name,
                                                                                      split=split_name,
                                                                                      gene=gene_allele,
                                                                                      bad=set(gene_sequence)))

    def write_new_alleles(self, sample, hash_str, gene_seq):
        print(hash_str, gene_seq, sample.general.closestrefseqgenus, sample[self.analysistype].targetpath)
        print(self.reportpath)
        report_alleles = os.path.join(self.reportpath, 'new_cgmlst_alleles_{genus}.fasta'
                                      .format(genus=sample.general.closestrefseqgenus))
        db_alleles = os.path.join(sample[self.analysistype].targetpath, 'novel_alleles.fna')
        record = SeqRecord(Seq(gene_seq),
                           id=hash_str,
                           description=str())
        if not os.path.isfile(report_alleles):
            with open(report_alleles, 'w') as report_allele:
                SeqIO.write(record, report_allele, 'fasta')
        else:
            record_dict = SeqIO.to_dict(SeqIO.parse(report_alleles, 'fasta'))
            if record.id not in record_dict:
                with open(report_alleles, 'a+') as report_allele:
                    SeqIO.write(record, report_allele, 'fasta')

        if not os.path.isfile(db_alleles):
            with open(db_alleles, 'w') as db_allele:
                SeqIO.write(record, db_allele, 'fasta')
        else:
            record_dict = SeqIO.to_dict(SeqIO.parse(db_alleles, 'fasta'))
            if record.id not in record_dict:
                with open(db_alleles, 'a+') as db_allele:
                    SeqIO.write(record, db_allele, 'fasta')


