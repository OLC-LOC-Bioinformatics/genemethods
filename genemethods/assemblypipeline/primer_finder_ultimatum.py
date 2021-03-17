#!/usr/bin/env python
from olctools.accessoryFunctions.accessoryFunctions import GenObject, make_path, MetadataObject, SetupLogging
from genemethods.assemblypipeline.primer_finder_ipcress import best_hit, best_hits_singletons, \
    blast_primer_mismatch_details, create_fasta_format, fasta_primers, make_blastdb, parse_blast, run_blast
from genemethods.assemblypipeline.legacy_vtyper import epcr_primers, Filer
from genemethods.assemblypipeline.primer_finder_bbduk import ampliconfile
from argparse import ArgumentParser
import multiprocessing
import logging
import shutil
import os


class Ultimatum(object):

    def main(self):
        if self.primer_format == 'fasta':
            self.forward_dict, self.reverse_dict = fasta_primers(primerfile=self.primerfile,
                                                                 forward_dict=self.forward_dict,
                                                                 reverse_dict=self.reverse_dict)
            self.fastaprimers = self.primerfile
        else:
            self.forward_dict, self.reverse_dict = epcr_primers(primerfile=self.primerfile,
                                                                forward_dict=self.forward_dict,
                                                                reverse_dict=self.reverse_dict,
                                                                formattedprimers=self.formattedprimers)
            create_fasta_format(forward_dict=self.forward_dict,
                                reverse_dict=self.reverse_dict,
                                fastaprimerfile=self.fastaprimers)
        # Get the metadata object ready to use the methods from primer_finder_ipcress
        self.metadata_prep()
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
                                    mismatches=self.mismatches + 1,
                                    min_hits=1)
        self.metadata = blast_primer_mismatch_details(metadata=self.metadata,
                                                      analysistype=self.analysistype,
                                                      iupac=self.iupac)

        self.report_prep()
        self.metadata = best_hit(metadata=self.metadata,
                                 analysistype=self.analysistype)
        self.metadata = best_hits_singletons(metadata=self.metadata,
                                             analysistype=self.analysistype)
        self.amplicon_singleton_report()
        if self.export_amplicons:
            self.amplicon_prep()

    def metadata_prep(self):
        """
        Prepare the metadata objects to be compatible with the methods from primer_finder_ipcress
        :return:
        """
        for sample in self.metadata:
            setattr(sample, self.analysistype, GenObject())
            # Initialise the results GenObject
            sample[self.analysistype].results = GenObject()
            # Make the output path
            make_path(sample.general.outputdirectory)

    def report_prep(self):
        """
        Prepare the metadata objects for the report method
        """
        for sample in self.metadata:
            for experiment in sample[self.analysistype].results.datastore:
                for contig in sample[self.analysistype].results[experiment].datastore:
                    # Ensure that all the necessary GenObjects exist, and create them as required
                    if not GenObject.isattr(sample[self.analysistype].results[experiment][contig],
                                            'forward_mismatch_details'):
                        sample[self.analysistype].results[experiment][contig].forward_mismatch_details = str()
                    if not GenObject.isattr(sample[self.analysistype].results[experiment][contig], 'forward_query'):
                        sample[self.analysistype].results[experiment][contig].forward_query = str()
                    if not GenObject.isattr(sample[self.analysistype].results[experiment][contig], 'forward_mismatch'):
                        sample[self.analysistype].results[experiment][contig].forward_mismatch = str()
                    if not GenObject.isattr(sample[self.analysistype].results[experiment][contig], 'forward_ref'):
                        sample[self.analysistype].results[experiment][contig].forward_ref = str()
                    if not GenObject.isattr(sample[self.analysistype].results[experiment][contig],
                                            'reverse_mismatch_details'):
                        sample[self.analysistype].results[experiment][contig].reverse_mismatch_details = str()
                    if not GenObject.isattr(sample[self.analysistype].results[experiment][contig], 'reverse_query'):
                        sample[self.analysistype].results[experiment][contig].reverse_query = str()
                    if not GenObject.isattr(sample[self.analysistype].results[experiment][contig], 'reverse_mismatch'):
                        sample[self.analysistype].results[experiment][contig].reverse_mismatch = str()
                    if not GenObject.isattr(sample[self.analysistype].results[experiment][contig], 'reverse_ref'):
                        sample[self.analysistype].results[experiment][contig].reverse_ref = str()
                    # Determine the amplicon location - ensure that both the forward and reverse primer have results
                    if GenObject.isattr(sample[self.analysistype].results[experiment][contig],
                                        'forward_range') and GenObject.\
                            isattr(sample[self.analysistype].results[experiment][contig], 'reverse_range'):
                        # Find the minimum (start) and maximum (end) positions in the range variables
                        minimum = min(sample[self.analysistype].results[experiment][contig].forward_range[0],
                                      sample[self.analysistype].results[experiment][contig].forward_range[-1],
                                      sample[self.analysistype].results[experiment][contig].reverse_range[0],
                                      sample[self.analysistype].results[experiment][contig].reverse_range[-1])
                        maximum = max(sample[self.analysistype].results[experiment][contig].forward_range[0],
                                      sample[self.analysistype].results[experiment][contig].forward_range[-1],
                                      sample[self.analysistype].results[experiment][contig].reverse_range[0],
                                      sample[self.analysistype].results[experiment][contig].reverse_range[-1])
                        # Set the amplicon as min->max
                        sample[self.analysistype].results[experiment][contig].location = f'{minimum}-{maximum}'
                        sample[self.analysistype].results[experiment][contig].amplicon_range = range(minimum, maximum)
                        # Set the amplicon size as max-min
                        sample[self.analysistype].results[experiment][contig].amplicon_length = maximum - minimum
                        # The method for calculating orientation is only for singleton primer direction; use the forward
                        # and reverse ranges
                        sample[self.analysistype].results[experiment][contig].direction = 'forward' if sample[
                            self.analysistype].results[experiment][contig].forward_range[0] < sample[
                            self.analysistype].results[experiment][contig].reverse_range[0] else 'reverse'
                        # Split the contig attribute;the BLAST parser assumes that the forward and reverse primers will
                        # hit on separate contigs, but since they are on the same contig, we only need the reverse info
                        sample[self.analysistype].results[experiment][contig].contig = \
                            sample[self.analysistype].results[experiment][contig].contig.split('reverse_')[-1]
                        # Calculate the total number of mismatches
                        sample[self.analysistype].results[experiment][contig].total_mismatch = \
                            int(sample[self.analysistype].results[experiment][contig].forward_mismatch) + \
                            int(sample[self.analysistype].results[experiment][contig].reverse_mismatch)

                    else:
                        # Since no amplicons exist, create empty attributes
                        if not GenObject.isattr(sample[self.analysistype].results[experiment][contig], 'location'):
                            sample[self.analysistype].results[experiment][contig].location = str()
                        if not GenObject.isattr(sample[self.analysistype].results[experiment][contig],
                                                'amplicon_length'):
                            sample[self.analysistype].results[experiment][contig].amplicon_length = str()
                        # Set the .direction attribute
                        if GenObject.isattr(sample[self.analysistype].results[experiment][contig], 'forward_direction'):
                            sample[self.analysistype].results[experiment][contig].direction = \
                                sample[self.analysistype].results[experiment][contig].forward_direction
                        else:
                            sample[self.analysistype].results[experiment][contig].direction = \
                                sample[self.analysistype].results[experiment][contig].reverse_direction

    def amplicon_singleton_report(self):
        """
        Create a custom report that summarises the outputs for both amplicon-creating and singleton primer sets
        :return:
        """
        logging.info('Creating summary report')
        with open(self.report, 'w') as report:
            data = 'Sample,Gene,Contig,GenomeLocation,AmpliconSize,Orientation,ForwardMismatches,' \
                   'ForwardMismatchDetails,ForwardLength,ReverseMismatches,ReverseMismatchDetails,ReverseLength,' \
                   'ForwardPrimer,ForwardQuery,ReversePrimer,ReverseQuery\n'
            for sample in self.metadata:
                # Initialise a variable to store whether results were found for a sample
                results = False
                for contig, amplicon_dict in sample[self.analysistype].best_hits.items():
                    # Use the best_hits dictionary to extract only primer sets considered to be the best hit in a range
                    for amplicon_range, hit_dict in amplicon_dict.items():
                        experiment = hit_dict['experiment']
                        # Set the results to true
                        results = True
                        data += \
                            '{sn},{gene},{contig},{loc},{a_sz},{ori},{f_mm},{f_md},{f_l},{r_mm},{r_md},{r_l},{fp},' \
                            '{fq},{rp},{rq}\n'\
                            .format(sn=sample.name,
                                    gene=sample[self.analysistype].results[experiment][contig].primer_set,
                                    contig=sample[self.analysistype].results[experiment][contig].contig,
                                    loc=sample[self.analysistype].results[experiment][contig].location,
                                    a_sz=sample[self.analysistype].results[experiment][contig].amplicon_length,
                                    ori=sample[self.analysistype].results[experiment][contig].direction,
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

    def amplicon_prep(self):
        """
        Get the metadata objects ready to be used in the amplicon creating method
        """
        logging.info('Creating amplicons')
        for sample in self.metadata:
            sample[self.analysistype].assemblyfile = sample.general.bestassemblyfile
            sample[self.analysistype].ampliconfile = os.path.join(self.reportpath, 'alleles.fasta')
            for experiment in sample[self.analysistype].results.datastore:
                for contig in sample[self.analysistype].results[experiment].datastore:
                    try:
                        amplicon_range = [
                            sample[self.analysistype].results[experiment][contig].amplicon_range[0],
                            sample[self.analysistype].results[experiment][contig].amplicon_range[-1],
                        ]
                        ampliconfile(sample=sample,
                                     analysistype=self.analysistype,
                                     contig=contig,
                                     amplicon_range=amplicon_range,
                                     forward_primer=[experiment + '-F'],
                                     reverse_primer=[experiment + '-R'])
                    except AttributeError:
                        pass

    def __init__(self, metadataobject, sequencepath, reportpath, primerfile, primer_format, mismatches=2,
                 export_amplicons=False):
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
        self.primer_format = primer_format
        self.mismatches = mismatches
        self.export_amplicons = export_amplicons
        self.formattedprimers = os.path.join(os.path.dirname(self.primerfile), 'epcr_formatted_primers',
                                             'formatted_primers.txt')
        self.fastaprimers = os.path.join(os.path.dirname(self.primerfile), 'formattedprimers.fa')
        make_path(os.path.dirname(self.formattedprimers))
        try:
            shutil.rmtree(self.reportpath)
        except FileNotFoundError:
            pass
        make_path(self.reportpath)
        self.forward_dict = dict()
        self.reverse_dict = dict()
        self.threads = multiprocessing.cpu_count() - 1
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
    parser = ArgumentParser(description='Use BLAST to search for primer hits in a sequence')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path to folder containing FASTA files')
    parser.add_argument('-m', '--mismatches',
                        default=2,
                        help='Number of mismatches to allow for ipcress searches. Default is 2')
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
    parser.add_argument('-e', '--export_amplicons',
                        action='store_true',
                        help='Export the sequence of the calculated amplicons. Default is False')
    # Get the arguments into an object
    arguments = parser.parse_args()
    SetupLogging(debug=True)
    arguments.reportpath = os.path.join(arguments.sequencepath, 'reports')
    arguments.runmetadata = MetadataObject()
    # Create metadata objects for the samples
    arguments.runmetadata.samples = Filer.filer(arguments)
    finder = Ultimatum(metadataobject=arguments.runmetadata,
                       sequencepath=arguments.sequencepath,
                       reportpath=arguments.reportpath,
                       primerfile=arguments.primerfile,
                       primer_format=arguments.primer_format,
                       mismatches=arguments.mismatches,
                       export_amplicons=arguments.export_amplicons)
    finder.main()
    logging.info('Primer Finder Ultimatum Complete!')


if __name__ == '__main__':
    cli()
