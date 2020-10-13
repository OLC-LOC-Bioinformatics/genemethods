#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import GenObject, make_path, run_subprocess, write_to_logfile
import logging
import os
__author__ = 'adamkoziol'


class SeqSero(object):

    def main(self):
        logging.info('Performing SeqSero2 analyses on samples')
        self.seqsero()
        self.parse_seqsero()
        self.seqsero_report()

    def seqsero(self):
        """
        Run seqsero2 on Salmonella samples
        """
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Create the GenObject, and necessary attributes
                setattr(sample, self.analysistype, GenObject())
                sample[self.analysistype].outputdir = os.path.join(sample.general.outputdirectory, self.analysistype)
                sample[self.analysistype].log = os.path.join(sample[self.analysistype].outputdir, 'log')
                sample[self.analysistype].err = os.path.join(sample[self.analysistype].outputdir, 'err')
                sample[self.analysistype].report = os.path.join(sample[self.analysistype].outputdir,
                                                                'SeqSero_result.txt')
                make_path(sample[self.analysistype].outputdir)
                # Extract the calculated genus
                try:
                    genus = sample.general.closestrefseqgenus
                except AttributeError:
                    try:
                        genus = sample.general.referencegenus
                    except AttributeError:
                        genus = 'ND'
                # Only run seqsero on Salmonella samples
                if genus == 'Salmonella':
                    # Make the system call - use the default allele micro-assembly workflow
                    sample[self.analysistype].seqsero_cmd = 'SeqSero2_package.py -p {cpus} -t {numreads} -i {fastq} ' \
                                                            '-d {outdir}'\
                        .format(cpus=self.cpus,
                                numreads=len(sample.general.trimmedcorrectedfastqfiles),
                                fastq=' '.join(sample.general.trimmedcorrectedfastqfiles),
                                outdir=sample[self.analysistype].outputdir)
                    # Run the system call if the seqsero report doesn't exist
                    if not os.path.isfile(os.path.join(sample[self.analysistype].report)):
                        out, err = run_subprocess(command=sample[self.analysistype].seqsero_cmd)
                        write_to_logfile(out='{cmd}\n{out}'.format(cmd=sample[self.analysistype].seqsero_cmd,
                                                                   out=out),
                                         err=err,
                                         logfile=self.logfile,
                                         samplelog=sample.general.logout,
                                         sampleerr=sample.general.logerr,
                                         analysislog=sample[self.analysistype].log,
                                         analysiserr=sample[self.analysistype].err)

    def parse_seqsero(self):
        """
        Parse the output file, and create attributes as necessary
        """
        for sample in self.metadata:
            if sample.general.bestassemblyfile != 'NA':
                # Non-Salmonella samples will not have the report
                try:
                    with open(sample[self.analysistype].report, 'r') as report:
                        for line in report:
                            # Lines corresponding to input files, and notes, yield list sizes of 3, and 1,
                            # respectively. These values can be skipped
                            try:
                                key, value = line.rstrip().split('\t')
                                # Clean the string by removing illegal/unnecessary characters/spaces
                                key = key.replace(':', '').replace(' ', '_').replace('(', '_').replace(')', '').lower()
                                # Add the cleaned key:value pair as attribute: value to the GenObject
                                setattr(sample[self.analysistype], key, value)
                            except ValueError:
                                # The 'Note' field is not tab-separated. If empty, it is "Note:", while a report with
                                # a note will look like this: "Note: This predicted.... "
                                if 'Note' in line:
                                    # Split the note on the Note:[SPACE]
                                    try:
                                        sample[self.analysistype].note = line.rstrip().split('Note: ')[1]
                                    except IndexError:
                                        sample[self.analysistype].note = str()
                except (AttributeError, FileNotFoundError):
                    pass

    def seqsero_report(self):
        """
        Create a summary report of the SeqSero outputs
        """
        with open(os.path.join(self.reportpath, 'seqsero.tsv'), 'w') as report:
            data = 'Strain\tO_antigen_prediction\tH1_antigen_prediction\tH2_antigen_prediction\t' \
                   'Predicted_subspecies\tPredicted_antigenic_profile\tPredicted_serotype\tNote\n'
            for sample in self.metadata:
                data += '{sn}\t'.format(sn=sample.name)
                # Samples that are either not Salmonella, or on which seqsero was unable to complete successfully will
                # not have all the required attributes
                try:
                    data += '{O_ant}\t{H1}\t{H2}\t{sub}\t{profile}\t{sero}\t{note}\n'\
                        .format(O_ant=sample[self.analysistype].o_antigen_prediction,
                                H1=sample[self.analysistype].h1_antigen_prediction_flic,
                                H2=sample[self.analysistype].h2_antigen_prediction_fljb,
                                sub=sample[self.analysistype].predicted_subspecies,
                                profile=sample[self.analysistype].predicted_antigenic_profile,
                                sero=sample[self.analysistype].predicted_serotype,
                                note=sample[self.analysistype].note)
                except AttributeError:
                    data += '\n'
            # Write the report
            report.write(data)

    def __init__(self, inputobject):
        self.metadata = inputobject.runmetadata.samples
        self.cpus = inputobject.cpus
        self.analysistype = 'seqsero'
        self.logfile = inputobject.logfile
        self.reportpath = inputobject.reportpath
        try:
            self.threads = int(self.cpus / len(self.metadata)) if self.cpus / len(self.metadata) > 1 else 1
        except TypeError:
            self.threads = self.cpus
