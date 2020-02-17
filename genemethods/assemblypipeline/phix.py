#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import run_subprocess
from argparse import ArgumentParser
import logging
import os

__author__ = 'adamkoziol'


class PhiX(object):

    def main(self):
        logging.info('Attempting to extract PhiX mapping data')
        interop_folder = os.path.join(self.path, 'InterOp')
        # Determine if the InterOp folder is present
        if os.path.isdir(interop_folder):
            # Try to extract the relevant information
            try:
                self.run_interop_summary()
                if self.pipeline:
                    self.interop_parse()
            # interop.py_interop_comm.file_not_found_exception: RunParameters.xml required for legacy run folders with
            # missing channel names
            except:
                if self.pipeline:
                    for sample in self.metadata:
                        sample.run.actual_yield = 'ND'
                        sample.run.projected_yield = 'ND'
                        sample.run.phix_aligned = 'ND'
                        sample.run.error_rate = 'ND'
                        sample.run.over_q30 = 'ND'
        else:
            if self.pipeline:
                # Create attributes reflecting the lack of the InterOp folder
                for sample in self.metadata:
                    sample.run.actual_yield = 'ND'
                    sample.run.projected_yield = 'ND'
                    sample.run.phix_aligned = 'ND'
                    sample.run.error_rate = 'ND'
                    sample.run.over_q30 = 'ND'

    def run_interop_summary(self):
        """
        Run the interop_summary script to create a .csv file summarising the run stats
        """
        interop_summary_cmd = 'interop_summary {path} > {report}'.format(path=self.path,
                                                                         report=self.interop_summary_report)
        if not os.path.isfile(self.interop_summary_report):
            run_subprocess(command=interop_summary_cmd)

    def interop_parse(self):
        """
        Parse the .csv summary file to extract the percent PhiX aligned and the error rate
        """
        if os.path.isfile(self.interop_summary_report):
            run_yield = str()
            projected_yield = str()
            phix_aligned = str()
            error_rate = str()
            over_q30 = str()
            cluster_density = str()
            with open(self.interop_summary_report, 'r') as summary:
                for line in summary:
                    # Only the line starting with 'Total' is required for most of the values
                    # Level        Yield    Projected Yield	  Aligned  	Error Rate  Intensity C1   	%>=Q30
                    # Total        17.23	17.23	          1.87	    2.75	    246	            65.19
                    if line.startswith('Total'):
                        total, run_yield, projected_yield, phix_aligned, error_rate, intensity_c1, over_q30 \
                            = line.replace(' ', '').rstrip().split(',')
                    # Extract the cluster density from the detailed outputs
                    # Lane ... Density, Cluster PF, ... Reads PF, %>=Q30, Yield, Cycles Error, Aligned, Error ...
                    if 'Density' in line:
                        for sub_line in summary:
                            # Extract the cleaned value for the cluster density from the line
                            cluster_density = sub_line.rstrip().split(',')[3].replace('    ', '')
                            # Only need to extract the value for cluster density from the first line
                            break
                        # Only extract the value once
                        break
                for sample in self.metadata:
                    sample.run.actual_yield = run_yield
                    sample.run.projected_yield = projected_yield
                    sample.run.phix_aligned = phix_aligned
                    sample.run.error_rate = error_rate
                    sample.run.over_q30 = over_q30
                    sample.run.cluster_density = cluster_density
        else:
            for sample in self.metadata:
                sample.run.actual_yield = 'ND'
                sample.run.projected_yield = 'ND'
                sample.run.phix_aligned = 'ND'
                sample.run.error_rate = 'ND'
                sample.run.over_q30 = 'ND'
                sample.cluster_density = 'ND'

    def __init__(self, inputobject, pipeline=True):
        self.path = inputobject.path
        self.reportpath = os.path.join(self.path, 'reports')
        self.interop_summary_report = os.path.join(self.reportpath, 'interop_summary_report.csv')
        self.pipeline = pipeline
        if self.pipeline:
            self.metadata = inputobject.runmetadata.samples
        else:
            self.metadata = list()


def cli():
    # Parser for arguments
    parser = ArgumentParser(description='Extract run information from InterOp folder')
    parser.add_argument('-p', '--path',
                        required=True,
                        help='Path to folder containing InterOp folder')
    # Get the arguments into an object
    arguments = parser.parse_args()
    interop = PhiX(inputobject=arguments,
                   pipeline=False)
    interop.main()


if __name__ == '__main__':
    cli()
