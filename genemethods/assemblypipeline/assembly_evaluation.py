#!/usr/bin/env python3

"""
Assembly evaluation using various tools including Bowtie2, Samtools, Quast,
and Qualimap.
"""

# Standard imports
from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from glob import glob
from io import StringIO
import logging
import os
import multiprocessing
from queue import Queue
import re
import shutil

# Third-party imports
from Bio import SeqIO
from Bio.Application import ApplicationError
from Bio.Sequencing.Applications import SamtoolsIndexCommandline
from olctools.accessoryFunctions.accessoryFunctions import (
    log_str,
    run_subprocess,
    write_to_log_file
)
from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'


class AssemblyEvaluation:
    """
    Class to evaluate genome assemblies using various tools.
    """

    def __init__(self, log_file, metadata, sequence_path, threads) -> None:
        """
        Initialize the AssemblyEvaluation class with the provided input object.

        Args:
            inputobject (Any): Input object containing metadata and
            configuration.
        """
        self.metadata = metadata
        self.threads = threads
        self.log_file = log_file
        self.sequence_path = sequence_path

        # Initialize queues
        self.qualimap_queue = Queue(maxsize=self.threads)
        self.pilon_queue = Queue(maxsize=self.threads)
        self.index_queue = Queue(maxsize=self.threads)
        self.filter_queue = Queue(maxsize=self.threads)

    def main(self) -> None:
        """
        Run the methods in the correct order.
        """
        self.bowtie_build()
        self.bowtie_run()
        self.indexing()
        self.quast()
        self.parse_quast_report()
        self.qualimapper()
        self.parse_qualimap_report()
        self.clean_quast()
        self.extract_insert_size()
        self.calculate_weighted_insert_size()
        self.pilon()
        self.filter()
        self.clear()

        return self.metadata

    def bowtie_build(self) -> None:
        """
        Use bowtie2-build to index each target file.

        This method iterates over the metadata samples, constructs the
        bowtie2-build command, and executes it if the index files do not
        already exist. The output and errors are logged appropriately.

        Preconditions:
        - `self.metadata` should contain the sample data with necessary
          attributes.
        - `self.log_file` should be the path to the log file.

        Raises:
        - FileNotFoundError: If any required file paths are not found.
        - RuntimeError: If the subprocess command fails.
        """
        logging.info('Preparing targets for reference mapping')

        for sample in self.metadata:
            logging.debug("Processing sample: %s", sample.name)

            # Create and populate the quast GenObject
            sample.quast = CustomBox()
            sample.quast.base_name = os.path.splitext(
                sample.general.assembly_file
            )[0]
            sample.quast.build_command = (
                f'bowtie2-build {sample.general.assembly_file} '
                f'{sample.quast.base_name}'
            )
            logging.debug(
                "Constructed build command: %s", sample.quast.build_command
            )

            # Ensure the outputs don't already exist before running bowtie2
            if not os.path.isfile(f'{sample.quast.base_name}.1.bt2'):
                logging.info(
                    "Index files not found for sample %s. Running "
                    "bowtie2-build", sample.name
                )
                try:
                    out, err = run_subprocess(sample.quast.build_command)
                    logging.debug("Command output: %s", out)
                    logging.debug("Command error: %s", err)

                    # Write the appropriate information to the log_file
                    write_to_log_file(
                        out=f'{sample.quast.build_command}\n{out}',
                        err=err,
                        log_file=self.log_file,
                        sample_log=sample.general.log_out,
                        sample_err=sample.general.log_err
                    )
                except RuntimeError as e:
                    logging.error(
                        "Failed to run build command for sample %s: %s",
                        sample.name, e
                    )
            else:
                logging.info(
                    "Index files already exist for sample %s, skipping build",
                    sample.name
                )

    def bowtie_run(self) -> None:
        """
        Map the FASTQ reads against the appropriate target file using Bowtie2
        and Samtools.

        This method iterates over the metadata samples, constructs the Bowtie2
        mapping command, and executes it if the sorted BAM file does not
        already exist. The output and errors are logged appropriately.

        Preconditions:
        - `self.metadata` should contain the sample data with necessary
          attributes.
        - `self.threads` should be set to the number of CPUs
          and threads to use.
        - `self.log_file` should be the path to the log file.

        Raises:
        - FileNotFoundError: If any required file paths are not found.
        - RuntimeError: If the subprocess command fails.
        """
        logging.info('Performing reference mapping for quality evaluation')

        for sample in self.metadata:
            logging.debug("Processing sample: %s", sample.name)

            # Ensure the output directory exists
            sample.quast.output_dir = os.path.join(
                sample.general.output_directory, 'quast'
            )
            if sample.general.best_assembly_file == "NA":
                logging.warning(
                    "Skipping sample %s as best assembly file is 'NA'",
                    sample.name
                )
                continue

            os.makedirs(sample.quast.output_dir, exist_ok=True)
            sample.quast.sorted_bam = os.path.join(
                sample.quast.output_dir, f'{sample.name}_sorted.bam'
            )

            # Construct the Bowtie2 mapping command
            sample.quast.map_command = f'bowtie2 -x {sample.quast.base_name}'

            # Add FASTQ files to the command
            if len(sample.general.trimmed_corrected_fastq_files) == 1:
                sample.quast.map_command += (
                    f' -U {sample.general.trimmed_corrected_fastq_files[0]}'
                )
            else:
                sample.quast.map_command += (
                    f' -1 {sample.general.trimmed_corrected_fastq_files[0]} '
                    f'-2 {sample.general.trimmed_corrected_fastq_files[1]}'
                )

            # Add Samtools commands to convert and sort the BAM file
            sample.quast.map_command += (
                f' -p {self.threads} -X 1000 | '
                f'samtools view -@ {self.threads} -h -F 4 -bT '
                f'{sample.general.assembly_file} - | '
                f'samtools sort - -@ {self.threads} '
                f'-o {sample.quast.sorted_bam}'
            )
            logging.debug(
                "Constructed map command: %s", sample.quast.map_command
            )

            # Run the mapping command if the sorted BAM file doesn't exist
            if not os.path.isfile(sample.quast.sorted_bam):
                logging.info(
                    "Sorted BAM file not found for sample %s, running command",
                    sample.name
                )
                try:
                    out, err = run_subprocess(sample.quast.map_command)
                    logging.debug("Command output: %s", out)
                    logging.debug("Command error: %s", err)

                    # Write the appropriate information to the log_file
                    write_to_log_file(
                        out=f'{sample.quast.map_command}\n{out}',
                        err=err,
                        log_file=self.log_file,
                        sample_log=sample.general.log_out,
                        sample_err=sample.general.log_err,
                    )
                except RuntimeError as e:
                    logging.error(
                        "Failed to run mapping command for sample %s: %s",
                        sample.name, e
                    )
            else:
                logging.info(
                    "Sorted BAM file already exists for sample %s, skipping",
                    sample.name
                )

    def indexing(self) -> None:
        """
        Use samtools index to index the sorted BAM files.

        This method uses a ThreadPoolExecutor to manage threads and submit the
        indexing tasks. It enqueues tasks for each sample that needs indexing
        and waits for all tasks to be completed.

        Preconditions:
        - `self.metadata` should contain the sample data with necessary
          attributes.
        - `self.threads` should be set to the number of threads to use.
        - `self.index_queue` should be a Queue instance.

        Raises:
        - FileNotFoundError: If any required file paths are not found.
        - RuntimeError: If the subprocess command fails.
        """
        logging.info('Indexing sorted BAM files')

        # Use ThreadPoolExecutor to manage threads
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            # Submit the index method to the executor
            for _ in range(self.threads):
                executor.submit(self.index)

            # Enqueue tasks
            for sample in self.metadata:
                if sample.general.best_assembly_file != 'NA':
                    bam_index = SamtoolsIndexCommandline(
                        input=sample.quast.sorted_bam
                    )
                    sample.quast.sorted_bai = sample.quast.sorted_bam + '.bai'
                    sample.quast.bam_index = str(bam_index)
                    logging.debug(
                        "Enqueuing indexing task for sample: %s", sample.name
                    )
                    self.index_queue.put((sample, bam_index))

            # Wait for all tasks to be completed
            self.index_queue.join()

    def index(self) -> None:
        """
        Worker method to process the indexing tasks.

        This method runs in a loop, processing tasks from the index_queue.
        It uses samtools to index the BAM files and logs any errors.
        """
        while True:
            try:
                sample, bam_index = self.index_queue.get()
                logging.debug(
                    "Processing indexing task for sample: %s",
                    sample.name)

                # Only make the call if the .bai file doesn't already exist
                if not os.path.isfile(sample.quast.sorted_bai):
                    logging.info(
                        "Index file not found for sample %s, running "
                        "samtools index",
                        sample.name
                    )
                    # Use StringIO streams to handle output
                    stdout, stderr = map(
                        StringIO, bam_index(cwd=sample.quast.output_dir)
                    )
                    if stderr:
                        logging.error(
                            "Error indexing BAM file for sample %s: %s",
                            sample.name, stderr.getvalue()
                        )

                        # Set the name of the samtools log file
                        samtools_log = os.path.join(
                            sample.quast.output_dir,
                            'indexing_samtools_bam_index.log'
                        )

                        # Write the standard error to log
                        with open(samtools_log, 'a+', encoding='utf-8') as log:
                            log.writelines(
                                log_str(
                                    bam_index,
                                    stderr.getvalue(),
                                    stdout.getvalue()
                                )
                            )
                    stderr.close()
            except ApplicationError as exc:
                logging.error(
                    "Application error while indexing BAM file for "
                    "sample %s: %s",
                    sample.name, exc
                )
            finally:
                self.index_queue.task_done()

    def quast(self) -> None:
        """
        Run quast on the samples.

        This method iterates over the metadata samples, constructs the quast
        command, and executes it if the final quast report does not already
        exist. The output and errors are logged appropriately.

        Preconditions:
        - `self.metadata` should contain the sample data with necessary
          attributes.
        - `self.threads` should be set to the number of threads to use.
        - `self.log_file` should be the path to the log file.

        Raises:
        - FileNotFoundError: If any required file paths are not found.
        - RuntimeError: If the subprocess command fails.
        """
        logging.info('Running Quast on assemblies')

        for sample in self.metadata:
            logging.debug("Processing sample: %s", sample.name)

            # Populate necessary attributes
            sample.quast.report = os.path.join(
                sample.quast.output_dir, 'report.tsv'
            )
            if sample.general.best_assembly_file == "NA":
                logging.warning(
                    "Skipping sample %s as best assembly file is 'NA'",
                    sample.name
                )
                continue

            # Allow for non-paired samples
            if len(sample.general.trimmed_corrected_fastq_files) == 2:
                sample.quast.cmd = (
                    f'quast.py --pe1 '
                    f'{sample.general.trimmed_corrected_fastq_files[0]} '
                    f'--pe2 '
                    f'{sample.general.trimmed_corrected_fastq_files[1]}'
                )
            else:
                sample.quast.cmd = (
                    f'quast.py --single '
                    f'{sample.general.trimmed_corrected_fastq_files[0]}'
                )

            # Both paired and unpaired samples share the rest of the system
            # call
            sample.quast.cmd += (
                f' --ref-bam {sample.quast.sorted_bam} -t {self.threads} '
                f'--k-mer-stats --circos --rna-finding '
                f'--conserved-genes-finding -o {sample.quast.output_dir} '
                f'--debug {sample.general.assembly_file} '
                f'--threads {self.threads}'
            )
            logging.debug("Constructed quast command: %s", sample.quast.cmd)

            # Run the quast system call if the final quast report doesn't exist
            if not os.path.isfile(sample.quast.report):
                logging.info(
                    "Quast report not found for sample %s, running command",
                    sample.name
                )
                try:
                    out, err = run_subprocess(sample.quast.cmd)
                    logging.debug("Command output: %s", out)
                    logging.debug("Command error: %s", err)

                    # Write the appropriate information to the log_file
                    write_to_log_file(
                        out=f'{sample.quast.cmd}\n{out}',
                        err=err,
                        log_file=self.log_file,
                        sample_log=sample.general.log_out,
                        sample_err=sample.general.log_err,
                    )
                except RuntimeError as exc:
                    logging.error(
                        "Failed to run quast command for sample %s: %s",
                        sample.name, exc
                    )
            else:
                logging.info(
                    "Quast report already exists for sample %s, skipping",
                    sample.name
                )

    def parse_quast_report(self) -> None:
        """
        Parse the quast report, and populate the metadata object with the
        extracted key: value pairs.

        This method iterates over the metadata samples, reads the quast report,
        and populates the metadata object with the extracted key-value pairs.

        Preconditions:
        - `self.metadata` should contain the sample data with necessary
          attributes.

        Raises:
        - FileNotFoundError: If the quast report file is not found.
        """
        logging.info('Parsing Quast reports')

        for sample in self.metadata:
            logging.debug("Processing sample: %s", sample.name)

            if not os.path.isfile(sample.quast.report):
                logging.warning(
                    "Quast report not found for sample %s, skipping",
                    sample.name
                )
                continue

            logging.info(
                "Quast report found for sample %s, parsing report",
                sample.name
            )
            self._parse_report(sample)

    def _parse_report(self, sample) -> None:
        """
        Parse the quast report file and extract key-value pairs.

        Args:
            sample: The sample object to parse the report for.
        """
        try:
            with open(sample.quast.report, 'r', encoding='utf-8') as report:
                for line in report:
                    key, value = self.analyze(line)
                    logging.debug(
                        "Extracted key-value pair: %s: %s", key, value
                    )
                    setattr(sample.quast, key, value)
        except FileNotFoundError as exc:
            logging.error(
                "Quast report file not found for sample %s: %s",
                sample.name, exc
            )

    def clean_quast(self) -> None:
        """
        Remove all the unnecessary temporary files created by quast.

        This method iterates over the metadata samples, finds all files in the
        quast directory, and removes large files that do not have .err or .html
        extensions.

        Preconditions:
        - `self.metadata` should contain the sample data with necessary
          attributes.

        Raises:
        - FileNotFoundError: If any required file paths are not found.
        """
        logging.info('Cleaning Quast temporary files')

        for sample in self.metadata:
            logging.debug("Processing sample: %s", sample.name)

            if not os.path.isdir(sample.quast.output_dir):
                logging.warning(
                    "Quast output directory not found for sample %s, skipping",
                    sample.name
                )
                continue

            self._clean_sample_quast_files(sample)

    def _clean_sample_quast_files(self, sample) -> None:
        """
        Clean unnecessary temporary files for a sample.

        Args:
            sample: The sample object to clean files for.
        """
        for path, _, files in os.walk(sample.quast.output_dir):
            for quast_file in files:
                file_path = os.path.join(path, quast_file)
                self._remove_large_files(file_path, quast_file, sample)

    def _remove_large_files(self, file_path, quast_file, sample) -> None:
        """
        Remove large files that do not have .err or .html extensions.

        Args:
            file_path: The absolute path of the file.
            quast_file: The name of the file.
            sample: The sample object to log information for.
        """
        try:
            if (
                os.path.getsize(file_path) > 100000
                and '_sorted.bam' not in quast_file
                and '.err' not in quast_file
                and '.html' not in quast_file
            ):
                logging.info(
                    "Removing file %s for sample %s",
                    file_path, sample.name
                )
                os.remove(file_path)
        except FileNotFoundError as exc:
            logging.error(
                "File not found %s for sample %s: %s",
                file_path, sample.name, exc
            )
        except OSError as exc:
            logging.error(
                "Error removing file %s for sample %s: %s",
                file_path, sample.name, exc
            )

    def extract_insert_size(self) -> None:
        """
        Parse the bwa index log information to extract insert size estimations.

        This method iterates over the metadata samples, reads the bwa index
        log, and extracts insert size estimations. The extracted data is
        stored in the metadata object.

        Preconditions:
        - `self.metadata` should contain the sample data with necessary
          attributes.

        Raises:
        - FileNotFoundError: If the reads stats file is not found.
        """
        logging.info('Calculating insert size')

        for sample in self.metadata:
            logging.debug("Processing sample: %s", sample.name)

            sample.quast.reads_stats_file = os.path.join(
                sample.quast.output_dir, 'reads_stats', 'reads_stats.err'
            )

            if os.path.isfile(sample.quast.reads_stats_file):
                logging.info(
                    "Reads stats file found for sample %s, parsing file",
                    sample.name
                )

                self._initialize_sample_attributes(sample)
                self._parse_reads_stats_file(sample)
            else:
                logging.warning(
                    "Reads stats file not found for sample %s, skipping",
                    sample.name
                )
                sample.quast.insert_mean = 'ND'
                sample.quast.insert_std = 'ND'

    def _initialize_sample_attributes(self, sample) -> None:
        """
        Initialize attributes for the insert size estimation.

        Args:
            sample: The sample object to initialize attributes for.
        """
        sample.quast.total_reads = 0
        sample.quast.insert_mean = []
        sample.quast.insert_std = []
        sample.quast.read_blocks = []

    def _parse_reads_stats_file(self, sample) -> None:
        """
        Parse the reads stats file to extract insert size estimations.

        Args:
            sample: The sample object to populate with extracted data.
        """
        current_reads = 0

        # Create a variable to store the long attribute name
        stats_file = sample.quast.reads_stats_file

        with open(stats_file, 'r', encoding='utf-8') as read_stats:
            for line in read_stats:

                # Create a variable to store the long FR section line for
                # easier parsing
                size = 'analyzing insert size distribution for orientation FR'

                # BWA estimates the insert size distribution per 256*1024
                # read pairs. Extract the number of reads present in the
                # current block being processed e.g. # candidate unique pairs
                # for (FF, FR, RF, RR): (46,226102, 14, 28)
                if '# candidate unique pairs for' in line:
                    current_reads = int(
                        line.rstrip().replace(',', '').replace('(', '')
                        .replace(')', '').split()[-3]
                    )
                    sample.quast.total_reads += current_reads

                # Continue parsing to find the FR section of the current block
                elif size in line:
                    self._extract_insert_size_distribution(
                        read_stats, sample, current_reads
                    )

    def _extract_insert_size_distribution(
        self, read_stats, sample, current_reads
    ) -> None:
        """
        Extract the mean and standard deviation of the insert size for a block.

        Args:
            read_stats: The file object for the reads stats file.
            sample: The sample object to populate with extracted data.
            current_reads: The number of reads in the current block.
        """
        for sub_line in read_stats:

            # Extract the mean and standard deviation of the insert size for
            # this block [M::mem_pestat] mean and std.dev: (487.88, 246.14)
            if '[M::mem_pestat] mean and std.dev:' in sub_line:
                split_line = (
                    sub_line.rstrip().replace(',', '').replace('(', '')
                    .replace(')', '').split()
                )
                mean = float(split_line[-2])
                std = float(split_line[-1])
                sample.quast.insert_mean.append(mean)
                sample.quast.insert_std.append(std)
                sample.quast.read_blocks.append(current_reads)
                break

    def calculate_weighted_insert_size(self) -> None:
        """
        Calculate the weighted mean and standard deviation of the insert size
        from the extracted bwa blocks.

        This method iterates over the metadata samples, calculates the weighted
        mean and standard deviation of the insert size, and stores the results
        in the metadata object.

        Preconditions:
        - `self.metadata` should contain the sample data with necessary
          attributes.
        """
        logging.info('Calculating weighted insert size')

        for sample in self.metadata:
            logging.debug("Processing sample: %s", sample.name)

            # Initialize attributes to store the calculated weighted mean and
            # standard deviation of insert sizes
            sample.quast.mean_insert = float()
            sample.quast.std_insert = float()

            # Guard statement to check if insert_mean is 'ND'
            if sample.quast.insert_mean == 'ND':
                logging.warning(
                    "Insert mean is 'ND' for sample %s, skipping",
                    sample.name
                )
                continue

            # Guard statement to check if read_blocks is empty
            if not sample.quast.read_blocks:
                logging.warning(
                    "No read blocks found for sample %s, skipping",
                    sample.name
                )
                continue

            # Iterate through all the read blocks present in the sample
            for i, read_block in enumerate(sample.quast.read_blocks):

                # Calculate the weight of the current block by dividing it
                # (current number of reads) by the total number of reads in
                # the sample
                weight = read_block / sample.quast.total_reads

                # Multiply the mean for this block to obtain the weighted
                # mean, and add it to the total mean
                sample.quast.mean_insert += (
                    sample.quast.insert_mean[i] * weight
                )

                # Same calculation, but for standard deviation
                sample.quast.std_insert += (
                    sample.quast.insert_std[i] * weight
                )

            # Set the attributes to floats with two decimal places
            sample.quast.mean_insert = float(
                f'{sample.quast.mean_insert:.2f}'
            )
            sample.quast.std_insert = float(
                f'{sample.quast.std_insert:.2f}'
            )

    def qualimapper(self) -> None:
        """
        Create threads and commands for performing reference mapping for
        qualimap analyses.

        This method initializes threads and enqueues tasks for qualimap
        analysis. It waits for all tasks to be completed.

        Preconditions:
        - `self.metadata` should contain the sample data with necessary
          attributes.
        - `self.threads` should be set to the number of threads to use.
        - `self.log_file` should be the path to the log file.
        """
        logging.info('Running qualimap on samples')

        # Use ThreadPoolExecutor to manage threads
        with ThreadPoolExecutor(max_workers=self.threads) as executor:

            # Enqueue tasks for each sample
            for sample in self.metadata:
                self._initialize_sample_qualimap(sample)
                if sample.general.best_assembly_file != "NA":
                    executor.submit(self._run_qualimap, sample)

    def _initialize_sample_qualimap(self, sample) -> None:
        """
        Initialize the qualimap attributes for a sample.

        Args:
            sample: The sample object to initialize attributes for.
        """
        # Create and populate the qualimap attribute
        sample.qualimap = CustomBox()
        sample.qualimap.output_dir = os.path.join(
            sample.general.output_directory, 'qualimap'
        )
        os.makedirs(sample.qualimap.output_dir, exist_ok=True)
        sample.qualimap.report_file = os.path.join(
            sample.qualimap.output_dir, 'genome_results.txt'
        )

        # Initialize dictionaries to store qualimap results
        sample.qualimap.length = {}
        sample.qualimap.bases = {}
        sample.qualimap.coverage = {}
        sample.qualimap.stddev = {}

    def _run_qualimap(self, sample) -> None:
        """
        Run qualimap for a sample.

        Args:
            sample: The sample object to run qualimap on.
        """
        if not self._is_valid_sample(sample):
            logging.warning(
                "Invalid sample %s, skipping qualimap run", sample.name
            )
            return

        qualimap_call = self._construct_qualimap_command(sample)
        sample.commands.qualimap = qualimap_call

        if not os.path.isfile(sample.qualimap.report_file):
            self._execute_qualimap_command(sample)

    def _is_valid_sample(self, sample) -> bool:
        """
        Check if the sample is valid for running qualimap.

        Args:
            sample: The sample object to check.

        Returns:
            bool: True if the sample is valid, False otherwise.
        """
        if sample.general.best_assembly_file == "NA":
            return False
        if not hasattr(
                sample,
                'quast') or not hasattr(
                sample.quast,
                'sorted_bam'):
            return False
        if not hasattr(
                sample,
                'qualimap') or not hasattr(
                sample.qualimap,
                'output_dir'):
            return False
        return True

    def _construct_qualimap_command(self, sample) -> str:
        """
        Construct the qualimap command for a sample.

        Args:
            sample: The sample object to construct the command for.

        Returns:
            str: The constructed qualimap command.
        """
        return (
            f'qualimap bamqc -bam {sample.quast.sorted_bam} '
            f'-out_dir {sample.qualimap.output_dir}'
        )

    def _execute_qualimap_command(self, sample) -> None:
        """
        Execute the qualimap command for a sample.

        Args:
            sample: The sample object to run the command on.
        """
        try:
            out, err = run_subprocess(sample.commands.qualimap)
            self._log_qualimap_output(sample, out, err)
        except RuntimeError as exc:
            logging.error(
                "Failed to run qualimap command for sample %s: %s",
                sample.name, exc
            )

    def _log_qualimap_output(self, sample, out, err) -> None:
        """
        Log the output and errors from the qualimap command.

        Args:
            sample: The sample object to log information for.
            out: The standard output from the qualimap command.
            err: The standard error from the qualimap command.
        """
        write_to_log_file(
            out=f'{sample.commands.qualimap}\n{out}',
            err=err,
            log_file=self.log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err,
        )

    def parse_qualimap_report(self) -> None:
        """
        Parse the qualimap report.

        This method iterates over the metadata samples, reads the qualimap
        report, and populates the metadata object with the extracted key-value
        pairs.

        Preconditions:
        - `self.metadata` should contain the sample data with necessary
          attributes.
        """
        logging.info('Parsing Qualimap reports')

        for sample in self.metadata:
            logging.debug("Processing sample: %s", sample.name)

            # Check if the sample is valid for parsing
            if not self._is_valid_qualimap_sample(sample):
                logging.warning(
                    "Skipping sample %s as best assembly file is 'NA' or "
                    "report file not found", sample.name
                )
                continue

            # Parse the report file and extract key-value pairs
            qualimap_dict = self._parse_report_file(sample)
            # Update the sample object with the extracted key-value pairs
            self._update_sample_with_qualimap_dict(sample, qualimap_dict)

    def _is_valid_qualimap_sample(self, sample) -> bool:
        """
        Check if the sample is valid for parsing qualimap report.

        Args:
            sample: The sample object to check.

        Returns:
            bool: True if the sample is valid, False otherwise.
        """
        # Check if the best assembly file is not 'NA'
        if sample.general.best_assembly_file == "NA":
            return False
        # Check if the qualimap report file exists
        if not os.path.isfile(sample.qualimap.report_file):
            return False
        return True

    def _parse_report_file(self, sample) -> dict:
        """
        Parse the qualimap report file and extract key-value pairs.

        Args:
            sample: The sample object to parse the report for.

        Returns:
            A dictionary containing the extracted key-value pairs.
        """
        qualimap_dict = {}
        try:
            with open(sample.qualimap.report_file, encoding='utf-8') as report:
                for line in report:
                    # Sanitize the keys and values using self.qualimap_analyze
                    key, value = self.qualimap_analyze(line)

                    # If the keys and values exist, enter them into the dict
                    if (key, value) != (None, None):
                        # Only keep two decimal places for float values
                        if isinstance(value, float):
                            value = float(f'{value:.2f}')
                        qualimap_dict[key] = value

                    # Parse the coverage per contig section
                    if 'Coverage per contig' in line:
                        self._parse_contig_coverage(report, sample)
        except (IOError, FileNotFoundError) as exc:
            logging.error(
                "Error reading qualimap report file for sample %s: %s",
                sample.name, exc
            )
        return qualimap_dict

    def _parse_contig_coverage(self, report, sample) -> None:
        """
        Parse the coverage per contig section of the qualimap report.

        Args:
            report: The file object for the qualimap report.
            sample: The sample object to populate with extracted data.
        """
        for contig_line in report:
            try:
                # Extract the contig coverage details
                _, name, length, bases, coverage, stddev = (
                    contig_line.rstrip().split('\t')
                )

                # Update the sample object with the extracted details
                sample.qualimap.length[name] = length
                sample.qualimap.bases[name] = bases
                sample.qualimap.coverage[name] = coverage
                sample.qualimap.stddev[name] = stddev
            except ValueError:
                logging.debug(
                    "Skipping malformed line in contig coverage section: %s",
                    contig_line
                )

    def _update_sample_with_qualimap_dict(self, sample, qualimap_dict) -> None:
        """
        Update the sample object with the extracted key-value pairs.

        Args:
            sample: The sample object to update.
            qualimap_dict: The dictionary containing the extracted key-value
                pairs.
        """
        if qualimap_dict:
            for attribute, value in qualimap_dict.items():
                # Remove the 'X' from the depth values e.g. 40.238X
                setattr(sample.qualimap, attribute, value.rstrip('X'))

    def pilon(self) -> None:
        """
        Run pilon to fix any misassemblies in the contigs - will look for SNPs
        and indels.

        This method initializes threads and enqueues tasks for pilon analysis.
        It waits for all tasks to be completed.

        Preconditions:
        - `self.metadata` should contain the sample data with necessary
          attributes.
        - `self.threads` should be set to the number of threads to use.
        - `self.log_file` should be the path to the log file.
        """
        logging.info('Improving quality of assembly with pilon')

        # Use ThreadPoolExecutor to manage threads
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            # Enqueue tasks for each sample
            for sample in self.metadata:
                if self._is_valid_pilon_sample(sample):
                    logging.debug(
                        "Initializing pilon for sample: %s", sample.name)
                    self._initialize_sample_pilon(sample)
                    logging.debug(
                        "Submitting pilon task for sample: %s", sample.name)
                    executor.submit(self._run_pilon, sample)

    def _is_valid_pilon_sample(self, sample) -> bool:
        """
        Check if the sample is valid for running pilon.

        Args:
            sample: The sample object to check.

        Returns:
            bool: True if the sample is valid, False otherwise.
        """
        logging.debug("Checking if sample is valid for pilon: %s", sample.name)

        # Check if the best assembly file is not 'NA'
        if sample.general.best_assembly_file == 'NA':
            logging.debug(
                "Sample %s is not valid for pilon: best_assembly_file is 'NA'",
                sample.name)
            return False

        # Check if the necessary attributes are present
        if not hasattr(
                sample,
                'quast') or not hasattr(
                sample.quast,
                'sorted_bam'):
            logging.debug(
                "Sample %s is not valid for pilon: missing quast outputs or "
                "sorted_bam file", sample.name
            )
            return False
        if not hasattr(
                sample,
                'general') or not hasattr(
                sample.general,
                'assembly_file'):
            logging.debug(
                "Sample %s is not valid for pilon: missing general or "
                "assembly_file", sample.name
            )
            return False
        return True

    def _initialize_sample_pilon(self, sample) -> None:
        """
        Initialize the pilon attributes for a sample.

        Args:
            sample: The sample object to initialize attributes for.
        """
        logging.debug(
            "Initializing pilon attributes for sample: %s",
            sample.name)
        # Initialize the pilon attribute
        sample.pilon = CustomBox()

        # Set the contigs file to the assembly file
        sample.general.contigs_file = sample.general.assembly_file

        # Set the output directory for pilon
        sample.pilon.out_dir = os.path.join(sample.quast.output_dir, 'pilon')
        os.makedirs(sample.pilon.out_dir, exist_ok=True)

        # Construct the pilon command using f-strings
        sample.pilon.cmd = (
            f'pilon --genome {sample.general.contigs_file} '
            f'--bam {sample.quast.sorted_bam} --fix bases '
            f'--threads {self.threads} --out_dir {sample.pilon.out_dir} '
            f'--changes --mindepth 0.25'
        )
        logging.debug(
            "Pilon command for sample %s: %s",
            sample.name,
            sample.pilon.cmd)

    def _run_pilon(self, sample) -> None:
        """
        Run pilon for a sample.

        Args:
            sample: The sample object to run pilon on.
        """
        logging.debug("Running pilon for sample: %s", sample.name)
        # Check if the contigs file already exists
        if not os.path.isfile(sample.general.contigs_file):
            logging.debug(
                "Contigs file does not exist for sample: %s",
                sample.name)
            self._execute_pilon_command(sample)

    def _execute_pilon_command(self, sample) -> None:
        """
        Execute the pilon command for a sample.

        Args:
            sample: The sample object to run the command on.
        """
        try:
            logging.debug(
                "Executing pilon command for sample: %s",
                sample.name)
            # Run the pilon command
            out, err = run_subprocess(sample.pilon.cmd)
            # Log the output and errors
            self._log_pilon_output(sample, out, err)
        except RuntimeError as exc:
            logging.error(
                "Failed to run pilon command for sample %s: %s",
                sample.name, exc
            )

    def _log_pilon_output(self, sample, out, err) -> None:
        """
        Log the output and errors from the pilon command.

        Args:
            sample: The sample object to log information for.
            out: The standard output from the pilon command.
            err: The standard error from the pilon command.
        """
        logging.debug("Logging pilon output for sample: %s", sample.name)
        write_to_log_file(
            out=f'{sample.pilon.cmd}\n{out}',
            err=err,
            log_file=self.log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err,
        )

    def filter(self) -> None:
        """
        Filter contigs based on depth and length.

        This method initializes threads and enqueues tasks for filtering
        contigs. It waits for all tasks to be completed.

        Preconditions:
        - `self.metadata` should contain the sample data with necessary
          attributes.
        - `self.threads` should be set to the number of threads to use.
        """
        logging.info('Filtering contigs')

        # Use ThreadPoolExecutor to manage threads
        with ThreadPoolExecutor(max_workers=self.threads) as executor:
            # Enqueue tasks for each sample
            for sample in self.metadata:
                if self._is_valid_filter_sample(sample):
                    logging.debug(
                        "Initializing filter for sample: %s", sample.name)
                    self._initialize_sample_filter(sample)
                    logging.debug(
                        "Submitting filter task for sample: %s", sample.name)
                    executor.submit(self._run_filter, sample)

    def _is_valid_filter_sample(self, sample) -> bool:
        """
        Check if the sample is valid for filtering.

        Args:
            sample: The sample object to check.

        Returns:
            bool: True if the sample is valid, False otherwise.
        """
        logging.debug(
            "Checking if sample is valid for filtering: %s",
            sample.name)

        # Check if the best assembly file is not 'NA'
        if sample.general.best_assembly_file == 'NA':
            logging.debug(
                "Sample %s is not valid for filtering: best_assembly_file is "
                "'NA'", sample.name
            )
            return False

        # Check if the necessary attributes are present
        if not hasattr(
                sample,
                'general') or not hasattr(
                sample.general,
                'assembly_file'):
            logging.debug(
                "Sample %s is not valid for filtering: missing assembly_file",
                sample.name
            )
            return False
        return True

    def _initialize_sample_filter(self, sample) -> None:
        """
        Initialize the filter attributes for a sample.

        Args:
            sample: The sample object to initialize attributes for.
        """
        logging.debug(
            "Initializing filter attributes for sample: %s",
            sample.name)
        # Set the name of the unfiltered assembly output file
        sample.general.contigs_file = sample.general.assembly_file

    def _run_filter(self, sample) -> None:
        """
        Run the filter process for a sample.

        Args:
            sample: The sample object to run the filter on.
        """
        logging.debug("Running filter for sample: %s", sample.name)
        # Check if the contigs file exists
        if os.path.isfile(sample.general.contigs_file):
            logging.debug("Contigs file exists for sample: %s", sample.name)
            # Filter the contigs based on depth and length
            self._filter_contigs(sample)
            # Copy the filtered file to the BestAssemblies folder
            self._copy_filtered_file(sample)
        else:
            logging.debug(
                "Contigs file does not exist for sample: %s",
                sample.name)
            # Set the best assembly file to 'NA' if contigs file doesn't exist
            sample.general.best_assembly_file = 'NA'

    def _filter_contigs(self, sample) -> None:
        """
        Filter contigs based on depth and length for a sample.

        Args:
            sample: The sample object to filter contigs for.
        """
        logging.debug("Filtering contigs for sample: %s", sample.name)
        pass_depth = []

        # Open the contigs file for reading
        with open(sample.general.contigs_file, 'r', encoding='utf-8') as file:
            # Parse the contigs file in FASTA format
            for record in SeqIO.parse(file, "fasta"):
                # Extract the contig name without '_pilon'
                contig = record.id.split('_pilon')[0]
                # Get the mean coverage and standard deviation
                coverage_mean = float(sample.qualimap.coverage[contig])
                coverage_std = float(sample.qualimap.stddev[contig])

                # Check if the contig passes the depth and length filters
                if (float(sample.qualimap.coverage[contig]) >
                    (coverage_mean - coverage_std * 1.5) and
                        len(record.seq) > 500):

                    # Replace 'Contig' in the record ID with the sample name
                    new_id = re.sub("Contig", sample.name, record.id)
                    record.id = new_id
                    # Clear the name and description attributes
                    record.name = ''
                    record.description = ''
                    # Add the record to the list of passing contigs
                    pass_depth.append(record)

        # Check if there are any contigs that passed the filters
        if pass_depth:
            logging.debug(
                "Contigs passed the filters for sample: %s",
                sample.name)
            # Set the name of the filtered file
            filtered = sample.general.filtered_file

            # Open the filtered file for writing
            with open(filtered, 'w', encoding='utf-8') as formatted:
                # Write the passing contigs to the filtered file
                SeqIO.write(pass_depth, formatted, 'fasta')

    def _copy_filtered_file(self, sample) -> None:
        """
        Copy the filtered file to the BestAssemblies folder.

        Args:
            sample: The sample object to copy the filtered file for.
        """
        logging.debug("Copying filtered file for sample: %s", sample.name)
        if os.path.isfile(sample.general.filtered_file):
            sample.general.best_assemblies_path = os.path.join(
                self.sequence_path, 'BestAssemblies'
            )
            best_assembly_file = os.path.join(
                sample.general.best_assemblies_path, f'{sample.name}.fasta'
            )
            sample.general.best_assembly_file = best_assembly_file

            if not os.path.isfile(best_assembly_file):
                logging.debug(
                    "Copying file to BestAssemblies for sample: %s",
                    sample.name
                )
                shutil.copyfile(
                    sample.general.filtered_file,
                    best_assembly_file)
        else:
            logging.debug(
                "Filtered file does not exist for sample: %s",
                sample.name
            )
            sample.general.best_assembly_file = 'NA'

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
        """
        Analyze a line from the qualimap report.

        Args:
            line: A line from the qualimap report.

        Returns:
            A tuple containing the sanitized key and value.
        """
        # Check if the line contains ' = '
        if ' = ' in line:
            key, value = line.split(' = ', 1)

            # Sanitize the key
            key = (
                key.replace('number of ', "")
                .replace("'", "")
                .title()
                .replace(" ", "")
            )

            # Sanitize the value
            value = value.replace(",", "").replace(" ", "").rstrip()
        else:
            # Set the keys and values to None if ' = ' is not found
            key, value = None, None

        return key, value

    @staticmethod
    def analyze(line):
        """
        Analyze a line from the report.

        Args:
            line: A line from the report.

        Returns:
            A tuple containing the sanitized key and value.
        """
        # Split the line on tab character
        key, value = line.rstrip().split('\t', 1)

        # Sanitize the key
        key = (
            key.replace(' (%)', '')
            .replace(' ', '_')
            .replace('#', 'num')
            .replace('(>=', 'greater_than')
            .replace(')', '')
            .replace('.', '')
            .replace('\'', '')
        )

        return key, value


class Parser:

    """
    Parser class to allow running the assembly evaluation from the command line
    """
    def __init__(self):
        """
        Initialize the parser and set up the arguments.
        """
        parser = ArgumentParser(
            description='Calculates coverage depth by mapping FASTQ reads '
                        'against assemblies'
        )
        parser.add_argument(
            '-p', '--path',
            default=os.getcwd(),
            help='Specify the path of the folder that either contains the '
                 'files of interest, or will be used to store the outputs'
        )
        parser.add_argument(
            '-a', '--assemblies',
            help='Path to a folder of assemblies. If not provided, the script '
                 'will look for .fa or .fasta files in the path'
        )
        parser.add_argument(
            '-f', '--fastq',
            help='Path to a folder of fastq files. If not provided, the '
            'script will look for fastq or .fastq.gz files in the path'
        )
        parser.add_argument(
            '-t', '--threads',
            type=int,
            default=multiprocessing.cpu_count(),
            help='Number of threads. Default is the number of cores in the '
                 'system'
        )

        # Get the arguments into an object
        args = parser.parse_args()

        # Define variables from the arguments
        self.sequence_path = os.path.join(args.path, '')
        self.assembly_path = os.path.join(
            args.assemblies, '') if args.assemblies else self.sequence_path
        self.fastq_path = os.path.join(
            args.fastq, '') if args.fastq else self.sequence_path
        self.threads = args.threads

        # Initialize variables
        self.strains = []
        self.metadata = []
        self.log_file = os.path.join(self.sequence_path, 'log_file.txt')

        # Associate the assemblies and fastq files in a metadata object
        self.associate()

        # Evaluate the assemblies
        AssemblyEvaluation(
            log_file=self.log_file,
            metadata=self.metadata,
            sequence_path=self.sequence_path,
            threads=self.threads
        )

    def associate(self):
        """
        Associate assemblies and fastq files in a metadata object.
        """
        # Get the sequences in the sequences folder into a list
        self.strains = [
            fasta
            for fasta in sorted(
                glob(os.path.join(self.assembly_path, '*.fa*')))
            if '.fastq' not in fasta]

        for strain in self.strains:
            # Extract the name of the strain from the path and file extension
            strain_name = os.path.split(strain)[1].split('.')[0]

            # Find the corresponding fastq files for each strain
            fastq_files = sorted(
                glob(
                    os.path.join(
                        self.fastq_path,
                        f'{strain_name}*fastq*')))

            # Ensure that fastq files are present for each assembly
            assert fastq_files, f'Cannot find fastq files for strain {
                strain_name} '

            # Create the metadata object
            metadata = CustomBox()
            metadata.name = strain_name
            metadata.general = CustomBox()
            metadata.general.best_assemblies_path = self.assembly_path
            metadata.general.trimmed_fastq_files = fastq_files
            metadata.general.output_directory = os.path.join(
                self.sequence_path, strain_name)
            metadata.mapping = CustomBox()

            # Append the metadata for each sample to the list of samples
            self.metadata.append(metadata)


if __name__ == '__main__':
    Parser()
