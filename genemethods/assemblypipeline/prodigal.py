#!/usr/bin/env python3
"""
Runs Prodigal for gene prediction on sequence data.
"""

# Standard imports
from concurrent.futures import ThreadPoolExecutor
import logging
import os
from queue import Queue
from threading import Lock
from typing import List

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    run_subprocess,
    write_to_log_file
)
from olctools.accessoryFunctions.metadata import CustomBox

__author__ = 'adamkoziol'

# Initialise a thread lock for thread-safe logging
thread_lock = Lock()


class Prodigal:
    """
    Class to perform gene predictions using Prodigal.
    """

    def __init__(
        self,
        log_file,
        metadata
    ):
        """
        Initialize the Prodigal class.

        Args:
            input_object: An object containing metadata and other parameters.

        Preconditions:

        """
        logging.debug("Initializing Prodigal class")

        # Extract necessary attributes from the input object
        self.metadata = metadata
        self.log_file = log_file

        # Initialize queues for threading
        self.predict_queue = Queue()
        self.parse_queue = Queue()

    def main(self) -> List:
        """
        Run the prediction and parsing processes
        """
        # Start the prediction and parsing processes
        self.predict_threads()
        self.prodigal_parse()

        # Return the updated metadata
        return self.metadata

    def predict_threads(self) -> None:
        """
        Create threads for gene prediction.

        Preconditions:
        - self.metadata must be a list of sample objects.
        """
        logging.info('Performing gene predictions')
        logging.debug("Creating threads for gene prediction")

        # Use ThreadPoolExecutor to manage threads
        with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
            for sample in self.metadata:
                if sample.general.best_assembly_file != 'NA':
                    logging.debug(
                        "Submitting predict task for sample: %s", sample.name
                    )
                    # Submit the predict task to the executor
                    executor.submit(self.predict, sample)

        # Iterate over the samples
        for sample in self.metadata:
            # Initialize the Prodigal attribute for each sample
            sample.prodigal = CustomBox()
            if sample.general.best_assembly_file != 'NA':
                logging.debug(
                    "Adding sample to predict queue: %s", sample.name
                )
                # Add the sample to the prediction queue
                self.predict_queue.put(sample)

        # Wait for all tasks in the queue to be completed
        self.predict_queue.join()

    def predict(self, sample) -> None:
        """
        Perform gene prediction on a sample.

        Args:
            sample: A sample object containing metadata and file paths.

        Preconditions:
        - sample must have attributes: general.best_assembly_file, name,
          general.output_directory.
        """
        logging.debug("Starting gene prediction for sample: %s", sample.name)

        while True:
            sample = self.predict_queue.get()

            # Populate attributes for Prodigal results
            self._populate_prodigal_attributes(sample)

            # Create the folder to store the reports
            self._create_report_directory(sample)

            # Check if the report already exists and is not empty
            if self._is_report_needed(sample):
                # Run the Prodigal command
                self._run_prodigal_command(sample)

            logging.debug(
                "Completed gene prediction for sample: %s", sample.name
            )

            # Mark the task as done
            self.predict_queue.task_done()

    def _populate_prodigal_attributes(self, sample) -> None:
        """
        Populate attributes for Prodigal results.

        Args:
            sample: A sample object containing metadata and file paths.
        """
        sample.prodigal.report_dir = os.path.join(
            sample.general.output_directory, 'prodigal'
        )
        sample.prodigal.results_file = os.path.join(
            sample.prodigal.report_dir,
            f'{sample.name}_prodigalresults.sco'
        )
        sample.prodigal.results = sample.prodigal.results_file

        # Create a variable to store the genes.fa path
        fasta_path = os.path.join(
            sample.prodigal.report_dir, f"{sample.name}_genes.fa"
        )

        sample.commands.prodigal = (
            f'prodigal -i {sample.general.best_assembly_file} '
            f'-o {sample.prodigal.results_file} -f sco -d '
            f'{fasta_path}'
        )

        logging.debug(
            "Prodigal command for sample %s: %s",
            sample.name, sample.commands.prodigal
        )

    def _create_report_directory(self, sample) -> None:
        """
        Create the folder to store the reports.

        Args:
            sample: A sample object containing metadata and file paths.
        """
        os.makedirs(sample.prodigal.report_dir, exist_ok=True)

    def _is_report_needed(self, sample) -> bool:
        """
        Check if the report already exists and is not empty.

        Args:
            sample: A sample object containing metadata and file paths.

        Returns:
            bool: True if the report is needed, False otherwise.
        """
        return not os.path.isfile(sample.prodigal.results_file) or \
            os.stat(sample.prodigal.results_file).st_size == 0

    def _run_prodigal_command(self, sample) -> None:
        """
        Run the Prodigal command and log the output.

        Args:
            sample: A sample object containing metadata and file paths.
        """
        out, err = run_subprocess(sample.commands.prodigal)
        with thread_lock:
            # Log the command and its output
            write_to_log_file(
                out=sample.commands.prodigal,
                err=sample.commands.prodigal,
                log_file=self.log_file,
                sample_log=sample.general.log_out,
                sample_err=sample.general.log_err
            )
            write_to_log_file(
                out=out,
                err=err,
                log_file=self.log_file,
                sample_log=sample.general.log_out,
                sample_err=sample.general.log_err
            )

    def prodigal_parse(self) -> None:
        """
        Parse the results of gene predictions.

        Preconditions:
        - self.metadata must be a list of sample objects.
        """
        logging.info('Parsing gene predictions')

        for sample in self.metadata:
            # Initialize Prodigal attributes for the sample
            self._initialize_prodigal_attributes(sample)

            if sample.general.best_assembly_file != 'NA':
                logging.debug(
                    "Parsing Prodigal results for sample: %s", sample.name
                )
                # Parse the Prodigal results for the sample
                self._parse_prodigal_results(sample)

    def _initialize_prodigal_attributes(self, sample) -> None:
        """
        Initialize Prodigal attributes for a sample.

        Args:
            sample: A sample object containing metadata and file paths.

        Preconditions:
        - sample must have a prodigal attribute.
        """

        # Initialize counters for predicted genes
        sample.prodigal.predicted_genes_total = 0
        sample.prodigal.predicted_genes_over_3000bp = 0
        sample.prodigal.predicted_genes_over_1000bp = 0
        sample.prodigal.predicted_genes_over_500bp = 0
        sample.prodigal.predicted_genes_under_500bp = 0

    def _parse_prodigal_results(self, sample) -> None:
        """
        Parse the Prodigal results for a sample.

        Args:
            sample: A sample object containing metadata and file paths.

        Preconditions:
        - sample.prodigal.results must be a valid file path.
        """
        logging.debug(
            "Opening Prodigal results file for sample: %s", sample.name
        )

        # Open the Prodigal results file for reading
        with open(sample.prodigal.results, 'r', encoding='utf-8') as results:
            for line in results:
                if line.startswith('>'):
                    # Extract the start and end positions of the gene
                    start, end = map(int, line.split('_')[1:3])
                    length = abs(start - end)
                    sample.prodigal.predicted_genes_total += 1

                    # Categorize the gene based on its length
                    if length > 3000:
                        sample.prodigal.predicted_genes_over_3000bp += 1
                    elif length > 1000:
                        sample.prodigal.predicted_genes_over_1000bp += 1
                    elif length > 500:
                        sample.prodigal.predicted_genes_over_500bp += 1
                    else:
                        sample.prodigal.predicted_genes_under_500bp += 1
