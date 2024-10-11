#!/usr/bin/env python3

"""
FastQC functions
"""

# Standard imports
from concurrent.futures import as_completed, ThreadPoolExecutor
from glob import glob
import logging
import os
import shutil
from queue import Queue
from typing import List, Tuple

# Third-party imports
from olctools.accessoryFunctions.accessoryFunctions import (
    run_subprocess,
    write_to_log_file,
)
from olctools.accessoryFunctions.metadata import CustomBox

qc_queue = Queue()
__author__ = 'adamkoziol'


def fastqc_threader(
        level: str, log_file: str, metadata: List[CustomBox], threads: int
) -> List[CustomBox]:
    """
    Run quality control on FASTQ files using FastQC.

    Args:
        level (str): The level of processing (e.g., 'Trimmed', 'merged').
        log_file (str): Path to the log file.
        metadata (List[CustomBox]): List of metadata sample objects.
        threads (int): Number of threads to use.

    Returns:
        List[CustomBox]: Updated metadata sample objects.
    """
    logging.info('Running quality control on %s fastq files', level)

    # Create and start threads for each FASTQ file in the list
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [
            executor.submit(
                fastqc,
                sample=sample,
                level=level,
                log_file=log_file,
                threads=threads) for sample in metadata if isinstance(
                sample.general.fastq_files,
                list)]

    # Wait on the queue until everything has been processed
    qc_queue.join()

    # Collect results from futures
    updated_metadata = []
    for future in as_completed(futures):
        updated_metadata.append(future.result())

    return updated_metadata


def fastqc(
        sample: CustomBox, level: str, log_file: str, threads: int
) -> CustomBox:
    """
    Run FastQC on the given sample.

    Args:
        sample (CustomBox): Metadata sample object.
        level (str): The level of processing (e.g., 'Trimmed', 'merged').
        log_file (str): Path to the log file.
        threads (int): Number of threads to use.

    Returns:
        CustomBox: Updated metadata sample object.
    """
    logging.debug('Preparing FastQC call for sample: %s', sample.name)

    # Prepare the FastQC system call and reads call
    fastqc_call, fastqc_reads = prepare_fastqc_call(
        sample=sample, level=level, threads=threads
    )

    # If a FastQC call is prepared, add it to the queue
    if fastqc_call:
        sample.commands.fastqc = fastqc_call
        setattr(sample.commands, f'fastqc_{level.lower()}', fastqc_call)
        qc_queue.put((sample, fastqc_call, fastqc_reads, log_file, level))
        logging.debug('Added FastQC call to queue for sample: %s', sample.name)

    # Process the FastQC queue
    while not qc_queue.empty():
        process_fastqc_queue()

    return sample


def prepare_fastqc_call(
        sample: CustomBox, level: str, threads: int
) -> Tuple[str, str]:
    """
    Prepare the FastQC system call based on the processing level.

    Args:
        sample (CustomBox): Metadata sample object.
        level (str): The level of processing (e.g., 'trimmed', 'merged').
        threads (int): Number of threads to use.

    Returns:
        tuple: FastQC system call and reads call.
    """
    logging.debug(
        'Preparing FastQC call for level: %s, sample: %s',
        level,
        sample.name)

    reader = 'cat'
    fastq_files = None

    # Determine the reader and FASTQ files based on the processing level
    if level == 'trimmed':
        reader, fastq_files = get_reader_and_files(
            sample=sample, attr='trimmed_fastq_files'
        )
    elif level == 'trimmed_corrected':
        reader, fastq_files = get_reader_and_files(
            sample=sample, attr='trimmed_corrected_fastq_files'
        )
    elif level == 'normalised':
        reader, fastq_files = get_reader_and_files(
            sample=sample, attr='normalised_reads'
        )
    elif level == 'merged':
        reader, fastq_files = get_reader_and_files(
            sample=sample, attr='merged_reads', single_file=True
        )
    else:
        reader, fastq_files = get_reader_and_files(
            sample=sample, attr='fastq_files'
        )

    # If no valid FASTQ files are found, return empty strings
    if not isinstance(fastq_files, list):
        logging.debug(
            'No valid %s FASTQ files found for sample: %s',
            level, sample.name
        )
        return '', ''

    # Set the output directory for FastQC results
    out_dir = os.path.join(sample.general.output_directory, 'fastqc', level)
    os.makedirs(out_dir, exist_ok=True)

    # Prepare the FastQC system call and reads call based on the number of
    # FASTQ files
    if len(fastq_files) == 2:
        fastqc_call = (
            f'{reader} {fastq_files[0]} {fastq_files[1]} | fastqc -q -t '
            f'{threads} stdin -o {out_dir}'
        )
        fastqc_reads = (
            f"fastqc {fastq_files[0]} {fastq_files[1]} -q -o {out_dir} -t "
            f"{threads}"
        )
    elif len(fastq_files) == 1:
        fastqc_call = (
            f'{reader} {fastq_files[0]} | fastqc -q -t {threads} stdin -o '
            f'{out_dir}'
        )
        fastqc_reads = (
            f"fastqc {fastq_files[0]} -q -o {out_dir} -t {threads}"
        )
    else:
        fastqc_call = ''
        fastqc_reads = ''

    logging.debug(
        'Prepared FastQC call for %s reads for sample: %s',
        level, sample.name
    )
    return fastqc_call, fastqc_reads


def get_reader_and_files(
        sample: CustomBox, attr: str, single_file: bool = False
) -> Tuple[str, List[str]]:
    """
    Get the appropriate reader and FASTQ files based on the attribute.

    Args:
        sample (CustomBox): Metadata sample object.
        attr (str): Attribute name to get the FASTQ files.
        single_file (bool): Whether to expect a single file or a list of files.

    Returns:
        tuple: Reader command and FASTQ files.
    """
    try:
        # Get the FASTQ files from the specified attribute
        fastq_files = getattr(sample.general, attr)
        if single_file:
            fastq_files = [fastq_files]

        # Determine the appropriate reader based on the file extension
        if '.gz' in fastq_files[0]:
            reader = 'gunzip --to-stdout'
        elif '.bz2' in fastq_files[0]:
            reader = 'bunzip2 --stdout'
        else:
            reader = 'cat'
    except AttributeError:
        logging.debug(
            'AttributeError: %s not found in sample: %s',
            attr,
            sample.name
        )
        reader = 'cat'
        fastq_files = []

    logging.debug(
        'Reader: %s, FASTQ files: %s for sample: %s',
        reader,
        fastq_files,
        sample.name
    )
    return reader, fastq_files


def process_fastqc_queue() -> None:
    """
    Process the FastQC queue.
    """
    while not qc_queue.empty():
        # Get the next item from the queue
        sample, system_call, fastqc_reads, log_file, level = qc_queue.get()
        output_dir = os.path.join(sample.general.output_directory, 'fastqc')

        try:
            # Check if the FastQC output HTML file already exists
            _ = glob(os.path.join(output_dir, '*.html'))[0]
            logging.debug(
                'FastQC output already exists for sample: %s',
                sample.name)
        except IndexError:
            # If the output file does not exist, create the output directory
            os.makedirs(output_dir, exist_ok=True)
            logging.debug(
                'Running FastQC for %s reads for sample: %s',
                level, sample.name
            )

            # Run the FastQC system calls and log the output
            run_fastqc(
                system_call=system_call,
                fastqc_reads=fastqc_reads,
                log_file=log_file,
                level=level,
                output_dir=output_dir,
                sample=sample
            )
            # Rename the FastQC output files
            rename_fastqc_outputs(
                level=level,
                output_dir=output_dir,
                sample=sample
            )

        # Mark the task as done
        qc_queue.task_done()


def run_fastqc(
        system_call: str,
        fastqc_reads: str,
        log_file: str,
        level: str,
        output_dir: str,
        sample: CustomBox
) -> Tuple[str, str]:
    """
    Run the FastQC system calls and log the output.

    Args:
        system_call (str): FastQC system call.
        fastqc_reads (str): FastQC reads call.
        log_file (str): Path to the log file.
        level (str): The level of processing (e.g., 'trimmed', 'merged').
        output_dir (str): Output directory.
        sample (CustomBox): Metadata sample object.

    Returns:
        tuple: Standard output and error strings.
    """
    logging.debug(
        'Running FastQC system call for %s reads for sample: %s',
        level, sample.name
    )

    # Set the name of the output file
    fastqc_html = os.path.join(output_dir, f'{sample.name}_fastqc.html')

    # Only run the system call if the FastQC outputs don't already exist
    if os.path.isfile(fastqc_html):
        logging.debug(
            'FastQC output already exists for sample: %s, skipping FastQC',
            sample.name
        )
        return "", ""

    out_str, err_str = '', ''

    # Run the first FastQC system call
    out, err = _run_command(system_call)
    out_str += out
    err_str += err

    logging.debug(
        'Running FastQC reads call for %s reads for sample: %s',
        level, sample.name
    )

    # Run the second FastQC reads call
    out, err = _run_command(fastqc_reads)
    out_str += out
    err_str += err

    # Write the system call and reads call to the log file
    _log_fastqc_output(
        system_call,
        fastqc_reads,
        out_str,
        err_str,
        log_file,
        sample)

    logging.debug(
        'Completed FastQC for %s reads for sample: %s',
        level, sample.name
    )
    return out_str, err_str


def _run_command(command: str) -> Tuple[str, str]:
    """
    Run a subprocess command and return the output and error.

    Args:
        command (str): The command to run.

    Returns:
        tuple: Standard output and error strings.
    """
    logging.debug('Running command: %s', command)
    out, err = run_subprocess(command=command)
    return out, err


def _log_fastqc_output(
        system_call: str,
        fastqc_reads: str,
        out_str: str,
        err_str: str,
        log_file: str,
        sample: CustomBox
) -> None:
    """
    Log the output and errors from the FastQC commands.

    Args:
        system_call (str): FastQC system call.
        fastqc_reads (str): FastQC reads call.
        out_str (str): Standard output string.
        err_str (str): Standard error string.
        log_file (str): Path to the log file.
        sample (CustomBox): Metadata sample object.
    """
    logging.debug('Logging FastQC output for sample: %s', sample.name)
    write_to_log_file(
        out=system_call,
        err=system_call,
        log_file=log_file,
        sample_log=sample.general.log_out,
        sample_err=sample.general.log_err
    )
    write_to_log_file(
        out=fastqc_reads,
        err=fastqc_reads,
        log_file=log_file,
        sample_log=sample.general.log_out,
        sample_err=sample.general.log_err
    )
    write_to_log_file(
        out=out_str,
        err=err_str,
        log_file=log_file,
        sample_log=sample.general.log_out,
        sample_err=sample.general.log_err
    )


def rename_fastqc_outputs(
    level: str,
    output_dir: str,
    sample: CustomBox
) -> None:
    """
    Rename the FastQC output files.

    Args:
        level (str): The level of processing (e.g., 'trimmed', 'merged').
        output_dir (str): Output directory.
        sample (CustomBox): Metadata sample object.
    """
    try:
        # Rename the FastQC HTML output file
        shutil.move(
            src=os.path.join(output_dir, 'stdin_fastqc.html'),
            dst=os.path.join(output_dir, f'{sample.name}_fastqc.html')
        )
        # Rename the FastQC ZIP output file
        shutil.move(
            src=os.path.join(output_dir, 'stdin_fastqc.zip'),
            dst=os.path.join(output_dir, f'{sample.name}_fastqc.zip')
        )
        logging.debug(
            'Renamed FastQC output files for %s reads for sample: %s',
            level, sample.name
        )
    except IOError:
        logging.error(
            'Failed to rename FastQC output files for %s reads for sample: %s',
            level, sample.name
        )
