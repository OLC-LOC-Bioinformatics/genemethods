#!/usr/bin/env python3

"""
FASTQ validation, reformatting, and repair functions
"""

# Standard imports
import logging
import os
from subprocess import CalledProcessError
import traceback
from typing import List

# Third-party imports
from genewrappers.biotools import bbtools
from olctools.accessoryFunctions.accessoryFunctions import (
    write_to_log_file,
    write_metadata_to_file
)
from olctools.accessoryFunctions.metadata import CustomBox


def validate_fastq(
        log_file: str,
        metadata: List[CustomBox]) -> List[CustomBox]:
    """
    Runs reformat.sh on the FASTQ files. If a CalledProcessError arises, do
    not proceed with the assembly of these files.

    Args:
        log_file (str): Path to the log file.
        metadata (List[CustomBox]): List of metadata sample objects.

    Returns:
        List[CustomBox]: Updated metadata sample objects.
    """
    logging.info('Validating FASTQ files')

    # Iterate over each sample in the metadata
    for sample in metadata:
        logging.debug('Validating sample: %s', sample.name)

        # Check if the FASTQ file size is valid
        if not validate_file_size(sample):
            # Log an error if the file size is too small
            logging.error('File size too small for sample: %s', sample.name)
            error(sample, 'files_too_small')
            continue

        try:
            # Validate the reads using reformat.sh
            logging.debug('Running reformat.sh for sample: %s', sample.name)
            validate_reads(log_file, sample)
            logging.debug('Validation successful for sample: %s', sample.name)
        except CalledProcessError as e:
            # Handle any validation errors by attempting to repair the reads
            logging.error(
                'Validation failed for sample: %s with error: %s',
                sample.name, str(e)
            )
            handle_validation_error(log_file, sample)

    # Write the updated metadata to a file
    logging.debug('Writing updated metadata to file')
    write_metadata_to_file(metadata)

    logging.info('FASTQ validation completed')
    return metadata


def validate_file_size(sample: CustomBox) -> bool:
    """
    Validate the size of the FASTQ file.

    Args:
        sample (CustomBox): Metadata sample object.

    Returns:
        bool: True if the file size is valid, False otherwise.
    """
    logging.debug('Validating file size for sample: %s', sample.name)

    # Get the size of the first FASTQ file
    fastq_file = sample.general.fastq_files[0]
    size = os.path.getsize(fastq_file)
    logging.debug('File size for %s: %d bytes', fastq_file, size)

    # Check if the file size is greater than or equal to 1,000,000 bytes
    is_valid = size >= 1000000
    if not is_valid:
        logging.warning('File size too small for sample: %s', sample.name)

    return is_valid


def validate_reads(log_file: str, sample: CustomBox) -> None:
    """
    Validate the reads using reformat.sh.

    Args:
        log_file (str): Path to the log file.
        sample (CustomBox): Metadata sample object.

    Raises:
        CalledProcessError: If reformat.sh fails.
    """
    logging.debug('Validating reads for sample: %s', sample.name)

    try:
        # Run reformat.sh to validate the reads
        logging.debug('Running reformat.sh for sample: %s', sample.name)
        out, err, _ = bbtools.validate_reads(
            forward_in=sample.general.fastq_files[0],
            returncmd=True
        )
        logging.debug('reformat.sh completed for sample: %s', sample.name)

        # Write the output and error messages to the log file
        logging.debug(
            'Writing reformat.sh output to log file for sample: %s',
            sample.name)
        write_to_log_file(
            out=out,
            err=err,
            log_file=log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err,
        )
    except CalledProcessError as exc:
        logging.error(
                'reformat.sh failed for sample: %s with error: %s',
                sample.name, str(exc)
            )
        logging.error(traceback.format_exc())
        raise


def handle_validation_error(log_file: str, sample: CustomBox) -> None:
    """
    Handle validation errors by attempting to repair the reads.

    Args:
        log_file (str): Path to the log file.
        sample (CustomBox): Metadata sample object.
    """
    # Log a warning about detected errors in the FASTQ files
    logging.warning(
        'Errors detected in FASTQ files for sample %s. '
        'Please check the following files for details %s %s %s. '
        'The pipeline will use reformat.sh to attempt to repair issues',
        sample.name, log_file, sample.general.log_out,
        sample.general.log_err
    )

    # Get the paths for reformatted and repaired files
    reformatted_forward, reformatted_reverse, repair_forward, \
        repair_reverse = get_reformatted_and_repair_paths(sample)
    logging.debug(
        'Reformatted and repair paths for sample %s: %s, %s, %s, %s',
        sample.name, reformatted_forward, reformatted_reverse,
        repair_forward, repair_reverse
    )

    # If the reformatted forward file already exists, return early
    if os.path.isfile(reformatted_forward):
        logging.debug(
            'Reformatted forward file already exists for sample %s: %s',
            sample.name, reformatted_forward
        )
        return

    try:
        # Attempt to reformat the reads
        logging.debug(
            'Attempting to reformat reads for sample %s',
            sample.name)
        reformat_reads(
            log_file=log_file,
            sample=sample,
            reformatted_forward=reformatted_forward
        )
        logging.debug('Reformatting completed for sample %s', sample.name)

        # Attempt to repair the reads
        logging.debug('Attempting to repair reads for sample %s', sample.name)
        repair_reads(
            log_file=log_file,
            sample=sample,
            reformatted_forward=reformatted_forward,
            reformatted_reverse=reformatted_reverse,
            repair_forward=repair_forward,
            repair_reverse=repair_reverse
        )
        logging.debug('Repairing completed for sample %s', sample.name)
    except CalledProcessError as exc:
        # Handle any errors that occur during the repair process
        logging.error(
            'Repair process failed for sample %s with error: %s',
            sample.name, str(exc)
        )
        logging.error(traceback.format_exc())
        handle_repair_error(
            log_file=log_file,
            sample=sample,
            reformatted_forward=reformatted_forward,
            reformatted_reverse=reformatted_reverse,
            repair_forward=repair_forward,
            repair_reverse=repair_reverse
        )


def get_reformatted_and_repair_paths(sample: CustomBox) -> tuple:
    """
    Get the paths for reformatted and repaired files.

    Args:
        sample (CustomBox): Metadata sample object.

    Returns:
        tuple: Paths for reformatted and repaired files.
    """
    logging.debug(
        'Getting reformatted and repair paths for sample: %s',
        sample.name)

    # Define the path for the reformatted forward file
    reformatted_forward = os.path.join(
        sample.general.output_directory,
        f'{sample.name}_reformatted_R1.fastq.gz'
    )
    logging.debug('Reformatted forward path: %s', reformatted_forward)

    # Define the path for the repaired forward file
    repair_forward = os.path.join(
        sample.general.output_directory,
        f'{sample.name}_repaired_R1.fastq.gz'
    )
    logging.debug('Repair forward path: %s', repair_forward)

    # Check if there are two FASTQ files
    if len(sample.general.fastq_files) == 2:
        # Define the paths for the reformatted and repaired reverse files
        reformatted_reverse = os.path.join(
            sample.general.output_directory,
            f'{sample.name}_reformatted_R2.fastq.gz'
        )
        repair_reverse = os.path.join(
            sample.general.output_directory,
            f'{sample.name}_repaired_R2.fastq.gz'
        )
        logging.debug('Reformatted reverse path: %s', reformatted_reverse)
        logging.debug('Repair reverse path: %s', repair_reverse)
    else:
        # If there is only one FASTQ file, set the reverse paths to empty
        # strings
        reformatted_reverse = ''
        repair_reverse = ''
        logging.debug('Single FASTQ file detected, no reverse paths needed.')

    # Return the paths for the reformatted and repaired files
    return reformatted_forward, reformatted_reverse, repair_forward, \
        repair_reverse


def reformat_reads(
    log_file: str, sample: CustomBox, reformatted_forward: str
) -> None:
    """
    Reformat the reads using reformat.sh.

    Args:
        log_file (str): Path to the log file.
        sample (CustomBox): Metadata sample object.
        reformatted_forward (str): Path to the reformatted forward file.

    Raises:
        CalledProcessError: If reformat.sh fails.
    """
    logging.debug('Reformatting reads for sample: %s', sample.name)

    try:
        # Run reformat.sh to reformat the reads
        logging.debug(
            'Running reformat.sh for sample: %s, output: %s',
            sample.name, reformatted_forward
        )
        out, err, _ = bbtools.reformat_reads(
            forward_in=sample.general.fastq_files[0],
            forward_out=reformatted_forward,
            returncmd=True
        )
        logging.debug('reformat.sh completed for sample: %s', sample.name)

        # Write the output and error messages to the log file
        logging.debug(
            'Writing reformat.sh output to log file for sample: %s',
            sample.name
        )
        write_to_log_file(
            out=out,
            err=err,
            log_file=log_file,
            sample_log=sample.general.log_out,
            sample_err=sample.general.log_err,
        )
    except CalledProcessError as exc:
        logging.error(
            'reformat.sh failed for sample: %s with error: %s',
            sample.name, str(exc)
        )
        logging.error(traceback.format_exc())
        raise


def repair_reads(
    log_file: str, sample: CustomBox, reformatted_forward: str,
    reformatted_reverse: str, repair_forward: str, repair_reverse: str
) -> None:
    """
    Repair the reads using repair.sh.

    Args:
        log_file (str): Path to the log file.
        sample (CustomBox): Metadata sample object.
        reformatted_forward (str): Path to the reformatted forward file.
        reformatted_reverse (str): Path to the reformatted reverse file.
        repair_forward (str): Path to the repaired forward file.
        repair_reverse (str): Path to the repaired reverse file.

    Raises:
        CalledProcessError: If repair.sh fails.
    """
    logging.debug('Repairing reads for sample: %s', sample.name)

    try:
        # Check if there is a reformatted reverse file
        if reformatted_reverse:
            logging.debug(
                'Running repair.sh for sample: %s, forward: %s, reverse: %s',
                sample.name, reformatted_forward, reformatted_reverse
            )
            # Run repair.sh to repair the reads
            out, err, _ = bbtools.repair_reads(
                forward_in=reformatted_forward,
                reverse_in=reformatted_reverse,
                forward_out=repair_forward,
                reverse_out=repair_reverse,
                returncmd=True
            )
            logging.debug('repair.sh completed for sample: %s', sample.name)

            # Write the output and error messages to the log file
            logging.debug(
                'Writing repair.sh output to log file for sample: %s',
                sample.name
            )
            write_to_log_file(
                out=out,
                err=err,
                log_file=log_file,
                sample_log=sample.general.log_out,
                sample_err=sample.general.log_err
            )

        # Check if the reformatted forward file exists
        if os.path.isfile(reformatted_forward):
            logging.debug(
                'Updating fastq_files attribute for sample: %s', sample.name
            )
            # Update the fastq_files attribute to point to the repaired files
            sample.general.fastq_files = (
                [repair_forward, repair_reverse]
                if repair_reverse else [reformatted_forward]
            )
    except CalledProcessError as exc:
        logging.error(
            'repair.sh failed for sample: %s with error: %s',
            sample.name, str(exc)
        )
        logging.error(traceback.format_exc())
        raise


def handle_repair_error(
    log_file: str, sample: CustomBox, reformatted_forward: str,
    reformatted_reverse: str, repair_forward: str, repair_reverse: str
) -> None:
    """
    Handle errors during the repair process.

    Args:
        log_file (str): Path to the log file.
        sample (CustomBox): Metadata sample object.
        reformatted_forward (str): Path to the reformatted forward file.
        reformatted_reverse (str): Path to the reformatted reverse file.
        repair_forward (str): Path to the repaired forward file.
        repair_reverse (str): Path to the repaired reverse file.
    """
    logging.debug('Handling repair error for sample: %s', sample.name)

    # Check if both reformatted forward and reverse files exist
    if os.path.isfile(reformatted_forward) and os.path.isfile(
            reformatted_reverse):
        try:
            logging.debug(
                'Running repair.sh for sample: %s, forward: %s, reverse: %s',
                sample.name, reformatted_forward, reformatted_reverse
            )
            # Run repair.sh to repair the reads
            out, err, _ = bbtools.repair_reads(
                forward_in=reformatted_forward,
                reverse_in=reformatted_reverse,
                forward_out=repair_forward,
                reverse_out=repair_reverse,
                returncmd=True
            )
            logging.debug('repair.sh completed for sample: %s', sample.name)

            # Write the output and error messages to the log file
            logging.debug(
                'Writing repair.sh output to log file for sample: %s',
                sample.name
            )
            write_to_log_file(
                out=out,
                err=err,
                log_file=log_file,
                sample_log=sample.general.log_out,
                sample_err=sample.general.log_err
            )

            # Update the fastq_files attribute to point to the repaired files
            logging.debug(
                'Updating fastq_files attribute for sample: %s', sample.name
            )
            sample.general.fastq_files = (
                [repair_forward, repair_reverse]
                if repair_reverse else [repair_forward]
            )
        except CalledProcessError as exc:
            # Log an error if the repair process fails
            logging.error(
                'repair.sh failed for sample: %s with error: %s',
                sample.name, str(exc)
            )
            logging.error(traceback.format_exc())
            log_fastq_error(log_file, sample)
            error(sample, 'fastq_error')
    else:
        # Log an error if the reformatted files do not exist
        logging.error(
            'Reformatted files do not exist for sample: %s', sample.name
        )
        log_fastq_error(log_file, sample)
        error(sample, 'fastq_error')


def log_fastq_error(log_file: str, sample: CustomBox) -> None:
    """
    Log an error message for FASTQ file issues.

    Args:
        log_file (str): Path to the log file.
        sample (CustomBox): Metadata sample object.
    """
    logging.debug('Logging FASTQ error for sample: %s', sample.name)

    # Create an error message indicating issues with the FASTQ files
    message = (
        f'An error was detected in the FASTQ files for sample {sample.name}. '
        'These files will not be processed further'
    )
    logging.error(message)

    # Write the error message to the log file and sample-specific logs
    logging.debug(
        'Writing FASTQ error message to log file for sample: %s', sample.name
    )
    write_to_log_file(
        out=message,
        err=message,
        log_file=log_file,
        sample_log=sample.general.log_out,
        sample_err=sample.general.log_err
    )


def error(sample: CustomBox, message: str) -> None:
    """
    Check to see if the run CustomBox exists. If so, update the
    run.status to reflect the error.

    Args:
        sample (CustomBox): Metadata sample object.
        message (str): Error message to add to the sample.run.status attribute.
    """
    logging.debug('Setting error status for sample: %s', sample.name)

    # Set the .fastq_files attribute to an empty list
    sample.general.fastq_files = []
    logging.debug('Cleared fastq_files for sample: %s', sample.name)

    # Ensure that the run attribute exists and is a CustomBox
    if not hasattr(sample, 'run'):
        sample.run = CustomBox()
        logging.debug('Created run attribute for sample: %s', sample.name)

    # Set the status attribute in the run CustomBox
    sample.run.status = message
    logging.debug(
        'Set run.status to "%s" for sample: %s', message, sample.name
    )
