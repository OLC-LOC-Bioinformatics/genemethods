#!/usr/bin/env python3
from olctools.accessoryFunctions.accessoryFunctions import make_path, SetupLogging
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio import SeqIO
from argparse import ArgumentParser
from glob import glob
import logging
import os


class SequenceExtractor(object):

    def main(self):
        self.parse_details()
        self.sequence_extract()

    def parse_details(self):
        """
        Parse sequence extraction details from the supplied file
        """
        logging.info('Parsing details from file')
        with open(self.fastadetails, 'r') as fastadetails:
            for line in fastadetails:
                # Split the line on semi-colons, and remove newline characters
                seqid, contig, start, stop = line.rstrip().split(';')
                logging.debug(f'SEQID: {seqid}, contig: {contig}, start: {start}, stop: {stop}')
                # Initialise the list to store details for each supplied SEQID
                if seqid not in self.details_dict:
                    self.details_dict[seqid] = list()
                # Ensure that the min and max are in order, and subtract one from the start due to zero-indexing
                start_int = min([int(start), int(stop)]) - 1
                # If the start position is provided as 0, ensure that it stays at 0, rather than being -1
                if int(start) != 0:
                    logging.debug(f'The provided start position, {start} is negative, and has been adjusted to 1')
                start_int = start_int if start_int >= 0 else 0
                stop_int = max([int(start), int(stop)])
                # Initialise a dictionary to store the individual details
                info_dict = {
                    'contig': contig,
                    'start': start_int,
                    'stop': stop_int
                }
                # Update the list with the populated dictionary
                self.details_dict[seqid].append(info_dict)

    def sequence_extract(self):
        """
        Extract the relevant sequence(s) from the FASTA files. Write the extracted sequences with appropriate headers
        to an output file
        """
        logging.info('Extracting sequences')
        # Create the output path as required
        output_path = os.path.join(self.sequencepath, 'output')
        make_path(output_path)
        # Iterate through all the supplied SEQIDs in the FASTA details file
        for seqid, fasta_details_list in self.details_dict.items():
            # Try to find the file that matches the supplied SEQID
            try:
                fasta_file = glob(os.path.join(self.sequencepath, f'*{seqid}*.*'))[0]
            except IndexError:
                logging.warning(f'Could not locate a file containing the string \"{seqid}\" in {self.sequencepath}')
                continue
            for fasta_details in fasta_details_list:
                # Boolean to track whether a contig is present in the file
                record_present = False
                # Iterate through all the records in the FASTA file
                for record in SeqIO.parse(fasta_file, 'fasta'):
                    if record.id != fasta_details['contig']:
                        continue
                    record_present = True
                    # Ensure that both the specified start and stop positions are within the contig
                    if fasta_details['start'] not in range(len(record.seq)) and fasta_details['stop'] not in \
                            range(len(record.seq)):
                        logging.warning(f'Could not find the specified start position: {fasta_details["start"]} or '
                                        f'the specified stop position: {fasta_details["stop"]} in '
                                        f'contig {fasta_details["contig"]} in file {fasta_file}')
                        continue
                    # Only the start position is invalid
                    elif fasta_details['start'] not in range(len(record.seq)):
                        logging.warning(f'Could not find the specified start position: {fasta_details["start"]} in '
                                        f'contig {fasta_details["contig"]} in file {fasta_file}')
                        continue
                    # Only the stop position is invalid
                    elif fasta_details['stop'] not in range(len(record.seq)):
                        logging.warning(f'Could not find the specified stop position: {fasta_details["stop"]} in '
                                        f'contig {fasta_details["contig"]} in file {fasta_file}')
                        continue
                    # Modify the fasta details dictionary to extract the entire contig if start and stop are 0
                    if fasta_details['start'] == 0 and fasta_details['stop'] == 0:
                        fasta_details['stop'] = len(str(record.seq))
                        logging.debug(f'Setting stop coordinate for SEQID: {seqid}, contig: {fasta_details["contig"]} '
                                      f'to fasta_details["stop"]')
                    # Create a new header for the extracted sequence. SEQID_contig_start(plus 1 to match the
                    # supplied position)_stop
                    newid = f'{seqid}_{fasta_details["contig"]}_{fasta_details["start"] + 1}_' \
                            f'{fasta_details["stop"]}'
                    # Create a record for the extracted sequence
                    newrecord = SeqRecord(seq=Seq(record.seq[fasta_details['start']:fasta_details['stop']]),
                                          id=newid,
                                          name=str(),
                                          description=str())
                    logging.debug(f'SEQID: {seqid}, contig: {fasta_details["contig"]}, '
                                  f'start: {fasta_details["start"] + 1}, stop: {fasta_details["stop"]}, '
                                  f'sequence: {record.seq[fasta_details["start"]:fasta_details["stop"]]}')
                    # Add the record to the list of all records
                    self.newrecords.append(newrecord)
                if not record_present:
                    logging.warning(f'Could not find {fasta_details["contig"]} in file {fasta_file}')
        # Write the records to the output file
        output_file = os.path.join(output_path, 'extracted_sequences.fasta')
        logging.info(f'Writing extracted sequences to {output_file}')
        with open(output_file, 'w') as output:
            SeqIO.write(self.newrecords, output, 'fasta')

    def __init__(self, fastadetails, sequencepath, debug=False):
        # Setup the logging
        try:
            SetupLogging(debug=debug,
                         filehandle=os.path.join(sequencepath, 'log'),
                         logfile_level=logging.WARNING)
        except FileNotFoundError:
            SetupLogging(debug=debug,
                         logfile_level=logging.WARNING)
        logging.info('Welcome to the SequenceExtractor')
        # Determine the path in which the sequence files are located. Allow for ~ expansion
        if sequencepath.startswith('~'):
            self.sequencepath = os.path.abspath(os.path.expanduser(os.path.join(sequencepath)))
        else:
            self.sequencepath = os.path.abspath(os.path.join(sequencepath))
        # Ensure that the path exists
        assert os.path.isdir(self.sequencepath), 'Invalid path specified: {path}'.format(path=self.sequencepath)
        logging.debug('Supplied sequence path: \n{path}'.format(path=self.sequencepath))
        # Determine the path in which the fasta details file is located. Allow for ~ expansion
        if fastadetails.startswith('~'):
            self.fastadetails = os.path.abspath(os.path.expanduser(os.path.join(fastadetails)))
        else:
            self.fastadetails = os.path.abspath(os.path.join(fastadetails))
        # Ensure that the file exists
        assert os.path.isfile(self.fastadetails), 'Invalid file specified: {file}'.format(file=self.fastadetails)
        logging.debug('Supplied fasta details file: \n{file}'.format(file=self.fastadetails))
        self.details_dict = dict()
        self.newrecords = list()


def main():
    parser = ArgumentParser(description='Extracts sequences from FASTA files from a given header name, as well as '
                                        'start and stop coordinates')
    parser.add_argument('-f', '--fastadetails',
                        required=True,
                        help='Name and path of the file containing the details of the sequence(s) to extract. The '
                             'format must be (one entry per line): '
                             'sequence file name(no extension);contig name;start coordinates; stop coordinates '
                             'e.g. 2019-SEQ-0848;Contig_1_149.079_Circ;1;50 '
                             'NOTE the semi-colon separated values. If you want the whole contig, specify the start '
                             'and stop positions to be 0')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Name and path to folder containing FASTA files from which sequences are to be extracted')
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        help='Enable this if you want debug level messages printed to the console')
    args = parser.parse_args()
    sequence_extractor = SequenceExtractor(fastadetails=args.fastadetails,
                                           sequencepath=args.sequencepath,
                                           debug=args.debug)
    sequence_extractor.main()
    logging.info('Sequence extraction complete')


if __name__ == '__main__':
    main()
