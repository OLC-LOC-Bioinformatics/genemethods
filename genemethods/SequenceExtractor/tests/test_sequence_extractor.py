#!/usr/bin/env python3
from genemethods.SequenceExtractor.src.sequenceExtractor import SequenceExtractor
import shutil
import logging
import pytest
import os

test_path = os.path.abspath(os.path.dirname(__file__))
sequencepath = os.path.join(test_path, 'sequences')
outputpath = os.path.join(sequencepath, 'output')
outputfile = os.path.join(outputpath, 'extracted_sequences.fasta')
detailspath = os.path.join(test_path, 'fastadetails')
logfile = os.path.join(sequencepath, 'log')


def test_foo(caplog):
    caplog.set_level(logging.INFO)
    pass


def test_invalid_details_path():
    with pytest.raises(AssertionError):
        SequenceExtractor(fastadetails=os.path.join(detailspath, 'not_a_file.txt'),
                          sequencepath=sequencepath)


def test_invalid_sequence_path():
    with pytest.raises(AssertionError):
        SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_vstop.txt'),
                          sequencepath='not_a_real_path')


def test_invalid_delimiter():
    with pytest.raises(ValueError):
        se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_vstop_idelimiter.txt'),
                               sequencepath=sequencepath)
        se.parse_details()


def test_invalid_number_of_details():
    with pytest.raises(ValueError):
        se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_mstop.txt'),
                               sequencepath=sequencepath)
        se.parse_details()


def test_invalid_start():
    with pytest.raises(ValueError):
        se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_istart_vstop.txt'),
                               sequencepath=sequencepath)
        se.parse_details()


def test_valid_basic_fasta_details():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_vstop.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    assert se.details_dict['2019-SEQ-0848'][0]['contig'] == 'Contig_1_149.079_Circ'


def test_valid_two_different_contig_fasta_details():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_vstop_2_separate_contigs.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    assert se.details_dict['2019-SEQ-0848'][0]['contig'] == 'Contig_1_149.079_Circ'
    assert se.details_dict['2019-SEQ-0848'][1]['start'] == 1


def test_valid_two_same_contig_fasta_details():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_vstop_2_same_contig.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    assert se.details_dict['2019-SEQ-0848'][0]['contig'] == 'Contig_1_149.079_Circ'
    assert se.details_dict['2019-SEQ-0848'][1]['stop'] == 77


def test_valid_two_same_contig_fasta_details_reversed():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath,
                                                     'vseqid_vcontig_vstart_vstop_2_same_contig_reversed.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    assert se.details_dict['2019-SEQ-0848'][0]['contig'] == 'Contig_1_149.079_Circ'
    assert se.details_dict['2019-SEQ-0848'][1]['stop'] == 50


def test_valid_two_seqids_contig_fasta_details():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_vstop_2_seqids.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    assert se.details_dict['2019-SEQ-0848'][0]['contig'] == 'Contig_1_149.079_Circ'
    assert se.details_dict['2019-SEQ-1019'][0]['stop'] == 77


def test_valid_two_seqids_multiple_contig_fasta_details():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_vstop_2_seqids_multiple.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    assert se.details_dict['2019-SEQ-0848'][0]['contig'] == 'Contig_1_149.079_Circ'
    assert se.details_dict['2019-SEQ-1019'][0]['stop'] == 77
    assert se.details_dict['2019-SEQ-0848'][1]['contig'] == 'Contig_2_392.879_Circ'


def test_invalid_seqid():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'iseqid_vcontig_vstart_vstop.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    se.sequence_extract()
    assert not se.newrecords


def test_invalid_contig():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_icontig_vstart_vstop.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    assert se.details_dict['2019-SEQ-0848'][0]['contig'] == 'Contig_11_149.079_Circ'
    se.sequence_extract()
    assert not se.newrecords


def test_negative_start():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_nstart_vstop.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    se.sequence_extract()
    assert not se.newrecords


def test_out_of_range_end():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_bstop.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    se.sequence_extract()
    assert not se.newrecords


def test_basic():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_vstop.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    se.sequence_extract()
    assert se.newrecords[0].seq == 'AAAAAAAACAAATATATACTTTGATGATAACTTTCTAAATATCTACAAAA'


def test_basic_reverse():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_rstart_rstop.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    se.sequence_extract()
    assert se.newrecords[0].seq == 'AAAAAAAACAAATATATACTTTGATGATAACTTTCTAAATATCTACAAAA'


def test_two_contigs():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_vstop_2_separate_contigs.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    se.sequence_extract()
    assert se.newrecords[0].seq == 'AAAAAAAACAAATATATACTTTGATGATAACTTTCTAAATATCTACAAAA'
    assert se.newrecords[1].seq == 'AAAAAAACAATAAAAAACACCGCAAAAATGGATTGTTA'


def test_extract_two_same_contigs():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_vstop_2_same_contig.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    se.sequence_extract()
    assert se.newrecords[0].seq == 'AAAAAAAACAAATATATACTTTGATGATAACTTTCTAAATATCTACAAAA'
    assert se.newrecords[1].id == '2019-SEQ-0848_Contig_1_149.079_Circ_14_77'


def test_extract_two_same_contigs_overlap():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath,
                                                     'vseqid_vcontig_vstart_vstop_2_same_contig_overlap.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    se.sequence_extract()
    assert se.newrecords[0].seq == 'ACTTTGATGATAACTTTCTAAATATCTACAAAA'
    assert se.newrecords[1].id == '2019-SEQ-0848_Contig_1_149.079_Circ_14_77'


def test_multiple_extractions():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath, 'vseqid_vcontig_vstart_vstop_2_seqids_multiple.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    se.sequence_extract()
    assert se.newrecords[0].seq == 'AAAAAAAACAAATATATACTTTGATGATAACTTTCTAAATATCTACAAAA'
    assert se.newrecords[1].id == '2019-SEQ-0848_Contig_2_392.879_Circ_2_39'
    assert se.newrecords[2].id == '2019-SEQ-1019_Contig_1_388.862_Circ_14_77'
    assert se.newrecords[3].id == '2019-SEQ-1019_Contig_3_52.4575_5_22'


def test_multiple_extractions_one_invalid_contig():
    se = SequenceExtractor(fastadetails=os.path.join(detailspath,
                                                     'vseqid_vcontig_vstart_vstop_2_seqids_multiple_invalid.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    se.sequence_extract()
    assert se.newrecords[0].seq == 'AAAAAAAACAAATATATACTTTGATGATAACTTTCTAAATATCTACAAAA'
    assert se.newrecords[1].id == '2019-SEQ-0848_Contig_2_392.879_Circ_2_39'
    assert se.newrecords[2].id == '2019-SEQ-0848_Contig_3_847.01_Circ_3_43'
    assert se.newrecords[3].id == '2019-SEQ-0848_Contig_1_149.079_Circ_19_47'
    assert se.newrecords[4].id == '2019-SEQ-1019_Contig_3_52.4575_5_22'


def test_remove_output():
    shutil.rmtree(outputpath)
    assert not os.path.isdir(outputpath)


def test_log(caplog):
    se = SequenceExtractor(fastadetails=os.path.join(detailspath,
                                                     'vseqid_vcontig_vstart_vstop_2_seqids_multiple_invalid.txt'),
                           sequencepath=sequencepath)
    se.parse_details()
    se.sequence_extract()
    assert 'Welcome to the SequenceExtractor' in caplog.records[0].message
    assert caplog.records[1].levelname == 'INFO'
    assert 'Could not find Contig_11_388.862_Circ' in caplog.records[3].message
