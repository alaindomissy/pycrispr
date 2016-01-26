import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from Bio.SeqRecord import SeqRecord

from buffet.prime import *


def test_mask_sequence():
    assert( mask_sequence(SeqRecord('acgtACGTacgtACGT')).seq == 'NNNNACGTNNNNACGT' )
