import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))


from Bio.SeqRecord import SeqRecord

from crispr.amplicon import Amplicon, mask_sequence


def test_mask_sequence():
    assert(mask_sequence(SeqRecord('acgtACGTacgtACGT')).seq == 'NNNNACGTNNNNACGT' )




def test_amplicon():
    run = []
    guidelist = []
    a = Amplicon
    # a1 = Amplicon(run, guidelist, 'hg38', 'chr21')
    assert(True)

