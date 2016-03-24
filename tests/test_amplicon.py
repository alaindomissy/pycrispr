from __future__ import absolute_import, division, print_function
import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))


from Bio.SeqRecord import SeqRecord

from crispr.amplicon import *


def test_amplicon():
    run = []
    guidelist = []
    # _ = Amplicon
    # a1 = Amplicon(run, guidelist, 'hg38', 'chr21')
    assert(True)
