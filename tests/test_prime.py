from __future__ import absolute_import, division, print_function
import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from Bio.SeqRecord import SeqRecord

from crispr.prime import *

def test_prime():
    assert(
        prime(["ampl1", "ampl2", "ampl3"])
        ==
        "Amplicon 0 : ampl1\nAmplicon 1 : ampl2\nAmplicon 2 : ampl3"
    )

