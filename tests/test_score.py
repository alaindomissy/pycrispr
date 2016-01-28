import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from buffet.score import *


def test_is_valid_pam():
    assert(is_valid_pam('AGG'))
    assert(is_valid_pam('agg'))
    assert(is_valid_pam('GGG'))
    assert(not is_valid_pam('AGC'))
    assert(not is_valid_pam('ACG'))
    assert(not is_valid_pam('GGC'))
    assert(not is_valid_pam('GTT'))

def test_load_genome_dict():
    pass