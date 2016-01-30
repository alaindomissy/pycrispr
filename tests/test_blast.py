import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from crispr.blast import *


def test_grouper_longest():
    assert(
        list(grouper_longest(range(10),3,99)) == [(0, 1, 2), (3, 4, 5), (6, 7, 8), (9, 99, 99)]
    )
