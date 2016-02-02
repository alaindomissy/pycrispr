import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from crispr.bedtuple import *


def test_bedtuple():

    bt1 = Bedtuple('chr1',1000,1100)
    assert(Bedtuple(scaffold='chr1', beg=1000, end=1100))
    assert(bt1.length == 100)
    assert(bt1.beg1 == '1001')
    assert(bt1.filename == 'chr1_1000-1100_100')

    # has_dash
    assert(Bedtuple.from_coord('chr6:136640001-136680000').bedlist == [('chr6', '136640000', '136680000')])
    assert(Bedtuple.from_coord('chr6:136640001-136680000').filename == 'chr6_136640000-136680000_40000')

    # has_under
    assert(Bedtuple.from_coord('chr6:136640001_40000').bedlist  == [('chr6', '136640000', '136680000')])
    assert(Bedtuple.from_coord('chr6:136640001_40000').filename == 'chr6_136640000-136680000_40000')
