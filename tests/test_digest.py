import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

import buffet.digest as digest



def test_coord_to_bedtuple_filename():
    assert(digest.coord_to_bedtuple_filename('chr6:136640001-136680000') ==
        ([('chr6', '136640000', '136680000')], 'chr6+136640001-136680000_40000')
    )


def digesttest_stretch_to_bedtuple_filename():
    assert(stretch_to_bedtuple_filename('chr6:136640001_40000') ==
        ([('chr6', '136640000', '136680000')], 'chr6+136640001-136680000_40000')
    )

