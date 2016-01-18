import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

import Bio


import Bio.Restriction
from buffet.analysis import *


def test_create_enzyme():
    enz = create_enzyme('BfaI')
    assert(str(enz)) == 'BfaI'
    assert(type(enz)) == Bio.Restriction.Restriction.RestrictionType


def test_create_batch():
    batch1 = create_enzyme("ScrFI") + create_enzyme("HpaII") + create_enzyme("BfaI")
    batch2 = create_batch(['BfaI', 'HpaII', 'ScrFI'])
    batch3 = Bio.Restriction.RestrictionBatch(['BfaI', 'HpaII', 'ScrFI'])
