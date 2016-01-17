import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))


from Bio.Restriction import RestrictionBatch
from buffet.analysis import *


# def test_create_enzyme():
#     assert(str(create_enzyme("ScrFI") + create_enzyme("HpaII") + create_enzyme("BfaI")) ==
#            RestrictionBatch(['BfaI', 'HpaII', 'ScrFI']))