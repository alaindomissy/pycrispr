import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

import Bio
from Bio.Seq import Seq

import Bio.Restriction
from crispr.analyse import create_enzyme, create_batch, create_analysis, analyse


def test_create_enzyme():
    enz = create_enzyme('BfaI')
    assert(str(enz)) == 'BfaI'
    assert(type(enz)) == Bio.Restriction.Restriction.RestrictionType


def test_create_batch():
    batch1 = create_enzyme("ScrFI") + create_enzyme("HpaII") + create_enzyme("BfaI")
    batch2 = create_batch(['BfaI', 'HpaII', 'ScrFI'])
    batch3 = Bio.Restriction.RestrictionBatch(['BfaI', 'HpaII', 'ScrFI'])



# BfAI
# 5'...C TA G ...3'
# 3'...G AT C ...5'

# HpaII
# 5'...C CG G ...3'
# 3'...G GC C ...5'

# ScrFI
# 5'...CC N GG ...3'
# 3'...GG N CC ...5'




# 42nt long section of of Scaffold02974 around the first cut from enzyme BfaI
#############################################################################
# with recognition site C^TA_G at one-based coords: 20,21,22,23
# 01 GCGCTGGCCA
# 11 GAACGTTCTC^TA_GGAATCGTGG
# 31 AGAAGACATT

#
# def test_analyse_42():
#     bioseq1 = Seq(u'''GCGCTGGCCAGAACGTTCTCTAGGAATCGTGGAGAAGACATT''')
#     expected_dict_default = '{BfaI: [21], HpaII: [], ScrFI: []}'
#     expected_dict_BfaI = '{BfaI: [21]}'
#     expected_dict_HpaII = '{HpaII: []}'
#     expected_dict_ScrFI = '{ScrFI: []}'
#     assert(str(analyse(bioseq1)) == expected_dict_default)
#     assert(str(analyse(bioseq1,['BfaI', 'HpaII', 'ScrFI'])) == expected_dict_default)
#     assert(str(analyse(bioseq1,['BfaI'])) == expected_dict_BfaI)
#     assert(str(analyse(bioseq1,['HpaII'])) == expected_dict_HpaII)
#     assert(str(analyse(bioseq1,['ScrFI'])) == expected_dict_ScrFI)
#



# FIRST 1500bps from xenopus-laevis scaffold used in crispr-eating paper
########################################################################
# 0001 TCTCAACAAACCCATAAGTCACTGGTACATTAACAGATAACGCTAATAATAAAATGGCAGAATGACTCAGATACAATTCCCATATAAAAAGCCATTAATA
# 0101 TGAACAGAACTTACTCGTAAAAATCCAGTGAAAATGAACAGGGTGATGATTAGTGACCTCATGGTGGCATTCAGTCACTAATCTGACAATCCCAACTGTA
# 0201 AACAGTCCTTATACCAGCTTGTGTCAGGGGCACCATTATGAGGTCATCACTAACAATCATGGCGACCATTGTGATGTCGCTACTAACAATCATGTGACAT
# 0301 CAACTGCTACAGGCCAAGTCAGCTTTCTGTTACGGCAAATATAGAATGACTGTTAGGGAGACAATATAAGAGCACCTTACAGTAATAAAACTTTATACAT
# 0401 AAATTCTACTTCATAAATTCATATTTTTTTATTGAAAATCACGATTTATTTAGATTACATTCAACAGCTAATCTTAAATTGACTCAACCTTAGTTATTCT
# 0501 TATTGCCTGTAATCAATAACAGGTTTTGGGTTGGACATAGTTACATATATTACTTTAAATTCCAAAAGGACAAAAGTCTGTCAAGTTCAACCCCTCCAAT
# 0601TTATCCCCAGCATATGTGTGTACATATATATATACTGATGCACACTGGACGTCCACAGAAATGTTGCTACCTATGTCGGATCACAGATAAGGAGATTATT
# 0701 TAATGCAATTTAATAAAGTGCAACAGATGGGGGGAGGGACAGCAGCTGGG  CC T GG  GGGGGGGGCAAGAAAGATAAATTTGGTCGTTGCGCTGGCCAGAAC
# 0801 GTTCT  C TA G  GAATCGTGGAGAAGACATTGGCCAAAAGAGGAAGGTCTCACCAATGTGCCTTGATGAATAGAACAGTTAAAATATTTCAGCTGTAGCCCTC
# 0901 CAGCTCCAGAACTATAATTACCAGAATCCACCTCAGTTGGAGGCTGAAGATGCCTTAATTCCAAAATGTATTTATTCTTATGCCTTCTTTTCTTATTTAG
# 1001 TTTCAGTTTTTTTCTGAACTTTTCTTCTATTAATTTTTACCCTTCCCCAACCAACCCCCATCTTTCCTCACATTTCTCAGTCCCTTATGACCCCCCTTCA
# 1101 TCTCTTTTTCTGCTCCTTGTTCCTACTCATCTTCTTCTCTGTTCTGCATATTTTGCTCTCTCTCTATGTGTCATGGTGGCTTCAGCTTTTTTCTGAATAG
# 1201 AGCTCAGCTGAGTAGGGATATTGATTGTGATTGGGAGGTGAACTTCAACCTTAAATCTCCTATTTAGTGCCATGTGAGTGCAGCAATACAACTGTTGGGA
# 1301 AGGGACTTTGATGGGAGGAATTTGGATGGAAGGGTCTAAAAGAGAGACAGGCCAAGGAAAGGAATCATCTGATTGGCTTGTTTTGATTACTGGCGTAAGA
# 1401 CAGTTAGAATGCTGAAAGCAACTACTACTGGTCTGTGCTGCTTAAAGGAGAAGGAAACCCC  C TA G  GCACAAAAATCCCTCCCCTCTCCCCTGTGTTGTC


# def test_analyse_1500():
#     bioseq1 = Bio.Seq.Seq(u'''TCTCAACAAACCCATAAGTCACTGGTACATTAACAGATAACGCTAATAATAAAATGGCAGAATGACTCAGATACAATTCCCATATAAAAAGCCATTAATA
# TGAACAGAACTTACTCGTAAAAATCCAGTGAAAATGAACAGGGTGATGATTAGTGACCTCATGGTGGCATTCAGTCACTAATCTGACAATCCCAACTGTA
# AACAGTCCTTATACCAGCTTGTGTCAGGGGCACCATTATGAGGTCATCACTAACAATCATGGCGACCATTGTGATGTCGCTACTAACAATCATGTGACAT
# CAACTGCTACAGGCCAAGTCAGCTTTCTGTTACGGCAAATATAGAATGACTGTTAGGGAGACAATATAAGAGCACCTTACAGTAATAAAACTTTATACAT
# AAATTCTACTTCATAAATTCATATTTTTTTATTGAAAATCACGATTTATTTAGATTACATTCAACAGCTAATCTTAAATTGACTCAACCTTAGTTATTCT
# TATTGCCTGTAATCAATAACAGGTTTTGGGTTGGACATAGTTACATATATTACTTTAAATTCCAAAAGGACAAAAGTCTGTCAAGTTCAACCCCTCCAAT
# TTATCCCCAGCATATGTGTGTACATATATATATACTGATGCACACTGGACGTCCACAGAAATGTTGCTACCTATGTCGGATCACAGATAAGGAGATTATT
# TAATGCAATTTAATAAAGTGCAACAGATGGGGGGAGGGACAGCAGCTGGGCCTGGGGGGGGGGCAAGAAAGATAAATTTGGTCGTTGCGCTGGCCAGAAC
# GTTCTCTAGGAATCGTGGAGAAGACATTGGCCAAAAGAGGAAGGTCTCACCAATGTGCCTTGATGAATAGAACAGTTAAAATATTTCAGCTGTAGCCCTC
# CAGCTCCAGAACTATAATTACCAGAATCCACCTCAGTTGGAGGCTGAAGATGCCTTAATTCCAAAATGTATTTATTCTTATGCCTTCTTTTCTTATTTAG
# TTTCAGTTTTTTTCTGAACTTTTCTTCTATTAATTTTTACCCTTCCCCAACCAACCCCCATCTTTCCTCACATTTCTCAGTCCCTTATGACCCCCCTTCA
# TCTCTTTTTCTGCTCCTTGTTCCTACTCATCTTCTTCTCTGTTCTGCATATTTTGCTCTCTCTCTATGTGTCATGGTGGCTTCAGCTTTTTTCTGAATAG
# AGCTCAGCTGAGTAGGGATATTGATTGTGATTGGGAGGTGAACTTCAACCTTAAATCTCCTATTTAGTGCCATGTGAGTGCAGCAATACAACTGTTGGGA
# AGGGACTTTGATGGGAGGAATTTGGATGGAAGGGTCTAAAAGAGAGACAGGCCAAGGAAAGGAATCATCTGATTGGCTTGTTTTGATTACTGGCGTAAGA
# CAGTTAGAATGCTGAAAGCAACTACTACTGGTCTGTGCTGCTTAAAGGAGAAGGAAACCCCCTAGGCACAAAAATCCCTCCCCTCTCCCCTGTGTTGTC
# ''')
#     expected_dict_default = '{BfaI: [807, 1463], HpaII: [], ScrFI: [753]}'
#     expected_dict_BfaI = '{BfaI: [807, 1463]}'
#     expected_dict_HpaII = '{HpaII: []}'
#     expected_dict_ScrFI = '{ ScrFI: [753]}'
#     assert(str(analyse(bioseq1)) == expected_dict_default)
#     assert(str(analyse(bioseq1,['BfaI', 'HpaII', 'ScrFI'])) == expected_dict_default)
#     assert(str(analyse(bioseq1,['BfaI'])) == expected_dict_BfaI)
#     assert(str(analyse(bioseq1,['HpaII'])) == expected_dict_HpaII)
#     assert(str(analyse(bioseq1,['ScrFI'])) == expected_dict_ScrFI)




# problem with print or str function on the anlysis object ? similar to other Bio.Restriction bug already fixed ?
#
# ________________________________________________________________________________________________ test_analyse_42 _________________________________________________________________________________________________
#
#     def test_analyse_42():
#         bioseq1 = Seq(u'''GCGCTGGCCAGAACGTTCTCTAGGAATCGTGGAGAAGACATT''')
#         expected_dict_default = '{BfaI: [21], HpaII: [], ScrFI: []}'
#         expected_dict_BfaI = '{BfaI: [21]}'
#         expected_dict_HpaII = '{HpaII: []}'
#         expected_dict_ScrFI = '{ScrFI: []}'
#         assert(str(analyse(bioseq1)) == expected_dict_default)
#         assert(str(analyse(bioseq1,['BfaI', 'HpaII', 'ScrFI'])) == expected_dict_default)
# >       assert(str(analyse(bioseq1,['BfaI'])) == expected_dict_BfaI)
#
# tests/test_analyse.py:55:
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
# crispr/analyse.py:77: in analyse
#     return create_analysis(seq, enzyme_names).full()     # .only_between(20, -20)  would get rid of both sides, not good
# crispr/analyse.py:62: in create_analysis
#     return Analysis(batch, seq, linear= True)
# /root/anaconda3/envs/pycrispr/lib/python2.7/site-packages/Bio/Restriction/Restriction.py:2092: in __init__
#     RestrictionBatch.__init__(self, restrictionbatch)
# _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _ _
#
# self = <[AttributeError("'Analysis' object has no attribute 'rb'") raised in repr()] SafeRepr object at 0x7f5d849cc878>, first = BfaI, suppliers = []
#
#     def __init__(self, first=[], suppliers=[]):
#         """RestrictionBatch([sequence]) -> new RestrictionBatch."""
# >       first = [self.format(x) for x in first]
# E       TypeError: 'RestrictionType' object is not iterable
#
# /root/anaconda3/envs/pycrispr/lib/python2.7/site-packages/Bio/Restriction/Restriction.py:1857: TypeError

