import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from buffet.cut import cut_fastafile, cut_unicodestring





# 42nt long section of of Scaffold02974 around the first cut from enzyme BfaI
#############################################################################
# with recognition site C^TA_G at one-based coords: 21,22,23,24
# GCGCTGGCCAGAACGTTCTC^TA_GGAATCGTGGAGAAGACATT
def test_cut_fastafile_42():

    assert (
        cut_fastafile('tests/data/fourtytwobp.fasta') == [
            'fourtytwobps\t0\t20\tBfaI\t1000\t+\n',
            'fourtytwobps\t22\t42\tBfaI\t1000\t-\n']
    )


def test_cut_unicodestring_42():
    unicodestring = u'''>fourtytwobps
GCGCTGGCCAGAACGTTCTCTAGGAATCGTGGAGAAGACATT
'''
    expected_list_of_bedtules = ['fourtytwobps\t0\t20\tBfaI\t1000\t+\n', 'fourtytwobps\t22\t42\tBfaI\t1000\t-\n']
    assert(cut_unicodestring(unicodestring)==expected_list_of_bedtules)




# From crisper-eating paper, we should get:
#
# [seqio.SeqRecord(seq=Seq('GCGCTGGCCAGAACGTTCTC', SingleLetterAlphabet()),
#                  id='1_F', name='788', description='BfaI', dbxrefs=['Scaffold102974']),
#  seqio.SeqRecord(seq=Seq('AATGTCTTCTCCACGATTCC', SingleLetterAlphabet()),
#                  id='2_R', name='810', description='BfaI', dbxrefs=['Scaffold102974']),
#  seqio.SeqRecord(seq=Seq('TAAAGGAGAAGGAAACCCCC', SingleLetterAlphabet()),
#                  id='2_F', name='1444', description='BfaI', dbxrefs=['Scaffold102974']),
#  seqio.SeqRecord(seq=Seq('AGGGGAGGGATTTTTGTGCC', SingleLetterAlphabet()),
#                  id='3_R', name='1465', description='BfaI', dbxrefs=['Scaffold102974']),
#  seqio.SeqRecord(seq=Seq('GGAGGGACAGCAGCTGGGCC', SingleLetterAlphabet()),
#                  id='4_F', name='734', description='ScrFI', dbxrefs=['Scaffold102974']),
#  seqio.SeqRecord(seq=Seq('ATCTTTCTTGCCCCCCCCCC', SingleLetterAlphabet()),
#                  id='5_R', name='754', description='ScrFI', dbxrefs=['Scaffold102974'])
#  ]
def test_cut_fastafile_1500():
    cuts1500 = cut_fastafile('tests/data/Scaffold102974:1-1500_1500.fasta')
    expected_list_of_bedtules = ['Scaffold102974:1-1500()\t786\t806\tBfaI\t1000\t+\n',
                                 'Scaffold102974:1-1500()\t808\t828\tBfaI\t1000\t-\n',
                                 'Scaffold102974:1-1500()\t1442\t1462\tBfaI\t1000\t+\n',
                                 'Scaffold102974:1-1500()\t1464\t1484\tBfaI\t1000\t-\n',
                                 'Scaffold102974:1-1500()\t732\t752\tScrFI\t1000\t+\n',
                                 'Scaffold102974:1-1500()\t753\t773\tScrFI\t1000\t-\n']
    assert(cuts1500==expected_list_of_bedtules)


def test_cut_unicodestring_1500():
    unicodestring = u'''>Scaffold102974:1-1500()
TCTCAACAAACCCATAAGTCACTGGTACATTAACAGATAACGCTAATAATAAAATGGCAGAATGACTCAGATACAATTCCCATATAAAAAGCCATTAATA
TGAACAGAACTTACTCGTAAAAATCCAGTGAAAATGAACAGGGTGATGATTAGTGACCTCATGGTGGCATTCAGTCACTAATCTGACAATCCCAACTGTA
AACAGTCCTTATACCAGCTTGTGTCAGGGGCACCATTATGAGGTCATCACTAACAATCATGGCGACCATTGTGATGTCGCTACTAACAATCATGTGACAT
CAACTGCTACAGGCCAAGTCAGCTTTCTGTTACGGCAAATATAGAATGACTGTTAGGGAGACAATATAAGAGCACCTTACAGTAATAAAACTTTATACAT
AAATTCTACTTCATAAATTCATATTTTTTTATTGAAAATCACGATTTATTTAGATTACATTCAACAGCTAATCTTAAATTGACTCAACCTTAGTTATTCT
TATTGCCTGTAATCAATAACAGGTTTTGGGTTGGACATAGTTACATATATTACTTTAAATTCCAAAAGGACAAAAGTCTGTCAAGTTCAACCCCTCCAAT
TTATCCCCAGCATATGTGTGTACATATATATATACTGATGCACACTGGACGTCCACAGAAATGTTGCTACCTATGTCGGATCACAGATAAGGAGATTATT
TAATGCAATTTAATAAAGTGCAACAGATGGGGGGAGGGACAGCAGCTGGGCCTGGGGGGGGGGCAAGAAAGATAAATTTGGTCGTTGCGCTGGCCAGAAC
GTTCTCTAGGAATCGTGGAGAAGACATTGGCCAAAAGAGGAAGGTCTCACCAATGTGCCTTGATGAATAGAACAGTTAAAATATTTCAGCTGTAGCCCTC
CAGCTCCAGAACTATAATTACCAGAATCCACCTCAGTTGGAGGCTGAAGATGCCTTAATTCCAAAATGTATTTATTCTTATGCCTTCTTTTCTTATTTAG
TTTCAGTTTTTTTCTGAACTTTTCTTCTATTAATTTTTACCCTTCCCCAACCAACCCCCATCTTTCCTCACATTTCTCAGTCCCTTATGACCCCCCTTCA
TCTCTTTTTCTGCTCCTTGTTCCTACTCATCTTCTTCTCTGTTCTGCATATTTTGCTCTCTCTCTATGTGTCATGGTGGCTTCAGCTTTTTTCTGAATAG
AGCTCAGCTGAGTAGGGATATTGATTGTGATTGGGAGGTGAACTTCAACCTTAAATCTCCTATTTAGTGCCATGTGAGTGCAGCAATACAACTGTTGGGA
AGGGACTTTGATGGGAGGAATTTGGATGGAAGGGTCTAAAAGAGAGACAGGCCAAGGAAAGGAATCATCTGATTGGCTTGTTTTGATTACTGGCGTAAGA
CAGTTAGAATGCTGAAAGCAACTACTACTGGTCTGTGCTGCTTAAAGGAGAAGGAAACCCCCTAGGCACAAAAATCCCTCCCCTCTCCCCTGTGTTGTC
'''
    expected_list_of_bedtules = ['Scaffold102974:1-1500()\t786\t806\tBfaI\t1000\t+\n',
         'Scaffold102974:1-1500()\t808\t828\tBfaI\t1000\t-\n',
          'Scaffold102974:1-1500()\t1442\t1462\tBfaI\t1000\t+\n',
          'Scaffold102974:1-1500()\t1464\t1484\tBfaI\t1000\t-\n',
          'Scaffold102974:1-1500()\t732\t752\tScrFI\t1000\t+\n',
          'Scaffold102974:1-1500()\t753\t773\tScrFI\t1000\t-\n']
    assert(cut_unicodestring(unicodestring)==expected_list_of_bedtules)

