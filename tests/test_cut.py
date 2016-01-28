import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from buffet.cut import cut_file



def test_cut_file_42():

    # 42nt long section of of Scaffold02974 around the first cut from enzyme BfaI
    #############################################################################
    # with recognition site C^TA_G at one-based coords: 21,22,23,24
    # GCGCTGGCCAGAACGTTCTC^TA_GGAATCGTGGAGAAGACATT

    assert (
        cut_file('tests/data/fourtytwobp.fasta') == ['fourtytwobps\t0\t20\tBfaI\t1000\t+\n',
                                                     'fourtytwobps\t22\t42\tBfaI\t1000\t-\n']
    )

def test_cut_file_1500():

    # FIRST 1500bps from xenopus-laevis scaffold used in crispr-eating paper
    ########################################################################

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

    cuts500 = cut_file('tests/data/Scaffold102974:1-1500_1500.fasta')

    assert(cuts500 == ['Scaffold102974:1-1500()\t786\t806\tBfaI\t1000\t+\n',
         'Scaffold102974:1-1500()\t808\t828\tBfaI\t1000\t-\n',
          'Scaffold102974:1-1500()\t1442\t1462\tBfaI\t1000\t+\n',
          'Scaffold102974:1-1500()\t1464\t1484\tBfaI\t1000\t-\n',
          'Scaffold102974:1-1500()\t732\t752\tScrFI\t1000\t+\n',
          'Scaffold102974:1-1500()\t753\t773\tScrFI\t1000\t-\n']
    )

# test_cut_file_42()
#
# test_cut_file_1500()