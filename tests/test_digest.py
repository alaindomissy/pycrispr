import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from crispr.digest import coord_to_bedtuple_and_filename, digest_genome


def test_coord_to_bedtuple_filename():

    # has_dash
    assert(coord_to_bedtuple_and_filename('chr6:136640001-136680000') ==
           ([('chr6', '136640000', '136680000')], 'chr6_136640001-136680000_40000.bed')
           )

    # has_under()
    assert(coord_to_bedtuple_and_filename('chr6:136640001_40000') ==
           ([('chr6', '136640000', '136680000')], 'chr6_136640001-136680000_40000.bed')
           )



# 42nt long section of of Scaffold02974 around the first cut from enzyme BfaI
#############################################################################
# with recognition site C^TA_G at one-based coords: 21,22,23,24
# GCGCTGGCCAGAACGTTCTC^TA_GGAATCGTGGAGAAGACATT

def test_digest_fastafile_42():
    pass

    # this requires bedtools being installed
    #
    # digest_genome('tests/data/fourtytwobp.fasta')
    #
    # with open('tests/data/fourtytwobp.fasta.prsp.bed') as prspbed:
    #     str(prspbed.read())=='fourtytwobps\t0\t20\tBfaI\t1000\t+\nfourtytwobps\t22\t42\tBfaI\t1000\t-\n'
    #
    # with open('tests/data/fourtytwobp.fasta.prsp.fasta') as prspfasta:
    #      str(prspfasta.read())=='fourtytwobps\t0\t20\tBfaI\t1000\t+\nfourtytwobps\t22\t42\tBfaI\t1000\t-\n'



"""
In [101]: cat phix_1-4000_4000.prsp.bed
phix	3116	3136	BfaI	1000	+
phix	3138	3158	BfaI	1000	-
phix	3887	3907	BfaI	1000	+
phix	3909	3929	BfaI	1000	-
phix	709	729	HpaII	1000	+
phix	731	751	HpaII	1000	-
phix	1083	1103	HpaII	1000	+
phix	1105	1125	HpaII	1000	-
phix	2780	2800	HpaII	1000	+
phix	2802	2822	HpaII	1000	-
phix	2999	3019	HpaII	1000	+
phix	3021	3041	HpaII	1000	-
phix	3347	3367	HpaII	1000	+
phix	3369	3389	HpaII	1000	-
phix	862	882	ScrFI	1000	+
phix	883	903	ScrFI	1000	-
phix	2781	2801	ScrFI	1000	+
phix	2802	2822	ScrFI	1000	-
phix	3481	3501	ScrFI	1000	+
phix	3502	3522	ScrFI	1000	-

In [102]: cat phix_1-4000_4000.prsp.fasta
>phix:3116-3136(+)
AACCCTGATGAGGCCGCCCC
>phix:3138-3158(-)
CATAGCACCAGAAACAAAAC
>phix:3887-3907(+)
GGTATTGGCTCTAATTTGTC
>phix:3909-3929(-)
CAATCCTGACGGTTATTTCC
>phix:709-729(+)
TATTAATGGCGTCGAGCGTC
>phix:731-751(-)
AACAATTCAGCGGCTTTAAC
>phix:1083-1103(+)
TATTACCATTTCAACTACTC
>phix:1105-1125(-)
AAGGAGTCGCCAGCGATAAC
>phix:2780-2800(+)
TGACGTCCTTCCTCGTACGC
>phix:2802-2822(-)
CCAACATAAACATTATTGCC
>phix:2999-3019(+)
GGCGGTCAAAAAGCCGCCTC
>phix:3021-3041(-)
CACATCACCTTGAATGCCAC
>phix:3347-3367(+)
TCTGCTGGTATGGTTGACGC
>phix:3369-3389(-)
CTCTTTTTGATTCTCAAATC
>phix:862-882(+)
AAACGTTCTGGCGCTCGCCC
>phix:883-903(-)
CGCAACGGCTGCGGACGACC
>phix:2781-2801(+)
GACGTCCTTCCTCGTACGCC
>phix:2802-2822(-)
CCAACATAAACATTATTGCC
>phix:3481-3501(+)
ACGCCAGAATACGAAAGACC
>phix:3502-3522(-)
TCTCATTTTGTGCATATACC
"""