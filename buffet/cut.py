#######################################################################################################################
#  Digestion Cuts Generator
#
# Given a sequence faile in fasta/fastq format, returms a list of bed-formatted cuts
# for the crispr-eatimg dgestion process
#
# API functions: create_cutbedlines_from_seq_file
#
#######################################################################################################################

from os.path import splitext

import gzip
from string import lstrip

from Bio import SeqIO as seqio

from analysis import create_analysis


#######################################################################################################################
# Core formating function - depends on create_analysis (cutomized Bio.Restriction functionality)

def tabbed_string_from_list(list):
    return '\t'.join(map(str, list)) + '\n'


# TODO somehow retain the info about enzyme name (goes in seqrecord.description)
# TODO store a serial id number for each cut, concatenated with with _F or _R (goes in seqrecord.id)
# TODO store cutpos (goes in seqrecords.name)
def cut(enzyme, cutpos, substrate_id, substrate_end, sense='+'):
    """

    :param enzyme:
    :param cutpos:
    :param substrate_id:
    :param substrate_end:
    :param sense:
    :return:
    """
    if sense == '+':
        bedend = cutpos - 1
        bedstart = max(bedend - 20, 0)          # TODO do we still need the max?
    else:
        assert sense=='-'
        bedstart = cutpos - enzyme.ovhg - 1     # TODO do we still need the min?
        bedend = min(bedstart + 20, substrate_end)
    # cutinfo = [substrate_id, bedstart, bedend, str(enzyme), '1000', sense, bedstart, bedend, '255,0,0']
    cutinfo = [substrate_id, bedstart, bedend, str(enzyme), '1000', sense]
    cutbedline = tabbed_string_from_list(cutinfo)
    return cutbedline


############
# ITERATIONS
############

def cut_cutdict(cutdict, substrate_id, substrate_end):
    return [cut(enzyme, cutpos, substrate_id, substrate_end, sense)
            for enzyme in cutdict.keys()
            for index, cutpos in enumerate(cutdict[enzyme])
            for sense in ['+', '-']
            ]

def cut_seqrecord(seqrecord, enzyme_names=['BfaI', 'HpaII', 'ScrFI'] ):
    # cut_dict = create_analysis(seqrecord.seq, enzyme_names).full()        # TODO we should only loose 1 side, not both
    cut_dict = create_analysis(seqrecord.seq, enzyme_names).only_between(20, len(seqrecord)-20)   # TODO is 20 the best?
    return cut_cutdict(cut_dict, seqrecord.id, len(seqrecord))


def cut_seqrecords(seqrecords, enzyme_names=['BfaI', 'HpaII', 'ScrFI'] ):
    cuts = []
    for seqrecord in seqrecords:
        cuts.extend(cut_seqrecord(seqrecord, enzyme_names))
    return cuts

# IOs
######
# TODO can seqio handle the decompression itself?
# TODO write similar utlity to create a dict using SeqIO.to_dict
def create_seqrecords_from_file(filepath, fileformat):
    if fileformat[-3:] == '.gz':
        string_or_handle = gzip.open(filepath, 'r')
        parseformat = fileformat[:-3]
    else:
        string_or_handle = filepath
        parseformat = fileformat
    if parseformat=='fa':
        parseformat = 'fasta'
    return list(seqio.parse(string_or_handle, parseformat))


###################
# MAIN API FUNCTION
###################

def cut_file(filepath):
    """

    :param filepath: a sequence file name to load, extension should be one of: fasta, fastq, fastq.gz or fastq.gz
    :return: a list of strings, each a bed-formatted line for a found protospacer in
    """
    #
    fileformat = lstrip(splitext(filepath)[1], '.')
    seqrecords = create_seqrecords_from_file(filepath, fileformat)
    return cut_seqrecords(seqrecords)



# cut.py executed as script
###########################

if __name__ == '__main__':


    # 42nt long section of of Scaffold02974 around the first cut from enzyme BfaI
    # with recognition site C^TA_G at one-based coords: 21,22,23,24
    # GCGCTGGCCAGAACGTTCTC^TA_GGAATCGTGGAGAAGACATT
    cuts42 = cut_file('tests/data/fourtytwobp.fasta')

    assert(
        cuts42[0] == 'fourtytwobps\t0\t20\tBfaI\t1000\t+\n'
    )
    assert(
        cuts42[1] == 'fourtytwobps\t22\t42\tBfaI\t1000\t-\n'
    )

    cutbedlines1500 = cut_file('tests/data/Scaffold102974:1-1500_1500.fasta')

    assert(
        cutbedlines1500
        ==
        ['Scaffold102974:1-1500()\t786\t806\tBfaI\t1000\t+\n',
         'Scaffold102974:1-1500()\t808\t828\tBfaI\t1000\t-\n',
          'Scaffold102974:1-1500()\t1442\t1462\tBfaI\t1000\t+\n',
          'Scaffold102974:1-1500()\t1464\t1484\tBfaI\t1000\t-\n',
          'Scaffold102974:1-1500()\t732\t752\tScrFI\t1000\t+\n',
          'Scaffold102974:1-1500()\t753\t773\tScrFI\t1000\t-\n']
    )
