#######################################################################################################################
#  Digestion Cuts Generator
#
# Given a sequence faile in fasta/fastq format, returms a list of bed-formatted cuts
# for the crispr-eatimg dgestion process
#
# API functions: create_cutbedlines_from_seq_file
#
#######################################################################################################################

import gzip

from Bio import SeqIO as seqio

from analysis import create_analysis


#######################################################################################################################
# Core formating function - depends on create_analysis (cutomized Bio.Restriction functionality)

def tabbed_string_from_list(list):
    return '\t'.join(map(str, list)) + '\n'

# TODO somehow retain the info about enzyme name (goes in seqrecord.description)
# TODO store a serial id number for each cut, concatenated with with _F or _R (goes in seqrecord.id)
# TODO store cutpos (goes in seqrecords.name)
def create_cutbedline(enzyme, cutpos, substrate_id, substrate_end, sense='+'):
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

def create_cutbedlines(cutdict, substrate_id, substrate_end):
    return [create_cutbedline(enzyme, cutpos, substrate_id, substrate_end, sense)
            for enzyme in cutdict.keys()
            for index, cutpos in enumerate(cutdict[enzyme])
            for sense in ['+', '-']
            ]

def digest_seqrecord_to_cutbedlines(seqrecord, enzyme_names=['BfaI', 'HpaII', 'ScrFI'] ):
    # cut_dict = create_analysis(seqrecord.seq, enzyme_names).full()        # TODO we should only loose 1 side, not both
    cut_dict = create_analysis(seqrecord.seq, enzyme_names).only_between(20, len(seqrecord)-20)   # TODO is 20 the best?
    return create_cutbedlines(cut_dict, seqrecord.id, len(seqrecord))


def digest_seqrecords_to_cutbedlines(seqrecords, enzyme_names=['BfaI', 'HpaII', 'ScrFI'] ):
    cutinfos = []
    for seqrecord in seqrecords:
        cutinfos.extend(digest_seqrecord_to_cutbedlines(seqrecord, enzyme_names))
    return cutinfos


######
# IOs
######
# TODO can seqio handle the decompression itself?
# TODO write similar utlity to creat a dict using SeqIO.to_dict
def create_seqrecords_from_file(direct, fn_noext, fileformat):
    if fileformat[-3:] == '.gz':
        string_or_handle = gzip.open(direct + fn_noext, + '.' + fileformat, 'r')
        parseformat = fileformat[:-3]
    else:
        string_or_handle = direct + fn_noext + '.'+ fileformat
        parseformat = fileformat
    if parseformat=='fa':
        parseformat = 'fasta'
    return list(seqio.parse(string_or_handle, parseformat))


###################
# MAIN API FUNCTION
###################

def create_cutbedlines_from_seq_file(direct, fn_noext , fileformat='fasta'):
    """

    :param direct: folder to look in
    :param fn_noext: a sequence file name to load, no extension
    :param fileformat: fasta, fastq, fastq.gz or fastq.gz, format and extensiom of file to load
    :return: a list of strings, each a bed-formatted line for a found protospacer in
    """
    seqrecords = create_seqrecords_from_file(direct, fn_noext , fileformat)
    return digest_seqrecords_to_cutbedlines(seqrecords)

