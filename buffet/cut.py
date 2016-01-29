#######################################################################################################################
#  Digestion Cuts Generator
#
# Given a sequence file in fasta/fastq format, returms a list of bed-formatted cuts
# for the crispr-eating dgestion process
#
# API functions: create_cutbedlines_from_seq_file
#
#######################################################################################################################

from os.path import splitext
import gzip
try:
    from io import StringIO        # python3
except ImportError:
    from StringIO import StringIO  # python2

from Bio import SeqIO as seqio

from buffet.settings import RESTRICTION_ENZYMES_LIST
from buffet.analysis import create_analysis



#######################################################################################################################
# Core formatting function - depends on create_analysis (cutomized Bio.Restriction functionality)

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
    print('> parsing guides from digestion loci', cutdict)
    # print('> parsing guides from digestion loci', str(enzyme))    # TODO insert log within enzyme loop level
    return [cut(enzyme, cutpos, substrate_id, substrate_end, sense)
            for enzyme in cutdict.keys()
            for index, cutpos in enumerate(cutdict[enzyme])
            for sense in ['+', '-']
            ]

def cut_seqrecord(seqrecord, enzyme_names=RESTRICTION_ENZYMES_LIST):
    # cut_dict = create_analysis(seqrecord.seq, enzyme_names).full()        # TODO we should only loose 1 side, not both
    print('> digesting using enzymes', enzyme_names)
    cut_dict = create_analysis(seqrecord.seq, enzyme_names).only_between(20, len(seqrecord)-20)   # TODO is 20 the best?
    return cut_cutdict(cut_dict, seqrecord.id, len(seqrecord))


def cut_seqrecords(seqrecords, enzyme_names=RESTRICTION_ENZYMES_LIST):
    return sum([cut_seqrecord(seqrecord, enzyme_names) for seqrecord in seqrecords])

# IOs
######
# TODO can seqio handle the decompression itself?
# TODO write similar utlity to create a dict using SeqIO.to_dict
def create_seqrecords_from_fastafilepath(fastafilepath, fileformat):
    print('> loading sequence from ref genome for interval', fastafilepath)
    if fileformat[-3:] == '.gz':
        string_or_handle = gzip.open(fastafilepath, 'r')
        parseformat = fileformat[:-3]
    else:
        string_or_handle = fastafilepath
        parseformat = fileformat
    if parseformat=='fa':
        parseformat = 'fasta'
    return list(seqio.parse(string_or_handle, parseformat))


###################
# MAIN API FUNCTION
###################

def cut_fastafile(fastafilepath):
    """

    :param fastafilepath: a sequence file name to load, extension should be one of: fasta, fastq, fastq.gz or fastq.gz
    :return: a list of strings, each a bed-formatted line for a found protospacer in fastafilepath
    """
    #
    fileformat = splitext(fastafilepath)[1].strip('.')
    seqrecords = create_seqrecords_from_fastafilepath(fastafilepath, fileformat)
    return cut_seqrecords(seqrecords)

def cut_unicodestring(unicodestring):
    """

    :param unicodestring: a sequence as a unicode string
    :return:  a list of strings, each a bed-formatted line for a found protospacer in unicodestring
    """
    cut_file(StringIO(unicodestring))