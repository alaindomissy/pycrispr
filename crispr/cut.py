#######################################################################################################################
#  Digestion Cuts Generator
#
# Given a sequence file in fasta/fastq format, returms a list of bed-formatted cuts
# for the crispr-eating dgestion process
#
# API functions: create_cutbedlines_from_seq_file
#
#######################################################################################################################

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
# from __future__ import unicode_literals
try:
    from io import StringIO        # python3
except ImportError:
    from StringIO import StringIO  # python2
from os.path import splitext
import gzip
from Bio import SeqIO as seqio
from .config import RESTRICTION_ENZYMES_LIST
from .analyse import analyse



#######################################################################################################################
# Core formatting function - depends on create_analysis (cutomized Bio.Restriction functionality)

def tabbed_string_from_list(list):
    """
    :param list:
    :return:
    >>>tabbed_string_from_list('Scaffold102974:1-1500()', 786, 806, 'BfaI', 1000)
    'Scaffold102974:1-1500()\t786\t806\tBfaI\t1000\n'
    """
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
    print('> buid forward and reverse guide annotations for all restriction cut loci')
    # print cutdict
    # print('> parsing guides from digestion loci', str(enzyme))    # TODO insert log within enzyme loop level
    return [cut(enzyme, cutpos, substrate_id, substrate_end, sense)
            for enzyme in cutdict.keys()
            for index, cutpos in enumerate(cutdict[enzyme])
            for sense in ['+', '-']
            ]

# this is where we call analyse
def cut_seqrecord(seqrecord, enzyme_names=RESTRICTION_ENZYMES_LIST):
    # cut_dict = create_analysis(seqrecord.seq, enzyme_names).full()
    print('> analyse loaded sequences with restriction batch of following enzymes', enzyme_names)
    cut_dict = analyse(seqrecord.seq, enzyme_names)
    return cut_cutdict(cut_dict, seqrecord.id, len(seqrecord))


def cut_seqrecords(seqrecords, enzyme_names=RESTRICTION_ENZYMES_LIST):
    return sum([cut_seqrecord(seqrecord, enzyme_names) for seqrecord in seqrecords],[])

# IOs
######
# TODO make seqio handle decompression?
# TODO use SeqIO.to_dict to create a dict
def create_seqrecords_from_fastafilepath(fastafilepath, fileformat):
    print('> load sequence from fasta file', fastafilepath)
    if fileformat[-3:] == '.gz':
        string_or_handle = gzip.open(fastafilepath, 'r')
        parseformat = fileformat[:-3]
    else:
        string_or_handle = fastafilepath
        parseformat = fileformat
    if parseformat=='fa':
        parseformat = 'fasta'
    return list(seqio.parse(string_or_handle, parseformat))


def create_seqrecords_from_unicodestring(unicodestring):
    string_or_handle = StringIO(unicodestring)
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
    seqrecords = create_seqrecords_from_unicodestring(unicodestring)
    return cut_seqrecords(seqrecords)


