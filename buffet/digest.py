########################################################################################################################
# GENOMIC COORDOINATES INTERFACE for digestion
#
# - accept bed format string input to specify and retrieve fasta scaffold to process
# - export and save bed formatted data for found protospacers
# - refocus on input genomic interval (via bed genomic arithmetic) and export bed and fasta for focused protospacers
#
########################################################################################################################

from __future__ import print_function
from os.path import splitext
import pybedtools

from cut import cut_file


def bedlines_save(cutbedlines, filepath):
    """
    saves a list of bed formatted strings as lines into a bed formatted file
    :param cutbedlines: list of bed formatted strings
    :param filepath: where to save bed file
    :return:
    """
    with open(filepath, 'w') as handle:
        for cutbedline in cutbedlines:
            handle.write(cutbedline)


def bed_to_fasta(filepath, reference):
    """
    given a bedfile, extrat corresponding sequence from a reference fasta file and save as fasta file
    :param filepath: path to a bed file to be saved as fasta
    :param reference: file path to a fasta file used as sequence reference
    :return:
    """
    root, ext = splitext(filepath)
    assert(ext=='.bed')
    reference = pybedtools.BedTool(reference)
    bedtool = pybedtools.BedTool(filepath)
    bedtool.sequence(fi=reference, s=True).save_seqs(root + '.fasta')


# INITIALIZATION FUNCTION
#########################

def digest_file(filepath):
    """
    Only needed for initialization, the first time a genome is being worked on.... TBA
    :param filepath:
    :return:

    >>>digest_file('hg38.fa')
    "1000000 protospacers saved as hg38.fasta.prsp.bed and hg38.fasta.prsp.fasta"
    """
    cutbedlines = cut_file(filepath)
    saveasbedpath = filepath + '.prsp.bed'
    saveasfastapath = filepath + '.prsp.fasta'
    # TODO move this to right spot so we know how many prsps we have ahead of blasting
    print("%s protospacers saved as %s and %s" % (len(cutbedlines), saveasbedpath, saveasfastapath))
    bedlines_save(cutbedlines, saveasbedpath)
    bed_to_fasta(filepath, saveasfastapath )


# THE WORKHORSE FUNCTION
########################

def digest_focused(focusfn, reference):
    """
    The core inner function handling digest. Saves 4 files
    :param focusfn:
    :param reference:
    :return:
    """
    focus_bedtool = pybedtools.BedTool(focusfn +'.bed')

    whole_bedtool = pybedtools.BedTool(reference + '.prsp.bed')
    whole_bedtool.intersect(focus_bedtool).moveto(focusfn + ".prsp.bed")

    bed_to_fasta(focusfn, reference)
    bed_to_fasta(focusfn + '.prsp', reference)


# INPUT HANDLING
################

# TODO merge coord and stretch input options into a singke function
def coord_to_bedtuple_filename(coord):
    chrom = coord.split(':')[0]
    start, end = coord.split(':')[1].split('-')
    # bed coords are zero-based
    start0 = str(int(start) - 1)
    length = str(int(end) - int(start) + 1)
    filename = '%s:%s-%s_%s' % (chrom, start, end, length)
    return [(chrom, start0, end)], filename

assert(coord_to_bedtuple_filename('chr6:136640001-136680000')
    == ([('chr6', '136640000', '136680000')], 'chr6:136640001-136680000_40000'))


def stretch_to_bedtuple_filename(stretch):
    chrom = stretch.split(':')[0]
    start, length = stretch.split(':')[1].split('_')
    # bed coords are zero-based
    start0 = str(int(start) - 1)
    end = str(int(start) + int(length) - 1)
    filename =  '%s:%s-%s_%s' % (chrom, start, end, length)
    return [(chrom, start0, end)], filename


assert( stretch_to_bedtuple_filename('chr6:136640001_40000')
    == ([('chr6', '136640000', '136680000')], 'chr6:136640001-136680000_40000'))


# MAIN API FUNCTIONS
####################

def digest_coord(direct, coord, reference):
    bedtuplelist, focusfn = coord_to_bedtuple_filename(coord)
    print('digesting ',  focusfn, end='')
    pybedtools.BedTool(bedtuplelist).moveto(direct + focusfn + ".bed")
    digest_focused(direct + '/' + focusfn, reference)
    print('Done')
    return focusfn


def digest_stretch(direct, stretch, reference):
    bedtuplelist, focusfn = stretch_to_bedtuple_filename(stretch)
    print('digesting ',  focusfn, end='')
    pybedtools.BedTool(bedtuplelist).moveto(direct + focusfn + ".bed")
    digest_focused(direct + '/' + focusfn, reference)
    print('Done')
    return focusfn
