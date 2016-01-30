########################################################################################################################
#
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
from cut import cut_fastafile


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
    print("> saving %s protospacers" % (len(cutbedlines),))
    print("> saved as bed file %s" % (filepath,))


def bed_to_fasta(bedfilepath, referencefastafilepath):
    """
    given a bedfile, extrat corresponding sequence from a reference fasta file and save as fasta file
    :param filepath: path to a bed file to be saved as fasta
    :param referencefastafilepath: file path to a fasta file used as sequence reference
    :return:
    """
    root, ext = splitext(bedfilepath)
    assert(ext=='.bed')
    saveasfastapath = root + '.fasta'
    reference = pybedtools.BedTool(referencefastafilepath)
    bedtool = pybedtools.BedTool(bedfilepath)
    bedtool.sequence(fi=reference, s=True).save_seqs(saveasfastapath)
    print("fasta file %s" % (saveasfastapath))


# INITIALIZATION FUNCTION
#########################

def digest_fastafile(fastafilepath):
    """
    Only needed for initialization, the first time a genome is being worked on.... TBA
    Creates files: filepath.prsp.bed and filepath.prsp.fasta
    Also creates index file filepath.fai if it was not there before
    This is the only  place that a call to cut_fastafile is made

    :param filepath:
    :return:
    >>>digest_fastafile('hg38.fa')
    "1000000 protospacers saved as hg38.fasta.prsp.bed and hg38.fasta.prsp.fasta"
        ('> loading sequence from ref genome for interval', 'mm8.fasta')
        ('> digesting using enzymes', ['BfaI', 'HpaII', 'ScrFI'])
        ('> parsing guides from digestion loci', .....
        ('> digesting using enzymes', ['BfaI', 'HpaII', 'ScrFI'])
        > saved as bed file mm8.fasta.prsp.bed
        > saved as fasta file mm8.fasta.prsp.fasta
        index file mm8.fasta.fai not found, generating...
    """
    cutbedlines = cut_fastafile(fastafilepath)
    bedpath = fastafilepath + '.prsp.bed'
    bedlines_save(cutbedlines, bedpath)
    bed_to_fasta(bedpath, fastafilepath)   # input fasta filepath serves as its own reference
    return(cutbedlines)

# THE WORKHORSE FUNCTION
########################

def digest_focused(focusfn, referencefastafilepath):
    """
    The core inner function handling digest. Saves 4 files
    :param focusfn:
    :param reference:
    :return:
    """
    focus_bedtool = pybedtools.BedTool(focusfn +'.bed')
    referenceprspbedfilepath = referencefastafilepath + ".prsp.bed"
    print("> load protospacers from reference bed file %s" % referenceprspbedfilepath)
    whole_bedtool = pybedtools.BedTool(referencefastafilepath + '.prsp.bed')
    print("> intersect target with reference bed file %s" % referenceprspbedfilepath)
    whole_bedtool.intersect(focus_bedtool).moveto(focusfn + ".prsp.bed")
    print("> save in-target protospacers as bed file %s" % (focusfn + ".prsp.bed",))

    print("> save target as ", end='')
    bed_to_fasta(focusfn + '.bed', referencefastafilepath)
    print("> save in-target protospacers as ", end='')
    bed_to_fasta(focusfn + '.prsp.bed', referencefastafilepath)


# INPUT HANDLING
################

# TODO merge coord and stretch input options into a singke function
def coord_to_bedtuple_filename(coord):
    """

    :param coord:
    :return:
    >>> coord_to_bedtuple_filename('chr6:136640001-136680000')
    ([('chr6', '136640000', '136680000')], 'chr6+136640001-136680000_40000')
    """
    chrom = coord.split(':')[0]
    start, end = coord.split(':')[1].split('-')
    # bed coords are zero-based
    start0 = str(int(start) - 1)
    length = str(int(end) - int(start) + 1)
    filename = '%s+%s-%s_%s' % (chrom, start, end, length)
    return [(chrom, start0, end)], filename

assert(coord_to_bedtuple_filename('chr6:136640001-136680000')
    == ([('chr6', '136640000', '136680000')], 'chr6+136640001-136680000_40000'))


def stretch_to_bedtuple_filename(stretch):
    """

    :param stretch:
    :return:
    >>> stretch_to_bedtuple_filename('chr6:136640001_40000')
    ([('chr6', '136640000', '136680000')], 'chr6+136640001-136680000_40000')
    """
    chrom = stretch.split(':')[0]
    start, length = stretch.split(':')[1].split('_')
    # bed coords are zero-based
    start0 = str(int(start) - 1)
    end = str(int(start) + int(length) - 1)
    filename =  '%s+%s-%s_%s' % (chrom, start, end, length)
    return [(chrom, start0, end)], filename

assert( stretch_to_bedtuple_filename('chr6:136640001_40000')
    == ([('chr6', '136640000', '136680000')], 'chr6+136640001-136680000_40000'))


# MAIN API FUNCTIONS
####################

def digest_coord(direct, coord, reference, restriction_enzymes=[u'BfaI', u'ScrFI', u'HpaII']):
    """
    reference must be a path to a fastafile.
    You mast have run digest_fastafile on that file already,
    thereby creatang protospacers files .prsp.bed' and .prsp.fasta'

    :param direct:
    :param coord:
    :param reference:
    :return:
    >>> digest_coord('.', 'chr6:136640001-136680000', 'chr6.fasta')
    """
    bedtuplelist, focusfn = coord_to_bedtuple_filename(coord)
    print('bedtuplelist:', bedtuplelist, '\t', 'focusfn:', focusfn)
    print('\nDIGEST GENOMIC INTERVAL', focusfn)
    pybedtools.BedTool(bedtuplelist).moveto(direct + focusfn + ".bed")
    print("> save target as bed file %s" % (direct + focusfn + ".bed",))
    digest_focused(direct + focusfn, reference)
    return focusfn


def digest_stretch(direct, stretch, reference, restriction_enzymes=[u'BfaI', u'ScrFI', u'HpaII']):
    """

    :param direct:
    :param stretch:
    :param reference:
    :return:
    """
    bedtuplelist, focusfn = stretch_to_bedtuple_filename(stretch)
    print('\nDIGEST GENOMIC INTERVAL\n', focusfn)
    pybedtools.BedTool(bedtuplelist).moveto(direct + focusfn + ".bed")
    print("> save target as bed file %s" % (focusfn + ".bed",))
    digest_focused(direct + '/' + focusfn, reference)
    print('...done')
    return focusfn




# interface to prime
####################


def digest_target(target):
    guidelist = []
    return guidelist


def count_non_overlapping_guides(guidelist, binding_interference_spacing=20):
    '''
    Sequence position information must be encoded in sequence name attribute,
    e.g. name="100" indicates that left edge of guide (regardless of strand)
    starts 100nt along the target scaffold.

    Optional argument binding_interference_spacing specifies the number of
    nucleotides apart that a guides must be to contribute to unique
    labeling events (for example, of two guides only 10nt apart, it's impossible
    for both to bind at once)

    From the spCas9 crystal structure, it looks like guides will need to be
    at least 24-25 nucleotides apart to bind simultaneously, and possibly
    more.
    '''
    prev_position = 0
    guidecount = 0
    for item in guidelist:
        distance_from_last = int(item.name) - prev_position
        if distance_from_last > binding_interference_spacing:
            guidecount = guidecount + 1
        prev_position = int(item.name)
    return guidecount


def nonoverlapping_guidecount(target):
    count_non_overlapping_guides(digest_target(target))