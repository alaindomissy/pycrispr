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
from config import digestlog

def referencebedlines_save(cutbedlines, filepath):
    """
    saves a list of bed formatted strings as lines into a bed formatted file
    :param cutbedlines: list of bed formatted strings
    :param filepath: where to save bed file
    :return:
    """
    with open(filepath, 'w') as handle:
        for cutbedline in cutbedlines:
            handle.write(cutbedline)
    print("> save %s reference protospacers as bed file %s" % (len(cutbedlines), filepath))


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
    # test for empty list of prsp in the intersection and avoid or handle bedtool.seqence error due to empty lisy
    bedtool.sequence(fi=reference, s=True)
    bedtool.save_seqs(saveasfastapath)
    print("fasta file %s" % (saveasfastapath))


# INITIALIZATION FUNCTION
#########################

def digest_referencefastafile(fastafilepath, restriction_enzymes=[u'BfaI', u'ScrFI', u'HpaII']):
    """
    Only needed for initialization, the first time a genome is being worked on.... TBA
    Creates files: filepath.prsp.bed and filepath.prsp.fasta
    Also creates index file filepath.fai if it was not there before
    This is the only  place that a call to cut_fastafile is made
    :param filepath:
    :param restriction_enzymes:
    :return:

>>>crispr.digest.digest_referencefastafile('phix.fasta')
('> load sequence from fasta file', 'phix.fasta')
('> analyse loaded sequences with restriction batch of following enzymes', ['BfaI', 'HpaII', 'ScrFI'])
> buid forward and reverse guide annotations for all restriction cut loci
> save 22 refrence protospacers as bed file phix.fasta.prsp.bed
> save refrence protospacers as fasta file phix.fasta.prsp.fasta
['phix\t3116\t3136\tBfaI\t1000\t+\n',
 'phix\t3138\t3158\tBfaI\t1000\t-\n',
 'phix\t3887\t3907\tBfaI\t1000\t+\n',
 'phix\t3909\t3929\tBfaI\t1000\t-\n',
 'phix\t5041\t5061\tBfaI\t1000\t+\n',
 'phix\t5063\t5083\tBfaI\t1000\t-\n',
 'phix\t709\t729\tHpaII\t1000\t+\n',
 'phix\t731\t751\tHpaII\t1000\t-\n',
 'phix\t1083\t1103\tHpaII\t1000\t+\n',
 'phix\t1105\t1125\tHpaII\t1000\t-\n',
 'phix\t2780\t2800\tHpaII\t1000\t+\n',
 'phix\t2802\t2822\tHpaII\t1000\t-\n',
 'phix\t2999\t3019\tHpaII\t1000\t+\n',
 'phix\t3021\t3041\tHpaII\t1000\t-\n',
 'phix\t3347\t3367\tHpaII\t1000\t+\n',
 'phix\t3369\t3389\tHpaII\t1000\t-\n',
 'phix\t862\t882\tScrFI\t1000\t+\n',
 'phix\t883\t903\tScrFI\t1000\t-\n',
 'phix\t2781\t2801\tScrFI\t1000\t+\n',
 'phix\t2802\t2822\tScrFI\t1000\t-\n',
 'phix\t3481\t3501\tScrFI\t1000\t+\n',
 'phix\t3502\t3522\tScrFI\t1000\t-\n']
    """
    cutbedlines = cut_fastafile(fastafilepath)
    bedpath = fastafilepath + '.prsp.bed'
    referencebedlines_save(cutbedlines, bedpath)
    # NOT NEEDED !?
    # print("> save reference protospacers as ", end='')
    # bed_to_fasta(bedpath, fastafilepath)   # input fasta filepath serves as its own reference
    return(cutbedlines)


# THE WORKHORSE FUNCTION
########################
def digest_focused(focusfn, referencefastafilepath, restriction_enzymes=[u'BfaI', u'ScrFI', u'HpaII']):
    """
    The core inner function handling digest. Saves 4 files
    :param focusfn:
    :param reference:
    :param restriction_enzymes:
    :return:
    """
    targetbedfilepath = focusfn +'.bed'
    referenceprspbedfilepath = referencefastafilepath + '.prsp.bed'
    intargetprspbedfilepath = focusfn + ".prsp.bed"

    focus_bedtool = pybedtools.BedTool(targetbedfilepath)

    digestlog("> load protospacers from reference bed file %s" % referenceprspbedfilepath)
    whole_bedtool = pybedtools.BedTool(referenceprspbedfilepath)

    digestlog("> intersect target with reference bed file %s" % referenceprspbedfilepath)
    digestlog("> save in-target protospacers as bed file %s" % (intargetprspbedfilepath,))
    whole_bedtool.intersect(focus_bedtool).moveto(intargetprspbedfilepath)

    digestlog("> save target as ", end='')
    bed_to_fasta(targetbedfilepath, referencefastafilepath)

    digestlog("> save in-target protospacers as ", end='')
    bed_to_fasta(intargetprspbedfilepath, referencefastafilepath)


# INPUT HANDLING
################

def coord_to_bedtuple_and_filename(coord):
    """
    :param coord: scaffold:start-end or scaffold:start_length
    :return: a tuple of: a list of 1 3cols-bed-tuple, and a filename

    works with start and end:

    >>> coord_to_bedtuple_and_filename('chr6:136640001-136680000')
    ([('chr6', '136640000', '136680000')], 'chr6+136640001-136680000_40000')

    and also works with start and length:

    >>> coord_to_bedtuple_and_filename('chr6:136640001_40000')
    ([('chr6', '136640000', '136680000')], 'chr6+136640001-136680000_40000')
    """
    has_dash = '-' in coord
    has_under = '_' in coord
    assert(has_dash ^ has_under)
    if has_dash:
        start, end = coord.split(':')[1].split('-')
        length = str(int(end) - int(start) + 1)
    if has_under:
        start, length = stretch.split(':')[1].split('_')
        end = str(int(start) + int(length) - 1)
    chrom = coord.split(':')[0]
    start0 = str(int(start) - 1)    # bed coords are zero-based
    filename = '%s_%s-%s_%s' % (chrom, start, end, length)
    return [(chrom, start0, end)], filename



# MAIN API FUNCTION
###################

def digest_coord(direct, coord, reference, restriction_enzymes=[u'BfaI', u'ScrFI', u'HpaII']):
    """
    reference must be a path to a fastafile.
    You mast have run digest_referencefastafile on that file already,
    thereby creating protospacers files .prsp.bed' and .prsp.fasta'
    :param direct:
    :param coord: scaffold:start-end or scaffold:start_length
    :param reference: reference fasta file
    :return:cris
    >>> digest_coord('.', 'chr6:136640001-136680000', 'chr6.fasta')

    >>> digest_coord('.', 'chr6:136640001_40000', 'chr6.fasta')

    >>> crispr.digest.digest_coord('./', 'phix:1-4000', './phix.fasta')

bedtuplelist: [('phix', '0', '4000')] 	 focusfn: phix_1-4000_4000
DIGEST GENOMIC INTERVAL phix_1-4000_4000

> save target as bed file ./phix_1-4000_4000.bed
> load protospacers from reference bed file ./phix.fasta.prsp.bed
> intersect target with reference bed file ./phix.fasta.prsp.bed
> save in-target protospacers as bed file ./phix_1-4000_4000.prsp.bed
> save target as fasta file ./phix_1-4000_4000.fasta
> save in-target protospacers as fasta file ./phix_1-4000_4000.prsp.fasta

'phix_1-4000_4000'
 """
    bedtuplelist, focusfn = coord_to_bedtuple_and_filename(coord)
    # digestlog('bedtuplelist:', bedtuplelist, '\t', 'focusfn:', focusfn)
    digestlog('\nDIGEST GENOMIC INTERVAL', focusfn, '\n')
    pybedtools.BedTool(bedtuplelist).moveto(direct + focusfn + ".bed")
    digestlog("> save target as bed file %s" % (direct + focusfn + ".bed",))
    digest_focused(direct + focusfn, reference, restriction_enzymes)
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
