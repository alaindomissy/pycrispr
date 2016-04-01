########################################################################################################################
#
# GENOMIC COORDINATES INTERFACE for digestion
#
# - accept bed format string input to specify and retrieve fasta scaffold to process
# - export and save bed formatted data for found protospacers
# - refocus on input genomic interval (via bed genomic arithmetic) and export bed and fasta for focused protospacers
#
########################################################################################################################

from __future__ import absolute_import, division, print_function
# from __future__ import unicode_literals
from os.path import splitext
from pybedtools import BedTool
from .cut import cut_fastafile
from .config import digestlog, genomes_path, protosp_path
from .bedtuple import Bedtuple

def genome_bedlines_save(cutbedlines,filepath):
    """
    saves a list of bed formatted strings as lines into a bed formatted file
    :param cutbedlines: list of bed formatted strings
    :param filepath: where to save bed file
    :return:
    """
    with open(filepath, 'w') as handle:
        for cutbedline in cutbedlines:
            handle.write(cutbedline)
    print("> save %s reference protospacers as %s" % (len(cutbedlines), filepath))


def bed_to_fasta(bedpath, referencefastapath):
    """
    given a bedfile, extract corresponding sequence from a reference fasta file and save as fasta file
    :param filepath: path to a bed file to be saved as fasta
    :param referencefastafilepath: file path to a fasta file used as sequence reference
    :return:
    """
    root, ext = splitext(bedpath)
    assert(ext=='.bed')
    fastapath = root + '.fasta'

    # TODO test for empty list of prsp in the intersection and avoid or handle bedtool.seqence error due to empty list
    BedTool(bedpath).sequence(fi=BedTool(referencefastapath), s=True).save_seqs(fastapath)

    return fastapath


# INITIALIZATION FUNCTION
#########################

def digest_genome(genome, restriction_enzymes=[u'BfaI', u'ScrFI', u'HpaII']):
    """
    Only needed for initialization, the first time a genome is being worked on.... TBA
    Creates files: filepath.prsp.bed and filepath.prsp.fasta
    Also creates index file filepath.fai if it was not there before
    This is the only  place that a call to cut_fastafile is made
    :param genome:
    :param restriction_enzymes:
    :return:

>>>digest_genome('phix')
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
    cutbedlines = cut_fastafile(genomes_path(genome))
    genome_bedlines_save(cutbedlines, protosp_path(genome))
    # NOT NEEDED
    # print("> save reference protospacers as ", end='')
    # bed_to_fasta(bedpath, fastafilepath)   # input fasta filepath serves as its own reference
    return(cutbedlines)


# THE WORKHORSE FUNCTION
########################
def digest_bedfile(filename, genome='mm8', directory='./', restriction_enzymes=(u'BfaI', u'ScrFI', u'HpaII')):
    """
    The core inner function handling digest. Saves 4 files
    :param bedfile:
    :param genome:
    :param restriction_enzymes:
    :return:
    """
    # root, ext = splitext(bedfile)
    # assert(ext == '.bed')
    # prspbefile = root + '.prsp.bed'
    bedpath = directory + filename + '.bed'
    prspbedpath = directory + filename + '.prsp.bed'
    digestlog("> load reference prsps at %s" % protosp_path(genome))
    digestlog("> intersect target with reference")
    digestlog("> save protospacers as", prspbedpath)
    BedTool(protosp_path(genome)).intersect(BedTool(bedpath)).moveto(prspbedpath)

    fastapath = bed_to_fasta(prspbedpath, genomes_path(genome))
    digestlog("> save protospacers as", fastapath)


# MAIN API FUNCTION
###################

def digest_coord(coord, genome='mm8', directory='./', restriction_enzymes=(u'BfaI', u'ScrFI', u'HpaII')):
    """
    You must have run digest_genome genome already,
    thereby creating protospacers files .prsp.bed' and .prsp.fasta'
    :param directory:
    :param coord: scaffold:start-end or scaffold:start_length
    :param genome:
    :return:cris
    >>>digest_coord('chr6:136640001-136680000', 'mm8')
    >>>digest_coord('chr6:136640001_40000', 'mm8')
    >>>digest_coord('phix:1-4000', 'phix')
    >>>digest_coord('chr:101001-105000', 'ecoli')
    """

    bedt = Bedtuple.from_coord(coord)
    bedfile = directory + bedt.filename + '.bed'

    digestlog('\n*******************************************************************')
    digestlog('DIGEST GENOMIC INTERVAL', bedt.filename)
    digestlog('*******************************************************************')

    BedTool(bedt.bedlist).moveto(bedfile)
    digestlog("> save target as", bedfile)

    fastapath = bed_to_fasta(bedfile, genomes_path(genome))
    digestlog("> save target as", fastapath)

    digest_bedfile(bedt.filename, genome, directory, restriction_enzymes)

    # TODO handle case when no prsp found, for ex if chosen coord is NNNN... only
    return bedt.filename


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
