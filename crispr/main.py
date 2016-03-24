########################################################################################################################
#
# MAIN
#
# API functions: digest_and_blast_and_score_coord
#
########################################################################################################################

from __future__ import absolute_import, division, print_function
# from __future__ import unicode_literals
from .bedtuple import filename_from_coord
from .digest import digest_coord
from .blast import blast
from .score import score
from .stretch import stretch, run
from .amplicon import amplicon
from .prime import prime


def blast_coord(coord, genome, directory, max_hsps, chunk_size):
    """
    Prerequisite is to have digested coord, therby created the corresponding file
    :param coord:
    :param genome:
    :param directory: default is ./
    :param max_hsps: default is 10
    :param chunk_size: default is 10
    :return:
    """
    return blast(filename_from_coord(coord), genome, directory, max_hsps, chunk_size)


def score_coord(nbrofchunks, coord, genome, directory, chunk_size, load_genome=False):
    """
    Prerequisite is to have digested coord, therby created the corresponding file
    :param nbrofchunks:
    :param coord:
    :param genome:
    :param directory:
    :param chunk_size:
    :param load_genome:
    :return:
    """
    return score(nbrofchunks, filename_from_coord(coord), genome, directory, chunk_size, load_genome=False)


def stretch_coord(coord, directory, low, high, howmany):
    """
    Prerequisite is to have digested coord, therby created the corresponding file
    :param coord:
    :param directory:
    :param low:
    :param high:
    :param howmany:
    :return:
    """
    return stretch(filename_from_coord(coord), directory, low, high, howmany)

def run_coord(coord, directory, high):
    """
    Prerequisite is to have digested coord, therby created the corresponding file
    :param guides:
    :param coord:
    :param directory:
    :param high:
    :return:
    """
    return run(filename_from_coord(coord), directory, high)

def amplicon_coord(coord, directory, genome, substrate, threshold):
    """
    Prerequisite is to have digested coord, therby created the corresponding file
    :param guides:
    :param coord:
    :param directory:
    :param high:
    :return:
    """
    return amplicon(filename_from_coord(coord), directory, genome, substrate, threshold)

def prime_coord(coord, directory, genome, high):
    """
    Prerequisite is to have digested coord, therby created the corresponding file
    :param guides:
    :param coord:
    :param directory:
    :param high:
    :return:
    """
    return prime(filename_from_coord(coord), directory, genome, high)


# def digest_blast_score_cluster_prime_coord(coord, genome, directory, max_hsps, chunk_size, high):
#     """
#     Prerequisite is to have digested coord, therby created the corresponding file
#     :param guides:
#     :param coord:
#     :param directory:
#     :param high:
#     :return:
#     """
#     return digest_blast_score_cluster_prime(filename_from_coord(coord), genome, directory, max_hsps, chunk_size,  high)


##########
# MAIN API
##########
# requires:
# blastdb database
# dict      a genome fasta file against which pam lookup is optionnaly done for guides candidates
# genome    same, used when (correctly) blasting agaist whole geome
#

def digest_blast_score_cluster_prime(coord, genome, directory, max_hsps, chunk_size,
                                     low=75, high=94, load_genome=False, howmany=0,
                                     restriction_enzymes=(u'BfaI', u'ScrFI', u'HpaII')):

    filename = digest_coord(coord, genome, directory, restriction_enzymes)
    print("filename : " , filename)
    # nbr =
    nbr = blast_coord(coord, genome, directory, max_hsps, chunk_size)
    # guides =
    score_coord(nbr, coord, genome, directory, chunk_size, load_genome)

    guides, stretches = stretch_coord(coord, directory, low, high, howmany)

    # guides, runs = run_coord(coord, directory, high)
    # amplicons = amplicon_coord(coord, directory, genome, high)      # TODO de-hardcode substrate
    prime_coord(coord, directory, genome, high)

    return guides, stretches  #, amplicons
