########################################################################################################################
#
# MAIN
#
# API functions: digest_and_blast_and_score_coord
#
########################################################################################################################

from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
# from __future__ import unicode_literals
from .bedtuple import filename_from_coord
from .digest import digest_coord
from .blast import blast
from .score import score
from .cluster import cluster


def blast_coord(coord, genome, direct, max_hsps, chunk_size):
    """
    Prerequisite is to have digested coord, therby created the corresponding file
    :param coord:
    :param blastdb:
    :param direct:
    :param max_hsps:
    :param chunk_size:
    :return:
    """
    return blast(filename_from_coord(coord), genome, direct, max_hsps, chunk_size)


def score_coord(nbrofchunks, coord, genome, direct, chunk_size, load_genome=False):
    """
    Prerequisite is to have digested coord, therby created the corresponding file
    :param nbrofchunks:
    :param coord:
    :param blastdb:
    :param direct:
    :param chunk_size:
    :param reref_substrate_id:
    :param oad_genome:
    :return:
    """
    return score(nbrofchunks, filename_from_coord(coord), genome, direct, chunk_size, load_genome=False)


def cluster_coord(guides, coord, direct, low, high, howmany):
    """
    Prerequisite is to have digested coord, therby created the corresponding file
    :param guides:
    :param coord:
    :param direct:
    :param reref_substrate_id:
    :param low:
    :param high:
    :param howmany:
    :return:
    """
    return cluster(guides, filename_from_coord(coord), direct, low, high, howmany)


##########
# MAIN API
##########
# requires:
# blastdb database
# dict      a genome fasta file against which pam lookup is optionnaly done for guides candidates
# genome    same, used when (correctly) blasting agaist whole geome
#

def digest_blast_score_cluster_prime(coord, genome, direct, max_hsps, chunk_size,
                                     low=75, high=75, load_genome=False, howmany=None,
                                     restriction_enzymes=(u'BfaI', u'ScrFI', u'HpaII')):

    _ = digest_coord(coord, genome, direct, restriction_enzymes)

    nbr = blast_coord(coord, genome, direct, max_hsps, chunk_size)

    guides = score_coord(nbr, coord, genome, direct, chunk_size, load_genome=load_genome)

    guides, groups = cluster_coord(guides, coord, direct, low, high, howmany)

    print('\nDESIGN PRIMERS')

    return guides, groups
