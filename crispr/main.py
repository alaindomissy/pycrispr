########################################################################################################################
#
# MAIN
#
# API functions: digest_and_blast_and_score_coord
#
########################################################################################################################


from bedtuple import filename_from_coord
from digest import digest_coord
from blast import blast
from score import score
from cluster import cluster


def blast_coord(coord, blastdb='mm8', dir='./', max_hsps=10, chunk_size=50):
    """
    Prerequisite is to have digested coord, therby created the corresponding file
    :param coord:
    :param blastdb:
    :param dir:
    :param max_hsps:
    :param chunk_size:
    :return:
    """
    return blast(filename_from_coord(coord), blastdb='mm8', dir='./', max_hsps=10, chunk_size=50)


def score_coord(nbrofchunks, coord, blastdb='mm8', dir='./', chunk_size=20,
          reref_substrate_id=None, load_genome=False):
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
    score(nbrofchunks, filename_from_coord(coord), blastdb='mm8', direct='./', chunk_size=20,
          reref_substrate_id=None, load_genome=False)


def cluster_coord(guides, coord, direct='./', reref_substrate_id='mm8', low=75, high=75, howmany=12):
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

    cluster(guides, filename_from_coord(coord), direct='./', reref_substrate_id='mm8', low=75, high=75, howmany=12)


##########
# MAIN API
##########
# requires:
# blastdb database
# dict      a genome fasta file against which pam lookup is optionnaly done for guides candidates
# genome    same, used when (correctly) blasting agaist whole geome
#
# TODO get rid of reref_substrate_id
def digest_and_blast_and_score_coord(coord, genome='mm8', dir='./', max_hsps=10, chunk_size=20,
                                     reref_substrate_id=None, low=75, high=75, load_genome=False, howmany=None,
                                     restriction_enzymes=(u'BfaI', u'ScrFI', u'HpaII')):

    _ = digest_coord(coord, genome, dir, restriction_enzymes)

    nbr = blast_coord(coord, blastdb=genome, dir=dir, max_hsps=max_hsps, chunk_size=chunk_size)

    guides = score_coord(nbr, coord, blastdb=genome, dir=dir, chunk_size=chunk_size,
                   reref_substrate_id='chr6', load_genome=load_genome)

    guides, groups = cluster_coord(guides, coord, dir, reref_substrate_id, low, high, howmany)

    return guides, groups
