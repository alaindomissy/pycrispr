########################################################################################################################
#
# BLAST
#
#
#
########################################################################################################################

from digest import digest_coord, digest_stretch
from blast import blast
from score import score
from cluster import cluster


#########################################################
# DIGEST AND BLAST
#########################################################

def digest_and_blast_coord(direct, coord, blastdb_db, fileformat, chunksize=100, max_hsps=100):
    focusfn = digest_coord(direct, coord, blastdb_db, fileformat)
    nbr = blast(direct, focusfn + '.prsp', '/media/mis/BLASTDB/', blastdb_db, chunksize=chunksize, max_hsps=max_hsps)
    return nbr

def digest_and_blast_stretch(direct, stretch, blastdb_db, fileformat, chunksize=100, max_hsps=100):
    focusfn = digest_stretch(direct, stretch, blastdb_db, fileformat)
    nbr = blast(direct, focusfn + '.prsp', '/media/mis/BLASTDB/', blastdb_db, chunksize=chunksize, max_hsps=max_hsps)
    return nbr


#########################################################
# DIGEST AND BLAST AND SCORE
#########################################################

# requires:
# blastdb database
# dict         a genome fasta file against which pam lookup is optionnaly done for guides candidates
# genome       same, used when (correctly) blasting agaist whole geome
#
def digest_and_blast_and_score_coord(direct, coord, blastdb_db, fileformat, chunksize=100, max_hsps=100,
                                     reref_substrate_id='chr6', low=50, high=75, load_genome=True, howmany=None):
    """

    :param direct:
    :param coord:
    :param blastdb_db:
    :param fileformat:
    :param chunksize:
    :param max_hsps:
    :param reref_substrate_id:
    :param low:
    :param high:
    :param load_genome:
    :param howmany:
    :return:
    """
    focusfn = digest_coord(direct, coord, blastdb_db, fileformat)
    fn_noext = focusfn + '.prsp'
    blastdb_directory = '/media/mis/BLASTDB/'
    nbr = blast(direct, fn_noext ,blastdb_directory , blastdb_db, chunksize=chunksize, max_hsps=max_hsps)
    guides = score(direct, fn_noext, blastdb_directory, blastdb_db, chunksize=chunksize, nbrofchunks=nbr,
          reref_substrate_id='chr6', load_genome=load_genome)
    return cluster(guides, direct, reref_substrate_id, low=50, high=75, howmany=None)


def digest_and_blast_and_score_stretch(direct, stretch, blastdb_db, fileformat, chunksize=100, max_hsps=100,
                                     reref_substrate_id='chr6', low=50, high=75, load_genome=True, howmany=None):
    """

    :param direct:
    :param stretch:
    :param blastdb_db:
    :param fileformat:
    :param chunksize:
    :param max_hsps:
    :param reref_substrate_id:
    :param low:
    :param high:
    :param load_genome:
    :param howmany:
    :return:
    """
    focusfn = digest_stretch(direct, stretch, blastdb_db, fileformat)
    fn_noext = focusfn + '.prsp'
    blastdb_directory = '/media/mis/BLASTDB/'
    nbr = blast(direct, fn_noext , '/media/mis/BLASTDB/', blastdb_db, chunksize=chunksize, max_hsps=max_hsps)
    guides = score(direct, fn_noext, blastdb_directory, blastdb_db, chunksize=chunksize, nbrofchunks=nbr,
          reref_substrate_id='chr6', load_genome=load_genome)
    return cluster(guides, direct, reref_substrate_id, low=50, high=75, howmany=None)


#########################################################
# BLAST AND SCORE
#########################################################

def blast_and_score(direct, fn_noext, blastdb_directory, blastdb_db, chunksize=100, max_hsps=None,
                    reref_substrate_id='chr6', low=50, high=75, load_genome=True, howmany=24):

    nbr = blast(direct, fn_noext, blastdb_directory, blastdb_db, chunksize=chunksize, max_hsps=max_hsps)
    guides = score(direct, fn_noext, blastdb_directory, blastdb_db, chunksize=chunksize, nbrofchunks=nbr,
          reref_substrate_id='chr6', load_genome=load_genome)
    return cluster(guides, direct, reref_substrate_id, low=50, high=75, howmany=None)

