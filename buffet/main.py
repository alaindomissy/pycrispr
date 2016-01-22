########################################################################################################################
#
# MAIN
#
# API functions: digest_and_blast_and_score_coord
#
########################################################################################################################


##########
# TODO design primers
# TODO implement choice of enzymes
# TODO exclude redundant overlapping good guides from cluster yield count
# TODO option to batch together or separate each enzyme digestion
# TODO address truncated guides due to <20bp-distnat enzyme cuts (same or different): dicount goodones, watch bad ones
# TODO parameter(s) for priming: margin in/out from cluster end, possibly loose ending good guides
# TODO same MT constarint on all clusters primers, parameter for MT window
# TODO id the clusters
# TODO parameter chosen coord defined by length, common length menu
# TODO input bedfile as coord input
# TODO multi run (severeal coord, alternative enzyme choices, ...
# TODO cache all resuls to speed up furure use
# TODO logging utility function
# TODO adapt logs to basespace monitoring log: line length, emphasis
# TODO parameter logging options
# TODO pipeline upstream/ downstream other analises / apps in basespace
# TODO compare / integrate with usegalaxy, genebase ...
# TODO app intro section
# TODO app background literate code bacground: bio, algo, app tutorial
# TODO dockerfile and docker registry CI
# TODO use illumina docker registry, or private one
# TODO retry amazon docker machine
# TODO genomes data move to S3
# TODO rename github project
# TODO readme.md
# TODO conda package
# TODO reference Berkely
# TODO clean room
# TODO scikitbio inspiration for pycripr
# TODO testing framework
# TODO doctest



from buffet.digest import digest_coord, digest_stretch
from buffet.blast import blast
from buffet.score import score
from buffet.cluster import cluster



# DIGEST AND BLAST
#########################################################

def digest_and_blast_coord(direct, coord, reference, blastdb_db, chunk_size=50, max_hsps=50):
    focusfn = digest_coord(direct, coord, reference)
    nbr = blast(direct, focusfn + '.prsp', blastdb_db, chunk_size=chunk_size, max_hsps=max_hsps)
    return nbr

def digest_and_blast_stretch(direct, stretch, reference, blastdb_db, chunk_size=100, max_hsps=100):
    focusfn = digest_stretch(direct, stretch, reference)
    nbr = blast(direct, focusfn + '.prsp', blastdb_db, chunk_size=chunk_size, max_hsps=max_hsps)
    return nbr



# BLAST AND SCORE
#########################################################

def blast_and_score(direct, fn_noext, blastdb_db, chunk_size=50, max_hsps=50,
                    reref_substrate_id='chr6', low=75, high=75, load_genome=True, howmany=24):
    nbr = blast(direct, fn_noext, blastdb_db, chunk_size=chunk_size, max_hsps=max_hsps)
    guides = score(direct, fn_noext, blastdb_db, chunk_size=chunk_size, nbrofchunks=nbr,
          reref_substrate_id='chr6', load_genome=load_genome)
    return cluster(guides, direct, reref_substrate_id, low=75, high=75, howmany=None)


#########################################################
# MAIN API
#
# DIGEST AND BLAST AND SCORE
#########################################################

# requires:
# blastdb database
# dict      a genome fasta file against which pam lookup is optionnaly done for guides candidates
# genome    same, used when (correctly) blasting agaist whole geome
#
# TODO get rid of reref_substrate_id122
def digest_and_blast_and_score_coord(direct, coord, reference, blastdb_db, chunk_size=50, max_hsps=50,
                                     reref_substrate_id='chr6', low=75, high=75, load_genome=False, howmany=None):
    """

    :param direct:
    :param coord:
    :param reference:
    :param blastdb_db:
    :param chunk_size:
    :param max_hsps:
    :param reref_substrate_id:
    :param low:
    :param high:
    :param load_genome:
    :param howmany:
    :return:
    """
    focusfn = digest_coord(direct, coord, reference)
    fn_noext = focusfn + '.prsp'
    nbr = blast(direct, fn_noext, blastdb_db, chunk_size=chunk_size, max_hsps=max_hsps)
    guides = score(direct, fn_noext, blastdb_db, chunk_size=chunk_size, nbrofchunks=nbr,
          reref_substrate_id='chr6', load_genome=load_genome)
    return cluster(guides, direct, reref_substrate_id, low, high, howmany, fn_noext)


# TODO merge with the above onto one flexible input function
def digest_and_blast_and_score_stretch(direct, stretch, reference, blastdb_db, chunk_size=50, max_hsps=50,
                                       reref_substrate_id='chr6', low=75, high=75, load_genome=False, howmany=None):
    """

    :param direct:
    :param stretch:
    :param reference:
    :param blastdb_db:
    :param chunk_size:
    :param max_hsps:
    :param reref_substrate_id:
    :param low:
    :param high:
    :param load_genome:
    :param howmany:
    :return:
    """
    # TODO validation of input as either coord with a - or stretch with a _
    # assert type(stretch) is type(""), "requires a string"
    focusfn = digest_stretch(direct, stretch, reference)
    fn_noext = focusfn + '.prsp'
    nbr = blast(direct, fn_noext, blastdb_db, chunk_size=chunk_size, max_hsps=max_hsps)
    guides = score(direct, fn_noext, blastdb_db, chunk_size=chunk_size, nbrofchunks=nbr,
          reref_substrate_id='chr6', load_genome=load_genome)
    return cluster(guides, direct, reref_substrate_id, low, high, howmany, fn_noext)
