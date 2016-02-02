########################################################################################################################
#
# MAIN
#
# API functions: digest_and_blast_and_score_coord
#
########################################################################################################################
# TODO design primers
# TODO implement choice of enzymes
# TODO exclude redundant overlapping good guides from cluster yield count
# TODO option to batch together or separate each enzyme digestion
# TODO address truncated guides due to <20bp-distnat enzyme cuts (same or different): dicount goodones, watch bad ones
# TODO parameter(s) for priming: margin in/out from cluster end, possibly loose ending good guides
# TODO same MT constraint on all clusters primers, parameter for MT window
# TODO id the clusters
# TODO parameter chosen coord defined by length, common length menu
# TODO input bedfile as coord input
# TODO multi run (severeal coord, alternative enzyme choices, ...
# TODO cache all results to speed up furure use
# TODO logging utility function
# TODO adapt logs to basespace monitoring log: line length, emphasis
# TODO parameter logging options
# TODO pipeline upstream/ downstream other analises / apps in basespace
# TODO compare / integrate with usegalaxy, genebase ...
# TODO app intro section
# TODO app background literate code bacground: bio, algo, app tutorial
# TODO dockerfile and docker registry CI
# TODO new doker image lighter weight and less layers
# TODO use illumina docker registry, or private one
# TODO retry amazon docker machine
# TODO genomes data move to S3
# TODO rename github project
# TODO conda package
# TODO reference Berkely
# TODO clean room
# TODO scikitbio inspiration for pycripr
# TODO testing fixtures
# TODO doctest
# TODO gui to input number of good guides desired
# TODO max length of a cluster 2kb
# TODO figure out which ecoli K12 or ?
# TODO use human blast data from igenomes
# TODO check plastid, metaseq, bioconda, track
# TODO bigBed bigWig


from bedtuple import filename_from_coord
from digest import digest_coord
from blast import blast
from score import score
from cluster import cluster


def blast_coord(coord, blastdb='mm8', dir='./', max_hsps=10, chunk_size=50):
    return blast(filename_from_coord(coord), blastdb='mm8', dir='./', max_hsps=10, chunk_size=50)


# # DIGEST AND BLAST
# ##################
#
# def digest_and_blast_coord(coord, genome='mm8', direct='./', blastdb='mm8', chunk_size=50, max_hsps=50):
#
#     filenamenoext = digest_coord(coord, genome, direct)
#
#     nbr = blast(filenamenoext + '.bed.prsp', blastdb, direct, max_hsps=max_hsps, chunk_size=chunk_size)
#
#     return nbr
#
#
# # BLAST AND SCORE
# #################
#
# def blast_and_score(filenamenoext, direct='./', blastdb='mm8',
#                     chunk_size=50, max_hsps=50,
#                     reref_substrate_id=None, low=75, high=75, load_genome=True, howmany=None):
#
#     nbr = blast(filenamenoext, blastdb, direct, max_hsps=max_hsps, chunk_size=chunk_size)
#
#     guides = score(direct, filenamenoext, blastdb, chunk_size=chunk_size, nbrofchunks=nbr,
#                    reref_substrate_id=None, load_genome=load_genome)
#
#     return cluster(guides, direct, reref_substrate_id, low=low, high=high, howmany=howmany)


##########
# MAIN API
##########


# DIGEST AND BLAST AND SCORE
############################

# requires:
# blastdb database
# dict      a genome fasta file against which pam lookup is optionnaly done for guides candidates
# genome    same, used when (correctly) blasting agaist whole geome
#
# TODO get rid of reref_substrate_id
def digest_and_blast_and_score_coord(coord, genome='mm8', blastdb='mm8', direct='./', max_hsps=10, chunk_size=20,
                                     reref_substrate_id=None, low=75, high=75, load_genome=False, howmany=None,
                                     restriction_enzymes=[u'BfaI', u'ScrFI', u'HpaII']):

    filename = digest_coord(coord, genome, direct, restriction_enzymes)

    nbr = blast(filename + '.prsp', blastdb, direct, max_hsps=max_hsps, chunk_size=chunk_size)

    guides = score(nbrofchunks=nbr, filename=filename + '.prsp', blastdb=blastdb, direct=direct,
                   chunk_size=chunk_size,
                   reref_substrate_id='chr6', load_genome=load_genome)

    guides, groups = cluster(guides, filename + '.prsp', direct, reref_substrate_id, low, high, howmany)

    return guides, groups
