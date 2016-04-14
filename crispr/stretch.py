########################################################################################################################
# STRETCH
#
# API functions:  stretch, run
#
########################################################################################################################

from __future__ import absolute_import, division, print_function   # , unicode_literals
import pickle
import matplotlib as mpl
# Agg is a backend allowing to create plots without a running X server
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import mean
import collections
from .cut import tabbed_string_from_list, savelinestofile
from .config import stretchlog
from .score import get_score, get_score_s, get_substrate_position_score_tuples, load_guides


# PRINT SCORES
##############

def print_scores_info(guides):
    stretchlog()
    stretchlog('distribution of scores among all', len(guides), 'guides (score, nb of guides with that score)')
    stretchlog(sorted(collections.Counter(get_score_s(guides)).items()))
    return guides

def print_scores_histo(direct, guides, fn_noext):
    scores = [guide.annotations['score'] for guide in guides]
    bins = range(0,106,5)
    plt.hist(scores, bins, color="gray")
    plt.tick_params(axis=u'both', labelsize=18)
    plt.savefig(direct + fn_noext + '.guides.pdf', format="pdf")


# STRETCHES
###########

def regroup(guides, high):
    substrate_position_score_tuples = get_substrate_position_score_tuples(guides)
    stretches=[]
    goodones = 0
    tryfrom = 0
    end = 0
    substr = 'none'
    for tuple in substrate_position_score_tuples:
        #stretchlog(tuple)
        # if this is a good guide
        if tuple[2] >= high:
            # only consider a good guide if it starts after the tryfrom poition
            if tuple[1] >= tryfrom:
                goodones += 1
                end = tuple[1] + 20
                tryfrom = tuple[1] + 21
                # if this is the first good, reset begin to this position
                if goodones == 1:
                    begin = tuple[1]
                    substr = tuple[0]
        # if this is not a good guide
        else:
            tryfrom = tuple[1] + 21  # TODO we could actually start sooner with an overlap of the bad guide
            # TODO fix with bedtools arithmetic, start with set of good guides not overlapping bad guides, and loop
            end = min(end,  tuple[1])
            # if we were in a good stretch until now, reset goodones to 0 and save the just finished good interval
            if goodones > 0:
                # TODO why does end==begin happen ?
                if end > begin:
                    stretches.append((substr, begin, end, end-begin, goodones, 1000*goodones/(end-begin)))
                goodones = 0

    stretches.sort(key=lambda x : -x[4])
    return stretches


# STRETCHES BED UTILITIES
#########################

def stretchbedline(stretch, index, howmanythreshold):
    substrate_id = stretch[0]
    bedstart = stretch[1]
    bedend = stretch[2]
    score = stretch[4]
    # TODO use gradual shades reflecting toatl good guide yield
    if index < howmanythreshold:
        color = '0,255,0'      # green
    else:
        color = '0,127,0'       # pale green
    cutinfo = [substrate_id, bedstart, bedend, 'stretch', score, '.', bedstart, bedend, color]
    bedline = tabbed_string_from_list(cutinfo)
    return bedline


def stretchbedlines(stretches, howmanythreshold):
    if howmanythreshold == 0:
        howmanythreshold = len(stretches) // 2
    return [stretchbedline(stretch, index, howmanythreshold) for index, stretch in enumerate(stretches)]


# STRETCHES IO
##############

def save_stretches(stretches, filename, directory, howmanythreshold):
    """
    save stretches
    :param guides:
    :param filename:
    :param directory:
    :return:
    """
    with  open(directory + filename + ".stretches.pkl", "w") as output_pkl:
        pickle.dump(stretches, output_pkl)

    savelinestofile(directory + filename + '.stretch.bed',
                    stretchbedlines(stretches, howmanythreshold)
                    )

    return stretches


def load_stretches(filename, directory):
    """
    load stretches
    :param filename:
    :param directory:
    :return:
    """
    with  open(directory + filename + ".stretches.pkl", "r") as input_pkl:
        stretches = pickle.load(input_pkl)
    return stretches


# SRETCHES YIELD INFO
#####################
# TODO add cumulative guide yield, cumulatibve bp length coverage for top n stretches,
# TODO create tabular data and export as CSV file
# TODO view in ipython ? or in IGV?
def print_stretches_yield_info(stretches, howmanythreshold=None):
    if not howmanythreshold:
        howmanythreshold = len(stretches) // 2
    stretchlog('')
    stretchlog('%s stretches totaling %s bps yielding a total of %s guides'
            % (len(stretches), sum([t[3] for t in stretches]), sum([t[4] for t in stretches]))
            )
    stretchlog('with average group of %s bps yielding %s guides'
            % (mean([t[3] for t in stretches])//1, mean([t[4] for t in stretches])//1)
            )
    stretchlog('')
    stretchlog('top %s stretches totaling %s bps yielding a total of %s guides'
            % (howmanythreshold, sum([t[3] for t in stretches[0:howmanythreshold]]), sum([t[4] for t in stretches[0:howmanythreshold]]))
            )
    stretchlog('with average group of %s bps yielding %s guides'
            % (mean([t[3] for t in stretches[0:howmanythreshold]])//1, mean([t[4] for t in stretches[0:howmanythreshold]])//1)
            )
    stretchlog('')
    for stretch in stretches[0:howmanythreshold]:
        stretchlog("%s:%s-%s (%s bps) yields %s guides (%s guides per 1000bps )" % stretch)
    stretchlog('')
    for stretch in stretches[howmanythreshold:]:
        stretchlog("%s:%s-%s (%s bps) yields %s guides (%s guides per 1000bps )" % stretch)
    stretchlog('')
    return stretches


###################
# MAIN API FUNCTION
###################

def stretch(filename, directory, low=75, high=94, howmanythreshold=0):
    """
    given a directory and filename, pickle-loads a list of guides in a region.
    This list must be: SeqRecord objects with :

    Then: plots scores info and scores histogram, find good guides stretches, rank stretches by guides yield
    ...
    :param filename:
    :param directory:
    :param low:
    :param high:
    :param howmany:
    :return:
    """
    guides = load_guides(filename, directory)

    stretchlog('\n*******************************************************************')
    stretchlog('STRETCHES:  SCORES INFO AND HISTOGRAM')
    stretchlog('*******************************************************************')
    print_scores_info(guides)
    print_scores_histo(directory, guides, filename)

    stretchlog('\n*******************************************************************')
    stretchlog('FIND GOOD GUIDES STRETCHES')
    stretchlog('*******************************************************************')
    stretches = regroup(guides, high)
    stretchlog(("substrate, start, end, length, good guides, good guides per thousand nucl."))
    for idx, stretch in enumerate(stretches):
        stretchlog(idx, stretch)

    stretchlog('\nSAVE STRETCHES')
    save_stretches(stretches, filename, directory, howmanythreshold)

    stretchlog('\nRANK STRETCHES BY GUIDES YIELD')
    print_stretches_yield_info(stretches, howmanythreshold)
    return guides, stretches


def run(theguides, filename, directory, threshold=94):
    """
    given a directory and filename, pickle-loads a list of guides in a region.
    This list must be: SeqRecord objects with :
            - their position within the region in their "name" attribute
            - their specificity score in their "annotations["score"]" attribute
            - their length obtained with the len function
    :param filename: filename to pickle-load from
    :param directory: working directory
    :param threshold: (default=94) the score cutoff below which a run will be "broken"
            Raising this leads to shorter runs (shorter PCR products) leading
            to fewer guides produced on average. Lowering it produces more PCR
            products, but of lower mean specificity.
    :return: tuple of:
    guides_sorted: list of guides sorted by positions
    runs: list of (start, end) indexes within guides_sorted list, of high-scoring guides
    which if PCR-amplified as sub-regions,will generate a set of highly-specific sgRNAs.
    """

    # guides = load_guides(filename, directory)
    guides = theguides

    stretchlog('\nPRINT GUIDES from stretch')
    stretchlog("============")
    stretchlog("len(guides)", len(guides))
    stretchlog("guides", guides)

    starts=[]
    ends=[]
    prev_item_good = False
    for index, guide in enumerate(guides):

        this_item_good = guide.annotations["score"] >= threshold and len(guide) > 19
        stretchlog("\n\n===\nindex", index, "\nguide\n", guide, "\nthis_item_good", this_item_good)

        # starting a new good guides strech
        if not prev_item_good and this_item_good:
            # store as a run start the index of the first good guide in the stretch
            starts.append(index)
            prev_item_good = True
        # ending a current good guide strech
        if prev_item_good and not this_item_good:
            # store as a run end the index of the first bad guide after the stretch
            ends.append(index)
            prev_item_good = False

    # only retain stretches that include at least 3 good guides
    # TODO why this criteria of at least 3 good guides ?

    stretchlog("zip(starts, ends)", zip(starts, ends))
    runs = [(start, end) for (start, end) in zip(starts, ends) if end - start > 3]

    # stretchlog('\nPRINT GUIDES')
    # for idx, guide in enumerate(guides):
    #     stretchlog(idx, guide.id, get_score(guide), "GOOD" if get_score(guide) >= threshold else "BAD")

    stretchlog('\nPRINT RUNS from stretch')
    stretchlog("============")
    stretchlog("len(runs)", len(runs))
    stretchlog("runs", runs)

    return guides, runs


# def count_non_overlapping_guides(guidelist, binding_interference_spacing=20):
#     '''
#     Sequence position information must be encoded in sequence name attribute,
#     e.g. name="100" indicates that left edge of guide (regardless of strand)
#     starts 100nt along the target scaffold.
#
#     Optional argument binding_interference_spacing specifies the number of
#     nucleotides apart that a guides must be to contribute to unique
#     labeling events (for example, of two guides only 10nt apart, it's impossible
#     for both to bind at once)
#
#     From the spCas9 crystal structure, it looks like guides will need to be
#     at least 24-25 nucleotides apart to bind simultaneously, and possibly
#     more.
#     '''
#     prev_position = 0
#     guidecount = 0
#     for item in guidelist:
#         distance_from_last = int(item.name) - prev_position
#         if distance_from_last > binding_interference_spacing:
#             guidecount = guidecount + 1
#         prev_position = int(item.name)
#     return guidecount


# see: http://sebastianraschka.com/Articles/2014_multiprocessing_intro.html#An-introduction-to-parallel-programming-using-Python's-multiprocessing-module
# import multiprocessing as mp
#
# def multi_protospacers_score(dir, blastdb_directory, blastdb_db):
#     pool = mp.Pool(processes=2)
#     guides = [pool.apply(protospacers_score, args=(x,)) for x in cutslist]
#
# def multiscore_pool(x):
#    score = al_scoreguide(x, "xl71", xl71genomedict)
#    return score


# missing the sequences! now implemented with pybedtools
# def save_to_fastafile(bedtuples, outputdirpath=''):
#     with open(outputdirpath + 'cuts_nameonly.fasta', 'w') as outfp:
#         idx = 2
#         for tuple in bedtuples:
#             # (substrate_id, bedstart, bedend, enzyme_name, score, sense, bbedstart, bbedend, color) = tuple
#             (substrate_id, bedstart, bedend, enzyme_name, score, sense) = tuple
#             id = str(idx//2) + '_' + ('F' if sense=='+' else 'R')
#             idx += 10
#             outfp.write(">lcl|" + substrate_id + "|" + enzyme_name + "|" + str(bedstart) + "|" + id +"\n")
#             # outfp.write(str(item.seq) + "\n")
#             outfp.write("\n")
