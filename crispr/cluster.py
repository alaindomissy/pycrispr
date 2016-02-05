########################################################################################################################
# CLUSTER
#
# API functions:  cluster
#
########################################################################################################################

from __future__ import print_function
import matplotlib as mpl
# Agg is a backend allowing to create plots without a running X server
mpl.use('Agg')
import matplotlib.pyplot as plt
from numpy import mean
import collections
from cut import tabbed_string_from_list
from config import clusterlog

###################
# MAIN API FUNCTION
###################

def cluster(guides, filename, direct, low=75, high=75, howmany=12):

    filename = filename + '.prsp'

    print('\nSAVE SCORES')
    bedlines_s = scorebedlines(guides, low, high)
    savelinestofile(direct + filename + '.scores.bed', bedlines_s)
    # if reref_substrate_id:
    #     guides_scorebedfile(direct, filename + '.scores.abs', guides, reref_substrate_id, low, high)


    print('\nPLOT SCORES')
    print_scores_info(scores(guides))
    histo(direct, guides, filename)


    print('\nCLUSTER GOOD GUIDES')
    groups = regroup(substr_pos_score(guides), high)
    bedlines_c = groupbedlines(groups, howmany)
    savelinestofile(direct + filename + '.groups.bed', bedlines_c)

    print('\nRANK CLUSTERS')
    print_groups_info(groups, howmany)

    return guides, groups





# HISTOGRAM
############

def histo(direct, guides, fn_noext):
    scores = [guide.annotations['score'] for guide in guides]
    bins = range(0,106,5)
    plt.hist(scores, bins, color="gray")
    plt.tick_params(axis=u'both', labelsize=18)
    plt.savefig(direct + fn_noext + '.scores.pdf', format="pdf")





# BED FILES
############

def scorebedline(guide, low, high):
    substrate_id = guide.id.split(':')[0]
    bedstart = guide.id.split(':')[1].split('-')[0]
    # rename scaffold field and make coord 'absolute' (ie based on higher level)
    # if reref_substrate_id:
    #     substrate_id = reref_substrate_id
    #     bedstart = str(int(bedstart) + int(guide.id.split(':')[1].split('-')[0]))
    bedend = str(int(bedstart) + 20)
    sense = guide.id[-2]
    try:
        score = str(guide.annotations['score'])
    except:
        score ='0'
    # TODO use gradual shades reflecting scores
    # suddenly got this error invalid literal for int() with base 10:0.0 Now changed int() to float() But why ?
    if float(score) < low:
        color = '255,0,0'       # red
    elif int(score) >= high:
        color = '0,255,0'       # green
    else:
        color = '0,0,255'       # blue
    cutinfo = [substrate_id, bedstart, bedend, 'crRNA', score, sense, bedstart, bedend, color]
    # cutinfo = [substrate_id, bedstart, bedend, 'crspr_guide', 'score', sense]
    bedline = tabbed_string_from_list(cutinfo)
    return bedline


def scorebedlines(guides, low, high):
    return [scorebedline(guide, low, high) for guide in guides]


def savelinestofile(path, lines):
    with open(path, 'w') as handle:
       for line in lines:
          handle.write(line)






# GROUPS
########

def regroup(substr_pos_score_tuples, high=90):
    groups=[]
    goodones = 0
    tryfrom = 0
    end = 0
    substr = 'none'
    for tuple in substr_pos_score_tuples:
        #print(tuple)
        if tuple[2]>=high:
            # only consider a good guide if it starts after the tryfrom poition
            if tuple[1] >= tryfrom:
                goodones+=1
                end = tuple[1] + 20
                tryfrom = tuple[1] + 21
                # if this is the first good, reset begin to this position
                if goodones==1:
                    begin = tuple[1]
                    substr = tuple[0]
        else:
            tryfrom = tuple[1] + 21
            # TODO fix with bedtools arithmetic, start with set of good guides not overlapping bad guides, and loop
            end = min(end,  tuple[1])
            # if we were in a good stretch until now, reset goodones to 0 and save the just finished good interval
            if goodones>0:
                # TODO why does end==begin happen ?
                if end>begin:
                    groups.append((substr, begin, end, end-begin, goodones, 1000*goodones/(end-begin)))
                goodones = 0

    groups.sort(key=lambda x : -x[4])
    return groups


def groupbedline(group, howmany, index):
    substrate_id = group[0]
    bedstart = group[1]
    bedend = group[2]
    score = group[4]
    # TODO use gradual shades reflecting toatl good guide yield
    if index < howmany:
        color = '0,255,0'      # green
    else:
        color = '0,127,0'       # pale green
    cutinfo = [substrate_id, bedstart, bedend, 'group', score, '.', bedstart, bedend, color]
    bedline = tabbed_string_from_list(cutinfo)
    return bedline


def groupbedlines(groups, howmany=None):
    if not howmany:
        howmany = len(groups) // 2
    return [groupbedline(group, howmany, index) for index, group in enumerate(groups)]






# GUIDES GET INFO UTILITIES
###########################

def scores(guides):
    return [guide.annotations['score'] for guide in guides]

def pos_score(guides):
    return[(int(guide.name.split(':')[1].split('-')[0]) , guide.annotations['score'])
           for guide in guides]

def substr_pos_score(guides):
    return[(guide.id.split(':')[0], int(guide.name.split(':')[1].split('-')[0]) , guide.annotations['score'])
           for guide in guides]

# def index_score(guides):
#     return[(guide.annotations['blastindex'], guide.annotations['score'])
#            for guide in guides]




# GROUPS INFO
###############
# TODO add cumulative guide yield, cumulatibve bp length coverage for top n groups,
# TODO create tabular data and export as CSV file
# TODO view in ipython ? or in IGV?
def print_groups_info(groups, howmany=None):
    if not howmany:
        howmany = len(groups) // 2
    print('')
    print('%s groups totaling %s bps yielding a total of %s guides'
            % (len(groups), sum([t[3] for t in groups]), sum([t[4] for t in groups]))
            )
    print('with average group of %s bps yielding %s guides'
            % (mean([t[3] for t in groups])//1, mean([t[4] for t in groups])//1)
            )
    print('')
    print('top %s groups totaling %s bps yielding a total of %s guides'
            % (howmany, sum([t[3] for t in groups[0:howmany]]), sum([t[4] for t in groups[0:howmany]]))
            )
    print('with average group of %s bps yielding %s guides'
            % (mean([t[3] for t in groups[0:howmany]])//1, mean([t[4] for t in groups[0:howmany]])//1)
            )
    print('')
    for group in groups[0:howmany]:
        print("%s:%s-%s (%s bps) yields %s guides (%s guides per 1000bps )" % group)
    print('')
    for group in groups[howmany:]:
        print("%s:%s-%s (%s bps) yields %s guides (%s guides per 1000bps )" % group)
    print('')
    return groups


# SCORES INFO
#############
def print_scores_info(scores):
    print()
    print('distribution of scores among all', len(scores), 'guides (score, nb of guides with that score)')
    print(sorted(collections.Counter(scores).items()))
    return scores



#
#
#
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




###
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
###





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
30