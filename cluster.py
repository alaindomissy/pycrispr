########################################################################################################################
# CLUSTER
#
# API functions:  cluster
#
########################################################################################################################


import matplotlib.pyplot as plt
from numpy import mean
import collections

from cutbedlines import tabbed_string_from_list


###################
# MAIN API FUNCTION
###################

def cluster(guides, direct, reref_substrate_id, low=50, high=75, howmany=None):
    guides_scorebedfile(direct, 'scoreguides.rel', guides, None,low, high)
    if reref_substrate_id:
        guides_scorebedfile(direct, 'scoreguides.abs', guides, reref_substrate_id, low, high)

    # TODO guides not getting a score, how does this happen ? fix it better
    # kind of fix guides without a score
    for guide in guides:
        try:
            guide.annotations['score']
        except:
            # print 'guide %s does not have a score' % guide
            guide.annotations['score'] = 0.0

    ampls = amplicons(pos_score(guides), high)
    print_ampls_info(ampls, howmany)
    print_scores_info(scores(guides))
    histo(direct, guides)

    return guides, ampls



def histo(direct, guides):
    scores = [guide.annotations['score'] for guide in guides]
    bins = range(0,106,5)
    plt.hist(scores, bins, color="gray")
    plt.tick_params(axis=u'both', labelsize=18)
    plt.savefig(direct + 'scores.pdf', format="pdf")



def guides_scorebedfile(direct, fn_noext,
                        guides, reref_substrate_id, low, high):
    scorebedfile(direct, fn_noext,
                 scorebedlines(guides, reref_substrate_id, low, high)
                 )

def scorebedfile(outputdirpath, fn_noext, bedlines):
    with open(outputdirpath + fn_noext + '.bed', 'w') as bedfn:
       for bedline in bedlines:
          bedfn.write(bedline)

def scorebedlines(guides, reref_substrate_id, low, high):
    return [scorebedline(guide, reref_substrate_id, low, high) for guide in guides]

def scorebedline(guide, reref_substrate_id,low, high):
    substrate_id = guide.id.split('-')[0]
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
    if int(score) <= low:
        color = '255,0,0'       # red
    elif int(score) >= high:
        color = '0,255,0'       # green
    else:
        color = '0,0,255'       # blue
    cutinfo = [substrate_id, bedstart, bedend, 'crRNA', score, sense, bedstart, bedend, color]
    # cutinfo = [substrate_id, bedstart, bedend, 'crspr_guide', 'score', sense]
    bedline = tabbed_string_from_list(cutinfo)
    return bedline



###
# def index_score(guides):
#     return[(guide.annotations['blastindex'], guide.annotations['score'])
#            for guide in guides]

def scores(guides):
    return [guide.annotations['score'] for guide in guides]

def pos_score(guides):
    return[(int(guide.name.split(':')[1].split('-')[0]) , guide.annotations['score'])
           for guide in guides]

def amplicons(pos_score_tuples, high=90):
    ampls=[]
    goodones = 0
    for tuple in pos_score_tuples:
        #print tuple
        if tuple[1]>=high:
            goodones+=1
            # if this is the first good, reset begin to this position
            if goodones==1:
                begin = tuple[0]
        else:
            # if we were in a good stretch until now, reset goodones to 0
            if goodones>0:
                end = tuple[0] + 20
                ampls.append((begin, end, end-begin, goodones, 1000*goodones/(end-begin)))
                goodones = 0
    ampls.sort(key=lambda x : -x[3])
    return ampls


def print_scores_info(scores):
    print
    print 'distribution of scores among all', len(scores), 'guides (score, nb of guides with that score)'
    print sorted(collections.Counter(scores).items())


def print_ampls_info(ampls, howmany=None):
    if not howmany:
        howmany = len(ampls) // 2
    print('')
    print('%s amplicons totaling %s bps yielding a total of %s guides'
            % (len(ampls), sum([t[2] for t in ampls]), sum([t[3] for t in ampls]))
            )
    print('with average amplicon of %s bps yielding %s guides'
            % (mean([t[2] for t in ampls])//1, mean([t[3] for t in ampls])//1)
            )
    print('')
    print('top %s amplicons totaling %s bps yielding a total of %s guides'
            % (howmany, sum([t[2] for t in ampls[0:howmany]]), sum([t[3] for t in ampls[0:howmany]]))
            )
    print('with average amplicon of %s bps yielding %s guides'
            % (mean([t[2] for t in ampls[0:howmany]])//1, mean([t[3] for t in ampls[0:howmany]])//1)
            )
    print('')
    for ampl in ampls[0:howmany]:
        print("%s - %s (%s bps) yields %s guides (%s guides per 1000bps )" % ampl)
    print('')
    for ampl in ampls[howmany:2*howmany]:
        print("%s - %s (%s bps) yields %s guides (%s guides per 1000bps )" % ampl)
    print('')

    return ampls
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