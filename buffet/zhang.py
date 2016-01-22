##############################################################
# Zhang lab scoring formula

# see:
# DNA targeting specificity of RNA-guided Cas9 nucleases
# Patrick D Hsu,
# http://www.nature.com/nbt/journal/v31/n9/full/nbt.2647.html
###############################################################



# this module replaced the following code :
###########################################
# # make a list of mistmatching positions, 1based reverse indexed, starting from base next to Pam
# mmloc = []
# pos = 20
# for linestring in mmstr:
#     if linestring != "|":
#         mmloc.append(pos)
#         pos = pos - 1
# mmscore = [21 -x for x in mmloc]

# # Actually implement Zhang lab algorithm
# WEIGHT = [.000, .000, .014, .000, .000, .395, .317, .000,.389, .079, .445, .508, .613, .851, .732, .828, .615, .804, .685, .583]
# t1 = 1
# for mismatchlocation in mmscore:
#     t1 = t1 * (1.0 - WEIGHT[mismatchlocation - 1])
# if len(mmscore) > 1:
#     d = (float(max(mmscore)) - float(min(mmscore))) / float((len(mmscore) - 1))
# else:
#     d = 19
# t2 = 1 / ((((19.0 - d)/19.0) * 4) + 1)
# t3 = float(1)/ float(pow(len(mmscore), 2))
# match_score = 100.0 * t1 * t2 * t3


import numpy as np

off_target_likelyhoods  = np.array([.000, .000, .014, .000, .000, .395, .317, .000, .389, .079,
                                    .445, .508, .613, .851, .732, .828, .615, .804, .685, .583])
ontarget_likelyhoods_from_pam = 1 - off_target_likelyhoods[::-1]

# import matplotlib.pyplot as plt
# plt.plot(ontarget_likelyhoods_from_pam)
# plt.show()


def effect1(np_mismatches):
    """
    experimentally-determined effect of each mismatch position on targeting
    (relatively high disturbance of mismatches falling close to the PAM site)
    :param np_mismatches:
    :return:
    """
    return np.prod(ontarget_likelyhoods_from_pam ** np_mismatches)


def mean_pairwise_distance_between_mismatches(mismatches):
    """
    :param np_mismatches: 20-long list of 0 or 1 representing mismatches for a guide at a possible offtarget location
    :return: mean pairwise distance between mismatches
    """
    assert(len(mismatches)==20)
    d = len(mismatches) -1.0
    mmpos = np.where(mismatches==1.0)[0]
    nbr_of_pos = len(mmpos)
    if nbr_of_pos >=2:
        first_pos = mmpos[0]
        last_pos = mmpos[-1]
        d = float(last_pos - first_pos) / (nbr_of_pos - 1)
    return d


def effect2(np_mismatches):
    """
    effect of mean pairwise distance between mismatches (d), density of mistmatches reduce binding
    :param np_mismatches:
    :return:
    """
    d = mean_pairwise_distance_between_mismatches(np_mismatches)
    return 1 / ((((19.0 - d)/19.0) * 4) + 1)


def effect3(np_mismatches):
    """
    effect of total number of mismatches (dampening penalty for highly mismatched targets)
    :param mismatches:
    :return:
    """
    nbr_of_pos = len(np.where(np_mismatches==1.0)[0])
    return 1.0 if nbr_of_pos==0 else 1.0 / (nbr_of_pos ** 2)



def single_offtarget_score(mismatches):
    """
    :param mismatches: 20-long list of 0 or 1 representing mismatches for a guide at a possible offtarget location
    :return: a score between 0.0 and 100.0
    """
    np_mismatches = np.array(mismatches)
    t1 = effect1(np_mismatches)
    t2 = effect2(np_mismatches)
    t3 = effect3(np_mismatches)
    score = 100.0 * t1 * t2 * t3
    # print('t1:%5.1f t2:%5.1f t3:%5.1f score:%5.1f' % (t1*100, t2*100, t3*100, score))
    return 100.0 * t1 * t2 * t3


###################
# MAIN API FUNCTION
###################

def zhangscore(matchbars):
    """
    :param macthbars: a 20-long string of bars(|) and spaces ( ) for an alignement match resp. mismatch positions
    :return: a list of 0 and 1 representing the same match resp. mismatch positions
    """
    assert(len(matchbars)==20)
    mismatches = map(lambda char : 0 if char=='|' else 1, matchbars)
    return single_offtarget_score(mismatches)



