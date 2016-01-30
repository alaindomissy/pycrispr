import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from crispr.zhang import *


def test_mean_pairwise_distance_between_mismatches():
    assert(mean_pairwise_distance_between_mismatches(np.array([0,0,1,0,0] + 15*[0])) == 19.0)
    assert(mean_pairwise_distance_between_mismatches(np.array([1,0,1,0,1]+ 15*[0])) == 2.0)
    assert(mean_pairwise_distance_between_mismatches(np.array([1,1,0,0,1]+ 15*[0])) == 2.0)
    assert(mean_pairwise_distance_between_mismatches(np.array([0,1,1,1,1]+ 15*[0])) == 1.0)
    assert(mean_pairwise_distance_between_mismatches(np.array([1,0,0,0,0,0,0,0,0,0,0,0,1]+ 7*[0])) == 12.0)
    assert(mean_pairwise_distance_between_mismatches(np.array([1,1,0,0,0,0,0,0,0,1,0,0,1]+ 7*[0])) == 4.0)
    assert(mean_pairwise_distance_between_mismatches(np.zeros(20))==19.0)
    assert(mean_pairwise_distance_between_mismatches(np.ones(20))==1.0)


def test_zhangscore():
    assert( zhangscore(' '*20 ) == 8.6097000381855906e-08 )
    assert( zhangscore('|'*20 ) == 100.0 )

def test_ontarget_likelyhoods_from_pam():
    mm = np.zeros(20)
    mm[0:3]=1
    # assert( np.prod(ontarget_likelyhoods_from_pam ** mm) == np.prod(ontarget_likelyhoods_from_pam[:3]) )