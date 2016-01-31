from __future__ import print_function

import os

#BLASTDB = os.environ.get('BLASTDBD','/genomes/blast/')

SCAFFOLDS = os.environ.get('SCAFFOLDS','/genomes/scaffolds/')

SCRATCH = os.environ.get('SCRATCH','/data/scratch/')

APPSESSIONJSON = os.environ.get('APPSESSIONJSON', '/data/input/AppSession.json')

RESTRICTION_ENZYMES_LIST = ['BfaI', 'HpaII', 'ScrFI']

DIGEST_LOG_ON = True
BLAST_LOG_ON = True
SCORE_LOG_ON = False
CLUSTER_LOG_ON = True
PRIME_LOG_ON = True


def digestlog(*args, **kwargs):
    if DIGEST_LOG_ON:
        print(*args, **kwargs)

def blastlog(*args, **kwargs):
    if BLAST_LOG_ON:
        print(*args, **kwargs)

def scorelog(*args, **kwargs):
    if SCORE_LOG_ON:
        print(*args, **kwargs)

def clusterlog(*args, **kwargs):
    if CLUSTER_LOG_ON:
        print(*args, **kwargs)

def primelog(*args, **kwargs):
    if PRIME_LOG_ON:
        print(*args, **kwargs)


# Default parameters for primer3 primer design
# You ccan supply your own by copying this and modifyin
# These defaults were intended for NEB's Q5 HotStart polymerase
PRIMER3_PARAMETERS = {
    'PRIMER_NUM_RETURN': 20,
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_INTERNAL_MAX_SELF_END': 8,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 50.0, # was 56
    'PRIMER_MAX_TM': 75.0, # was 67
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_INTERNAL_MAX_POLY_X': 100,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_WT_SELF_ANY': 1, # added
    'PRIMER_PAIR_WT_COMPL_ANY':1, #added
    'PRIMER_PAIR_WT_COMPL_END':1, #added
    'PRIMER_PAIR_WT_SELF_END':1, #added
    'PRIMER_PAIR_WT_SELF_ANY':1, #added
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,
    'PRIMER_PRODUCT_SIZE_RANGE': [26, 18000],
    'PRIMER_MIN_THREE_PRIME_DISTANCE':3,
    'PRIMER_PAIR_MAX_DIFF_TM':5, #added
    'PRIMER_PAIR_WT_DIFF_TM':1, #added
    'PRIMER_PAIR_WT_PRODUCT_SIZE_LT': 0.5, #added
    }

