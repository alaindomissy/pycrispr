from __future__ import print_function

import os

#BLASTDB = os.environ.get('BLASTDBD','/genomes/blastdb/')

SCAFFOLDS = os.environ.get('SCAFFOLDS','/genomes/scaffolds/')

SCRATCH = os.environ.get('SCRATCH','/data/scratch/')

APPSESSIONJSON = os.environ.get('APPSESSIONJSON', '/data/input/AppSession.json')

RESTRICTION_ENZYMES_LIST = ['BfaI', 'HpaII', 'ScrFI']

SCORING_LOG_ON = False

def scorelog(*args, **kwargs):
    if SCORING_LOG_ON:
        print(*args, **kwargs)