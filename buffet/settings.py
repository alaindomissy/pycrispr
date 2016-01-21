import os

#BLASTDB = os.environ.get('BLASTDBD','/genomes/blastdb/')

SCAFFOLDS = os.environ.get('SCAFFOLDS','/genomes/scaffolds/')

SCRATCH = os.environ.get('SCRATCH','/data/scratch/')

APPSESSIONJSON = os.environ.get('APPSESSIONJSON', '/data/input/AppSession.json')

RESTRICTION_ENZYMES_LIST = ['BfaI', 'HpaII', 'ScrFI']