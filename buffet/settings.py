import os

BLASTDB = os.environ.get('BLASTDB','/genomes/blastdb/')

SCAFFOLDS = os.environ.get('SCAFFOLDS','/genomes/scaffolds/')

SCRATCH = os.environ.get('SCRATCH','/data/scratch/')

APPSESSIONJSON = os.environ.get('APPSESSION', '/data/input/AppSession.json')


RESTRICTION_ENZYMES_LIST = ['BfaI', 'HpaII', 'ScrFI']