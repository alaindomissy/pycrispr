########################################################################################################################
#
# PAYLOAD
#
# API functions: payload
#
########################################################################################################################

from __future__ import absolute_import, division, print_function   # , unicode_literals
from shutil import copytree
from datetime import datetime
from basespaceapp.app import main
from .config import SCRATCH, restore_genome, ARGUMENTS_WITH_CONTENT, ARGUMENTS_WITH_ITEMS
from .main import digest_blast_score_cluster_prime


###################
# MAIN API FUNCTION
###################

def payload(params_value, output_dir):

    coord = str(params_value.get('input.genomic_coord'))
    assert(coord)

    genome = str(params_value['input.refgenome_id'])
    restriction_enzymes =  params_value['input.restr_enzymes']
    chunk_size = int(params_value['input.blast_chunk_size'])
    max_hsps = int(params_value['input.blast_max_hsps'])
    # interference_gap = str(params_value['input.interference_gap'])


    restore_genome(genome)

    digest_blast_score_cluster_prime(coord, genome, direct=SCRATCH, max_hsps=max_hsps, chunk_size=chunk_size,
                                     low=75, high=75, load_genome=False, howmany=None,
                                     restriction_enzymes=restriction_enzymes)

    datetimenow = datetime.now().isoformat('_')
    # copytree(source, destination, ignore=_logpath)
    # copytree(SCRATCH, output_dir + '../sessiondetails_' + datetimenow)
    copytree(SCRATCH, output_dir + '../sessiondetails/')

    return "success _ results copied to putput directory at " + datetimenow


# to call from basespace
##############################
# commandLine: ["python", "-m", "crispr.payload.main", "/data/", "payload"],
# commandLine: ["python", "-m", "crispr.payload"],


# this file executed as script
##############################
if __name__ == '__main__':
    main('/data/', payload, ARGUMENTS_WITH_CONTENT, ARGUMENTS_WITH_ITEMS)
