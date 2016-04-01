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
from basespaceapp.app import main as appmain
from .config import SCRATCH, restore_genome, ARGUMENTS_WITH_CONTENT, ARGUMENTS_WITH_ITEMS  # ,\
    # DIGEST_LOG_ON, BLAST_LOG_ON, SCORE_LOG_ON, STRETCH_LOG_ON, PRIME_LOG_ON
from .main import digest_blast_score_cluster_prime


###################
# MAIN API FUNCTION
###################

def payload(params_value, output_dir):

    print('params_value', params_value)
    print('output_dir', output_dir)

    coord = str(params_value.get('input.genomic_coord'))
    genome = str(params_value['input.refgenome_id'])

    restriction_enzymes =  params_value['input.restr_enzymes']

    chunk_size = 10                  #int(params_value['input.blast_chunk_size'])
    max_hsps = 10                    # int(params_value['input.blast_max_hsps'])
    low_threshold = params_value['input.low_threshold']     #75
    high_threshold = params_value['input.high_threshold']   #94


    # interference_gap = str(params_value['input.interference_gap'])
    # primer_screening = params_value['input.primer_screening']

    checkbox_logging = params_value['input.checkbox_logging']
    DIGEST_LOG_ON = 'digest_log_on' in checkbox_logging
    BLAST_LOG_ON = 'blast_log_on' in checkbox_logging
    SCORE_LOG_ON = 'score_log_on' in checkbox_logging
    STRETCH_LOG_ON = 'stretch_log_on' in checkbox_logging
    PRIME_LOG_ON = 'prime_log_on' in checkbox_logging



    restore_genome(genome)
    digest_blast_score_cluster_prime(coord, genome, SCRATCH,
                                     max_hsps, chunk_size, low_threshold, high_threshold,
                                     load_genome=False, howmany=None,
                                     restriction_enzymes=restriction_enzymes)
    datetimenow = datetime.now().isoformat('_')
    # sessiondetails_dir = output_dir + '../sessiondetails/'
    sessiondetails_dir = output_dir + '../sessiondetails_' + datetimenow + '/'
    print('sessiondetails_dir', sessiondetails_dir)
    # copytree(source, destination, ignore=_logpath)
    # TODO why copytree is not working from here? currently it is called from basespaceapp.app.process_appsession
    # copytree(SCRATCH, sessiondetails_dir)

    return "success: results *NOT* copied to output directory: " + sessiondetails_dir + "at datetime: " + datetimenow


# to call from basespace
##############################
# commandLine: ["python", "-m", "crispr.payload.main", "/data/", "payload"],
# commandLine: ["python", "-m", "crispr.payload"],


# this file executed as script
##############################
if __name__ == '__main__':
    appmain('/data/', payload, ARGUMENTS_WITH_CONTENT, ARGUMENTS_WITH_ITEMS)
