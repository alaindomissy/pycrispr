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
from .config import SCRATCH, restore_genome
from .main import digest_blast_score_cluster_prime


ARGUMENTS_WITH_CONTENT = [
                            # 'input.interference_gap',
                            # 'input.blast_chunk_size',
                            # 'input.blast_max_hsps',
                            'input.refgenome_id',
                            'input.genomic_coord',
                            'input.low_threshold',
                            'input.high_threshold',
                            'input.primer_screening'
                         ]

ARGUMENTS_WITH_ITEMS = ['input.restr_enzymes', 'input.checkbox_logging']


###################
# MAIN API FUNCTION
###################

def payload(params_value, output_dir):

    # params_value = params_value or dict(
    #     {'input.samples': [],
    #       u'input.primer_screening': u'dumb',
    #       'input.projects': [{'href': u'v1pre3/projects/29085057', 'id': u'29085057', 'name': u'crispr-ecoli'}],
    #       u'input.restr_enzymes': [u'BfaI', u'ScrFI', u'HpaII'],
    #       u'input.checkbox_logging': [u'stretch_log_on', u'prime_log_on'],
    #       'output.projects': [{'href': u'v1pre3/projects/29085057', 'id': u'29085057', 'name': u'crispr-ecoli'}],
    #       u'input.refgenome_id': u'ecoli',
    #       u'input.genomic_coord': u'chr:102000-106999',
    #       u'input.high_threshold': u'99',
    #       u'input.low_threshold': u'75'
    #     })
    #
    # output_dir = output_dir or '/data/'

    print('params_value', params_value)
    print('output_dir', output_dir)

    coord = str(params_value.get('input.genomic_coord'))
    genome = str(params_value['input.refgenome_id'])

    chunk_size = 10                  #int(params_value['input.blast_chunk_size'])
    max_hsps = 10                    # int(params_value['input.blast_max_hsps'])
    low_threshold = params_value['input.low_threshold']     #75
    high_threshold = params_value['input.high_threshold']   #94
    # interference_gap = str(params_value['input.interference_gap'])
    # primer_screening = params_value['input.primer_screening']

    method = "dumb"
    tm = 40

    restriction_enzymes =  params_value['input.restr_enzymes']

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
                                     restriction_enzymes=restriction_enzymes,
                                     method=method, tm=tm)
    datetimenow = datetime.now().isoformat('_')
    # sessiondetails_dir = output_dir + '../sessiondetails/'
    sessiondetails_dir = output_dir + '../sessiondetails_' + datetimenow + '/'
    print('sessiondetails_dir', sessiondetails_dir)
    # copytree(source, destination, ignore=_logpath)    # TODO should we use the thrd parameterto copytree
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
