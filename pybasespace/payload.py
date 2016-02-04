########################################################################################################################
#
# PAYLOAD
#
# API functions: payload
#
########################################################################################################################

from shutil import copytree
from datetime import datetime
from crispr.config import SCRATCH, restore_genome
from crispr.main import digest_and_blast_and_score_coord


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
    # interference_spacing = str(params_value['input.interference_spacing'])


    restore_genome(genome)

    digest_and_blast_and_score_coord(coord, genome, direct=SCRATCH, max_hsps=max_hsps, chunk_size=chunk_size,
                                     low=75, high=75, load_genome=False, howmany=None,
                                     restriction_enzymes=restriction_enzymes)

    # copytree(source, destination, ignore=_logpath)
    copytree(SCRATCH, output_dir + '../sessiondetails_' + datetime.now().isoformat('_'))

    return
