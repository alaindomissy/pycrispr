import os

from buffet.cut import cut_file
from buffet.digest import digest_fastafile, digest_focused, digest_coord


BLASTDB = os.environ.get('BLASTDB','/genomes/blastdb/')

SCAFFOLDS = os.environ.get('SCAFFOLDS','/genomes/scaffolds/')

SCRATCH = os.environ.get('SCRATCH','/data/scratch/')




#
# main.blast(SCAFFOLDS, 'chr6:136640001-136641000_1000.prsp', '/BLASTDB/', 'mm8')
# from buffet import cut

# cut.cut_file(SCAFFOLDS + 'chr6-1-400' + 'fasta')
#
# from buffet import digest
# digest.digest_file(SCAFFOLDS + 'hg38' + 'fasta')
# digest.digest_coord(SCAFFOLDS +'chr21:42774475-42905100' + 'fasta' ,'mm8', SCAFFOLDS + 'mm8.fasta )


def payload(params_value):

    coord = params_value.get('input.genomic_coord')

    if not coord:
        return "no input.genomic_coord input"

    genome = params_value['input.genome_id']
    chunk_size = params_value['input.blast_chunk_size']
    max_hsps = params_value['input.blast_max_hsps']
    binding_interference_spacing = params_value['input.binding_interference_spacing']

    filepath = str(SCAFFOLDS + genome + '/' + coord + '.fasta')

    # cuts = cut_file(filepath)

    # digest_fastafile(filepath)

    # filename = 'chr6+47599949-47640339_40391'
    # digest_focused(str(SCAFFOLDS + genome + '/' + filename),
    #                str(SCAFFOLDS + genome + '/' + genome + '.fasta')
    #                )


    # coord = 'chr6:47599949-47640339'
    digest_coord(SCRATCH, str(coord),
                 str(SCAFFOLDS + genome + '/' + genome + '.fasta')
                 )

    # return cuts
    return
