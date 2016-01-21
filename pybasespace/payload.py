
from shutil import copytree
# import logging
#
# def _logpath(path, names):
#     logging.info('Working in %s' % path)
#     return []   # nothing will be ignored



from buffet.settings import SCAFFOLDS, SCRATCH
# from buffet.cut import cut_file
# from buffet.digest import digest_fastafile, digest_focused, digest_coord
# from buffet.blast import blast
from buffet.main import digest_and_blast_coord, digest_and_blast_and_score_coord


# main.blast(SCAFFOLDS, 'chr6:136640001-136641000_1000.prsp', 'mm8')
# from buffet import cut

# cut.cut_file(SCAFFOLDS + 'chr6-1-400' + 'fasta')
#
# from buffet import digest
# digest.digest_file(SCAFFOLDS + 'hg38' + 'fasta')
# digest.digest_coord(SCAFFOLDS +'chr21:42774475-42905100' + 'fasta' ,'mm8', SCAFFOLDS + 'mm8.fasta )


def payload(params_value, output_dir):

    coord = str(params_value.get('input.genomic_coord'))

    if not coord:
        print("no input.genomic_coord input")
        return "no input.genomic_coord input"

    genome = str(params_value['input.genome_id'])
    chunk_size = int(params_value['input.blast_chunk_size'])
    max_hsps = int(params_value['input.blast_max_hsps'])
    binding_interference_spacing = str(params_value['input.binding_interference_spacing'])

    filepath = str(SCAFFOLDS + genome + '/' + coord + '.fasta')

    # cuts = cut_file(filepath)
    # return cuts

    # digest_fastafile(filepath)

    # filename_noext = 'chr6+47599949-47640339_40391'
    # digest_focused(str(SCAFFOLDS + genome + '/' + filename_noext),
    #                str(SCAFFOLDS + genome + '/' + genome + '.fasta')
    #                )


    # coord = 'chr6:47599949-47640339'

    reference = str(SCAFFOLDS + genome + '/' + genome + '.fasta')

    # focusfn =  digest_coord(SCRATCH, str(coord), reference)
    # nbr = blast(SCRATCH, focusfn + '.prsp', genome, chunk_size=chunk_size, max_hsps=max_hsps)

    # digest_and_blast_coord(SCRATCH, coord, reference, genome, chunk_size=chunk_size, max_hsps=max_hsps)

    digest_and_blast_and_score_coord(SCRATCH, coord, reference, genome, chunk_size=chunk_size, max_hsps=max_hsps,
                                     reref_substrate_id='chr6',
                                     low=75, high=75, load_genome=False, howmany=None)


    # copytree(source, destination, ignore=_logpath)
    copytree(SCRATCH, output_dir)

    return