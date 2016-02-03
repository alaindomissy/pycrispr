########################################################################################################################
# BLAST
#
# API functions:  blast
#
########################################################################################################################

from __future__ import print_function
import os
try:                                     # Python 2
    from itertools import izip_longest
except ImportError:                      # Python 3
    from itertools import zip_longest as izip_longest
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Application import ApplicationError
from Bio import SeqIO as seqio    # TODO this is also imported from cut.py, ok?
from config import blastlog, blastdb_path


# SINGLE FILE BLASTING
######################

def blast1(prspfilename, genome, direct, max_hsps):

    blastdb = blastdb_path(genome)

    # TODO max)hsps stopped working, now replaced with max_hsps_per_subject , but what is going on ?
    blastn_cline = NcbiblastnCommandline(
        query=direct + prspfilename + '.fasta',
        out=direct + prspfilename + '.blast',
        outfmt=5,
        db=blastdb,
        max_target_seqs=25,
        num_threads=4,
        evalue=10,
        max_hsps_per_subject=max_hsps,
        task="blastn-short",
        dust="no",
        )
    # print('\nblastn_cline: %s' % blastn_cline)
    try:
        stdout, stderr = blastn_cline()
    except ApplicationError as err:
        # print('Bio:Application:ApplicationError number {0} : {1}'.format(err.errno, err.strerror), end='')
        # print('Bio:Application:ApplicationError: ', err.args)
        print('*********')
        print(str(err).split('message ')[1].strip('\''))
        print('*********')
        #print(err.args)
        return True
    print("...done")
    return False


# ITERATION UTILITY
###################

# def grouper(iterable, chunk_size):
#     args = [iter(iterable)] * chunk_size
#     return izip(*args)
def grouper_longest(iterable, chunk_size, fillvalue=None):
    args = [iter(iterable)] * chunk_size
    return izip_longest(*args, fillvalue=fillvalue)


###################
# MAIN API FUNCTION
###################

def blast(filename, genome, direct, max_hsps, chunk_size):
    """
    :param direct:
    :param filename:
    :param blastdb:
    :param chunk_size:
    :param max_hsps:
    :return:
    """
    filename = filename + '.prsp'


    blastlog("\nBLAST PROTOSPACERS\n")
    blastlog("\nload protospacers from", filename, end=' ')
    seqs = list(seqio.parse(direct + filename + '.fasta', 'fasta'))
    blastlog('done')
    nbr_of_prsps = len(seqs)
    nbr_of_chunks = 1+ len(seqs) // chunk_size
    blastlog("%s protospacers in total - will blast %s chunks of %s prsps each\n" % (nbr_of_prsps, nbr_of_chunks, chunk_size))
    nbr = 0
    nbrwrong = 0
    nameformat = '%s.%sseqs.%s'     # TODO add '.max%sHSPs.' % max_hsps
    for seqs_chunk in grouper_longest(seqs, chunk_size, None):
        # turn seqs_chunk into a list, not including the None fill-values
        seqs_chunk = [seq for seq in seqs_chunk if seq]
        nbr += 1
        # fn_code_noext = filename + + 'hsps.'+  str(chunk_size) + 'x'  + str(nbr)

        fn_code_noext = nameformat % (filename, chunk_size, nbr)
        rightfasta = direct + fn_code_noext + '.fasta'
        rightblast = direct + fn_code_noext + '.blast'
        wrongfasta = direct + fn_code_noext + '.fasta.FAILED'
        wrongblast = direct + fn_code_noext + '.blast.FAILED'
        rescuedfasta = direct + fn_code_noext + '.fasta.FAILED.RESCUED'
        rescuedblast = direct + fn_code_noext + '.blast.FAILED.RESCUED'
        if os.path.isfile(rightfasta) and os.path.isfile(rightblast):
            blastlog('> skip chunk ', str(nbr).zfill(3), fn_code_noext)
            continue
        elif os.path.isfile(wrongfasta) or os.path.isfile(wrongblast):
            blastlog('> rescue chunk ', str(nbr).zfill(3), fn_code_noext, end=' ')
        else:
            blastlog('> blast chunk ', str(nbr).zfill(3), fn_code_noext, end=' ')
        with open(direct + fn_code_noext + '.fasta', "w") as temp_hndl:
            seqio.write(seqs_chunk, temp_hndl, "fasta")
        ########################################################
        failed = blast1(fn_code_noext, genome, direct, max_hsps)
        ########################################################
        if failed:
            nbrwrong +=1
            print('rightfasta',': ',rightfasta)
            print('wrongfasta',': ',wrongfasta)
            print('rightblast',': ',rightblast)
            print('wrongblast',': ',wrongblast)
            os.rename(rightfasta, wrongfasta)
            os.rename(rightblast, wrongblast)
        else:
            if os.path.isfile(wrongfasta):
                 os.rename(wrongfasta, rescuedfasta)
            if os.path.isfile(wrongblast):
                 os.rename(wrongblast, rescuedblast)
    fn_code_noext = nameformat % (filename, chunk_size, nbr)
    with open(direct + fn_code_noext + '@' + str(max_hsps) + 'maxHSPs' + '.txt', "w") as temp_hndl:
            temp_hndl.write('all done - %s chunks' % nbr)
    blastlog('all', nbr, 'chunks blasted','with', nbrwrong, 'failed chunks' )
    return nbr
