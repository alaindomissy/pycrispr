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

# from crispr.settings import BLASTDBDIR


# SINGLE FILE BLASTING
######################

def blast1(dir, fn_noext, blastdb_db, max_hsps=50):
    # TODO max)hsps stopped working, now replaced with max_hsps_per_subject , but what is going on ?
    blastn_cline = NcbiblastnCommandline(
        query=dir + fn_noext + '.fasta',
        out=dir + fn_noext + '.blast',
        outfmt=5,
        db= blastdb_db,
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
        print(str(err).split('message ')[1].strip('\''))
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


def blast(dir, fn_noext, blastdb_db, chunk_size=50, max_hsps=50):
    print("\nBLAST PROTOSPACERS\n")
    print("\nload protospacers from fasta file", fn_noext, end=' ')
    seqs = seqio.parse(dir + fn_noext + '.fasta','fasta')
    print('done')
    nbr = 0
    nbrwrong = 0
    nameformat = '%s.%sseqs.%s'     # TODO add '.max%sHSPs.' % max_hsps
    for seqs_chunk in grouper_longest(seqs, chunk_size, None):
        seqs_chunk = [seq for seq in seqs_chunk if seq]
        nbr += 1
        # fn_code_noext = fn_noext + + 'hsps.'+  str(chunk_size) + 'x'  + str(nbr)

        fn_code_noext = nameformat % (fn_noext, chunk_size, nbr)
        rightfasta = dir + fn_code_noext + '.fasta'
        rightblast = dir + fn_code_noext + '.blast'
        wrongfasta = dir + fn_code_noext + '.fasta.FAILED'
        wrongblast = dir + fn_code_noext + '.blast.FAILED'
        rescuedfasta = dir + fn_code_noext + '.fasta.FAILED.RESCUED'
        rescuedblast = dir + fn_code_noext + '.blast.FAILED.RESCUED'
        if os.path.isfile(rightfasta) and os.path.isfile(rightblast):
            print('> skip chunk ', str(nbr).zfill(3), fn_code_noext)
            continue
        elif os.path.isfile(wrongfasta) or os.path.isfile(wrongblast):
            print('> rescue chunk ', str(nbr).zfill(3), fn_code_noext, end=' ')
        else:
            print('> blast chunk ', str(nbr).zfill(3), fn_code_noext, end=' ')
        with open(dir + fn_code_noext + '.fasta', "w") as temp_hndl:
            seqio.write(seqs_chunk, temp_hndl, "fasta")
        failed = blast1(dir, fn_code_noext, blastdb_db, max_hsps)
        if failed:
            nbrwrong +=1
            os.rename(rightfasta, wrongfasta)
            os.rename(rightblast, wrongblast)
        else:
            if os.path.isfile(wrongfasta):
                 os.rename(wrongfasta, rescuedfasta)
            if os.path.isfile(wrongblast):
                 os.rename(wrongblast, rescuedblast)
    fn_code_noext = nameformat % (fn_noext, chunk_size, nbr)
    with open(dir + fn_code_noext + '@' + str(max_hsps) + 'maxHSPs' + '.txt', "w") as temp_hndl:
            temp_hndl.write('all done - %s chunks' % nbr)
    print('all', nbr, 'chunks blasted','with', nbrwrong, 'failed chunks' )
    return nbr
