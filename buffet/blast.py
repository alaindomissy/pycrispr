########################################################################################################################
# BLAST
#
# API functions:  blast
#
########################################################################################################################

from __future__ import print_function
import os
from itertools import izip_longest
try:                                     # Python 2
    from itertools import izip_longest
except ImportError:                      # Python 3
    from itertools import zip_longest as izip_longest
## we no longer need izip, but if we did, since it does not exist in python3, we would do:
# try:
#     from itertools import izip, izip_longest
# except ImportError:
#     izip = zip

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Application import ApplicationError
from Bio import SeqIO as seqio    # TODO this is also imported from cut.py, ok?

# from buffet.settings import BLASTDBDIR


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
    print("\nBLASTING PROTOSPACERS fn_noext", fn_noext, end=' ')
    seqs = seqio.parse(dir + fn_noext + '.fasta','fasta')
    print('done')
    nbr = 0
    nameformat = '%s.%sseqs.%s'     # TODO add '.max%sHSPs.' % max_hsps
    for seqs_chunk in grouper_longest(seqs, chunk_size, None):
        seqs_chunk = [seq for seq in seqs_chunk if seq]
        nbr += 1
        # fn_code_noext = fn_noext + + 'hsps.'+  str(chunk_size) + 'x'  + str(nbr)

        fn_code_noext = nameformat % (fn_noext, chunk_size, nbr)
        rightfasta = dir + fn_code_noext + '.fasta'
        rightblast = dir + fn_code_noext + '.blast'
        wrongfasta = dir + fn_code_noext + '.fasta.WRONG'
        wrongblast = dir + fn_code_noext + '.blast.WRONG'
        rescuedfasta = dir + fn_code_noext + '.fasta.WRONG.RESCUED'
        rescuedblast = dir + fn_code_noext + '.blast.WRONG.RESCUED'
        if os.path.isfile(rightfasta) and os.path.isfile(rightblast):
            print('> skipping chunk ', nbr, fn_code_noext)
            continue
        elif os.path.isfile(wrongfasta) or os.path.isfile(wrongblast):
            print('> rescuing chunk ', nbr, fn_code_noext, end=' ')
        else:
            print('> blasting chunk ', nbr, fn_code_noext, end=' ')
        with open(dir + fn_code_noext + '.fasta', "w") as temp_hndl:
            seqio.write(seqs_chunk, temp_hndl, "fasta")
        failed = blast1(dir, fn_code_noext, blastdb_db, max_hsps)
        if failed:
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
    print('blasting all done:', nbr, 'chunks')
    return nbr
