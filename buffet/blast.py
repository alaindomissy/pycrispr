########################################################################################################################
# BLAST
#
# API functions:  blast
#
########################################################################################################################

from __future__ import print_function
import os
from itertools import izip, izip_longest

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Application import ApplicationError

# TODO this is also imported from cut.py, ok?
from Bio import SeqIO as seqio


# SINGLE FILE BLASTING

def blast1(dir, fn_noext, blastdb_directory, blastdb_db, max_hsps=100):
    blastn_cline = NcbiblastnCommandline(
        query=dir + fn_noext + '.fasta',
        out=dir + fn_noext + '.blast',
        outfmt=5,
        db=blastdb_directory + blastdb_db,
        max_target_seqs=25,
        num_threads=4,
        evalue=10,
        max_hsps=max_hsps,
        task="blastn-short",
        dust="no",
        )
    # print('blastn_cline:  ', end='')
    # print(blastn_cline)
    try:
        # blastn_cline()
        stdout, stderr = blastn_cline()
    except ApplicationError as err:
        # print('Bio:Application:ApplicationError number {0} : {1}'.format(err.errno, err.strerror), end='')
        # print('Bio:Application:ApplicationError: ', err.args)
        print(str(err).split('message ')[1].strip('\''))
        #print(err.args)
        return True
    print("Done")
    return False



# ITERATION UTILITY

# def grouper(iterable, chunksize):
#     args = [iter(iterable)] * chunksize
#     return izip(*args)

def grouper_longest(iterable, chunksize, fillvalue=None):
    args = [iter(iterable)] * chunksize
    return izip_longest(*args, fillvalue=fillvalue)


###################
# MAIN API FUNCTION
###################


def blast(dir, fn_noext, blastdb_directory, blastdb_db, chunksize=100, max_hsps=100):
    # print("starting blast function with " , fn_noext)
    print('parsing protospacers', fn_noext, end='')
    seqs = seqio.parse(dir + fn_noext + '.fasta','fasta')
    print('Done')
    nbr = 0
    nameformat = '%s.%sseqs.%s'     # TODO add '.max%sHSPs.' % max_hsps
    for seqs_chunk in grouper_longest(seqs, chunksize, None):
        seqs_chunk = [seq for seq in seqs_chunk if seq]
        nbr += 1
        # fn_code_noext = fn_noext + + 'hsps.'+  str(chunksize) + 'x'  + str(nbr)

        fn_code_noext = nameformat % (fn_noext, chunksize, nbr)
        rightfasta = dir + fn_code_noext + '.fasta'
        rightblast = dir + fn_code_noext + '.blast'
        wrongfasta = dir + fn_code_noext + '.fasta.WRONG'
        wrongblast = dir + fn_code_noext + '.blast.WRONG'
        rescuedfasta = dir + fn_code_noext + '.fasta.WRONG.RESCUED'
        rescuedblast = dir + fn_code_noext + '.blast.WRONG.RESCUED'
        if os.path.isfile(rightfasta) and os.path.isfile(rightblast):
            # print('skipping chunk ', nbr, fn_code_noext)
            continue
        elif os.path.isfile(wrongfasta) or os.path.isfile(wrongblast):
            print('rescuing chunk ', nbr, fn_code_noext, end='')
        else:
            print('blasting chunk ', nbr, fn_code_noext, end='')
        with open(dir + fn_code_noext + '.fasta', "w") as temp_hndl:
            seqio.write(seqs_chunk, temp_hndl, "fasta")
        failed = blast1(dir, fn_code_noext, blastdb_directory, blastdb_db, max_hsps)
        if failed:
            os.rename(rightfasta, wrongfasta)
            os.rename(rightblast, wrongblast)
        else:
            if os.path.isfile(wrongfasta):
                 os.rename(wrongfasta, rescuedfasta)
            if os.path.isfile(wrongblast):
                 os.rename(wrongblast, rescuedblast)
    fn_code_noext = nameformat % (fn_noext, chunksize, nbr)
    with open(dir + fn_code_noext + '@' + str(max_hsps) + 'maxHSPs' + '.txt', "w") as temp_hndl:
            temp_hndl.write('all done - %s chunks' % nbr)
    print('blasting all done:', nbr, 'chunks')
    return nbr

