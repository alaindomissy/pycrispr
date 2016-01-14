########################################################################################################################
# BLAST
#
# API functions:  blast
#
########################################################################################################################

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
        num_threads=6,
        evalue=10,
        max_hsps=max_hsps,
        task="blastn-short",
        dust="no",
        )
    # print 'blastn_cline:  ',
    # print blastn_cline
    try:
        # blastn_cline()
        stdout, stderr = blastn_cline()
    except ApplicationError as err:
        # print 'Bio:Application:ApplicationError number {0} : {1}'.format(err.errno, err.strerror),
        print '*failed*'
        print 'Bio:Application:ApplicationError: ', err
        return True
    print "done"
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
    print "starting blast function with " , fn_noext
    seqs = seqio.parse(dir + fn_noext + '.fasta','fasta')
    nbr = 0
    for seqs_chunk in grouper_longest(seqs, chunksize, None):
        seqs_chunk = [seq for seq in seqs_chunk if seq]
        nbr += 1
        # fn_code_noext = fn_noext + '.max' + str(max_hsps) + 'hsps.'+  str(chunksize) + 'x'  + str(nbr)
        fn_code_noext = fn_noext + str(chunksize) + 'x' + str(nbr)
        print 'blasting chunk ', fn_code_noext,
        with open(dir + fn_code_noext + '.fasta', "w") as temp_hndl:
            seqio.write(seqs_chunk, temp_hndl, "fasta")
        failed = blast1(dir, fn_code_noext, blastdb_directory, blastdb_db, max_hsps)
        if failed:
            # TODO write functions to report and attempt reblast of failed chunks
            os.rename(dir + fn_code_noext + '.fasta', dir + fn_code_noext + '.FAILED.fasta' )
            os.rename(dir + fn_code_noext + '.blast', dir + fn_code_noext + '.FAILED.blast' )
    fn_code_noext = fn_noext + str(chunksize) + 'x' + str(nbr)
    with open(dir + fn_code_noext  + '@max' + str(max_hsps) + 'hsps' + '.txt', "w") as temp_hndl:
            temp_hndl.write('all done - %s chunks' % nbr)
    print 'all done with', nbr, 'chunks'
    return nbr

