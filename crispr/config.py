from __future__ import absolute_import, division, print_function   # , unicode_literals
import os
from shutil import rmtree
from subprocess import check_output

GENOMES = os.environ.get('GENOMES', '/GENOMES/')
BLASTDB = os.environ.get('BLASTDB', '/RESTORE/')
PROTOSP = os.environ.get('PROTOSP','/RESTORE/')

SCRATCH = os.environ.get('SCRATCH','/data/scratch/')

APPSESS = os.environ.get('APPSESS', '/data/input/AppSession.json')


RESTRICTION_ENZYMES_LIST = ['BfaI', 'HpaII', 'ScrFI']

DIGEST_LOG_ON = False
BLAST_LOG_ON = False
SCORE_LOG_ON = False
STRETCH_LOG_ON = True
AMPLICON_LOG_ON = True
PRIME_LOG_ON = True


def digestlog(*args, **kwargs):
    if DIGEST_LOG_ON:
        print(*args, **kwargs)

def blastlog(*args, **kwargs):
    if BLAST_LOG_ON:
        print(*args, **kwargs)

def scorelog(*args, **kwargs):
    if SCORE_LOG_ON:
        print(*args, **kwargs)

def stretchlog(*args, **kwargs):
    if STRETCH_LOG_ON:
        print(*args, **kwargs)

def ampliconlog(*args, **kwargs):
    if AMPLICON_LOG_ON:
        print(*args, **kwargs)

def primelog(*args, **kwargs):
    if PRIME_LOG_ON:
        print(*args, **kwargs)


# Default parameters for primer3 primer design, intended for NEB's Q5 HotStart polymerase
PRIMER3_GLOBAL_ARGS = {
    'PRIMER_NUM_RETURN': 20,
    'PRIMER_OPT_SIZE': 20,
    'PRIMER_PICK_INTERNAL_OLIGO': 0,
    'PRIMER_INTERNAL_MAX_SELF_END': 8,
    'PRIMER_MIN_SIZE': 18,
    'PRIMER_MAX_SIZE': 25,
    'PRIMER_OPT_TM': 60.0,
    'PRIMER_MIN_TM': 50.0, # was 56
    'PRIMER_MAX_TM': 75.0, # was 67
    'PRIMER_MIN_GC': 20.0,
    'PRIMER_MAX_GC': 80.0,
    'PRIMER_MAX_POLY_X': 100,
    'PRIMER_INTERNAL_MAX_POLY_X': 100,
    'PRIMER_SALT_MONOVALENT': 50.0,
    'PRIMER_DNA_CONC': 50.0,
    'PRIMER_MAX_NS_ACCEPTED': 0,
    'PRIMER_MAX_SELF_ANY': 12,
    'PRIMER_MAX_SELF_END': 8,
    'PRIMER_PAIR_MAX_COMPL_ANY': 12,
    'PRIMER_PAIR_MAX_COMPL_END': 8,

    # 'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],[150,175],[175,200],[200,225]],
    # 'PRIMER_PRODUCT_SIZE_RANGE': [[26, 180]],
    'PRIMER_PRODUCT_SIZE_RANGE': [[26, 18000]],

    # 'PRIMER_WT_SELF_ANY': 1, # added
    # 'PRIMER_PAIR_WT_COMPL_ANY':1, #added
    # 'PRIMER_PAIR_WT_COMPL_END':1, #added
    # 'PRIMER_PAIR_WT_SELF_END':1, #added
    # 'PRIMER_PAIR_WT_SELF_ANY':1, #added

    'PRIMER_MIN_THREE_PRIME_DISTANCE':3,
    'PRIMER_PAIR_MAX_DIFF_TM':5, #added
    'PRIMER_PAIR_WT_DIFF_TM':1, #added
    'PRIMER_PAIR_WT_PRODUCT_SIZE_LT': 0.5, #added
    }



def sorted_unique_firstdotsplit(filenammes):
    """
    :param filenammes:
    :return:
    >>>filenames = ['phix.fasta', 'phix.fasta.fai', 'hg38.fasta.fai', 'hg38.fasta','saccer3.fasta', 'saccer3.fasta.fai']
    >>>sorted_unique_firstdotsplit(filenammes)
    ['hg38', 'phix', 'saccer3']cd /GE
    """
    return sorted(list(set([file.split('.')[0] for file in filenammes])))

# GENOMESPATHS = {}

GENOMESPATHS = {
    # 'hg18'     : '/genomes/Homo_sapiens/UCSC/hg18/Sequence/WholeGenomeFasta/genome.fa',
    'hg18'     : '/RESTORE/hg18/hg18.fasta',
    'hg19'     : '/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa',
    'hg38'     : '/genomes/Homo_sapiens/UCSC/hg38/Sequence/WholeGenomeFasta/genome.fa',
    'mm8'      : '/RESTORE/mm8/mm8.fasta',
    'mm9'      : '/genomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa',
    'mm10'     : '/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa',
    'ecoli'    : '/genomes/Escherichia_coli_K_12_DH10B/NCBI/2008-03-17/Sequence/WholeGenomeFasta/genome.fa',
    #'ecoli'    : '/RESTORE/ecoli/ecoli.fasta',
    'mycotube' : '/genomes/Mycobacterium_tuberculosis_H37RV/NCBI/2001-09-07/Sequence/WholeGenomeFasta/genome.fa',
    'phix'     : '/genomes/PhiX/Illumina/RTA/Sequence/WholeGenomeFasta/genome.fa',
    'saccer3'  : '/genomes/Saccharomyces_cerevisiae/UCSC/sacCer3/Sequence/WholeGenomeFasta/genome.fa'
}


def genomes_list():
    """
    :return:
    """
    return sorted_unique_firstdotsplit(os.listdir(GENOMES))

def blastdb_list():
    """
    :return:
    """
    return sorted_unique_firstdotsplit(os.listdir(BLASTDB))

def protosp_list():
    """
    :return:
    """
    return sorted_unique_firstdotsplit(os.listdir(PROTOSP))


def genomes_path(genome):
    """
    >>>genomes_path('mm8')
    /RESTORE/mm8/mm8.fasta
    :param genome:
    :return:
    """
    # print("> load reference genome at", GENOMESPATHS[genome] + '.fasta')
    # print("> load reference genome at", GENOMESPATHS[genome])
    # return GENOMES + genome + '/' + genome + '.fasta'
    return GENOMESPATHS[genome]


def blastdb_path(genome):
    """
    >>>blastdb_path('mm8')
    /RESTORE/mm8/mm8.fasta
    :param genome:
    :return:
    """
    return BLASTDB + genome + '/' + genome

def protosp_path(genome):
    """
    >>>protosp_path('mm8')
    /RESTORE/mm8/mm8.prsp.bed
    :param genome:
    :return:
    """
    return(PROTOSP + genome + '/' + genome + '.prsp.bed')


def status_genome(genome):
    print(check_output(['duply', genome, 'status']))


def backup_genome(genome):
    print(check_output(['duply', genome, 'backup']))


def restore_genome(genome, force=False):
    is_available = os.path.isdir('/RESTORE/' + genome)
    if not is_available or force:
        print('\nREMOVE/RESTORE REFERENCE GENOME', genome)
        # print(check_output(['rm', '/RESTORE/' + genome + '/*']))
        # print(check_output(['rmdir', '/RESTORE/' + genome]))
        # os.remove('/RESTORE/' + genome + '/' + genome + '.*')
        # os.removedirs('/RESTORE/' + genome)
        rmtree('/RESTORE/' + genome,  ignore_errors=True)
        output = check_output(['duply', genome, 'restore', '/RESTORE/' + genome])
        # print(output)
