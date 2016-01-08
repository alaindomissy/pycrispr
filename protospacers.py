
########################################################################################################################

# bug fix and customized Analysis maker that works

# error: http://nbviewer.ipython.org/gist/cfriedline/7790932
# fix: http://stackoverflow.com/questions/20381912/type-object-restrictiontype-has-no-attribute-size
# in ipython the following brings up error "type object 'RestrictionType' has no attribute 'size'"
# >>>from Bio import Restriction as rst
# >>>rst.EcoRI

import gzip

from Bio import SeqIO as seqio
from Bio.Restriction import Analysis
from Bio.Restriction import Restriction

from Bio.Restriction.Restriction_Dictionary import rest_dict, typedict


def keywords_tuple(enzyme_name):
    """
    :param enzyme_name: string name of an enzyme
    :return: tuple of enzyme type names for given enzyme name
    >>>keywords_tuple('BfaI')
    ('Palindromic', 'OneCut', 'Ov5', 'Defined', 'Meth_Undep', 'Commercially_available', 'AbstractCut', 'RestrictionType')
    """
    keywords_tuple_list = [keywords_tuple
            for re_type, (keywords_tuple, re_names_list) in typedict.items()
            if enzyme_name in re_names_list]
    # this list length should be 1, as each enzyme belongs to one and only one type entry of typedict
    return keywords_tuple_list[0]

def create_enzyme(enzyme_name):
    """
    >>>type(create_enzyme('BfaI'))
    Bio.Restriction.Restriction.RestrictionType)
    """
    enzyme_types = tuple(getattr(Restriction, x) for x in keywords_tuple(enzyme_name))
    return Restriction.RestrictionType(enzyme_name, enzyme_types, rest_dict[enzyme_name])

def create_batch(enzyme_names):
    return reduce(lambda x, y: x + y, map(create_enzyme, enzyme_names))

def create_analysis(seq, enzyme_names=['BfaI', 'HpaII', 'ScrFI']):
    batch = create_batch(enzyme_names)
    return Analysis(batch, seq, linear= True)

# ### this should now work in ipython notebook
# >>>from Bio import Restriction as rst
# >>>rst.EcoRI
#
# >>>restbatch = create_enzyme("ScrFI") + create_enzyme("HpaII") + create_enzyme("BfaI")
# >>>restbatch
# RestrictionBatch(['BfaI', 'HpaII', 'ScrFI'])

########################################################################################################################

# another possibly useful trick, not really needed

# http://stackoverflow.com/questions/30561459/user-input-to-check-a-dna-sequence-for-restriction-sites-with-biopython

# Get RestrictionType by name
# def getrestbyname(enzyme_name):
#     batch = Re.RestrictionBatch()
#     batch.add(enzyme_name)
#     return batch.get(enzyme_name)
# getrestbyname('BfaI')

########################################################################################################################

# TODO: subclass the class Analysis
# per : http://biopython.org/DIST/docs/cookbook/Restriction.html#mozTocId968583
# might allow combining the digest work with the followimg bed amd fasta making (?)


def tabbed_string_from_list(list):
    return '\t'.join(map(str, list)) + '\n'

# TODO somehow retain the info about enzyme name (goes in seqrecord.description)
# TODO store a serial id number for each cut, concatenated with with _F or _R (goes in seqrecord.id)
# TODO store cutpos (goes in seqrecords.name)
def create_cutbedline(enzyme, cutpos, substrate_id, substrate_end, sense='+'):
    if sense == '+':
        bedend = cutpos - 1
        bedstart = max(bedend - 20, 0)          # TODO do we still need the max?
    else:
        assert sense=='-'
        bedstart = cutpos - enzyme.ovhg - 1     # TODO do we still need the min?
        bedend = min(bedstart + 20, substrate_end)
    # cutinfo = [substrate_id, bedstart, bedend, str(enzyme), '1000', sense, bedstart, bedend, '255,0,0']
    cutinfo = [substrate_id, bedstart, bedend, str(enzyme), '1000', sense]
    cutbedline = tabbed_string_from_list(cutinfo)
    return cutbedline

def create_cutbedlines(cutdict, substrate_id, substrate_end):
    return [create_cutbedline(enzyme, cutpos, substrate_id, substrate_end, sense)
            for enzyme in cutdict.keys()
            for index, cutpos in enumerate(cutdict[enzyme])
            for sense in ['+', '-']
            ]


def digest_seqrecord_to_cutbedlines(seqrecord, enzyme_names=['BfaI', 'HpaII', 'ScrFI'] ):
    # cut_dict = create_analysis(seqrecord.seq, enzyme_names).full()        # TODO we should only loose 1 side, not both
    cut_dict = create_analysis(seqrecord.seq, enzyme_names).only_between(20, len(seqrecord)-20)   # TODO is 20 the best?
    return create_cutbedlines(cut_dict, seqrecord.id, len(seqrecord))

def digest_seqrecords_to_cutbedlines(seqrecords, enzyme_names=['BfaI', 'HpaII', 'ScrFI'] ):
    cutinfos = []
    for seqrecord in seqrecords:
        cutinfos.extend(digest_seqrecord_to_cutbedlines(seqrecord, enzyme_names))
    return cutinfos



def create_seqrecords_from_directory(dir, fileformat):
    if fileformat[-3:] == '.gz':
        string_or_handle = gzip.open(dir + "/" + "cat." + fileformat, 'r')
        parseformat = fileformat[:-3]
    else:
        string_or_handle = dir + "/" + "cat." + fileformat
        parseformat = fileformat
    return list(seqio.parse(string_or_handle, parseformat))

def digest_directory(dir, fileformat='fasta'):
    seqrecords = create_seqrecords_from_directory(dir, fileformat)
    return digest_seqrecords_to_cutbedlines(seqrecords)


#################
import pybedtools


def cutbedlines_to_bedfile(cutbedlines, outputdirpath=''):
    with open(outputdirpath + 'protospacers.bed', 'w') as bedfn:
       for cutbedline in cutbedlines:
          bedfn.write(cutbedline)


def bedfile_to_fastafile(dir, filename, fileformat):
    scaffoldfasta = pybedtools.BedTool(dir + 'cat.' + fileformat)
    bedtool = pybedtools.BedTool(dir + filename + '.bed').sequence(fi=scaffoldfasta)
    bedtoolstranded = pybedtools.BedTool(dir + filename + '.bed').sequence(fi=scaffoldfasta, s=True)
    bedtool.save_seqs(dir + filename + '.fasta')
    bedtoolstranded.save_seqs(dir + filename + '.stranded.fasta')


def protospacers_bedfile_and_fastafile_from_directory(dir, fileformat):
    cutbedlines = digest_directory(dir, fileformat)
    cutbedlines_to_bedfile(cutbedlines, dir)
    bedfile_to_fastafile(dir, 'protospacers', fileformat)


def protospacers_focused(focus_input, dir, fileformat):
    focus_bedtool = pybedtools.BedTool(focus_input)
    focus_protos = pybedtools.BedTool(dir + "protospacers.bed")
    focus_protos.intersect(focus_bedtool).moveto(dir + "protospacers_focused.bed")
    bedfile_to_fastafile(dir, 'protospacers_focused', fileformat )


#########################################################
# DIGEST AND BLAST AND SCORE
#########################################################

# requires:
# cat.fasta    the scaffold to digest and select guides in
# dict         a genome fasta file against which pam lookup is optionnaly done for guides candidates
# genome       same, used when (correctly) blasting agaist whole geome
#
def dgnblnsc(dir, fileformat,
             fn_noext, blastdb_directory, blastdb_db, chunksize=1000,
             max_hsps=100,
             reref_substrate_id=None, low=50, high=90, load_genome=False, howmany=24):

    protospacers_bedfile_and_fastafile_from_directory(dir, fileformat)

    blnsc(dir, fn_noext, blastdb_directory, blastdb_db, chunksize=1000,
          max_hsps=100,
          reref_substrate_id=None, low=50, high=90, load_genome=False, howmany=24):


#########################################################
# BLAST AND SCORE
#########################################################

def blnsc(dir, fn_noext, blastdb_directory, blastdb_db, chunksize=1000,
          max_hsps=100,
          reref_substrate_id=None, low=50, high=90, load_genome=False, howmany=24):

    nbr = blast(dir, fn_noext, blastdb_directory, blastdb_db, chunksize=chunksize,
          max_hsps=max_hsps)

    score(dir, fn_noext, blastdb_directory, blastdb_db, chunksize=chunksize,
          nbrofchunks = nbr,
          reref_substrate_id=reref_substrate_id, low=low, high=high, load_genome=load_genome, howmany=howmany)


#########################################################
# BLAST
#########################################################
import os
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Application import ApplicationError
from Bio.Blast import NCBIXML

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

###
from itertools import izip, izip_longest

# def grouper(iterable, chunksize):
#     args = [iter(iterable)] * chunksize
#     return izip(*args)

def grouper_longest(iterable, chunksize, fillvalue=None):
    args = [iter(iterable)] * chunksize
    return izip_longest(*args, fillvalue=fillvalue)

def blast(dir, fn_noext, blastdb_directory, blastdb_db, chunksize=1000, max_hsps=100):
    seqs = seqio.parse(dir + fn_noext + '.fasta','fasta')
    nbr = 0
    for seqs_chunk in grouper_longest(seqs, chunksize, None):
        seqs_chunk = [seq for seq in seqs_chunk if seq]
        nbr += 1
        fn_code_noext = fn_noext + str(nbr) + 'x' + str(chunksize)
        print 'blasting chunk ', fn_code_noext,
        with open(dir + fn_code_noext + '.fasta', "w") as temp_hndl:
            seqio.write(seqs_chunk, temp_hndl, "fasta")
        failed = blast1(dir, fn_code_noext, blastdb_directory, blastdb_db)
        if failed:
            os.rename(dir + fn_code_noext + '.fasta', dir + fn_code_noext + '.FAILED.fasta' )
            os.rename(dir + fn_code_noext + '.blast', dir + fn_code_noext + '.FAILED.blast' )
    with open(dir + fn_code_noext + '.chunks', "w") as temp_hndl:
            temp_hndl.write('all done - %s chunks' % nbr)
    print 'all done with', nbr, 'chunks'
    return nbr



####################################################
# SCORE
####################################################
import subprocess
import cStringIO
import collections
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

WEIGHT = [0, 0, 0.014, 0, 0, 0.395, 0.317, 0 ,0.389, 0.079, 0.445, 0.508, 0.613, 0.851, 0.732, 0.828, 0.615, 0.804, 0.685, 0.583]

def is_valid_pam(pam):
    return True
    # return (len(pam)==3 and (pam[1] == "G" or pam[1] == "A" or pam[1] == "g" or pam[1] == "a") and (pam[2] == "G" or pam[2] == "g"))


def score(dir, fn_noext, blastdb_directory, blastdb_db,
          chunksize, nbrofchunks, reref_substrate_id=None, low=50, high=90,
          load_genome=False, howmany=24):

    # genomedict = a dict made from a FASTA file identical to the one used to make the BLAST DB.
    # dict keys should be BLAST db hit_def.
    if load_genome:
        print 'reading genome',
        genomeseq = seqio.parse(open(dir+'dict.fasta', 'rb'), "fasta", alphabet=IUPACAmbiguousDNA())
        genomedict = {}
        for seq in genomeseq:
            genomedict[seq.id] = seq
        print 'done'

    with open(dir + fn_noext + '.fasta') as guidesfn:
        guides = list(seqio.parse(guidesfn, "fasta"))

    blastrecords = []
    for chunknbr in range(1,nbrofchunks+1):
        blastfn = dir  + fn_noext + str(chunknbr) + 'x' + str(chunksize)+ '.blast'
        try:
            with open(blastfn) as blasthndl:
                print 'parsing', blastfn
                print '',
                blastrecords.extend(list(NCBIXML.parse(blasthndl)))
                print 'done'
        except IOError as ioe:
            print 'skipping ', blastfn
            print '',

        # print blastrecords

    for blastindex, blastitem in enumerate(blastrecords):

        scorelist = []
        pamlist=[]

        for alignment in blastitem.alignments:    # alignment corresponds to a hit in the blast file
                                                  # a hit is a whole seq from  blastdb, many hsps can exist for 1 hit
            for hsp in alignment.hsps:

                foundfullmatch = False

                # getting the pam adjacent to the hsp's subject

                hit_threeprime_offset = len(guides[0]) - hsp.query_end
                pamstart = hsp.sbjct_end + hit_threeprime_offset

                # if sbjct_is_on_forward_strand, equivalent to (hsp.sbjct_end > hsp.sbjct_start)
                if (hsp.frame[1] > 0):
                    pam_start = hsp.sbjct_end + hit_threeprime_offset
                    pam_onebased_range = (pam_start + 1, pam_start + 3)
                    pam_zerobased_range = (pam_start , pam_start + 3)
                else:
                    pam_start = hsp.sbjct_end - hit_threeprime_offset
                    pam_onebased_range = (pam_start -3, pam_start -1)
                    pam_zerobased_range = (pam_start -4, pam_start -1)

                # TODO *** get the correct substrate id for subject (not same as hit!) - should be alignment.hit_def
                # TODO *** use betools and fai instead of loading full genome
                if load_genome:
                    # print 'using the genome dict to lookup pam'
                    lookup_context = genomedict[reref_substrate_id]
                    pam = lookup_context[pam_zerobased_range[0]:pam_zerobased_range[1]]
                else:
                    print 'using blastdbcmd to lookup pam'
                    context_lookup_command = "blastdbcmd -db " + blastdb_directory + blastdb_db \
                                             + " -dbtype nucl -entry " + alignment.accession \
                                             + " -range %s-%s" % pam_onebased_range
                    context_lookup_process = subprocess.Popen(context_lookup_command, stdout=subprocess.PIPE, shell = True)
                    fstring = context_lookup_process.communicate()
                    fstring = cStringIO.StringIO(fstring[0])
                    pam = seqio.read(fstring, "fasta")

                if not (hsp.frame[1] > 0):
                    pam = pam.reverse_complement()

                # make match string padded to query length, where bar(|) is match and space( ) is non-match
                mmloc = []
                mmstr = list(hsp.match)
                if hsp.query_start > 1:
                    mmstr = list(" " * (hsp.query_start - 1)) + mmstr
                if hsp.query_end < 20:
                    mmstr = mmstr + list(" " * (20 - hsp.query_end))
                mmstr = "".join(mmstr)

                # Test for PAM adjacency
                if len(pam) == 3:
                    if not is_valid_pam(pam):
                        scorelist.append({"match_score": float(0),
                                          "alignment": alignment.hit_def,
                                          "hit_location":hsp.sbjct_start,
                                          "hit_sequence": hsp.sbjct,
                                          "pam":str(pam.seq),
                                          "match_bars":mmstr})
                        pamlist.append(pam)
                    else:
                        if (hsp.positives >16 and hsp.positives < 20):
                            pos = 20
                            for linestring in mmstr:
                                if linestring != "|":
                                    mmloc.append(pos)
                                pos = pos - 1

                            # Actually implement Zhang lab algorithm
                            mmscore = [21 -x for x in mmloc]
                            t1 = 1
                            for mismatchlocation in mmscore:
                                t1 = t1 * (1.0 - float(WEIGHT[mismatchlocation - 1]))
                            if len(mmscore) > 1:
                                d = (float(max(mmscore)) - float(min(mmscore))) / float((len(mmscore) - 1))
                            else:
                                d = 19
                            t2 = 1 / ( (((19.0 - d)/19.0) * 4) + 1)
                            t3 = float(1)/ float(pow(len(mmloc), 2))
                            scorelist.append({"match_score": float(t1 * t2 * t3 * 100),
                                              "alignment": alignment.hit_def,
                                              "hit_location":hsp.sbjct_start,
                                              "hit_sequence": hsp.sbjct,
                                              "pam":str(pam.seq),
                                              "match_bars":mmstr})
                            pamlist.append(pam)
                        # Zhang lab algorithm doesn't handle perfect matches: give it a 50 if it's perfect
                        # we dont want more them one of these per hsps ?
                        if hsp.positives >= 20: #changed from == 20; might be worth keeping in mind for bugs
                            if not foundfullmatch:
                                scorelist.append({"match_score": float(50),
                                                  "alignment": alignment.hit_def,
                                                  "hit_location":hsp.sbjct_start,
                                                  "hit_sequence": hsp.sbjct,
                                                  "pam":str(pam.seq),
                                                  "match_bars":mmstr})
                                pamlist.append(pam)
                            foundfullmatch = True



        finalscore = int(10000.000 / (100.000 + float(sum(item["match_score"] for item in scorelist))))
        guides[blastindex].annotations['score'] = finalscore
        guides[blastindex].annotations["blastindex"] = blastindex

        # print "seq: ", guides[blastindex].seq, "score: ", finalscore, "id: ", guides[blastindex].id

        # print guides[blastindex].seq
        # print alignment.hit_def
        # print hsp.match
        # print hsp.sbjct
        # print pam.seq

    # print [(guide.seq, guide.annotations.get('score')) for guide in guides]

    # sort guides by position, using the seqrecord id
    guides.sort(key=lambda x: int(x.id.split(':')[2].split('-')[0]))

    guides_scorebedfile(dir, 'scoreguides.rel', guides, low, high)
    if reref_substrate_id:
        guides_scorebedfile(dir, 'scoreguides.abs',
                            guides, reref_substrate_id, low, high)

    # TODO guides not getting a score, how does this happen ? fix it better
    # kind of fix guides without a score
    for guide in guides:
        try:
            guide.annotations['score']
        except:
            # print 'guide %s does not have a score' % guide
            guide.annotations['score'] = 0.0

    ampls = amplicons(pos_score(guides), high)
    print_ampls_info(ampls, howmany)
    print_scores_info(scores(guides))
    histo(dir, guides)

    return guides, ampls



####
import matplotlib.pyplot as plt
from numpy import mean


def histo(dir, guides):
    scores = [guide.annotations['score'] for guide in guides]
    bins = range(0,106,5)
    plt.hist(scores, bins, color="gray")
    plt.tick_params(axis=u'both', labelsize=18)
    plt.savefig(dir + 'scores.pdf', format="pdf")



def guides_scorebedfile(dir, fn_noext,
                        guides, reref_substrate_id=None, low=30, high=30):
    scorebedfile(dir, fn_noext,
                 scorebedlines(guides, reref_substrate_id, low, high)
                 )

def scorebedfile(outputdirpath, fn_noext, bedlines):
    with open(outputdirpath + fn_noext + '.bed', 'w') as bedfn:
       for bedline in bedlines:
          bedfn.write(bedline)

def scorebedlines(guides, reref_substrate_id=None, thlow=50, high=90):
    return [scorebedline(guide, reref_substrate_id) for guide in guides]

def scorebedline(guide, reref_substrate_id=None,low=30, high=30):
    substrate_id = guide.id.split('-')[0]
    bedstart = guide.id.split(':')[2].split('-')[0]
    if reref_substrate_id:
        substrate_id = reref_substrate_id
        bedstart = str(int(bedstart) + int(guide.id.split(':')[1].split('-')[0]))
    bedend = str(int(bedstart) + 20)
    sense = guide.id[-2]
    try:
        score = str(guide.annotations['score'])
    except:
        score ='0'
    if int(score) <= low:
        color = '255,0,0'
    elif int(score) >= high:
        color = '0,255,0'
    else:
        color = '0,0,255'
    cutinfo = [substrate_id, bedstart, bedend, 'crRNA', score, sense, bedstart, bedend, color]
    # cutinfo = [substrate_id, bedstart, bedend, 'crspr_guide', 'score', sense]
    bedline = tabbed_string_from_list(cutinfo)
    return bedline



###
# def index_score(guides):
#     return[(guide.annotations['blastindex'], guide.annotations['score'])
#            for guide in guides]

def scores(guides):
    return [guide.annotations['score'] for guide in guides]

def pos_score(guides):
    return[(int(guide.name.split(':')[2].split('-')[0]) , guide.annotations['score'])
           for guide in guides]

def amplicons(pos_score_tuples, high=90):
    ampls=[]
    goodones = 0
    for tuple in pos_score_tuples:
        #print tuple
        if tuple[1]>=high:
            goodones+=1
            # if this is the first good, reset begin to this position
            if goodones==1:
                begin = tuple[0]
        else:
            # if we were in a good stretch until now, reset goodones to 0
            if goodones>0:
                end = tuple[0] + 20
                ampls.append((begin, end, end-begin, goodones, 1000*goodones/(end-begin)))
                goodones = 0
    ampls.sort(key=lambda x : -x[3])
    return ampls


def print_scores_info(scores):
    print
    print 'distribution of scores among all', len(scores), 'guides (score, nb of guides with that score)'
    print sorted(collections.Counter(scores).items())


def print_ampls_info(ampls, howmany=24):
    print('')
    print('%s amplicons totaling %s bps yielding a total of %s guides'
            % (len(ampls), sum([t[2] for t in ampls]), sum([t[3] for t in ampls]))
            )
    print('with average amplicon of %s bps yielding %s guides'
            % (mean([t[2] for t in ampls])//1, mean([t[3] for t in ampls])//1)
            )
    print('')
    print('top %s amplicons totaling %s bps yielding a total of %s guides'
            % (howmany, sum([t[2] for t in ampls[0:howmany]]), sum([t[3] for t in ampls[0:howmany]]))
            )
    print('with average amplicon of %s bps yielding %s guides'
            % (mean([t[2] for t in ampls[0:howmany]])//1, mean([t[3] for t in ampls[0:howmany]])//1)
            )
    print('')
    for ampl in ampls[0:howmany]:
        print("%s - %s (%s bps) yields %s guides (%s guides per 1000bps )" % ampl)
    print('')
    for ampl in ampls[howmany:2*howmany]:
        print("%s - %s (%s bps) yields %s guides (%s guides per 1000bps )" % ampl)
    print('')

    return ampls



###
# see: http://sebastianraschka.com/Articles/2014_multiprocessing_intro.html#An-introduction-to-parallel-programming-using-Python's-multiprocessing-module
# import multiprocessing as mp
#
# def multi_protospacers_score(dir, blastdb_directory, blastdb_db):
#     pool = mp.Pool(processes=2)
#     guides = [pool.apply(protospacers_score, args=(x,)) for x in cutslist]
#
# def multiscore_pool(x):
#    score = al_scoreguide(x, "xl71", xl71genomedict)
#    return score
###








 # missing the sequences! now implemented with pybedtools
# def save_to_fastafile(bedtuples, outputdirpath=''):
#     with open(outputdirpath + 'cuts_nameonly.fasta', 'w') as outfp:
#         idx = 2
#         for tuple in bedtuples:
#             # (substrate_id, bedstart, bedend, enzyme_name, score, sense, bbedstart, bbedend, color) = tuple
#             (substrate_id, bedstart, bedend, enzyme_name, score, sense) = tuple
#             id = str(idx//2) + '_' + ('F' if sense=='+' else 'R')
#             idx += 10
#             outfp.write(">lcl|" + substrate_id + "|" + enzyme_name + "|" + str(bedstart) + "|" + id +"\n")
#             # outfp.write(str(item.seq) + "\n")
#             outfp.write("\n")
30