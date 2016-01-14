########################################################################################################################
# GENOMIC COORDOINATES INTERFACE for digestion
#
# - accept bed format string input to specify and retrieve fasta scaffold to process
# - export and save bed formatted data for found protospacers
# - refocus on input genomic interval (via bed genomic arithmetic) and export bed and fasta for focused protospacers
#
########################################################################################################################


import pybedtools


from cut import create_cutbedlines_from_seq_file


def cutbedlines_to_bedfile(cutbedlines, direct, fn_noext):
    with open(direct + fn_noext + '.prsp.bed', 'w') as bedfn:
       for cutbedline in cutbedlines:
          bedfn.write(cutbedline)


def bed_to_fasta(direct, fastaref, fn_noext):
    fastarefbedtool = pybedtools.BedTool(direct + fastaref)
    bedtoolin = pybedtools.BedTool(direct + fn_noext + '.bed')
    bedtoolin.sequence(fi=fastarefbedtool, s=True).save_seqs(direct + fn_noext + '.fasta')


# INITIALIZATION FUNCTION
# TODO check if still working conssitent with the focused version
# TODO this one already has a fasta file! ok?
# digest is the new name for the protospacers_bedfile_and_fastafile_from_directory function
def digest_file(direct, fn_noext, fileformat):
    """
    Only needed for initialization, the frst time a genome is being worked on.... TBA
    :param dir:
    :param fileformat:
    :return:
    """
    cutbedlines = create_cutbedlines_from_seq_file(direct, fn_noext, fileformat)
    # TODO move this to right spot so we know how many prsps we have ahead of blasting
    print "digesting all done: %s protospacers found" % len(cutbedlines)
    cutbedlines_to_bedfile(cutbedlines, direct, fn_noext)
    bed_to_fasta(direct, fn_noext + '.'+ fileformat, fn_noext + '.prsp')


# THE WORKHORSE FUNCTION
def digest_focused(direct, focusfn, wholefn, fileformat):
    """
    The core inner function handling digest. Saves 4 files
    :param direct:
    :param focusfn:
    :param wholefn:
    :param fileformat:
    :return:
    """
    focus_bedtool = pybedtools.BedTool(direct + focusfn +'.bed')
    whole_bedtool = pybedtools.BedTool(direct + wholefn + '.prsp.bed')
    whole_bedtool.intersect(focus_bedtool).moveto(direct + focusfn + ".prsp.bed")
    bed_to_fasta(direct, wholefn + '.' + fileformat, focusfn )
    bed_to_fasta(direct, wholefn + '.' + fileformat, focusfn+ '.prsp' )


################
# INPUT HANDLING
################

# TODO merge coord and stretch input options into a singke function
def coord_to_bedtuple_filename(coord):
    chrom = coord.split(':')[0]
    start, end = coord.split(':')[1].split('-')
    # bed coords are zero-based
    start0 = str(int(start) - 1)
    length = str(int(end) - int(start) + 1)
    filename = '%s:%s-%s_%s' % (chrom, start, end, length)
    return [(chrom, start0, end)], filename

assert(coord_to_bedtuple_filename('chr6:136640001-136680000')
    == ([('chr6', '136640000', '136680000')], 'chr6:136640001-136680000_40000'))

def stretch_to_bedtuple_filename(stretch):
    chrom = stretch.split(':')[0]
    start, length = stretch.split(':')[1].split('_')
    # bed coords are zero-based
    start0 = str(int(start) - 1)
    end = str(int(start) + int(length) - 1)
    filename =  '%s:%s-%s_%s' % (chrom, start, end, length)
    return [(chrom, start0, end)], filename

assert( stretch_to_bedtuple_filename('chr6:136640001_40000')
    == ([('chr6', '136640000', '136680000')], 'chr6:136640001-136680000_40000'))


####################
# MAIN API FUNCTIONS
####################

def digest_coord(direct, coord, wholefn, fileformat):
    bedtuplelist, focusfn = coord_to_bedtuple_filename(coord)
    print 'digesting ',  focusfn,
    pybedtools.BedTool(bedtuplelist).moveto(direct + focusfn + ".bed")
    digest_focused(direct, focusfn, wholefn, fileformat)
    print 'Done'
    return focusfn

def digest_stretch(direct, stretch, wholefn, fileformat):
    bedtuplelist, focusfn = stretch_to_bedtuple_filename(stretch)
    print 'digesting ',  focusfn,
    pybedtools.BedTool(bedtuplelist).moveto(direct + focusfn + ".bed")
    digest_focused(direct, focusfn, wholefn, fileformat)
    print 'Done'
    return focusfn
