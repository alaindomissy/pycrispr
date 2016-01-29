########################################################################################################################
#
# BIOPYTHON RESTRICTION module corrections and additions
#
# bug fix and customized Analysis maker that works
#
# error: http://nbviewer.ipython.org/gist/cfriedline/7790932
# fix: http://stackoverflow.com/questions/20381912/type-object-restrictiontype-has-no-attribute-size
# in ipython the following brings up error "type object 'RestrictionType' has no attribute 'size'"
# >>>from Bio import Restriction as rst
# >>>rst.EcoRI
#
# API functions : create_analysis
#
########################################################################################################################

from functools import reduce
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


###################
# MAIN API FUNCTION
###################

def create_analysis(seq, enzyme_names=['BfaI', 'HpaII', 'ScrFI']):
    """
    bug-fixed analysis creator that is dedicated to the 3 specific enzymes for crisr-eating
    :param seq:
    :param enzyme_names: defaults to the names of the 3 specific enzymes needed by crispr-eating
    :return: a working Analysis instance that performs digestion for given list of enzyme names
    """
    batch = create_batch(enzyme_names)
    return Analysis(batch, seq, linear= True)

# TODO check and discard forward (resp reverse) side of prsps with a PAM less then 20bps from start (resp end)
def analyse(seq, enzyme_names=['BfaI', 'HpaII', 'ScrFI']):
    """
    returns 1-based positions of first nucleotide cut along the forward strand
    for example with enzyme
    # BfAI
    # 5'...C TA G ...3'
    # 3'...G AT C ...5'
    the position returned is the one-based position of the T nucleotide in forward strand
    :param seq:
    :param enzyme_names:
    :return: a dict, keys are object enzymes for given enzyme_names, values are a lsit of cut positions for that enzyme
    """
    return create_analysis(seq, enzyme_names).full()     # .only_between(20, -20)  would get rid of both sides, not good



# # this should now work in ipython notebook
# >>>from Bio import Restriction as rst
# >>>rst.EcoRI
#
# >>>restbatch = create_enzyme("ScrFI") + create_enzyme("HpaII") + create_enzyme("BfaI")
# >>>restbatch
# RestrictionBatch(['BfaI', 'HpaII', 'ScrFI'])


# another possibly useful trick, not really needed
####################################################
# http://stackoverflow.com/questions/30561459/user-input-to-check-a-dna-sequence-for-restriction-sites-with-biopython

# Get RestrictionType by name
# def getrestbyname(enzyme_name):
#     batch = Re.RestrictionBatch()
#     batch.add(enzyme_name)
#     return batch.get(enzyme_name)
# getrestbyname('BfaI')


# TODO: subclass the class Analysis
#########################################
# per : http://biopython.org/DIST/docs/cookbook/Restriction.html#mozTocId968583
# might allow combining the digest work with the followimg bed amd fasta making (?)
