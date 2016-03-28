from __future__ import absolute_import, division, print_function
# from __future__ import unicode_literals
from collections import namedtuple

class Bedtuple(namedtuple('Bedtuple', 'scaffold beg end')):
    """
    >>>bt1 = Bedtuple('chr1',1000,1100)
    Bedtuple(scaff='chr1', beg=1000, end=1100)
    >>>len(bt1)
    100
    bt1.beg1()
    1001
    bt1.filename()
    'chr1_1000-1100_100'
    """
    @property
    def length(self):
        return int(self.end) - int(self.beg)

    @property
    def beg1(self):
        return str(int(self.beg) + 1)

    @property
    def filename(self):
        return '%s_%s-%s_%s' % (self.scaffold, self.beg, self.end, self.length)

    @property
    def bedlist(self):
        return [(self.scaffold, self.beg, self.end)]

    @classmethod
    def from_coord(cls, coord):
        """
        :param coord: scaffold:start-end or scaffold:start_length (any commmas may be present and will be removed)
        :return: a tuple of: a list of 1 3cols-bed-tuple, and a filename
        works with start and end:
        >>> coord_to_bedtuple_and_filename('chr6:136640001-136680000')
        ([('chr6', '136640000', '136680000')], 'chr6+136640001-136680000_40000.bed')
        and also works with start and length:
        >>> coord_to_bedtuple_and_filename('chr6:136640001_40000')
        ([('chr6', '136640000', '136680000')], 'chr6+136640001-136680000_40000.bed')
        """
        coord = coord.replace(",", "")
        has_dash = '-' in coord
        has_under = '_' in coord
        assert(has_dash ^ has_under)
        if has_dash:
            start1, end = coord.split(':')[1].split('-')
            length = str(int(end) - int(start1) + 1)
        if has_under:
            start1, length = coord.split(':')[1].split('_')
            end = str(int(start1) + int(length) - 1)
        chrom = coord.split(':')[0]
        start0 = str(int(start1) - 1)    # bed coords are zero-based
        return cls._make([chrom, start0, end])



def filename_from_coord(coord):
    return Bedtuple.from_coord(coord).filename


def bedlist_from_coord(coord):
    return Bedtuple.from_coord(coord).bedlist

