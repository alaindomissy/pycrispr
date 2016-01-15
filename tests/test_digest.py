from buffet.digest import *


assert(
    coord_to_bedtuple_filename('chr6:136640001-136680000')
    ==
    ([('chr6', '136640000', '136680000')], 'chr6:136640001-136680000_40000')
)

assert(
    stretch_to_bedtuple_filename('chr6:136640001_40000')
    ==
    ([('chr6', '136640000', '136680000')], 'chr6:136640001-136680000_40000')
)