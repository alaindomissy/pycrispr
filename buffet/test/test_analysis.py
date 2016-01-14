from buffet.analysis import *


assert(
    str(create_enzyme("ScrFI") + create_enzyme("HpaII") + create_enzyme("BfaI"))
    ==
    "RestrictionBatch(['BfaI', 'HpaII', 'ScrFI'])"
)
