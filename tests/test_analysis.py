from unittest import TestCase


from buffet.analysis import *




class TestAnalysis(TestCase):

    def create_enzyme_does_not_crash(self):

        self.assertTrue(
            str(create_enzyme("ScrFI") + create_enzyme("HpaII") + create_enzyme("BfaI"))
            ==
            "RestrictionBatch(['BfaI', 'HpaII', 'ScrFI'])"
        )
