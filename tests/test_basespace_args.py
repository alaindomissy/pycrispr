import sys,os
sys.path.insert(0,os.path.abspath(__file__+"/../.."))

from pybasespace.basespace_args import *


def test_json_depcopy():
    assert(json_deepcopy(metadatatemplate) == json_deepcopy(metadatatemplate))
    assert(json_deepcopy(metadatatemplate) is not json_deepcopy(metadatatemplate))


