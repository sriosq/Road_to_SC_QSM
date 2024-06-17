
from .. import __dir_testing__
import nibabel as nib
import os

def test_helloworld():
    img = nib.load(os.path.join(__dir_testing__, "final_sc_seg.nii.gz"))
    #assert
    # Use assert to verify that the output is GOOD not only that the code runs
    