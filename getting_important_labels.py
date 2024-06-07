import os
import numpy as np
import nibabel as nib
import numpy as np
import os
from Volume import Volume

path_to_nii = "/home/soothsayer7/Downloads/expire_torso.nii"

img = nib.load(path_to_nii)

data = img.get_fdata()

ct_wb_seg = nib.load("/home/soothsayer7/Downloads/totalseg.nii.gz")

labels = ct_wb_seg.get_fdata()




