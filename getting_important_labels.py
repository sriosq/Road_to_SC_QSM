import os
import numpy as np
import nibabel as nib


path_to_nii = "/home/soothsayer7/Downloads/expire_torso.nii"

img = nib.load(path_to_nii)

data = img.get_fdata()

print(np.unique(data))


