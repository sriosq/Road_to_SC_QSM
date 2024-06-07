import nibabel as nib
import os
import numpy as np

from label import SegmentationLabel
from Volume import Volume


np.set_printoptions(precision=2, suppress = True)
 
#ct_wb = nib.load("data/ct_wb.nii.gz")
ct_wb_seg = nib.load("/home/soothsayer7/Downloads/totalseg.nii.gz")

# This variables are instances of the nibabel image
# To work with the data:
#ct_wb_data = ct_wb.get_fdata()
labels = ct_wb_seg.get_fdata()


a = np.arange(22,48)
b = np.arange(89,114)
c = np.append(a,[9,10,11,12,13])
d = np.append(b, [68,69,70,71,74,75,76,88])

final = set(np.append(c,d))

all_labels = set(np.unique(labels))
diff_list = sorted(list(all_labels - final))
#print(diff_list)
#print(len(all_labels))
#print(len(final))
#print(len(diff_list))

# Segmentating the fat

