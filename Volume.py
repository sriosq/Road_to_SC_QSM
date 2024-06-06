#Dependencies
import numpy as np
from label import SegmentationLabel

class Volume:
    #establcer relacion volumen label
    
    def __init__(self, volume):
        # we will recieve a 3D so we will follow the convention of labeling for full body scans
        # The 3D volume data will be a np array? or maybe it should be a nifti file and then I can manually
        # get_fdata?

        self.volume = volume
        self.uniq_labels = np.unique(volume)
        self.segmentation_labels = {}
        self.create_segmentation_labels()

    #def set_label_name(self, label_id, name):
        #for label in self.segmentation_labels:
            #if label.label_id == label_id:
                #label.set_name(name)
                
    def create_segmentation_labels(self):
        for label_id in self.uniq_labels:
            self.segmentation_labels[label_id] = SegmentationLabel(label_id)

    def set_label_name(self, label_id, name):
        if label_id in self.segmentation_labels:
            self.segmentation_labels[label_id].set_name(name)
        else:
            print(f"Label ID {label_id} not found.")

    def set_label_susceptibility(self, label_id, susceptibility):
        if label_id in self.segmentation_labels:
            self.segmentation_labels[label_id].set_susceptibility(susceptibility)
        else:
            print(f"Label ID {label_id} not found.")
 
    def display(self):
        pass

    def manual_labeling(self):
        #To iterate over self unique labels we need its total length
        for i in range(len(self.uniq_labels)):
            name = input(f"Enter name for label #{i}: ") 
            label = SegmentationLabel(i, name)
            self.segmentation_labels.append(label)
    
    def set_susceptibility(self):
        for i in self.segmentation_labels:
            name = self.segmentation_labels[i].name
            sus = input(f"Enter the Susceptibility for {name}")
            i.set_susceptibility(sus)


    def show_labels(self):
        for i in range(len(self.uniq_labels)):
            print(i)

    def __repr__(self):
        return f"SegmentationLabelManager == Volume"






