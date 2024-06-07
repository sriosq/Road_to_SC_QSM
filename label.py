#Dependencies
import numpy as np

class SegmentationLabel:
    def __init__(self, label_id, name=None, susceptibility=None, ct_number = None):
        self.label_id = label_id
        self.name = name
        self.susceptibility = susceptibility
        self.ct_number = ct_number
    
    def set_name(self, name):
        self.name = name
    
    def set_susceptibility(self, susceptibility):
        self.susceptibility = susceptibility
    
    def __repr__(self):
        return f"SegmentationLabel(label_id={self.label_id}, name={self.name}, susceptibility={self.susceptibility})"