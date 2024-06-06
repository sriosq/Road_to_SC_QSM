class SegmentationLabel:
    
    def __init__(self, label_id, name=None, susceptibility=None):
        # we will recieve a 3D so we will follow the convention of labeling for full body scans

        self.spinal_cord 

        self.label_id = label_id
        self.name = name
        self.susceptibility = susceptibility
        
    
    def set_name(self, name):
        self.name = name
    
    def set_susceptibility(self, susceptibility):
        self.susceptibility = susceptibility
    
    def __repr__(self):
        return f"SegmentationLabel(label_id={self.label_id}, name={self.name}, susceptibility={self.susceptibility})"
