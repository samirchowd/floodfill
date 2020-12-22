import numpy as np

class FloodFill():

    def __init__(self, data, adj, bandpass=(300,7500)):
        self.data = data
        self.adj = adj
        self.bandpass = bandpass

    def detectSpikes(self, weak=0, strong=0):
        # Initializing objkect varibles 
        spk = [] 
        if not weak and not strong
            sigma = np.median(np.absolute(self.data))/0.6745
            weak = 2*sigma
            strong = 2*sigma
        
        # Initializing Local Variables 
        res_bin = np.zeros(data.shape)
        
        
