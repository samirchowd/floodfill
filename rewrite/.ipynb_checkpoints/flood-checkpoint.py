import numpy as np
# Remove namespace np 

class FloodFill():

    def __init__(self, data, adj, bandpass=(300,7500)):
        self.data = data
        self.adj = adj
        self.bandpass = bandpass

    def detect_spikes(self, weakMul=2, strongMul=4):
        # Initializing objkect varibles 
        #BC: Change these to object vars, update rest of script 
        self.spk = [] 
        self.sigma = np.median(np.absolute(self.data))/0.6745
        self.weak = weakMul*sigma
        self.strong = strongMul*sigma
        
        # Initializing Local Variable 
        res_bin = np.zeros(data.shape)
        
        # Finding positive crossings 
        ix = cross_detect(self, strong)
        
        # Iterate through each positive crossing 
        for i in range(len(ix)):
            for j in range(len(ix[i])):
                t = ix[i][j] # Time point
                if not res_bin(t+1, i):
                    # Initialize result variables
                    spk_wt = []
                    spk_loc = []
                    # Call flood_fill on strong crossing 
                    spk_wt, spk_loc, res_bin = floodfill(self, t+1, i, weak, strong, spk_wt, spk_loc, res_bin)
                    
                    # Calculate Voltage Weighted Spk Center 
                    spk_wt =  spk_center(spk_wt)
                    
                    # Append results to spk object 
                    spk.append((spk_wt, spk_loc))
                    
        return spk 
        
    def flood_fill(self, t, c, weak, strong, spk_wt, spk_loc, res_bin):
        
        # Validating data 
        if validate_data(self, t, c, weak, res_bin):
            return
        
        # Calulating psi value and appending it and the (t,c) to matricies 
        # BC: Check append memory usage 
        spk_wt.append((t, psi(self, t, c, weak, strong)))
        
        
    def psi(self, t, c, weak, strong):
        return min(((-self.data[t,c] - weak) / (strong - weak)), 1);
    
    def spk_center(spk_wt):
        x = []
        y = [] 
        for tc 
    
        
    def cross_detect(self, strong):
        index = np.diff(np.sign(self.data-strong), axis=0) > 0
        ix = []
        # BC: Replace for loops with list-comprehension / lambda 
        for i in range(index.shape[1]):
            ix.append(np.where(index[:,i] == True))

        for i in range(len(ix)):
            ix[i] = ix[i][0]
        
        return ix
    
    # BC: Move this back into flood_fill, try using an assert 
    def validate_data(self, t, c, weak, res_bin):
        chanBound = self.data.shape[0]
        timeBound = self.data.shape[1] 
        
        if t < 0 or t > timeBound or c < 0 or c > chanBound:
            return False 
        
        if self.data[t,c] <= weak or res_bin[t,c]:
            return False
        
        return True 
        