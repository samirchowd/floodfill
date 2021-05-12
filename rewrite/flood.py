import numpy as np
import matplotlib.pyplot as plt 

class FloodFill():

    def __init__(self, data, adj, debug=False):
        self.data = data
        self.adj = adj
        self.debug = debug
    
    def preprocess(bandpass=(300,7500)):
        self.bandpass = bandpass
        

    def detect_spikes(self, weakMul=2, strongMul=4, refr=30, p=1):
        '''
        Returns the spike times and locations of a given set of EMG data.

                Parameters:
                        self (FloodFill): A flood fill object with data, an adjacency matrix, and a p-value
                        weakMul (int): An integer that indicates the multiplier for weak threshold crossings
                        strongMul (int): An integer that indicates the multiplier for strong threshold crossings
                        refr (int): An integer that indicates the refractory period between two spikes
                        p (int): #TODO 

                Returns:
                        spk (numpy array): An array where the first element is an array of spike times and the 
                                           second element is an array of arrays of [spike time, spike channel]                                       
        '''
        # Initializing object varibles 
        self.spk = [] 
        self.sigma = np.median(np.absolute(self.data))/0.6745
        self.weak = weakMul*self.sigma
        self.strong = strongMul*self.sigma
        self.refr = refr
        self.p = p 
        
        # DEBUG: for manually setting thresholds 
        if self.debug:
            self.weak = float(input("Please enter a value for the weak threshold: "))
            self.strong = float(input("Please enter a value for the strong threshold: "))
        
        # Initializing local variables 
        res_bin = np.zeros(self.data.shape)
        
        # Finding positive crossings 
        ix = self.cross_detect()
        
        # Iterate through each positive crossing 
        for i in range(len(ix)):
            for j in range(len(ix[i])):
                t = ix[i][j] # Time point
                if not res_bin[t+1][i]:
                    # Initialize result variables
                    spk_wt = []
                    spk_loc = []

                    # Call flood_fill on strong crossing 
                    spk_wt, spk_loc, res_bin = self.flood_fill(t+1, i, spk_wt, spk_loc, res_bin)

                    # Calculate Voltage Weighted Spk Center 
                    spk_wt =  self.spk_center(spk_wt)

                    # Append results to spk object 
                    self.spk.append((spk_wt, spk_loc))
        
        # Breaking up spks into spk times and spk locations 
        spkTimes = np.asarray([spkT[0] for spkT in self.spk])
        spkLocs = np.asarray([spkL[1] for spkL in self.spk])
        
        # Reformatting spk to be a numpy array 
        tmp = np.asarray([np.asarray(x) for x in spkLocs])
        self.spk = np.vstack((spkTimes, tmp))
        
        return self.spk
        
    def flood_fill(self, t, c, spk_wt, spk_loc, res_bin, recursive=True):
        '''
        Returns an updated listed of spike weights, locations, and visited points via the flood fill approach.

                Parameters:
                        self (FloodFill): A flood fill object with data, an adjacency matrix, and a p-value
                        t (int): the time index of the data to check 
                        c (int): the channel index of the data to check 
                        spk_wt (array): an array of psi-values calculated from other (t,c) pairs in the call stack
                        spk_loc (array): an array of arrays of all (t,c) pairs visited from the current call stack
                        res_bin (array): an binary array of the same dimensions of the data marking visited points
                        recursive (boolean) [optional]: boolean value whether to recursively call itself or not.
                Returns:
                        spk_wt (array): an updated array of psi-values calculated from other (t,c) pairs in the call stack
                        spk_loc (array): an updated array of arrays of all (t,c) pairs visited from the current call stack
                        res_bin (array): an updated binary array of the same dimensions of the data marking visited points                     
        '''
        
        # Validating data  
        if not self.validate_data(t, c, res_bin):
            return spk_wt, spk_loc, res_bin 
        
        # Calulating psi value and appending it and the (t,c) to matricies 
        spk_wt.append((t, self.psi(t, c)))
        spk_loc.append((t,c))
        
        res_bin[t][c] = 1
        
        for tc in self.adj: 
            if c == (int(tc[0])) and recursive: 
                if not res_bin[t][tc[1]]:
                    spk_wt, spk_loc, res_bin = self.flood_fill(t, tc[1], spk_wt, spk_loc, res_bin, recursive = True)
        
        for i in range(t-self.refr,t+self.refr):
            if not res_bin[i][c] and recursive:
                spk_wt, spk_loc, res_bin = self.flood_fill(i, c, spk_wt, spk_loc, res_bin, recursive=True)
        
        return spk_wt, spk_loc, res_bin 
        
    def psi(self, t, c):
        '''
        Returns the psi value of a given time and channel value on the data [See Rossant 2016]

                Parameters:
                        self (FloodFill): A flood fill object with data, an adjacency matrix, and a p-value
                        t (int): the time index of the data to check 
                        c (int): the channel index of the data to check 
                        
                Returns:
                        psi (int): The psi value of the given time channel pair                            
        '''
        return min(((-self.data[t,c] - self.weak) / (self.strong - self.weak)), 1);
    
    def spk_center(self, spk_wt):
        '''
        Returns the volatage weighted spike center of a given spike [See Rossant 2016]

                Parameters:
                        self (FloodFill): A flood fill object with data, an adjacency matrix, and a p-value
                        spk_wt (array): a full array of psi-values calculated from (t,c) pairs in a given spike
                        
                Returns:
                        spk_center (int): The value of the voltage weighted spike center of a given spike
                           
        '''
        x = np.asarray([x[0] for x in spk_wt])
        y = np.asarray([y[1] for y in spk_wt])
        return np.sum(np.power(x*y, self.p)) / np.sum(np.power(y, self.p))
        
    def cross_detect(self):
        '''
        Returns the positive threshold crossings of the given data using the strong threshold 

                Parameters:
                        self (FloodFill): A flood fill object with data, an adjacency matrix, and a p-value
                        
                Returns:
                        cross_detect (array): An array containing the indicies of all the positive threshold crossings                       
        '''
        index = np.diff(np.sign(self.data-self.strong), axis=0) > 0
        
        return np.asarray([np.where(i == True)[0] for i in index.transpose()])
    
    def validate_data(self, t, c, res_bin):
        '''
        Returns a boolean value whether the time point crosses the weak threshold or has not been checked before and is within range

                Parameters:
                        self (FloodFill): A flood fill object with data, an adjacency matrix, and a p-value
                        t (int): the time index of the data to check 
                        c (int): the channel index of the data to check 
                        res_bin (array): an binary array of the same dimensions of the data marking visited points
                        
                Returns:
                        validate_data (boolean): A boolean value where true indicates a valid data point to check                         
        '''
        timeBound = self.data.shape[0]-1
        chanBound = self.data.shape[1]-1
        if t < 0 or t > timeBound or c < 0 or c > chanBound:
            return False 
        
        if abs(self.data[t,c]) <= self.weak or res_bin[t,c]:
            return False
        
        return True 
    
    
        
    