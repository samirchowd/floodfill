import numpy as np
import matplotlib.pyplot as plt 

class FloodFill():

    def __init__(self, data, adj, p = 1, bandpass=(300,7500)):
        self.data = data
        self.adj = adj
        self.p = 1
        self.bandpass = bandpass

    def detect_spikes(self, weakMul=2, strongMul=4):
        # Initializing object varibles 
        self.spk = [] 
        self.sigma = np.median(np.absolute(self.data))/0.6745
        self.weak = weakMul*self.sigma
        self.strong = strongMul*self.sigma
        
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
        
    def flood_fill(self, t, c, spk_wt, spk_loc, res_bin):
        
        # TODO: SWITCH TO ASSERTS 
        # Validating data  
        if not self.validate_data(t, c, res_bin):
            return spk_wt, spk_loc, res_bin 
        
        # Calulating psi value and appending it and the (t,c) to matricies 
        spk_wt.append((t, self.psi(t, c)))
        spk_loc.append((t,c))
        
        res_bin[t][c] = 1
        
        # BC: Recursion on hold until alg verified 
        # TODO: Recursive Calls then switch to iterative 
        # DIG DEEP INTO THIS 
        for tc in self.adj: 
            if c == (int(tc[0])-1): 
                if not res_bin[t][tc[1]-1]:
                    spk_wt, spk_loc, res_bin = self.flood_fill(t, tc[1]-1, spk_wt, spk_loc, res_bin)
        
        
        # TODO: Check these outputs for continuity '
        # TODO: Refractory period = 1ms (30 samples), pass percentage of refr period and loop through checking that period 
        # Sub weak from data points to check (forward and backward), check only weak crossers
        if not res_bin[t+1][c]:
            spk_wt, spk_loc, res_bin = self.flood_fill(t+1, c, spk_wt, spk_loc, res_bin)
        
        if not res_bin[t-1][c]:
            spk_wt, spk_loc, res_bin = self.flood_fill(t-1, c, spk_wt, spk_loc, res_bin)
            
        return spk_wt, spk_loc, res_bin 
        
    def psi(self, t, c):
        return min(((-self.data[t,c] - self.weak) / (self.strong - self.weak)), 1);
    
    def spk_center(self, spk_wt):
        x = np.asarray([x[0] for x in spk_wt])
        y = np.asarray([y[1] for y in spk_wt])
        return np.sum(np.power(x*y, self.p)) / np.sum(np.power(y, self.p))
        
    def cross_detect(self):
        index = np.diff(np.sign(self.data-self.strong), axis=0) > 0
        
        return np.asarray([np.where(i == True)[0] for i in index.transpose()])
    
    # BC: Move this back into flood_fill, try using an assert 
    # SC: Assert will throw an error, we can try-catch the error. Is this something we want to do? 
    def validate_data(self, t, c, res_bin):
        timeBound = self.data.shape[0]-1
        chanBound = self.data.shape[1]-1
        if t < 0 or t > timeBound or c < 0 or c > chanBound:
            return False 
        
        if self.data[t,c] <= self.weak or res_bin[t,c]:
            return False
        
        return True 
        
    def plotSpk(self, N, threshold=True, save=False, fname=""):
        # Setting up input
        if not fname: 
            fname = "spk_{}".format(N)
        data = self.data.transpose()
        spkTime, spkLoc = self.spk[0][N], self.spk[1][N]
        chanFound = spkLoc[0][1]
        spkTime = int(spkTime)
        spkChans = set([x[1] for x in spkLoc])
        
        # Setting up grid 
        fig, axs = plt.subplots(4,4)
        fig.set_figheight(20)
        fig.set_figwidth(20)
        
        # Loop variables 
        count, minY, maxY = 0, float('inf'), -float('inf') 
        
        # Looping through all the channels 
        for i in range(4):
            for j in range(4):
                # Black or Red if spiking or not 
                if count in spkChans:
                    axs[i][j].plot(data[count][spkTime-100:spkTime+100],'k')
                else: 
                    axs[i][j].plot(data[count][spkTime-100:spkTime+100],'r--')
                
                # TODO: Add channel # to each plot 
                # Plotting detected points 
                detected_points = spkLoc[np.where(count == spkLoc[:,1])][:,0]
                detected_data = data[count][detected_points]
                adjusted_xVals = detected_points - spkTime + 100 
                axs[i][j].scatter(adjusted_xVals, detected_data)
                #print("CHANNEL {}: ".format(count), detected_points - spkTime + 100)

                # Plotting thresholds 
                if threshold:
                    axs[i][j].plot(np.arange(200), [self.weak]*200)
                    axs[i][j].plot(np.arange(200), [self.strong]*200)

                # Updating min/max Y vals to keep uniform scale 
                minY = min(min(data[count][spkTime-100:spkTime+100]),minY)
                maxY = max(max(data[count][spkTime-100:spkTime+100]),maxY)
                count+=1
        
        # Updating Y scale 
        for i in range(4):
            for j in range(4):
                axs[i][j].set_ylim((minY-1, maxY+1))
        
        if save:
            plt.savefig(fname)
    
"""
THINGS TO CHECK

IF ALG IS WRONG
Check return vals
Check validation case 
* Print out all the values that are being passed
* By hand, compare them to what should not count, and what should. Check to see if there's any discrepancy 
* print(failed: xyz) or print(good: xyz)
Check if iterating through adjaceny is correct
* print(checking channel x from  channel y)
* Compare that against adj matrix 
Make sure spk_loc is correct 
* take every point in spk loc and compare it to weak, if all > weak good, if not bad 
Check spk_center function (compare by hand)
* do a set of points by hand, then plug it into spk_center, see if the same 
Check the plotting 
* 
Check Strong crossings 
Check PSI 


IF ALG IS CORRECT 
Compare to Rossant
Try going out multiple time poitns
Play with threshold values 
Try absolute value of data (if rossant doesn't specify)

Think about spikes separated in space but not in time, alg cannot seperate 

Write unit tests basically 

TRY:
plot a spk, check for other spks near by (in a refr period?) plot them as well 

"""