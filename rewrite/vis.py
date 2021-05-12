import numpy as np 
import matplotlib.pyplot as plt

class vis():
    
    def __init__(self, ff):
        self.ff = ff
    
    
    def plot_waveforms(self):
        wave = self.ff.data.transpose()
        fig, axs = plt.subplots(4,4)
        fig.set_figheight(20)
        fig.set_figwidth(20)
        plt.grid=True
    
        # Loop variables 
        count, minY, maxY = 0, float('inf'), -float('inf') 

        for i in range(4):
            for j in range(4):
                axs[i][j].plot(wave[count])
                minY = min(min(wave[count]),minY)
                maxY = max(max(wave[count]),maxY)
                count += 1

        # Updating Y scale 
        for i in range(4):
            for j in range(4):
                axs[i][j].set_ylim((minY-1, maxY+1))

        plt.show()
        return fig, axs


    def plot_data(self):
        data = self.ff.data.transpose()
        fig, axs = plt.subplots(4,4)
        fig.set_figheight(20)
        fig.set_figwidth(20)
        plt.grid=True

        # Loop variables 
        count, minY, maxY = 0, float('inf'), -float('inf') 

        for i in range(4):
            for j in range(4):
                axs[i][j].plot(data[count])
                axs[i][j].set_title("Channel {}".format(count))
                minY = min(min(data[count]),minY)
                maxY = max(max(data[count]),maxY)
                count+= 1 

        # Updating Y scale 
        for i in range(4):
            for j in range(4):
                axs[i][j].set_ylim((minY-1, maxY+1))

        plt.show()
        return fig, axs
    
    
    def plot_waves(self, N, wavDur):
        wf = self.ff.waves[N]
        curr = 0
        for i in range(self.ff.data.shape[1]):
            plt.plot(range(curr,curr+wavDur), wf[curr:curr+wavDur])
            curr += wavDur
        
        
    def plot_spk(self, N, threshold=True, save=False, wf = False, seed=False, fname=""):
        # Setting up input
        if seed:
            np.random.seed(1738)
        if not fname: 
            fname = "spk_{}".format(N)
        data = self.ff.data.transpose()
        spkTime, spkLoc = self.ff.spk[0][N], self.ff.spk[1][N]
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
                    if wf: 
                        axs[i][j].plot(np.random.randn(200), 'k')
                    else: 
                        axs[i][j].plot(data[count][spkTime-100:spkTime+100],'r--')
                
                # Plotting detected points 
                detected_points = spkLoc[np.where(count == spkLoc[:,1])][:,0]
                detected_data = data[count][detected_points]
                adjusted_xVals = detected_points - spkTime + 100 
                axs[i][j].scatter(adjusted_xVals, detected_data)
                #print("CHANNEL {}: ".format(count), detected_points - spkTime + 100)

                # Plotting thresholds 
                if threshold:
                    axs[i][j].plot(np.arange(200), [self.ff.weak]*200)
                    axs[i][j].plot(np.arange(200), [self.ff.strong]*200)
                    axs[i][j].plot(np.arange(200), [-1*self.ff.weak]*200)
                    axs[i][j].plot(np.arange(200), [-1*self.ff.strong]*200)

                # Updating min/max Y vals to keep uniform scale 
                minY = min(min(data[count][spkTime-100:spkTime+100]),minY)
                maxY = max(max(data[count][spkTime-100:spkTime+100]),maxY)
                
                # Labeling each axis with its channel name 
                axs[i][j].set_title("Channel {}".format(count))
                
                # Plotting Spike Center 
                if count == chanFound: 
                    axs[i][j].scatter(100,data[chanFound][spkTime],marker='x')
                
                # Updating counter
                count+=1
        
        # Updating Y scale 
        for i in range(4):
            for j in range(4):
                axs[i][j].set_ylim((minY-1, maxY+1))
        
        if save:
            plt.savefig(fname)
        
        return fig, axs