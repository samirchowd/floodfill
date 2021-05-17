import numpy as np 
import matplotlib.pyplot as plt

class vis():
    
    def __init__(self, ff, waves=None):
        self.ff = ff
        self.wf = waves

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
        wf = self.wf.waves[N]
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
    
    def pca_vis(self):
        from sklearn.decomposition import PCA
        waves = self.wf.waves
        pca = PCA(n_components=3, whiten= True)
        pca.fit(waves)
        pca_data = pca.transform(waves)
        plt.scatter(pca_data[:, 0], pca_data[:, 1], 10, alpha = 0.2)
        plt.show()
        
    def clus_vis(self, K):
        plt.figure()
        labels = self.wf.clus.labels_
        
        from sklearn.decomposition import PCA
        pca = PCA(n_components=3, whiten= True)
        pca.fit(self.wf.waves)
        pca_data = pca.transform(self.wf.waves)
        
        for i in range(K):
            ixs = labels == i 
            x = pca_data[ixs,0]
            y = pca_data[ixs,1]
            plt.scatter(x,y)
    
    def data_spk_vis(self):
        plt.figure() 
        rowHeight = 50
        rowOffset = 60 
        for i in range(self.ff.data.shape[1]):
            tmpData = self.ff.data[:, i]
            maxH, minH = max(tmpData), min(tmpData)
            b = (tmpData - np.min(tmpData))/np.ptp(tmpData)
            b *= rowHeight
            b += i * rowOffset
            plt.plot(b)
        
        labels = self.wf.clus.labels_
        for i in range(5):
            ixs = labels == i 
            spkOfI = self.ff.spk[0][ixs]
            plt.scatter(spkOfI, [self.ff.data.shape[1]*rowOffset*1.1+(i*25)]*len(spkOfI))

        x = [a[0][0] for a in self.ff.spk[1]]
        y = [a[0][1]*rowOffset+rowHeight/2 for a in self.ff.spk[1]]
        plt.scatter(x,y, marker ='+')
        
    def elbow_plot(self):
        waves = self.wf.waves
        distortions = []
        K = range(1,10)
        for k in K:
            kmeanModel = KMeans(n_clusters=k)
            kmeanModel.fit(waves)
            distortions.append(kmeanModel.inertia_)
 
        plt.plot(K, distortions, 'bx-')
        plt.xlabel('k')
        plt.ylabel('Distortion')
        plt.title('The Elbow Method showing the optimal k')
        plt.show()
    