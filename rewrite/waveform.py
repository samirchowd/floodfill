import numpy as np
class waveforms():
    
    def __init__(self, data, spk_times, spk_chans):
        # Python is Pass-By-Reference, redundant data is fine here
        self.data = data 
        self.spk_centers = spk_times
        self.spk_chans = spk_chans
    
    def get_templates(self):
        '''
        Returns the clustered and realigned waveforms 

                Parameters:
                        self (waveforms): A waveforms object with data, spike centers, and spike channels

                Returns:
                        templates (array): An array containing realigned and clustered templates                                     
        '''
        pass
    
    def get_waves(self, wavDur, seed=False):
        '''
        Returns the waveforms on the data given a wave duration and spike times

                Parameters:
                        self (waveforms): A waveforms object with data, spike centers, and spike channels
                        wavDur (int): An integer indicating the duration of each waveform
                        seed (boolean): A boolean value where true seeds the random generator 
                Returns:
                        waves (array): An array containing each waves for the given data and spike times 
        '''
        # Setting up input
        if seed:
            np.random.seed(1738)
        
        # Setting up waveform array
        waves = [] 
        for N in range(len(self.spk_centers)):
            tmp = []
            spkTime, spkLoc = self.spk_centers[N], self.spk_chans[N]
            spkTime = int(spkTime)
            spkChans = set([x for x in spkLoc])
            for chan in range(self.data.shape[1]):
                if chan in spkChans:
                    tmp.extend(self.data[spkTime-wavDur:spkTime+wavDur,chan])
                else:
                    tmp.extend(np.random.randn(2*wavDur))
            waves.append(tmp)
        self.waves = waves
        return waves
    
    def cluster(self, n):
        '''
        Returns the clustered waveforms via PCA

                Parameters:
                        self (waveforms): A waveforms object with data, spike centers, and spike channels
                        n (int): Integer indicating number of clusters for K-Means 
                        
                Returns:
                        clus (object): Cluster object from sklearn of clustered waveforms  
        '''
        # Performing PCA
        from sklearn.decomposition import PCA
        pca = PCA(n_components=3, whiten= True)
        pca.fit(self.waves)
        pca_data = pca.transform(self.waves)
        
        # Clustering 
        from sklearn.cluster import KMeans
        clus = KMeans(n_clusters=n, random_state=0).fit(pca_data)
        self.clus = clus
        
        return clus 

    def _waveshift(self, n):
        '''
        Returns the realigned waveforms 

                Parameters:
                        self (waveforms): A waveforms object with data, spike centers, and spike channels
                        n (int): Integer indicating number of clusters for K-Means 
                        
                Returns:
                        waveshift (array): array of shifted waveforms 
        '''
        # Better clustering
        # Detect spikes and concatenate channels w/ noise and waveforms 
        # # Double check concatenation, make sure noie replacement is actually happening 
        # Run PCA and K-Means 

        # Cross Correlation Alignment 
        # Find time shift between each pair 
        # Truncate second waveform by offset
        # Add gaussian noise to pad truncated second waveform

        # Checking for merges 
        # Create a template, average across timepoints: One averaged waveform as long as all the waveforms
        # Distance metric between clusters between templates (averaged waveform)
        # threshold, if sim > some threshold, decrease k by 1 (aka merge)

        # Checking for splits
        # Create a matrix 
        # Spk x Spk (where spks are in a given cluster) (# of clus = # of matricies)
        # ith, jth element = distance between ith spk and jth spk in cluster
        # hope that matrix would converge to zero across all values
        # If [metric] > threshold, increase k by 1 (aka split)

        # repeat above from step 2 [until k stays the same or max itr]
        # Keep track of number of changes in assignments 

        waves = self.waves

        clus = self.cluster(n)

        # Cross corr alignments 
        from scipy import signal 
        newWavs = []
        labels = clus.labels_ 

        if vis:
            for i in range(1,max(labels)+1):
                ixs = labels == i 
                x = pca_data[ixs,0]
                y = pca_data[ixs,1]
                plt.scatter(x,y)

        for i in range(max(labels)+1):
            # Getting all the waveforms to a specific cluster 
            ixs = labels == i 
            wavTmp = np.asarray(waves)[ixs]
            y1 = wavTmp[0]
            tmp = [y1]
            for y2 in wavTmp: 
                # Calculating time offset between two signals 
                corr = signal.correlate(y2, y1, mode='same') / np.sqrt(signal.correlate(y1, y1, mode='same')[int(n/2)] * signal.correlate(y2, y2, mode='same')[int(n/2)])
                delay_arr = np.linspace(-0.5*n/sr, 0.5*n/sr, n)
                delay = np.argmax(corr) - len(corr) / 2 
                pad = np.random.randn(np.abs(int(delay)))
                if delay > 0:
                    # Cutting off from the front
                    y2 = y2[int(delay):]
                    y2 = np.concatenate((pad, y2))
                elif delay < 0:
                    # Cutting off from the back
                    delay = np.abs(delay)
                    y2 = y2[:-int(delay)]
                    y2 = np.concatenate((y2, pad))
                tmp.append(y2)
            newWavs.append(tmp)
        return newWavs