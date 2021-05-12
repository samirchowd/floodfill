class data_gen():
    
    def __init__():
        pass
    
    # Sinusodal data
    def gen_dummy_data():
        # TODO: Paramterize this so it can change size, perhaps function o_O
        return np.asarray([np.sin(np.arange(0,100,.1))]*16).transpose()

    def gen_dummy_adj():
        # TODO: No clue, this is hard coded
        tmp = []
        for i in range(1,16,2):
            tmp.append([i-1, i])
            tmp.append([i, i-1])
        return np.asarray(tmp)
    
    # Generating spike times 
    def gen_spike_times(numSpks, lim, numTemplates, numChans):
        """numSpks: # of spikes, lim: length of data, numTemplates: # of templates, numChans: # of channels"""
        spkTimes = []
        spkChans = []
        for i in range(numTemplates):
            spkTimes.append(np.random.choice(range(lim), numSpks, replace = False))
            spkChans.append(np.random.choice(range(numChans), numSpks, replace=True))
        return np.asarray(spkTimes), np.asarray(spkChans)

    # Generating spike template
    def gen_spike_templates(numSpks,muWidth=90,sigmaWidth=10, muFac=0.1, sigmaFac=0.02, muVert=200, sigmaVert=40):
        """
        numSpks: # of Spikes, muWidth: mean width of spike, sigmaWidth: stdev of spike width, muFac: mean width fator, 
        sigmaFac: stdev of spike factor, muVert: mean vertical scaling factor, sigmaVert: stdev of vertical scaling factor
        """
        from scipy import signal 

        waveforms = []

        for i in range(numSpks):
            width = sigmaWidth * np.random.randn() + muWidth 
            factor = (np.abs(sigmaFac*np.random.randn()+muFac))*width
            polarity = np.random.choice([1,-1], 1)[0]
            vertical = sigmaVert*np.random.randn()+muVert
            # Add a linear equations on top | super impose multiple rickers
            waveforms.append(polarity*signal.ricker(width, factor)*vertical)
        return np.asarray(waveforms, dtype=object)
    
    def generate_hat_data(numSpks, lim, numTemplates, numChans, muWidth=90, sigmaWidth=10, muFac=0.1, sigmaFac=0.02, muVert=25, sigmaVert=5, seed = False):
        if seed:
            np.random.seed(1738)

        times, chans = gen_spike_times(numSpks, lim, numTemplates, numChans)

        waveforms = gen_spike_templates(numTemplates, muWidth, sigmaWidth, muFac, sigmaFac, muVert, sigmaVert)

        # Generating zeros matrix 
        tmp = np.zeros((numChans, lim))
        spikes_inserted = []

        for i in range(numTemplates):
            for j in range(numSpks):
                wave = waveforms[i]
                spkTime = times[i][j]
                spkChan = chans[i][j]
                radius = wave.shape[0]//2
                if spkTime - radius > 0 and spkTime + radius + 1 < lim: 
                    if wave.shape[0] % 2 == 0:
                        tmp[spkChan][spkTime-radius:spkTime+radius] += wave
                    else:
                        tmp[spkChan][spkTime-radius:spkTime+radius+1] += wave

                    spikes_inserted.append([i,spkTime,spkChan])
        return tmp, np.asarray(spikes_inserted)
    
    def gen_rel_data(N, adj, nChans, nSpks, nTemplates, offsetMu = 40, offsetSig = 5, decayMu = 0.8, decaySig = 0.05, noiseMu = 0, noiseSig = 1, seed=False, noise=True):
        if seed:
            np.random.seed(1738)

        # Setup Steps
        data = np.zeros(shape=(nChans, N))

        # Generate Waveforms 
        waveforms = gen_spike_templates(nTemplates)

        # Generate Spiketimes 
        # Creating lead - follow pairs
        chanAssign = np.random.choice(range(nChans), size=nChans, replace=False)
    #     for i in range(nTemplates):
    #         chanAssign.append(np.random.choice(range(nChans), replace=True))

        # Generating sets of spiking channels
        chans = []
        for x in chanAssign:
            tmp = [x]
            for pair in adj: 
                if pair[0] == x:
                    tmp.append(pair[1])
            chans.append(tmp)

        # Assigning times to those pairs
        times = []
        for pair in chans:
            leadTimes = np.random.choice(range(100, N-100), nSpks, replace = False)
            offset = int(offsetSig*np.random.randn()+offsetMu)
            followTimes = leadTimes + offset
            times.append((leadTimes, followTimes))

        # Generating Data
        for i in range(nTemplates): 
            wave = waveforms[i]
            radius = wave.shape[0] // 2
            leadChan = chans[i][0]
            followChans = chans[i][1:]
            # Lead
            leadTimes = times[i][0]
            for t in leadTimes:
                if wave.shape[0] % 2 == 0:
                    data[leadChan][t-radius:t+radius] += wave
                else:
                    data[leadChan][t-radius:t+radius+1] += wave
            # Follow
            followTimes = times[i][1]
            decay = decaySig*np.random.randn() + decayMu
            for t in followTimes:
                if wave.shape[0] % 2 == 0:
                    for followChan in followChans:
                        data[followChan][t-radius:t+radius] += wave*decay
                else:
                    for followChan in followChans:
                        data[followChan][t-radius:t+radius+1] += wave*decay

            if noise:
                for channel in data:
                    channel += noiseSig*np.random.randn(N,) + noiseMu
        return data
            