from flood import FloodFill
import numpy as np
import time 
from tqdm import trange 
import matplotlib.pyplot as plt 

def gen_dummy_data():
    # TODO: Paramterize this so it can change size, perhaps function o_O
    # Generate any # of spk templates (amplitude, width in samples/time, wavelength (trough to trough))
    # Generate spike times that occur for each template across different channels 
    # Scale amplitude across adjacent channels (small on one, big on the other) 
    # Figure out how to convolve spike templates in some order 
    # Think about interpolation
    # Drop the sinusoid, look into mexican hat waveform
    # Think up of edge cases 
    # When it runs well, add noise
    # Gauge performance f(snr)
    # Loss Function = |floodfill(spikes) - actual(spikes)|/actual(spikes)
    # Look at ISI and width of spikes, see if width < ISI 
    return np.asarray([np.sin(np.arange(0,100,.1))]*16).transpose()
    
def gen_dummy_adj():
    # TODO: No clue, this is hard coded
    tmp = []
    for i in range(1,16,2):
        tmp.append([i-1, i])
        tmp.append([i, i-1])
    return np.asarray(tmp)

class Test():
    
    def __init__(self, floodfill=None, v=0):
        self.ff = FloodFill(gen_dummy_data(), gen_dummy_adj(), debug = True)
        self.v = v 
    
    def gen_sin_plots(self):
        
        spk = self.ff.detect_spikes()
        for i in range(spk.shape[1]):
            chanFound = spk[1][i][0][1]
            spkTime = spk[0][i]
            if spkTime > 100 and spkTime < 900:
                fn = "Outputs/Test/Channel {}/spk{}-spkTime{}.png".format(chanFound, i, spkTime)
                fig = self.ff.plotSpk(i, save=True, fname = fn)
                plt.close()
        plt.close('all')
    
    def gen_inj_plots(self):
        pass
    
    