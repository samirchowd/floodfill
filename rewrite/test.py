from flood import FloodFill
import numpy as np

adj =  np.genfromtxt('adj.csv', delimiter=',',dtype=np.int32)
data = np.genfromtxt('data.csv', delimiter=',')

spk_detector = FloodFill(data, adj)
spk = spk_detector.detect_spikes()