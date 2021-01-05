from flood import FloodFill
import numpy as np
import time 
from tqdm import trange 

adj =  np.genfromtxt('adj.csv', delimiter=',',dtype=np.int32)
data = np.genfromtxt('data.csv', delimiter=',')


spk_detector = FloodFill(data, adj)
N = 100 
times = []

for i in (t := trange(N)):
    start = time.time()
    spk_detector.detect_spikes()
    end = time.time()
    times.append(end-start)
    t.set_description("Time took %f" % (end-start) )

plt.plot(times) 
plt.show()