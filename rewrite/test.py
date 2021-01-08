from flood import FloodFill
import numpy as np
import time 
from tqdm import trange 

adj =  np.genfromtxt('adj.csv', delimiter=',',dtype=np.int32)
data = np.genfromtxt('data.csv', delimiter=',')

class Test():
    
    def __init__(self):
        pass