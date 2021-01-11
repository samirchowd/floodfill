from flood import FloodFill
import numpy as np
import time 
from tqdm import trange 
import matplotlib.pyplot as plt 

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

class Test():
    
    def __init__(self, floodfill=None, v=0):
        self.ff = FloodFill(gen_dummy_data(), gen_dummy_adj(), debug = True)
        self.v = v 
    
    def check_detect_spikes(self):
        spk = self.ff.detect_spikes()
        fig = self.ff.plotSpk(100)
        plt.show()
    
    def check_flood_fill(self):
        pass
    
    def check_psi(self):
        # Checking < 1
        assert self.ff.psi(0, 0) == -0.5 
        
        # Checking > 1 
        assert self.ff.psi(50,0) == 1
    
    def check_spk_center(self):
        pass
    
    def check_cross_detect(self):
        pass
    
    def check_validate_data(self):
        pass        
    

if __name__ == "__main__":
    
    # Running test cases 
    print("Running test cases validating FloodFill output")
    
    # Getting verbosity 
    v = 0
    v = int(input("Please enter 1 for verbose output: "))
    
    # Instantiating test object
    if v:
        print("Instantiating test object")
    test = Test(v = 1) 

    if v: 
        print("Running tests...")
    
    # check_detect_spikes()
    if v:
        print("Testing floodfill.detect_spikes...")
    try:
        test.check_detect_spikes()
    except AssertionError: 
        print("Failed on floodfill.detect_spikes...")
    
    # check_flood_fill()
    if v: 
        print("Testing floodfill.flood_fill...")
    try:
        pass
    except AssertionError:
        print("Failed on floodfill.flood_fill...")
    
    # check_psi()
    if v:
        print("Testing floodfill.psi...")
    try:
        test.check_psi()
    except AssertionError:
        print("Failed on floodfill.psi...")
    
    # check_spk_center()
    if v:
        print("Testing floodfill.spk_center...")
    try:
        pass
    except AssertionError:
        print("Failed on floodfill.spk_center...")
    
    # check_cross_detect()
    if v:
        print("Testing floodfill.cross_detect...")
    try:
        pass
    except AssertionError:
        print("Failed on floodfill.check_cross_detect...")
    
    # check_validate_data()
    if v:
        print("Testing floodfill.validate_data...")
    try:
        pass
    except AssertionError:
        print("Failed on testing floodfill.validate_data")