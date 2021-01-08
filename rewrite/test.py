from flood import FloodFill
import numpy as np
import time 
from tqdm import trange 

adj =  np.genfromtxt('adj.csv', delimiter=',',dtype=np.int32)
data = np.genfromtxt('data.csv', delimiter=',')

ff  = FloodFill(data, adj) 

class Test():
    
    def __init__(self, floodfill=None, v=0)
        if floodfill:
            self.ff = FloodFill(floodfill)
        else:
            self.ff = FloodFill(gen_dummy_data(), gen_dummy_adj)
        self.v = v 
    
    def check_detect_spikes(self):
        pass
    
    def check_flood_fill(self):
        pass
    
    def check_psi(self):
        pass
    
    def check_spk_center(self):
        pass
    
    def check_cross_detect(self):
        pass
    
    def check_validate_data(self):
        # Testing data is out of bounds 
        
    # Helper Methods
    def gen_dummy_data():
        # TODO: Paramterize this so it can change size, perhaps function o_O
        return np.asarray([np.cos(np.arange(0,100,.1))]*16).transpose()
    
    def gen_dummy_adj():
        # TODO: No clue, this is hard coded
        tmp = []
        for i in range(1,16,2):
            tmp.append([i-1, i])
            tmp.append([i, i-1])
        return np.asarray(tmp)
        
    
        
        
    

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
        # TODO: Add test params to this list 
        print("Test object initiated with following params: ")

    if v: 
        print("Running tests...")
    
    # check_detect_spiks()
    if v:
        print("Testing floodfill.detect_spikes...")
    try:
        pass
    except AssertionError: 
        pass 
    
    # check_flood_fill()
    if v: 
        print("Testing floodfill.flood_fill...")
    try:
        pass
    except AssertionError:
        pass
    
    # check_psi()
    if v:
        print("Testing floodfill.psi...")
    try:
        pass
    except AssertionError:
        pass
    
    # check_spk_center()
    if v:
        print("Testing floodfill.spk_center...")
    try:
        pass
    except AssertionError:
        pass 
    
    # check_cross_detect()
    if v:
        print("Testing floodfill.cross_detect...")
    try:
        pass
    except AssertionError:
        pass 
    
    # check_validate_data()
    if v:
        print("Testing floodfill.validate_data...")
    try:
        pass
    except AssertionError:
        pass 