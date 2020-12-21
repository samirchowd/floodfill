import numpy as np

class FloodFill():

    def __init__(self, data, adj, bandpass=(300,7500)):
        self.data = data
        self.adj = adj
        self.bandpass = bandpass

    def detectSpikes(self, weak=[], strong=[]):

