import flood 
import numpy as np

adj =  np.genfromtxt('adj.csv', delimiter=',')
data = np.genfromtxt('data.csv', delimiter=',')

print(type(data), type(adj))