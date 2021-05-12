#from sys import argv
import csv
import numpy as np
import pandas as pd
#import os
#from math import radians, cos, sin, asin, sqrt, exp
#import copy
#import random
import pyemd
from pyemd import emd_with_flow

# Load species seasonal abundance distributions (estimated from eBird data)
abundance_BR = pd.read_table('results/output/seasonalAbundance_BR.csv', sep=",")
abundance_NB = pd.read_table('results/output/seasonalAbundance_NB.csv', sep=",")
abundance_BR = abundance_BR.apply(lambda x: x/sum(x))
abundance_NB = abundance_NB.apply(lambda x: x/sum(x))

# Load matrix of pairwise distance between every hexagons on the grid
distanceMatrix = np.loadtxt('results/output/distanceMatrix.csv', delimiter=";")

# Compute optimal redistribution using the Earth Mover's Distance algorithm 
for s in abundance_BR.columns:
  EMD_results = emd_with_flow(np.array(abundance_BR[s]), np.array(abundance_NB[s]), distanceMatrix)
  EMD_results = np.array(EMD_results[1])[(np.array(abundance_BR[s]) > 0),:][:,(np.array(abundance_NB[s]) > 0)]
  np.savetxt("results/output/ORSIM_results_" + s + ".csv", EMD_results, delimiter=',') # Save simulated migratory connectivity 
