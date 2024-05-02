## CODE TO PLOT OTHER RESULTS, LIKE TIME

import numpy as np
import matplotlib.pyplot as plt

# Details... May not be needed directly
Nl = 1000
u = 20
c = 10
alp = 0
dt = 0.005
tf = 5
Nt = int(tf/dt) + 1

tSer = 123.6224
np = np.array([2,4,8,16])
tPar = np.array([80.4538, 69.8730, 55.4219, 59.8408])