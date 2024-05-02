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
nprocs = np.array([2,4,8,16])
tPar = np.array([80.4538,69.8730,55.4219,59.8408])

spdUp = tSer / tPar
print(spdUp)

eff = spdUp / nprocs
print(eff)

# e = (1/spdUp - 1/nprocs) / (1 - 1/nprocs)
e = (nprocs/spdUp - 1) / (nprocs - 1)
print(e)

plt.figure(figsize=(8,6))
plt.plot(nprocs,tPar,c='r',lw=1, marker='o', ms=5, mec='k', mfc='k')
plt.xlabel(r"Number of processors $p$",size=14)
plt.ylabel(r"Time taken $t$ (in s)",size=14)
plt.title("Parallel runtime",size=16)
plt.grid()
fig1 = plt.gcf()
plt.show()
# fig1.savefig("./Res/T1.png", dpi=300, bbox_inches='tight')


plt.figure(figsize=(8,6))
plt.plot(nprocs,spdUp,c='r',lw=1, marker='o', ms=5, mec='k', mfc='k')
plt.xlabel(r"Number of processors $p$",size=14)
plt.ylabel(r"Speed-up $\psi (n,p)$",size=14)
plt.title("Speed-up vs. Number of processors",size=16)
plt.grid()
fig2 = plt.gcf()
plt.show()
# fig2.savefig("./Res/T2.png", dpi=300, bbox_inches='tight')


plt.figure(figsize=(8,6))
plt.plot(nprocs,eff,c='r',lw=1, marker='o', ms=5, mec='k', mfc='k')
plt.xlabel(r"Number of processors $p$",size=14)
plt.ylabel(r"Efficiency $\epsilon (n,p)$",size=14)
plt.title("Efficiency vs. Number of processors",size=16)
plt.grid()
fig3 = plt.gcf()
plt.show()
# fig3.savefig("./Res/T3.png", dpi=300, bbox_inches='tight')


plt.figure(figsize=(8,6))
plt.plot(nprocs,spdUp,c='r',lw=1, marker='o', ms=5, mec='k', mfc='k')
plt.xlabel(r"Number of processors $p$",size=14)
plt.ylabel(r"Karp–Flatt Metric $e$",size=14)
plt.title("Karp–Flatt Metric vs. Number of processors",size=16)
plt.grid()
fig4 = plt.gcf()
plt.show()
# fig4.savefig("./Res/T4.png", dpi=300, bbox_inches='tight')