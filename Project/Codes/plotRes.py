## CODE FOR PLOTTING THE WAKE STRUCTURE

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

file1 = "./Res/Ser.txt"
dat = np.loadtxt(file1, max_rows=8)
Nl, Nt = map(int, dat[0:2])
dt, dx, rho, u, c, alp = dat[2:]

x0 = np.loadtxt(file1, max_rows=Nt, skiprows=9, usecols=0 ,delimiter=',')
y0 = np.loadtxt(file1, max_rows=Nt, skiprows=9, usecols=1 ,delimiter=',')

L = np.loadtxt(file1, max_rows=Nt, skiprows=(Nt+10), usecols=0 ,delimiter=',')
D = np.loadtxt(file1, max_rows=Nt, skiprows=(Nt+10), usecols=1 ,delimiter=',')

xwM = np.zeros((Nt,Nt))
ywM = np.zeros((Nt,Nt))

for i in range(Nt):
    sprw = (2 + i)*Nt + 11 + i # (2*Nt + 11) + i*(Nt+1)
    xwM[i][:] = np.loadtxt(file1, max_rows=Nt, skiprows=sprw, usecols=0 ,delimiter=',')
    ywM[i][:] = np.loadtxt(file1, max_rows=Nt, skiprows=sprw, usecols=1 ,delimiter=',')
# print(xwM[Nt-1][-10:],ywM[Nt-1][-10:])

fig1 = plt.figure(figsize=(8,6))
pdats = plt.scatter([])
pdat = pdats[0]

def animFunc(frame):
    pdat.set_data((xwM[frame][:], ywM[frame][:]))

anim = FuncAnimation(fig1, animFunc, frames=Nt, interval=100)
plt.show()