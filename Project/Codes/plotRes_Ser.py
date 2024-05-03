## CODE FOR PLOTTING THE WAKE STRUCTURE - SERIAL CODE

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

file1 = "./Res/Ser.txt"
print("Loading constants...")
dat = np.loadtxt(file1, max_rows=8)
Nl, Nt = map(int, dat[0:2])
dt, dx, rho, u, c, alp = dat[2:]
print("Loaded constants!!!")

print("Loading origin locations...")
x0 = np.loadtxt(file1, max_rows=Nt, skiprows=9, usecols=0 ,delimiter=',')
y0 = np.loadtxt(file1, max_rows=Nt, skiprows=9, usecols=1 ,delimiter=',')
print("Loaded origin locations!!!")

print("Loading forces...")
L = np.loadtxt(file1, max_rows=Nt, skiprows=(Nt+10), usecols=0 ,delimiter=',')
D = np.loadtxt(file1, max_rows=Nt, skiprows=(Nt+10), usecols=1 ,delimiter=',')
print("Loaded forces!!!")

xwM = np.zeros((Nt,Nt))
ywM = np.zeros((Nt,Nt))

print("Loading wake locations...")
fcnt = Nt
for i in range(fcnt):
    if i % 50 == 0:
        print("m = ",i,"/",fcnt)
    sprw = (2 + i)*Nt + 11 + i # (2*Nt + 11) + i*(Nt+1)
    xwM[i][:] = np.loadtxt(file1, max_rows=Nt, skiprows=sprw, usecols=0 ,delimiter=',')
    ywM[i][:] = np.loadtxt(file1, max_rows=Nt, skiprows=sprw, usecols=1 ,delimiter=',')
# print(xwM[Nt-1][-10:],ywM[Nt-1][-10:])
print("Loaded wake locations!!!")

# # fig1, ax1 = plt.subplots()
# # pdats = ax1.scatter([],[], s=3)
# fig1 = plt.figure(figsize=(8,6))
# pdats = plt.scatter([], [], s=4, c='r')
# plt.xlim(-110, 10*c)
# plt.ylim(-10, 10)
# plt.title("Serial Code", fontsize=16)
# # pdat = pdats[0]

# def animFunc(i):
#     pdats.set_offsets(np.c_[xwM[i][:i+1], ywM[i][:i+1]])
    
    
# # for i in range(fcnt):
# #     print(np.c_[xwM[i][:], ywM[i][:]])
# #     print("\n")

# anim = FuncAnimation(fig1, animFunc, frames=fcnt, interval=25)
# plt.show()

xte = x0[Nt-1] + c*np.cos(alp)
yte = y0[Nt-1] - c*np.sin(alp)

plt.figure(figsize=(8,6))
plt.scatter(x0[Nt-1],y0[Nt-1],c='k',s=8)
plt.plot([x0[Nt-1],xte],[y0[Nt-1],yte],c='k',lw=1)
plt.scatter(xwM[Nt-1][:],ywM[Nt-1][:],s=4, c='r')
plt.xlabel("x (in m)",size=14)
plt.ylabel("x (in m)",size=14)
plt.title("Wake structure (Serial)",size=16)
fig2 = plt.gcf()
plt.show()
fig2.savefig("./Res/Ser.png", dpi=300, bbox_inches='tight')