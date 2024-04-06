## CODE FOR PLOTTING THE RESULTS

import numpy as np
import matplotlib.pyplot as plt

# Q1 a) Serial

file1 = "./Res/Q1_Ser.txt"
act0 = np.loadtxt(file1, delimiter=",", usecols=1, max_rows=1001)
upw0 = np.loadtxt(file1, delimiter=",", usecols=2, max_rows=1001)
qck0 = np.loadtxt(file1, delimiter=",", usecols=3, max_rows=1001)

x = np.linspace(0, 2, num=act0.size)

plt.figure(figsize=(8,6))
ax = plt.subplot(111)
plt.plot(x, act0, 'r', lw=1.5, label="Actual")
plt.plot(x, upw0, 'k--', lw=0.75, label="Upwind")
plt.plot(x, qck0, 'b--', lw=1, label="QUICK")
plt.xlabel("x", size=14)
plt.ylabel("f'(x)", size=14)
plt.title("Wave at t = 0", size=14)
plt.grid()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
plt.legend(fontsize=14, loc="upper center", ncol=3, bbox_to_anchor=(0.5,-0.1))
fig1 = plt.gcf()
plt.show()
# fig1.savefig("./Res/LU.png", dpi=300, bbox_inches='tight')

act1 = np.loadtxt(file1, delimiter=",", usecols=1, max_rows=1001, skiprows=1002)
upw1 = np.loadtxt(file1, delimiter=",", usecols=2, max_rows=1001, skiprows=1002)
qck1 = np.loadtxt(file1, delimiter=",", usecols=3, max_rows=1001, skiprows=1002)

x = np.linspace(0, 2, num=act0.size)

plt.figure(figsize=(8,6))
ax = plt.subplot(111)
plt.plot(x, act1, 'r', lw=1.5, label="Actual")
plt.plot(x, upw1, 'k--', lw=0.75, label="Upwind")
plt.plot(x, qck1, 'b--', lw=1, label="QUICK")
plt.xlabel("x", size=14)
plt.ylabel("f'(x)", size=14)
plt.title("Wave at t = 0.5", size=14)
plt.grid()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
plt.legend(fontsize=14, loc="upper center", ncol=3, bbox_to_anchor=(0.5,-0.1))
fig2 = plt.gcf()
plt.show()
# fig2.savefig("./Res/LU.png", dpi=300, bbox_inches='tight')

act2 = np.loadtxt(file1, delimiter=",", usecols=1, max_rows=1001, skiprows=2004)
upw2 = np.loadtxt(file1, delimiter=",", usecols=2, max_rows=1001, skiprows=2004)
qck2 = np.loadtxt(file1, delimiter=",", usecols=3, max_rows=1001, skiprows=2004)

x = np.linspace(0, 2, num=act0.size)

plt.figure(figsize=(8,6))
ax = plt.subplot(111)
plt.plot(x, act2, 'r', lw=1.5, label="Actual")
plt.plot(x, upw2, 'k--', lw=0.75, label="Upwind")
plt.plot(x, qck2, 'b--', lw=1, label="QUICK")
plt.xlabel("x", size=14)
plt.ylabel("f'(x)", size=14)
plt.title("Wave at t = 1.0", size=14)
plt.grid()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
plt.legend(fontsize=14, loc="upper center", ncol=3, bbox_to_anchor=(0.5,-0.1))
fig3 = plt.gcf()
plt.show()
# fig3.savefig("./Res/LU.png", dpi=300, bbox_inches='tight')


# Q1 B) MPI

file1 = "./Res/Q1_MPI.txt"
act0 = np.loadtxt(file1, delimiter=",", usecols=1, max_rows=1001)
upw0 = np.loadtxt(file1, delimiter=",", usecols=2, max_rows=1001)
qck0 = np.loadtxt(file1, delimiter=",", usecols=3, max_rows=1001)

x = np.linspace(0, 2, num=act0.size)

plt.figure(figsize=(8,6))
ax = plt.subplot(111)
plt.plot(x, act0, 'r', lw=1.5, label="Actual")
plt.plot(x, upw0, 'k--', lw=0.75, label="Upwind")
plt.plot(x, qck0, 'b--', lw=1, label="QUICK")
plt.xlabel("x", size=14)
plt.ylabel("f'(x)", size=14)
plt.title("Wave at t = 0", size=14)
plt.grid()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
plt.legend(fontsize=14, loc="upper center", ncol=3, bbox_to_anchor=(0.5,-0.1))
fig4 = plt.gcf()
plt.show()
# fig4.savefig("./Res/LU.png", dpi=300, bbox_inches='tight')