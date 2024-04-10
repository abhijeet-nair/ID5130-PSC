## CODE FOR PLOTTING THE RESULTS for Q2

import numpy as np
import matplotlib.pyplot as plt

# Q2 a) Serial

file1 = "./Res/Q2_Ser.txt"
phivsx0 = np.loadtxt(file1, delimiter=",", usecols=0)
phivsy0 = np.loadtxt(file1, delimiter=",", usecols=1)

x = np.linspace(1, -1, num=phivsx0.size)

plt.figure(figsize=(8,6))
ax = plt.subplot(111)
plt.plot(x, phivsx0, 'r', lw=1.5, label="vs x")
plt.plot(x, phivsy0, 'k--', lw=1.5, label="vs y")
plt.xlabel("x", size=14)
plt.ylabel(r'$\phi$(x)', size=14)
plt.title("Solutions at x = 0 and y = 0", size=14)
plt.grid()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
plt.legend(fontsize=14, loc="upper center", ncol=3, bbox_to_anchor=(0.5,-0.1))
fig1 = plt.gcf()
plt.show()
# fig1.savefig("./Res/LU.png", dpi=300, bbox_inches='tight')