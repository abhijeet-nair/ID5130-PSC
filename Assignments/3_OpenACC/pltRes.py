## CODE FOR PLOTTING THE RESULTS

import numpy as np
import matplotlib.pyplot as plt


# Q1 a) Serial LU Plot

file1 = "./SerLU.txt"
num = np.loadtxt(file1, delimiter=",", usecols=0)
act = np.loadtxt(file1, delimiter=",", usecols=1)

print(np.linalg.norm(num - act))

x = np.linspace(0, 3, num=num.size)

plt.figure(figsize=(8,6))
ax = plt.subplot(111)
plt.plot(x, num, 'r', lw=1.5, label="Numerical")
plt.plot(x, act, 'k--', lw=0.75, label="Actual", marker='o', ms=5, mec='k', mfc='k')
plt.xlabel("x", size=14)
plt.ylabel("f'(x)", size=14)
plt.title("Derivative of f(x) - LU Method (Serial)", size=14)
plt.grid()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
fig1 = plt.gcf()
plt.show()
# fig1.savefig("./SerLU.png", dpi=300, bbox_inches='tight')


# Q1 b) Parallel LU Plot

file1 = "./ParLU.txt"
num = np.loadtxt(file1, delimiter=",", usecols=0)
act = np.loadtxt(file1, delimiter=",", usecols=1)

print(np.linalg.norm(num - act))

x = np.linspace(0, 3, num=num.size)

plt.figure(figsize=(8,6))
ax = plt.subplot(111)
plt.plot(x, num, 'r', lw=1.5, label="Numerical")
plt.plot(x, act, 'k--', lw=0.75, label="Actual", marker='o', ms=5, mec='k', mfc='k')
plt.xlabel("x", size=14)
plt.ylabel("f'(x)", size=14)
plt.title("Derivative of f(x) - LU Method (Parallel)", size=14)
plt.grid()
box = ax.get_position()
ax.set_position([box.x0, box.y0 + box.height * 0.1,
                 box.width, box.height * 0.9])
plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
fig1 = plt.gcf()
plt.show()
# fig1.savefig("./SerLU.png", dpi=300, bbox_inches='tight')