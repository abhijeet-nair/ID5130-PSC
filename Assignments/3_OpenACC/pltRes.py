## CODE FOR PLOTTING THE RESULTS

import numpy as np
import matplotlib.pyplot as plt


# # Q1 a) Serial LU Plot

# file1 = "./SerLU.txt"
# num = np.loadtxt(file1, delimiter=",", usecols=0)
# act = np.loadtxt(file1, delimiter=",", usecols=1)

# print(np.linalg.norm(num - act))

# x = np.linspace(0, 3, num=num.size)

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(x, num, 'r', lw=1.5, label="Numerical")
# plt.plot(x, act, 'k--', lw=0.75, label="Actual", marker='o', ms=5, mec='k', mfc='k')
# plt.xlabel("x", size=14)
# plt.ylabel("f'(x)", size=14)
# plt.title("Derivative of f(x) - LU Method (Serial)", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig1 = plt.gcf()
# plt.show()
# # fig1.savefig("./SerLU.png", dpi=300, bbox_inches='tight')


# # Q1 b) Parallel LU Plot

# file1 = "./ParLU.txt"
# num = np.loadtxt(file1, delimiter=",", usecols=0)
# act = np.loadtxt(file1, delimiter=",", usecols=1)

# print(np.linalg.norm(num - act))

# x = np.linspace(0, 3, num=num.size)

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(x, num, 'r', lw=1.5, label="Numerical")
# plt.plot(x, act, 'k--', lw=0.75, label="Actual", marker='o', ms=5, mec='k', mfc='k')
# plt.xlabel("x", size=14)
# plt.ylabel("f'(x)", size=14)
# plt.title("Derivative of f(x) - LU Method (Parallel)", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig2 = plt.gcf()
# plt.show()
# # fig2.savefig("./ParLU.png", dpi=300, bbox_inches='tight')


# # Q1 c) Time info
# tSer = 0.002691
# nG = np.array([10,100,1000])
# tCLK = np.array([2.464049, 2.386697, 2.414653])
# tPGI = np.array([3089, 2692, 2746])*1e-6

# plt.figure(figsize=(8,6))
# plt.semilogx(nG, tPGI*1e6, 'r', lw=1, marker='o', ms=5, mec='k', mfc='k')
# plt.xlabel("Number of gangs", size=14)
# plt.ylabel("Time taken (in us)", size=14)
# plt.title("Time taken vs. Number of gangs", size=16)
# plt.grid()
# plt.ylim(2500,3500)
# fig3 = plt.gcf()
# plt.show()
# # fig3.savefig("./Q1T.png", dpi=300, bbox_inches='tight')


# # Q2 a) Time info
# N = np.array([10,100,1000])
# tSer = np.array([0.000243, 0.02583, 2.627272])
# tPGI = np.array([0.000027, 0.003533, 0.799138])

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.semilogx(N, tSer, 'r', lw=1, label="Serial", marker='o', ms=5, mec='k', mfc='k')
# plt.semilogx(N, tPGI, 'b--', lw=1, label="Parallel", marker='o', ms=5, mec='k', mfc='k')
# plt.xlabel("Matrix Size", size=14)
# plt.ylabel("Time taken (in s)", size=14)
# plt.title("Time taken vs. Matrix Size", size=16)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.08,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig4 = plt.gcf()
# plt.show()
# # fig4.savefig("./Q2T.png", dpi=300, bbox_inches='tight')