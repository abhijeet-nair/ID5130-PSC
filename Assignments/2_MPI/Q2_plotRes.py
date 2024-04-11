## CODE FOR PLOTTING THE RESULTS for Q2

import numpy as np
import matplotlib.pyplot as plt


# Q2 a) Serial

# file1 = "./Res/Q2_Ser.txt"
# phivsx0s = np.loadtxt(file1, delimiter=",", usecols=0)
# phivsy0s = np.loadtxt(file1, delimiter=",", usecols=1)

# x = np.linspace(1, -1, num=phivsx0s.size)

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(np.flip(x), phivsx0s, 'r', lw=1.5, label="vs x", marker='o', ms=3, mec='k', mfc='k')
# plt.plot(x, phivsy0s, 'k--', lw=1.5, label="vs y", marker='o', ms=3, mec='k', mfc='k')
# plt.xlabel("x", size=14)
# plt.ylabel(r'$\phi$(x)', size=14)
# plt.title("Solutions at x = 0 and y = 0", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=3, bbox_to_anchor=(0.5,-0.1))
# fig1 = plt.gcf()
# plt.show()
# # fig1.savefig("./Res/Q2_Ser.png", dpi=300, bbox_inches='tight')


# Q2 b) Jacobi

# file1 = "./Res/Q2_MPI_J2.txt"
# phivsx01 = np.loadtxt(file1, delimiter=",", usecols=0)
# phivsy01 = np.loadtxt(file1, delimiter=",", usecols=1)

# file2 = "./Res/Q2_MPI_J4.txt"
# phivsx02 = np.loadtxt(file1, delimiter=",", usecols=0)
# phivsy02 = np.loadtxt(file1, delimiter=",", usecols=1)

# file3 = "./Res/Q2_MPI_J8.txt"
# phivsx03 = np.loadtxt(file1, delimiter=",", usecols=0)
# phivsy03 = np.loadtxt(file1, delimiter=",", usecols=1)

# x = np.linspace(1, -1, num=phivsx01.size)

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(np.flip(x), phivsx01, 'r', lw=1.5, label="p = 2")
# plt.plot(np.flip(x), phivsx02, 'b--', lw=1.5, label="p = 4")
# plt.plot(np.flip(x), phivsx03, 'k--', lw=1.5, label="p = 8")
# plt.plot(np.flip(x), phivsx0s, 'g--', lw=1.5, label="Serial")
# plt.xlabel("x", size=14)
# plt.ylabel(r'$\phi$(x)', size=14)
# plt.title("Solutions at x = 0 for different p", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig2 = plt.gcf()
# plt.show()
# # fig2.savefig("./Res/Q2_b_1.png", dpi=300, bbox_inches='tight')


# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(x, phivsy01, 'r', lw=1.5, label="p = 2")
# plt.plot(x, phivsy02, 'b--', lw=1.5, label="p = 4")
# plt.plot(x, phivsy03, 'k--', lw=1.5, label="p = 8")
# plt.plot(x, phivsy0s, 'g--', lw=1.5, label="Serial")
# plt.xlabel("x", size=14)
# plt.ylabel(r'$\phi$(x)', size=14)
# plt.title("Solutions at y = 0 for different p", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig3 = plt.gcf()
# plt.show()
# # fig3.savefig("./Res/Q2_b_2.png", dpi=300, bbox_inches='tight')



# # Q2 c) Gauss-Seidel

# file1 = "./Res/Q2_MPI_GS.txt"
# phivsx0 = np.loadtxt(file1, delimiter=",", usecols=0)
# phivsy0 = np.loadtxt(file1, delimiter=",", usecols=1)

# x = np.linspace(1, -1, num=phivsx0.size)

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(np.flip(x), phivsx0, 'r', lw=1.5, label="vs x")
# plt.plot(x, phivsy0, 'k--', lw=1.5, label="vs y")
# plt.xlabel("x", size=14)
# plt.ylabel(r'$\phi$(x)', size=14)
# plt.title("Solutions at x = 0 and y = 0", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=3, bbox_to_anchor=(0.5,-0.1))
# fig3 = plt.gcf()
# plt.show()
# # fig3.savefig("./Res/LU.png", dpi=300, bbox_inches='tight')