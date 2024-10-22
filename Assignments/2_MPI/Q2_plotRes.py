## CODE FOR PLOTTING THE RESULTS for Q2

import numpy as np
import matplotlib.pyplot as plt


# Q2 a) Serial

# file1 = "./Res/Q2_Ser_J_e1.txt"
# phivsx0 = np.loadtxt(file1, delimiter=",", usecols=0)
# phivsy0 = np.loadtxt(file1, delimiter=",", usecols=1)

# x = np.linspace(1, -1, num=phivsx0.size)

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(np.flip(x), phivsx0, 'r', lw=1.5, label=r'$\phi$ vs x', marker='o', ms=3, mec='k', mfc='k')
# plt.plot(x, phivsy0, 'k--', lw=1.5, label=r'$\phi$ vs y', marker='o', ms=3, mec='k', mfc='k')
# plt.xlabel("x", size=14)
# plt.ylabel(r'$\phi$', size=14)
# plt.title("Solutions at x = 0 and y = 0", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=3, bbox_to_anchor=(0.5,-0.1))
# fig1 = plt.gcf()
# plt.show()
# fig1.savefig("./Res/Q2_a.png", dpi=300, bbox_inches='tight')



# # Q2 b) Jacobi
# print("Parallel Jacobi\n")

# file1 = "./Res/Q2_Ser_J_e2.txt"
# phivsx0s = np.loadtxt(file1, delimiter=",", usecols=0)
# phivsy0s = np.loadtxt(file1, delimiter=",", usecols=1)

# file2 = "./Res/Q2_MPI_J2.txt"
# phivsx01 = np.loadtxt(file2, delimiter=",", usecols=0)
# phivsy01 = np.loadtxt(file2, delimiter=",", usecols=1)

# file3 = "./Res/Q2_MPI_J4.txt"
# phivsx02 = np.loadtxt(file3, delimiter=",", usecols=0)
# phivsy02 = np.loadtxt(file3, delimiter=",", usecols=1)

# file4 = "./Res/Q2_MPI_J8.txt"
# phivsx03 = np.loadtxt(file4, delimiter=",", usecols=0)
# phivsy03 = np.loadtxt(file4, delimiter=",", usecols=1)

# x = np.linspace(1, -1, num=phivsx01.size)

# errx1 = max(abs(phivsx01 - phivsx0s))
# errx2 = max(abs(phivsx02 - phivsx0s))
# errx3 = max(abs(phivsx03 - phivsx0s))

# erry1 = max(abs(phivsy01 - phivsy0s))
# erry2 = max(abs(phivsy02 - phivsy0s))
# erry3 = max(abs(phivsy03 - phivsy0s))

# print(errx1,errx2,errx3,sep="\t",end="\n\n")
# print(erry1,erry2,erry3,sep="\t",end="\n\n")

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(np.flip(x), phivsx0s, 'r', lw=1.5, label="Serial")
# plt.plot(np.flip(x), phivsx01, 'b--', lw=1.5, label="p = 2")
# plt.plot(np.flip(x), phivsx02, 'k--', lw=1.5, label="p = 4")
# plt.plot(np.flip(x), phivsx03, 'g--', lw=1.5, label="p = 8")
# plt.xlabel("x", size=14)
# plt.ylabel(r'$\phi$', size=14)
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
# plt.plot(x, phivsy0s, 'r', lw=1.5, label="Serial")
# plt.plot(x, phivsy01, 'b--', lw=1.5, label="p = 2")
# plt.plot(x, phivsy02, 'k--', lw=1.5, label="p = 4")
# plt.plot(x, phivsy03, 'g--', lw=1.5, label="p = 8")
# plt.xlabel("x", size=14)
# plt.ylabel(r'$\phi$', size=14)
# plt.title("Solutions at y = 0 for different p", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig3 = plt.gcf()
# plt.show()
# # fig3.savefig("./Res/Q2_b_2.png", dpi=300, bbox_inches='tight')



# # Q2 c1) Gauss-Seidel with J Serial
# print("Parallel GS with J Serial\n")

# file1 = "./Res/Q2_Ser_J_e2.txt"
# phivsx0s = np.loadtxt(file1, delimiter=",", usecols=0)
# phivsy0s = np.loadtxt(file1, delimiter=",", usecols=1)

# file2 = "./Res/Q2_MPI_GS2.txt"
# phivsx01 = np.loadtxt(file2, delimiter=",", usecols=0)
# phivsy01 = np.loadtxt(file2, delimiter=",", usecols=1)

# file3 = "./Res/Q2_MPI_GS4.txt"
# phivsx02 = np.loadtxt(file3, delimiter=",", usecols=0)
# phivsy02 = np.loadtxt(file3, delimiter=",", usecols=1)

# file4 = "./Res/Q2_MPI_GS8.txt"
# phivsx03 = np.loadtxt(file4, delimiter=",", usecols=0)
# phivsy03 = np.loadtxt(file4, delimiter=",", usecols=1)

# x = np.linspace(1, -1, num=phivsx01.size)

# errx1 = max(abs(phivsx01 - phivsx0s))
# errx2 = max(abs(phivsx02 - phivsx0s))
# errx3 = max(abs(phivsx03 - phivsx0s))

# erry1 = max(abs(phivsy01 - phivsy0s))
# erry2 = max(abs(phivsy02 - phivsy0s))
# erry3 = max(abs(phivsy03 - phivsy0s))

# print(errx1,errx2,errx3,sep="\t",end="\n\n")
# print(erry1,erry2,erry3,sep="\t",end="\n\n")

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(np.flip(x), phivsx0s, 'r', lw=1.5, label="Serial")
# plt.plot(np.flip(x), phivsx01, 'b--', lw=1.5, label="p = 2")
# plt.plot(np.flip(x), phivsx02, 'k--', lw=1.5, label="p = 4")
# plt.plot(np.flip(x), phivsx03, 'g--', lw=1.5, label="p = 8")
# plt.xlabel("x", size=14)
# plt.ylabel(r'$\phi$', size=14)
# plt.title("Solutions at x = 0 for different p", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig4 = plt.gcf()
# plt.show()
# # fig4.savefig("./Res/Q2_c_1.png", dpi=300, bbox_inches='tight')


# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(x, phivsy0s, 'r', lw=1.5, label="Serial")
# plt.plot(x, phivsy01, 'b--', lw=1.5, label="p = 2")
# plt.plot(x, phivsy02, 'k--', lw=1.5, label="p = 4")
# plt.plot(x, phivsy03, 'g--', lw=1.5, label="p = 8")
# plt.xlabel("x", size=14)
# plt.ylabel(r'$\phi$', size=14)
# plt.title("Solutions at y = 0 for different p", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig5 = plt.gcf()
# plt.show()
# # fig5.savefig("./Res/Q2_c_2.png", dpi=300, bbox_inches='tight')



# # Q2 c2) Gauss-Seidel with GS Serial
# print("Parallel GS with GS Serial\n")

# file1 = "./Res/Q2_Ser_GS_e2.txt"
# phivsx0s = np.loadtxt(file1, delimiter=",", usecols=0)
# phivsy0s = np.loadtxt(file1, delimiter=",", usecols=1)

# file2 = "./Res/Q2_MPI_GS2.txt"
# phivsx01 = np.loadtxt(file2, delimiter=",", usecols=0)
# phivsy01 = np.loadtxt(file2, delimiter=",", usecols=1)

# file3 = "./Res/Q2_MPI_GS4.txt"
# phivsx02 = np.loadtxt(file3, delimiter=",", usecols=0)
# phivsy02 = np.loadtxt(file3, delimiter=",", usecols=1)

# file4 = "./Res/Q2_MPI_GS8.txt"
# phivsx03 = np.loadtxt(file4, delimiter=",", usecols=0)
# phivsy03 = np.loadtxt(file4, delimiter=",", usecols=1)

# x = np.linspace(1, -1, num=phivsx01.size)

# errx1 = max(abs(phivsx01 - phivsx0s))
# errx2 = max(abs(phivsx02 - phivsx0s))
# errx3 = max(abs(phivsx03 - phivsx0s))

# erry1 = max(abs(phivsy01 - phivsy0s))
# erry2 = max(abs(phivsy02 - phivsy0s))
# erry3 = max(abs(phivsy03 - phivsy0s))

# print(errx1,errx2,errx3,sep="\t",end="\n\n")
# print(erry1,erry2,erry3,sep="\t",end="\n\n")

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(np.flip(x), phivsx0s, 'r', lw=1.5, label="Serial")
# plt.plot(np.flip(x), phivsx01, 'b--', lw=1.5, label="p = 2")
# plt.plot(np.flip(x), phivsx02, 'k--', lw=1.5, label="p = 4")
# plt.plot(np.flip(x), phivsx03, 'g--', lw=1.5, label="p = 8")
# plt.xlabel("x", size=14)
# plt.ylabel(r'$\phi$', size=14)
# plt.title("Solutions at x = 0 for different p", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig6 = plt.gcf()
# plt.show()
# # fig6.savefig("./Res/Q2_c_3.png", dpi=300, bbox_inches='tight')


# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(x, phivsy0s, 'r', lw=1.5, label="Serial")
# plt.plot(x, phivsy01, 'b--', lw=1.5, label="p = 2")
# plt.plot(x, phivsy02, 'k--', lw=1.5, label="p = 4")
# plt.plot(x, phivsy03, 'g--', lw=1.5, label="p = 8")
# plt.xlabel("x", size=14)
# plt.ylabel(r'$\phi$', size=14)
# plt.title("Solutions at y = 0 for different p", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig7 = plt.gcf()
# plt.show()
# # fig7.savefig("./Res/Q2_c_4.png", dpi=300, bbox_inches='tight')



# # Q2 d) Speed-up Calculation

# tJ  = np.array([736.3647, 635.5497, 193.0720, 107.2237, 62.8907])
# tGS = np.array([473.8178, 438.8389, 133.3876, 81.5027, 43.8397])
# p   = np.array([2, 4, 8, 16])

# psiJp  = tJ[0]/tJ[1:5]
# psiGSp = tGS[0]/tGS[1:5]

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(p, psiJp, 'r', lw=1.5, label='Jacobi', marker='o', ms=5, mec='k', mfc='k')
# plt.plot(p, psiGSp, 'b--', lw=1.5, label='Gauss-Seidel', marker='o', ms=5, mec='k', mfc='k')
# plt.xlabel(r'No. of processors $p$', size=14)
# plt.ylabel(r'Speed-up $\psi\,(n, p)$', size=14)
# plt.title('Speed-up as a function of processor count', size=16)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig8 = plt.gcf()
# plt.show()
# # fig8.savefig("./Res/Q2_d_1.png", dpi=300, bbox_inches='tight')

# delx = np.array([0.005, 0.01, 0.1])
# n    = 2/delx + 1
# print(n)

# tJ4  = np.array([193.0720, 16.3816, 0.0107])
# tJS  = np.array([736.3647, 60.7004, 0.0172])

# tGS4 = np.array([133.3876, 11.9735, 0.0101])
# tGSS = np.array([473.8178, 38.0178, 0.0163])

# psiJn  = tJS/tJ4
# psiGSn = tGSS/tGS4

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(n, psiJn, 'r', lw=1.5, label='Jacobi', marker='o', ms=5, mec='k', mfc='k')
# plt.plot(n, psiGSn, 'b--', lw=1.5, label='Gauss-Seidel', marker='o', ms=5, mec='k', mfc='k')
# plt.xlim(0, 420)
# plt.ylim(1.5, 4)
# plt.xlabel(r'No. of grid points $n$', size=14)
# plt.ylabel(r'Speed-up $\psi\,(n, p)$', size=14)
# plt.title('Speed-up as a function of problem size', size=16)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig9 = plt.gcf()
# plt.show()
# # fig9.savefig("./Res/Q2_d_2.png", dpi=300, bbox_inches='tight')