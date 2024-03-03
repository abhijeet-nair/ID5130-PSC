import numpy as np
import matplotlib.pyplot as plt


# # Q2 a) LU Plot
# file1 = "./Res/LU.txt"
# num = np.loadtxt(file1, delimiter=",", usecols=0)
# act = np.loadtxt(file1, delimiter=",", usecols=1)

# x = np.linspace(0, 3, num=num.size)

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(x, num, 'r', lw=1.5, label="Numerical")
# plt.plot(x, act, 'k--', lw=0.75, label="Actual", marker='o', ms=5, mec='k', mfc='k')
# plt.xlabel("x", size=14)
# plt.ylabel("f'(x)", size=14)
# plt.title("Derivative of f(x) - LU Method", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig1 = plt.gcf()
# plt.show()
# fig1.savefig("./Res/LU.png", dpi=300)



# # Q2 b) RDA Plot
# file2 = "./Res/RDA.txt"
# num = np.loadtxt(file2, delimiter=",", usecols=0)
# act = np.loadtxt(file2, delimiter=",", usecols=1)

# x = np.linspace(0, 3, num=num.size)

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(x, num, 'r', lw=1.5, label="Numerical")
# plt.plot(x, act, 'k--', lw=1.5, label="Actual")
# plt.xlabel("x", size=14)
# plt.ylabel("f'(x)", size=14)
# plt.title("Derivative of f(x) - RDA", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig2 = plt.gcf()
# plt.show()
# fig2.savefig("./Res/RDA.png", dpi=300)


# # Q2 b) RDA Time 
# p = np.array([2,4,8])
# t = np.array([0.000221, 0.000132, 0.000113])

# plt.figure(figsize=(8,6))
# plt.plot(p, t*1e6, 'r', lw=1.5, marker='o', ms=5, mec='k', mfc='k')
# plt.xlabel("No. of threads", size=14)
# plt.ylabel(r'Time taken (in $\mu$s)', size=14)
# plt.title("Time taken for RDA vs no. of threads", size=14)
# plt.grid()
# fig3 = plt.gcf()
# plt.show()
# # fig3.savefig("./Res/RDA_Time.png", dpi=300)


# # Q3 a) GS Serial
# file3 = "./Res/GS_Ser.txt"
# num = np.loadtxt(file3, delimiter=",", usecols=0)
# act = np.loadtxt(file3, delimiter=",", usecols=1)

# delta = 0.1
# x = -1 + np.arange(0, 21)*delta

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(x, num, 'r', lw=1.5, label="Numerical")
# plt.plot(x, act, 'k--', lw=0.75, label="Actual", marker='o', ms=5, mec='k', mfc='k')
# plt.xlabel("x", size=14)
# plt.ylabel(r"$\phi$(x)", size=14)
# plt.title("Solution to Poisson equation for y = 0 - GS Serial", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig3 = plt.gcf()
# plt.show()
# fig3.savefig("./Res/GS_Ser.png", dpi=300)


# # Q3 c) GS Serial, Diagonal, Red-Black
# file3 = "./Res/GS_Ser.txt"
# ser = np.loadtxt(file3, delimiter=",", usecols=0)

# file4 = "./Res/GS_Dia.txt"
# dia = np.loadtxt(file4, delimiter=",", usecols=0)

# file5 = "./Res/GS_RBC.txt"
# rbc = np.loadtxt(file5, delimiter=",", usecols=0)

# delta = 0.1
# x = -1 + np.arange(0, 21)*delta

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(x, dia, 'r', lw=1.5, label="Diagonal")
# plt.plot(x, ser, 'k--', lw=0.75, label="Serial", marker='o', ms=5, mec='k', mfc='k')
# plt.xlabel("x", size=14)
# plt.ylabel(r"$\phi$(x)", size=14)
# plt.title("Solution for y = 0 - GS Serial vs Diagonal", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig4 = plt.gcf()
# plt.show()
# fig4.savefig("./Res/GS_Dia.png", dpi=300)

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(x, rbc, 'r', lw=1.5, label="Red-Black Coloring")
# plt.plot(x, ser, 'k--', lw=0.75, label="Serial", marker='o', ms=5, mec='k', mfc='k')
# plt.xlabel("x", size=14)
# plt.ylabel(r"$\phi$(x)", size=14)
# plt.title("Solution for y = 0 - GS Serial vs Red-Black Coloring", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig5 = plt.gcf()
# plt.show()
# fig5.savefig("./Res/GS_RBC.png", dpi=300)

# # Q3 c) Time comparison
# delta = np.array([0.1, 0.01, 0.005])
# tSer = np.array([0.008798, 53.900578, 934.409124])

# # tDia = np.array([0.012313, 42.988759, 805.271786]) # 4
# # tDia = np.array([0.017157, 61.307862, 1078.219453]) # 8
# tDia2 = np.array([0.009243, 29.999306, 523.061110]) # 2
# tDia4 = np.array([0.009048, 17.905741, 352.446755]) # 4
# tDia8 = np.array([0.009218, 15.954448, 336.355502]) # 8
# tDiaL = np.array([523.061110, 352.446755, 336.355502]) # 2, 4, 8

# # tRBC = np.array([0.007855, 39.185373, 762.295155]) # 4
# # tRBC = np.array([0.008133, 53.647345, 1006.334569]) # 8
# tRBC2 = np.array([0.006428, 27.041639, 487.444447]) # 2
# tRBC4 = np.array([0.005856, 13.946496, 310.227622]) # 4
# tRBC8 = np.array([0.005207, 11.504038, 226.176504]) # 8
# tRBCL = np.array([487.444447, 310.227622, 226.176504]) # 2, 4, 8

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.loglog(delta, tSer, 'r', lw=1.5, label="Serial", marker='o', ms=5)
# plt.loglog(delta, tDia8, 'b', lw=1.5, label="Diagonal", marker='o', ms=5)
# plt.loglog(delta, tRBC8, 'k', lw=1.5, label="Red-Black Color", marker='o', ms=5)
# plt.xlabel(r"$\Delta$", size=14)
# plt.ylabel("Time taken (in s)", size=14)
# plt.title("Time taken for different methods with 8 threads", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig6 = plt.gcf()
# plt.show()
# # fig6.savefig("./Res/Q3_c_Time_log.png", dpi=300)

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot(delta, tSer, 'r', lw=1.5, label="Serial", marker='o', ms=5)
# plt.plot(delta, tDia8, 'b', lw=1.5, label="Diagonal", marker='o', ms=5)
# plt.plot(delta, tRBC8, 'k', lw=1.5, label="Red-Black Color", marker='o', ms=5)
# plt.xlabel(r"$\Delta$", size=14)
# plt.ylabel("Time taken (in s)", size=14)
# plt.title("Time taken for different methods with 8 threads", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig7 = plt.gcf()
# plt.show()
# # fig7.savefig("./Res/Q3_c_Time_plot.png", dpi=300)

# plt.figure(figsize=(8,6))
# ax = plt.subplot(111)
# plt.plot([2, 4, 8], tDiaL, 'b', lw=1.5, label="Diagonal", marker='o', ms=5)
# plt.plot([2, 4, 8], tRBCL, 'k', lw=1.5, label="Red-Black Color", marker='o', ms=5)
# plt.xlabel("No. of threads", size=14)
# plt.ylabel("Time taken (in s)", size=14)
# plt.title("Time taken for the two methods with different no. of threads", size=14)
# plt.grid()
# box = ax.get_position()
# ax.set_position([box.x0, box.y0 + box.height * 0.1,
#                  box.width, box.height * 0.9])
# plt.legend(fontsize=14, loc="upper center", ncol=2, bbox_to_anchor=(0.5,-0.1))
# fig8 = plt.gcf()
# plt.show()
# # fig8.savefig("./Res/Q3_d_Time.png", dpi=300)