import numpy as np
import matplotlib.pyplot as plt


# # Q2 a) LU Plot
# file1 = "./Res/LU.txt"
# num = np.loadtxt(file1, delimiter=",", usecols=0)
# act = np.loadtxt(file1, delimiter=",", usecols=1)

# x = np.linspace(0, 3, num=num.size)

# plt.figure(figsize=(8,6))
# plt.plot(x, num, 'r', lw=1.5, label="Numerical")
# plt.plot(x, act, 'b--', lw=1.5, label="Actual")
# plt.xlabel("x", size=14)
# plt.ylabel("f'(x)", size=14)
# plt.title("Derivative of f(x) - LU Method", size=14)
# plt.grid()
# fig1 = plt.gcf()
# plt.show()
# # fig1.savefig("./Res/LU.png", dpi=300)



# # Q2 b) RDA Plot
# file2 = "./Res/RDA.txt"
# num = np.loadtxt(file2, delimiter=",", usecols=0)
# act = np.loadtxt(file2, delimiter=",", usecols=1)

# x = np.linspace(0, 3, num=num.size)

# plt.figure(figsize=(8,6))
# plt.plot(x, num, 'r', lw=1.5, label="Numerical")
# plt.plot(x, act, 'b--', lw=1.5, label="Actual")
# plt.xlabel("x", size=14)
# plt.ylabel("f'(x)", size=14)
# plt.title("Derivative of f(x) - RDA", size=14)
# plt.grid()
# fig2 = plt.gcf()
# plt.show()
# # fig2.savefig("./Res/RDA.png", dpi=300)

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


# Q3 a) GS Serial
file3 = "./Res/GS_Ser.txt"
num = np.loadtxt(file3, delimiter=",", usecols=0)
act = np.loadtxt(file3, delimiter=",", usecols=1)

delta = 0.1
x = -1 + np.arange(0, 21)*delta

plt.figure(figsize=(8,6))
plt.plot(x, num, 'r', lw=1.5, label="Numerical")
plt.plot(x, act, 'b--', lw=1.5, label="Actual")
plt.xlabel("x", size=14)
plt.ylabel(r"$\phi$(x)", size=14)
plt.title("Solution to Poisson equation for y = 0 - GS Serial", size=14)
plt.grid()
fig3 = plt.gcf()
plt.show()
# fig3.savefig("./Res/GS_Ser.png", dpi=300)