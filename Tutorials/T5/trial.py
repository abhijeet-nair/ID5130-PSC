import numpy as np
import matplotlib.pyplot as plt

x = np.array([1,2,3,4,5,6])
y = np.array([1,4,9,16,25,36])

plt.figure(figsize=(8,6))
plt.plot(x,y)
plt.xlabel("x")
plt.ylabel("y")
plt.show()