import numpy as np
import matplotlib.pyplot as plt

N = np.array([256, 512, 1024, 2048])
t = np.array([0.0003842, 0.0012136, 0.0040756 ,0.0155444])

plt.figure(figsize=(10,8))
plt.plot(N,t*1e3, 'r--', linewidth=2, marker='o', markersize=10, markerfacecolor='k', markeredgecolor='k')
plt.xlabel('N', fontsize=18)
plt.ylabel('Time (in ms)', fontsize=18)
plt.title('Time vs. N', fontsize=20)
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.grid()
plt.show()