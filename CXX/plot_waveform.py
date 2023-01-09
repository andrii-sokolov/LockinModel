import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt("export.csv", delimiter=",").T

plt.plot(data[0], data[1])
ax = plt.twinx()
ax.plot(data[0], data[3], color="tab:orange")
plt.show()
