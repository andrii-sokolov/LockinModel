import matplotlib.pyplot as plt
import numpy as np

data = np.genfromtxt("export.csv", delimiter=",").T

fig, axs = plt.subplots(nrows=2, figsize=(10, 10))

axs[0].plot(data[0] * 1e6, data[1], ".", label="$\\leftarrow$ X")
axs[0].plot(data[0] * 1e6, data[3], ".", color="tab:red", label="$\\leftarrow$ R")
axs[0].set_xlabel("Pulse width [us]")
axs[0].set_ylabel("X,R [$V^2$]")
axs[0].legend()

ax = axs[0].twinx()
ax.plot(data[0] * 1e6, data[2], ".", color="tab:orange", label="$\\rightarrow$ Y")
ax.set_ylabel("Y [$V^2$]")
ax.legend()

axs[1].plot(data[0] * 1e6, data[4], "k.")
axs[1].set_xlabel("Pulse width [us]")
axs[1].set_ylabel("$\\Theta$ [rad]")

plt.show()
