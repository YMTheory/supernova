import numpy as np
import matplotlib.pyplot as plt

dist = [1.0, 2.0, 5.0, 10.0]
effNO = [0.884,     0.8706875, 0.8071875, 0.6954375]
effIO = [ 0.996375, 0.9615,    0.8128125, 0.661]


fig, ax = plt.subplots()
ax.plot(dist, effNO, "^-", ms=6, lw=2, color="blue",    label="Normal Ordering")
ax.plot(dist, effIO, "v-", ms=6, lw=2, color="crimson", label="Inverted Ordering")

ax.hlines(0.89225, 0, 10, lw=2, linestyle="--", color="blue")
ax.hlines(0.99975, 0, 10, lw=2, linestyle="--", color="crimson")
ax.text(4, 0.98, "statistics = 10000", fontsize=13)

ax.legend(loc="lower left", prop={"size":13})
ax.set_xlabel("Distance [kpc]", fontsize=13)
ax.set_ylabel("Efficiency", fontsize=13)
ax.grid(True)

plt.tight_layout()
plt.savefig("KStest_sens2dist.pdf")
plt.show()



