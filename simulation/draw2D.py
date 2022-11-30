import matplotlib.pyplot as plt
import numpy as np

arr = np.loadtxt("CEvNS.txt")
print(arr.shape)

fig, ax = plt.subplots(figsize=(6, 4))
im = ax.imshow(arr, extent=[-20, 40, 0, 0.1], aspect="auto", )
plt.colorbar(im, ax=ax)
ax.set_xlabel("post-bounce time [ms]", fontsize=14)
ax.set_ylabel(r"$E_\nu$ [MeV]", fontsize=14)
ax.tick_params(axis="both", labelsize=13)
plt.tight_layout()
plt.show()
