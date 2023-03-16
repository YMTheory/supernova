import numpy as np
import matplotlib.pyplot as plt
import numpy as np

chi2NO, chi2IO = [], []
dTNO, dTIO = [], []
masses = np.arange(0.0 , 2.1, 0.1)
masses1 = [0.0, 0.1, 0.2, 0.3, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
for MO in ["NO", "IO"]:
    for nuMass in masses:
        filename = f"../results/Garching82703_10kpc_{MO}_eESonly_asimovFit2D_nuMass{nuMass:.1f}eV.txt"
        arr = np.loadtxt(filename)
        
        if MO == "NO":
            dTNO.append(arr[0])
            chi2NO.append(2*arr[1])
        else:
            dTIO.append(arr[0])
            chi2IO.append(2*arr[1])

fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(masses, chi2NO - np.min(chi2NO), "-", color="blue", lw=2, label="NO")
ax.plot(masses, chi2IO - np.min(chi2IO), ":", color="red" , lw=2, label="IO")
#ax.plot(masses, dTNO, "-", color="blue", lw=2, label="NO")
#ax.plot(masses, dTIO, ":", color="red" , lw=2, label="IO")

ax.set_xlabel(r"$m_\nu$ [eV]", fontsize=18)
ax.set_ylabel(r"$\Delta\chi^2_\mathrm{m}$", fontsize=18)
ax.tick_params(axis="x", labelsize=16)
ax.tick_params(axis="y", labelsize=16)
ax.grid(True, linestyle=":")
ax.legend(prop={"size":11}, loc="upper left", frameon=True, ncol=1)

plt.tight_layout()
plt.savefig("../plots/Garching82703_10kpc_eESonly2D_nuMassSens.pdf")
plt.show()





