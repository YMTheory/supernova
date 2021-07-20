import numpy as np
import matplotlib.pyplot as plt
import argparse
from scipy import interpolate

#Parse command line
parser = argparse.ArgumentParser(description="Plot the neutrino energetics provided by Fornax.")
parser.add_argument("--time",help="time to plot spectra [0.5s]",type=float,default=0.500)
parser.add_argument("--filename",help="filename to use as input [spectra_share_12_0.dat]",type=str,default='spectra_share_12_0.dat')
args = parser.parse_args()

ngroup = 20 #number of neutrino energy bins per species
filename = args.filename #input file containing spectra

eg = np.loadtxt('energies.dat') #Loads energy bins (MeV) in first column
#and bin widths (MeV) in second column

time = np.loadtxt(filename,usecols=(0,)) #load simulation dump times
spec = np.zeros((len(time),ngroup+1))

#spectra
for i in range(ngroup+1):
  spec[:,i] = np.loadtxt(filename,usecols=(i,))

#returns spectra as a matrix with dimensions #timesteps by #ngroups+1,
#where the first element in each row is the time in seconds and subsequent elements
#the spectral energy in erg/s/Mev

#energy bins
eg_1 = eg[0][0:20] #electron neutrino species energies at bin center
eg_2 = eg[0][20:40] #anti-electron neutrino species energies at bin center
eg_3 = eg[0][40:60] #"mu" neutrino species energies at bin center, including mu, tau neutrinos and
                  #their antis

#energy width
de_1 = eg[1][0:20] #electron neutrino species energy widths
de_2 = eg[1][20:40] #anti-electron neutrino species widths
de_3 = eg[1][40:60] #"mu" neutrino species widths

#electron neutrino luminosity
lum = np.sum(spec[:,1:]*de_1,axis=1)

spec_idx = np.where(time<args.time)[0][-1] #find desired time index to plot spectra

#interpolate over spectra
eg_int = np.linspace(eg_1[0],eg_1[-1],1000)
tck0 = interpolate.splrep(eg_1, spec[spec_idx][1:], s=0)
spec_int = interpolate.splev(eg_int, tck0, der=0)
for i in range(len(spec_int)):
    spec_int[i] = np.max(spec_int[i],0)

#plot electron-neutrino spectra at desired time
plt.plot(eg_int,spec_int/100)
savefile_spec=args.filename[:19]
plt.xlabel('Time after bounce [s]')
plt.ylabel(r'$\delta L_{\nu_i}/d\epsilon$ [$10^{52}$ erg/s/MeV]')
plt.savefig(savefile_spec+'pdf')
plt.close()

#plot electron-neutrino luminosity
plt.plot(time,lum/100)
savefile_lum=args.filename[14:19]
plt.xlabel('Time after bounce [s]')
plt.ylabel(r'$L_{\nu}$ [$10^{52}$ erg s$^{-1}$]')
plt.savefig('lum'+savefile_lum+'pdf')
plt.close()
