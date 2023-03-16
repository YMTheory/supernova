N_ham = 4997
N_hQE_NNVT = 9895 
N_lQE_NNVT = 2720
Y = 1400

DNR_ham = 8.94
DNR_hQE_NNVT = 29.05
DNR_lQE_NNVT = 29.06 #kHz


t = 300 #ns
N_tot = (N_ham * DNR_ham + N_hQE_NNVT * DNR_hQE_NNVT + N_lQE_NNVT * DNR_lQE_NNVT ) * 1e3 * t * 1e-9


print(f"Total DN number in trigger window : {N_tot:.1f} with visible energy shift {N_tot/Y:.1f} MeV")
