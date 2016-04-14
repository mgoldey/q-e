def get_coupling(fn):
	import os
	coup=float(os.popen("grep eV "+fn).read().split()[-2])*1e3
	return coup

def plot_rate_vs_T(coup,l,dG,label="Marcus rate"):
	from scipy import constants as con
	import numpy as np
	import matplotlib.pyplot as plt
	coup/=1e3
	T=np.arange(1e-6,501,5e-2)
	rate=(2*np.pi*con.eV/con.hbar)*coup**2*(4*np.pi*l*con.k*T/con.eV)**-.5*np.exp(-(l+dG)**2/(4*l*con.k*T/con.eV))
	plt.plot(T,rate,label=label)
	plt.legend()
	plt.xlabel("Temperature [K]")
	plt.ylabel("Rate [1/s]")

def plot_files(fl,l,dG):
	import matplotlib.pyplot as plt
	fig,ax1=plt.subplots(1,1,sharex=True,figsize=(18,10))
	for ifl in fl:
		plot_rate_vs_T(get_coupling(ifl),l,dG,ifl.replace("_coupling.out",""))
	ax1.set_yscale("symlog",linthreshy=1e-4)
	plt.show()
