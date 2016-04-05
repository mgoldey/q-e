# convert files to inputs

import ase
from ase import io
import numpy as np
import os
fl=os.popen("ls *.xyz").read().split()

top="""&CONTROL
calculation  = "scf"
outdir="/scratch1/scratchdirs/mgoldey/hab11p"
prefix="foo"
wf_collect=.true.
/
&SYSTEM
ibrav = 0,
nat  =  num_atoms ,
ntyp  =  num_types ,
nbnd  =  num_bands,
ecutwfc   = 80,
assume_isolated = 'mt',
tot_charge = 1,
nspin     = 2,
tot_magnetization = 1
nosym = .true.
/
&ELECTRONS
electron_maxstep = 300
conv_thr    = 1.D-8,
diago_full_acc = .true.
/
"""
for ifl in fl:
    mol=io.read(ifl)
    natoms=len(mol.numbers)
    ntypes=len(np.unique(mol.numbers))
    aids=np.unique(mol.get_chemical_symbols())
    itop=top
    itop=itop.replace("num_atoms",str(natoms))
    itop=itop.replace("num_types",str(ntypes))
    nbnd=(mol.numbers-(mol.numbers>10)*8-(mol.numbers>2)*2).sum()/2+10
    
    itop=itop.replace("num_bands",str(nbnd))
    mol.set_pbc(True)
    mol.center(vacuum=4)
    
    d1,d2=[1,natoms/2]
    a1,a2=[natoms/2+1,natoms]
    
    # left
    of=open(ifl.replace(".xyz","_left.in"),'w')
    of.write(itop.replace("foo",ifl.replace(".xyz","_left")))
    of.write("K_POINTS {Gamma}\nATOMIC_SPECIES\n")
    for iat in aids:
        of.write(iat+" 0.0 "+iat+".pbe-mt_fhi.UPF\n")
    of.write("\nCELL_PARAMETERS (angstrom)\n")
    for j in range(3): of.write("\t".join([str(i) for i in mol.cell[j]])+"\n")
    of.write("\nATOMIC_POSITIONS (angstrom)\n")
    for i in range(natoms):
        of.write(mol.get_chemical_symbols()[i]+"\t"+"\t".join([str(j) for j in mol.positions[i]])+"\n")
    of.write("\nEPCDFT\n1 5e-5 1e-2 40\ndelta_charge "+" ".join([str(i) for i in [d1,d2,a1,a2]])+" 1.0 0.15\n\n")
    of.close()
      
    #right
    of=open(ifl.replace(".xyz","_right.in"),'w')
    of.write(itop.replace("foo",ifl.replace(".xyz","_right")))
    of.write("K_POINTS {Gamma}\nATOMIC_SPECIES\n")
    for iat in aids:
        of.write(iat+" 0.0 "+iat+".pbe-mt_fhi.UPF\n")
    of.write("\nCELL_PARAMETERS (angstrom)\n")
    for j in range(3): of.write("\t".join([str(i) for i in mol.cell[j]])+"\n")
    of.write("\nATOMIC_POSITIONS (angstrom)\n")
    for i in range(natoms):
        of.write(mol.get_chemical_symbols()[i]+"\t"+"\t".join([str(j) for j in mol.positions[i]])+"\n")
    of.write("\nEPCDFT\n1 5e-5 1e-2 40\ndelta_charge "+" ".join([str(i) for i in [a1,a2,d1,d2]])+" 1.0 0.15\n\n")
    of.close()
