# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* EPCDFT - External Potential Constrained DFT
* Version 

### How do I get set up? ###

* Clone the repo
* download QE 5.1.1
* tar -zxvf the QE  inside the repo you just downloaded
* run "git checkout ." inside the repo and this will restore the EPCDFT files
* install QE

### INPUT Flags ###
* tefield - needs to be set to .true.
* edir - atom # to center potential well around 
* emaxpos - radius of potential well in alat
* eamp - strength of potential in Ry a.u.
* eopreg - number of electrons needed in well
* example below

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact

### Example (included in the repo) ###
* &CONTROL
*  calculation  = "scf",
*  prefix       = "si",
*  pseudo_dir   = "./",
*  outdir       = "./out",
* !  wf_collect   = .TRUE.
*  tefield = .true.
* /
* &SYSTEM
*  ibrav     = 1,
*  a = 7,
*  nat       = 5,
*  ntyp      = 2,
*  ecutwfc   = 18,
*  !ecutrho   = 300,
*  ! nspin     = 2,
*  ! tot_magnetization = 1
*  edir = 1,
*  emaxpos = 0.2,
*  eamp = -1,
*  eopreg = 6
* /
* &ELECTRONS
*  conv_thr    = 1.D-6,
*  mixing_beta = 0.30,
* /
* &IONS
* /
* K_POINTS {Gamma}
* ATOMIC_SPECIES
* Si 28.0855 Si.pbe-n-van.UPF
* H 1.00794 H.pbe-rrkjus.UPF
* ATOMIC_POSITIONS (angstrom)
* Si        -0.19448        2.43790        0.00000
* H          1.27605        2.43790       -0.00000
* H         -0.68466        3.14888        1.19025
* H         -0.68466        1.05163        0.02060
* H         -0.68466        3.11319       -1.21085