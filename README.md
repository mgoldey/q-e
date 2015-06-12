# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* EPCDFT - External Potential Constrained DFT
* Version 

### How do I get set up? ###

* Clone the repo
* cd espresso-5.1.1
* ./configure -enable-openmp=yes -enable-parallel=yes
* make pw

### INPUT Flags ###
* tefield - needs to be set to .true.
* do_epcdft - needs to be set to .true.
* fragment_atom1 - atom # to designate start of fragment
* fragment_atom2 - atom # to designate end of fragment (set to zero if only one atom)
* epcdft_width - radius of potential well in bohr (for only one atom)
* epcdft_amp - strength of potential in Ry a.u.
* epcdft_electrons - number of electrons needed in well
* epcdft_thr - threshold on number of electrons in well (default 1.D-4)
* examples in test - run all using make

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact


### What to do if it doesn't work? ###
* Try diagonalization='cg'
* Try fewer electrons - compare with the N-1 electron system?
* Reduce beta mixing parameter


