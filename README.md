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
* do_epcdft - needs to be set to .true.
* hirshfeld - if set to .true. then hirshfeld partitioning is used rather than voronoi cell (default .false.)
* fragment_atom1 - atom # to designate start of fragment
* fragment_atom2 - atom # to designate end of fragment (set to zero if only one atom)
* epcdft_width - radius of potential well in bohr (for only one atom)
* epcdft_amp - strength of potential in Ry a.u.
* epcdft_electrons - charge difference between inside the well/acceptor and outside. Example: epcdft_electrons = 1 will force one additional electron in well/on acceptor
* epcdft_thr - threshold on number of electrons in well (default 1.D-4)
* examples in test - run all using make

### PP Flags ###
* Using plot_num = 12 will plot the weight times the lagrange multiplier (read from successful run)

### Contribution guidelines ###

* Writing tests - as you go, make small working examples
* Code review
* Other guidelines

### TODO ###

* input for largest step in potential

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact


### What to do if it doesn't work? ###
* Try diagonalization='cg'
* Try fewer electrons - compare with the N-1 electron system?
* Reduce beta mixing parameter