# README #

This README would normally document whatever steps are necessary to get your application up and running.

### What is this repository for? ###

* EPCDFT - External Potential Constrained DFT
* Version 

### How do I get set up? ###

* Clone the repo
* cd espresso-5.1.1
* ./configure -enable-openmp=yes -enable-parallel=yes
* make pw pp

### INPUT Flags ###
* assume_isolated = 'mt' - for isolated systems
* do_epcdft - needs to be set to .true.
* examples in test - run all using make
* NEW EPCDFT CARD for multiple constraints (all hirshfield for now)

* Example CARD INPUT:

* EPCDFT
* 1 1e-4 1e-2 20
* delta_charge 13 13 1 12 1.0 0.21
*
* Number of constraints, tolerance, delta_fld, update the potential every this many steps and at self consistency (default=40)
* type of constraint, acceptor start, acceptor end, donor start, donor end, number of electrons, initial lagrange multiplier guess

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
* add unoccupied states
* check spin density difference  rho_up(r) - rho_down(r) using pp.x
* try starting from a smaller cdft field
* Try fewer electrons - compare with the N-1 electron system?
* increase epcdft_thr
* use smearing till you converge W then restart the calculation with the same W and rho and back off on the smearing