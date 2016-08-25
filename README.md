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
* examples in test - run all using make
* NEW EPCDFT CARD for multiple constraints (all hirshfield for now)

* Example CARD INPUT:

* EPCDFT
* 1 1e-4 1e-2 20
* delta_charge 13 13 1 12 1.0 0.21

* The two lines under the EPCDFT card control the following inputs:
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


### Help! I can't achieve SCF convergence! ###

* Systems studied with CDFT are often difficult to converge. This is not a problem with CDFT but a matter of SCF convergence. Here is what you should do to achieve convergence.
* Turn off MT. We found that MT for the most part doesn't change the coupling but can really cause convergence issues even for simple systems.
* Turn on smearing start with a big value like 0.1 if you really don't have a good guess of the constraint potential.
* Run your CDFT calculation with loose convergence parameters (e.g. scf convergence threshold) and smearing and the program will being finding a better starting constraint potential.
* As your error in the constrained charge decreases, back off on the smearing and tighten your convergence parameters and continue restarting the CDFT calculation with new constraint potential.
* and iterate.
* In some cases, espescially for charged systems the individual sites have degenerate homos and you need to use small smearing to converge regardless if you are running CDFT.
* For difficult systems you can expect to run 600 or so SCF steps before converging. 
* Some other convergence tips: 
* try starting from a smaller cdft field
* Try fewer electrons - compare with the N-1 electron system?
* increase epcdft_thr

### Other tips ###
* check spin density difference  rho_up(r) - rho_down(r) using pp.x. This is extremely helpful for understanding the orbitals that are participating in charge transfer.