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
* do_epcdft - needs to be set to .true.
* fragment_atom1 - atom # to designate start of fragment
* fragment_atom2 - atom # to designate end of fragment (set to zero if only one atom)
* epcdft_width - radius of potential well in bohr (for only one atom)
* epcdft_amp - strength of potential in Ry a.u.
* epcdft_electrons - number of electrons needed in well
* examples in test - run all using make

### Contribution guidelines ###

* Writing tests
* Code review
* Other guidelines

### Who do I talk to? ###

* Repo owner or admin
* Other community or team contact
