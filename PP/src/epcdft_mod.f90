!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors to this file:
!   Nicholas P. Brawand (nicholasbrawand@gmail.com)
!   Matthew B. Goldey (matthew.goldey@gmail.com)
!
MODULE epcdft_mod
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL :: debug = .true.     ! print extra info
  LOGICAL :: debug2 = .false.   ! dump Smat, cofactor mat, W mat into files
  LOGICAL :: s_spin             ! calculate S matrix for each spin separately
  LOGICAL :: eig_of_w = .FALSE. ! not implemented yet. use eigenstates of W to orthog H (default is false so use lowdin)
  INTEGER :: iunwfc2 = 3636     ! unit for 2nd set of wfcs
  CHARACTER (len=256) :: outdir, outdir2, prefix2
  !
  ! wave vec vars
  !
  COMPLEX(DP), ALLOCATABLE :: evc1(:,:,:) ! ks vecs for system 1 evc1(npwx, nbnd, spin)
  COMPLEX(DP), ALLOCATABLE :: evc2(:,:,:) ! ks vecs for system 2
  COMPLEX(DP), ALLOCATABLE :: smat(:,:,:) ! det of overlap matrix  smat(  aa ab , ba bb, up down ) 
  INTEGER :: occup1, occup2               ! occupied up states for system 1, 2
  INTEGER :: occdown1, occdown2           ! occupied down states for system 1, 2
  !
  ! weight function vars
  !
  INTEGER :: fragment1_atom1, fragment1_atom2, fragment2_atom1, fragment2_atom2 ! atoms in acceptor fragments for sys 1 and 2
  REAL(DP), ALLOCATABLE :: w(:,:,:)         ! weight functions w( r , system, up down ) with amplitude
  COMPLEX(DP), ALLOCATABLE :: wmat(:,:,:) !  weigth matrix ( aa ab, ba bb, up down) without amplitude.
  REAL(DP), ALLOCATABLE :: lm(:,:) ! lagrange multipliers  ( constraint i,  system a or b)
  !
  ! energy vars
  !
  REAL(DP) :: free1, free2 ! free energies for system 1 and 2, 
                           ! Free = < H_KS + W(r) >, see text below Eq. 11b in J. Chem. Phys. 133, 244105 (2010)
                           !
  REAL(DP) :: cor1, cor2   ! energy corrections to the free energy for systems 1 and 2
  COMPLEX(DP) :: hc(2,2)   ! coupling matrix hc(a,b)
  COMPLEX(DP) :: ohc(2,2)  ! orthogonal coupling matrix ohc(a,b)
  !
ENDMODULE epcdft_mod
