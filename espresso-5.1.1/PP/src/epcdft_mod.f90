MODULE epcdft_mod
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL :: debug = .true.     ! print extra info
  LOGICAL :: s_spin             ! calculate S matrix for each spin separately
  LOGICAL :: eig_of_w = .FALSE. ! use eigenstates of W to orthog H (default is false so use lowdin)
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
  REAL(DP), ALLOCATABLE :: w(:,:)         ! weight functions w( r , system ) with amplitude
  REAL(DP) :: wamp1, wamp2                ! amplitude of external pot applied in cdft runs 1 and 2
  COMPLEX(DP), ALLOCATABLE :: wmat(:,:,:) !  weigth matrix ( aa ab, ba bb, up down) without amplitude.
  !
  ! energy vars
  !
  REAL(DP) :: free1, free2 ! free energies for system 1 and 2 (no correction)
  REAL(DP) :: cor1, cor2   ! energy corrections to the free energy for systems 1 and 2
  COMPLEX(DP) :: hc(2,2)   ! coupling matrix hc(a,b)
  COMPLEX(DP) :: ohc(2,2)  ! orthogonal coupling matrix ohc(a,b)
  !
ENDMODULE epcdft_mod
