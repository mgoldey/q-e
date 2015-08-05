MODULE epcdft_mod
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  SAVE
  !
  LOGICAL :: debug ! print extra info
  LOGICAL :: s_spin ! calculate S matrix for each spin separately
  LOGICAL :: det_by_zgedi
  INTEGER :: iunwfc2 = 3636 ! unit for 2nd set of wfcs
  CHARACTER (len=256) :: outdir, outdir2, prefix2
  INTEGER :: occup1, occup2 ! occupied up states for system 1, 2
  INTEGER :: occdown1, occdown2 ! occupied down states for system 1, 2
  INTEGER :: fragment1_atom1, fragment1_atom2, fragment2_atom1, fragment2_atom2 ! atoms in acceptor fragments for sys 1 and 2
  REAL(DP) :: fragment1_amp, fragment2_amp ! amplitudes of weight functions (external pots)
  REAL(DP) :: freeen1, freeen2 ! free energies system 1 and 2
  COMPLEX(DP), ALLOCATABLE :: evc2(:,:) ! ks vecs for system 2
  COMPLEX(DP), ALLOCATABLE :: smat(:,:,:) ! overlap matrix  smat(  aa ab , ba bb, up down ) 
  COMPLEX(DP), ALLOCATABLE :: wmat(:,:,:) !  weigth matrix ( aa ab, ba bb, up down)
  REAL(DP), ALLOCATABLE :: w(:,:) ! weight functions w( r , system )
  REAL(DP) :: free1, free2 ! free energies for system 1 and 2 (no correction)
  COMPLEX(DP) :: hc(2,2,2) ! coupling matrix hc(a,b,spin)
  COMPLEX(DP) :: ohc(2,2,2) ! orthogonal coupling matrix ohc(a,b,spin)
  !
ENDMODULE epcdft_mod
