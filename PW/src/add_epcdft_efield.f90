!
! Copyright (C) 2003-2017 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors to this file:
!   Nicholas P. Brawand (nicholasbrawand@gmail.com)
!   Matthew B. Goldey (matthew.goldey@gmail.com)
!
!--------------------------------------------------------------------------
SUBROUTINE add_epcdft_efield(vpoten,iflag)
  !--------------------------------------------------------------------------
  !
  !   This routine adds the constraining potential for cdft. 
  !
  USE kinds,            ONLY : DP
  USE epcdft,           ONLY : do_epcdft
  USE io_global,        ONLY : stdout,ionode
  USE lsda_mod,         ONLY : nspin
  USE fft_base,         ONLY : dfftp 
  USE uspp,             ONLY : okvan
  USE paw_variables,    ONLY : okpaw
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP),INTENT(INOUT) :: vpoten(dfftp%nnr,nspin) ! ef is added to this potential
  LOGICAL,INTENT(IN)     :: iflag ! set to true to force recalculation of field
  !
  ! local variables
  !
  REAL(DP) :: vconstr (dfftp%nnr,nspin) ! constraining potential shape*weight
  REAL(DP) :: x0(3) ! center of charge of system
  REAL(DP) :: qq ! total charge
  REAL(DP) :: dipole(3)!, quadrupole(3) ! total dips
  !
  LOGICAL :: first=.TRUE.
  LOGICAL :: hirshfeld=.TRUE.
  !
  INTEGER :: ier=1
  !
  SAVE first
  !
  !---------------------
  !  Execution control
  !---------------------
  !
  IF (.NOT. do_epcdft) RETURN
  !
  IF ((.NOT.first) .AND..NOT. iflag) RETURN
  !
  ! Check CDFT requirements
  !
  IF(first)THEN
    !
    IF(.NOT. gamma_only) CALL errore('add_epcdft_efield', 'CDFT requires gamma_only.', ier)  
    IF(okvan) CALL errore('add_epcdft_efield', 'CDFT: ultrasoft not implemented.', ier)  
    IF(okpaw) CALL errore('add_epcdft_efield', 'CDFT: PAW not implemented.', ier)  
    IF(noncolin) CALL errore('add_epcdft_efield', 'CDFT noncolin not implemented.', ier)  
    !
  ENDIF
  !
  first=.FALSE.
  !
  ! Necessary for restart/pp runs
  !
  CALL init_at_1
  !
  ! calculate hirshfeld potential
  !
  IF (hirshfeld) THEN 
    !
    vconstr = 0.D0
    CALL calc_hirshfeld_v(vconstr,0)
    vpoten= vpoten + vconstr 
    !
  ENDIF
  !
END SUBROUTINE add_epcdft_efield
!
!--------------------------------------------------------------------------
SUBROUTINE calc_hirshfeld_v( v,iconstraint)
  !--------------------------------------------------------------------------
  ! 
  !  Gamma point only for now ! If you extend this to k-points, make sure 
  !  to clearly define how many excess charges and what sizes of supercells
  !  are needed
  !
  !  Calculates hirshfeld potential and put into "v" following :
  !
  !     J. Chem. Phys. 133, 244105 (2010); http://dx.doi.org/10.1063/1.3507878 
  !     eqs. 6 & 7
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE wvfct,                ONLY : npwx, g2kin
  USE klist,                ONLY : xk, nks, igk_k,ngk,init_igk
  USE gvecw,                ONLY : gcutw
  USE gvect,                ONLY : ngm, g
  USE cell_base,            ONLY : tpiba2
  USE uspp_param,           ONLY : upf
  USE lsda_mod,      ONLY : nspin
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE fft_base,             ONLY : dfftp, dffts
  USE gvect,                ONLY : nl, nlm
  USE gvecs,                ONLY : nls, nlsm
  USE fft_interfaces,       ONLY : invfft
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_bcast, mp_sum
  USE io_global,     ONLY : stdout,ionode, ionode_id
  USE klist,         ONLY : nelec
  USE cell_base,     ONLY : omega
  USE wvfct,                ONLY : nbnd, npwx
  USE epcdft
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: v(dfftp%nnr,nspin) ! the hirshfeld potential
  INTEGER, INTENT(IN) :: iconstraint
  ! 
  INTEGER :: na, nt, nb, l, m,n,ir, nwfcm ,nwfc,icon, ik ,npw
  COMPLEX(DP), ALLOCATABLE :: wfcatomg (:,:) ! atomic wfcs in g
  COMPLEX(DP), ALLOCATABLE :: wfcatomr (:) ! atomic wfcs in r
  COMPLEX(DP), ALLOCATABLE :: total_atom_rho_r (:) ! total density of an atom in r
  INTEGER :: orbi ! orbital index for wfcatom
  INTEGER :: nfuncs, lmax ! number of functions as derived from lmax
  REAL(DP) :: orboc ! occupation of orbital
  COMPLEX(DP) :: vtop(dfftp%nnr,nspin) ! top of the hirshfeld potential fraction
  COMPLEX(DP) :: vbot(dfftp%nnr) ! bottom of the hirshfeld potential fraction
  COMPLEX(DP) :: normfac
  REAL(DP) :: cutoff
  COMPLEX(DP) :: vbottot
  REAL(DP) :: dv
  CHARACTER(len=1024) :: filename
  !
  CALL init_igk ( npwx, ngm, g, gcutw ) 
  n=dfftp%nnr
  ik=1
  npw=ngk(ik)
  !
  ALLOCATE( wfcatomr(n) )
  ALLOCATE( total_atom_rho_r(n) )
  !
  wfcatomr = 0.D0
  vtop = 0.D0
  vbot = 0.D0
  v = 0.D0
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  cutoff = 1.D-6
  total_atom_rho_r = 0.D0
  !
  ! construct hirshfeld looping over all atomic states
  !
  DO na = 1, nat ! for each atom
    !
    nt = ityp (na) ! get atom type
    nwfc=sum(upf(nt)%oc(:))
    total_atom_rho_r=0.D0
    !
    ALLOCATE( wfcatomg(npwx, nwfc) ) ! was npwx
    wfcatomg = 0.D0
    !
    CALL one_atom_wfc (1, wfcatomg, na,nwfc) 
    !
    orbi = 0 
    DO nb = 1, upf(nt)%nwfc ! for each orbital
      !
      IF (upf(nt)%oc(nb) > 0.d0) THEN ! if occupied
        !
        l = upf(nt)%lchi(nb) ! get l 
        !
        DO m = 1, 2*l+1 ! mag num
          !
          ! add all orbitals, each state weighted by orboc
          orbi = orbi + 1 ! update orbital index
          !  
          orboc = REAL( upf(nt)%oc(nb) ,KIND=DP) / REAL( 2*l+1 ,KIND=DP)
          !
          ! bring to real space
          !
          wfcatomr = 0.D0
          wfcatomr(nl(igk_k(1:npw,1)))  = wfcatomg(1:npw,orbi)
          IF(gamma_only) wfcatomr(nlm(igk_k(1:npw,1))) = CONJG( wfcatomg(1:npw,orbi) )
          !
          ! convert atomic(g) -> |atomic(r)|^2
          !
          CALL invfft ('Dense', wfcatomr(:), dfftp)
          wfcatomr(:) = wfcatomr(:) * CONJG(wfcatomr(:))
          !
          ! add orbital density to total density of atom
          !
          total_atom_rho_r(:) = total_atom_rho_r(:) + orboc * wfcatomr(:)
          !
        ENDDO ! m
      ENDIF ! end if occupied
    ENDDO ! nb 
    !
    ! add this atoms density to the hirshfeld potential, as appropriate
    !
    IF (iconstraint .eq. 0) THEN 
      DO icon=1,nconstr_epcdft
        !
        SELECT CASE( epcdft_type(icon) )
        CASE('charge')
          !
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            !
          ENDIF
          !
        CASE('spin')
          !
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + epcdft_guess(icon)*total_atom_rho_r( : ) 
            !
          ENDIF
          !
        CASE('delta_charge')
          !
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            !
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,1) = vtop(:,1) + epcdft_guess(icon)*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + epcdft_guess(icon)*total_atom_rho_r( : ) 
            !
          ENDIF   
          !
        CASE('delta_spin')
          !
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + epcdft_guess(icon)*total_atom_rho_r( : ) 
            !
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,1) = vtop(:,1) + epcdft_guess(icon)*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            !
          ENDIF
          !
        CASE('delta_alpha')
          !
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            !
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,1) = vtop(:,1) + epcdft_guess(icon)*total_atom_rho_r( : ) 
            !
          ENDIF
          !
        CASE('delta_beta')
          !
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,2) = vtop(:,2) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            !
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,2) = vtop(:,2) + epcdft_guess(icon)*total_atom_rho_r( : ) 
            !
          ENDIF
          !
        END SELECT
        !
      END DO
      !
    ELSE 
      !
      icon=iconstraint
      SELECT CASE( epcdft_type(icon) )
        CASE('charge')
          IF( na >= epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - total_atom_rho_r( : ) 
            !
          ENDIF
          !
        CASE('spin')
          !
          IF( na >= epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + total_atom_rho_r( : ) 
            !
          ENDIF
          !
        CASE('delta_charge')
          !
          IF( na >= epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - total_atom_rho_r( : ) 
            !
          ELSE IF( na >= epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,1) = vtop(:,1) + total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + total_atom_rho_r( : ) 
            !
          ENDIF   
          !
        CASE('delta_spin')
          !
          IF( na >= epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + total_atom_rho_r( : ) 
            !
          ELSE IF( na >= epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,1) = vtop(:,1) + total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - total_atom_rho_r( : ) 
            !
          ENDIF
          !
        CASE('delta_alpha')
          !
          IF( na >= epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - total_atom_rho_r( : ) 
            !
          ELSE IF( na >= epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,1) = vtop(:,1) + total_atom_rho_r( : ) 
            !
          ENDIF
          !
        CASE('delta_beta')
          !
          IF( na >= epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,2) = vtop(:,2) - total_atom_rho_r( : ) 
            !
          ELSE IF( na >= epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,2) = vtop(:,2) + total_atom_rho_r( : ) 
            !
          ENDIF
          !
      END SELECT
    ENDIF
    !
    !
    vbot(:) = vbot(:) + total_atom_rho_r( : ) 
    !
    DEALLOCATE( wfcatomg )
    !
    ! reset total atom density for next atom
    !
    total_atom_rho_r = 0.D0
    !
  ENDDO ! atom
  !
  ! hirshfeld normalization
  !
  vbottot = SUM(vbot)
  CALL mp_sum( vbottot, intra_bgrp_comm )
  normfac=REAL(nelec,KIND=DP)/(REAL(vbottot,KIND=DP) *dv)
  vtop = normfac * vtop
  vbot = normfac * vbot
  !
  ! Check if mag of vbot is less than cutoff
  !
  DO ir = 1, n
    !
    IF ( ABS(REAL( vbot(ir) )) .lt. cutoff ) THEN
      !
      ! vbot is less than cutoff, zero the output potential
      !
      vtop(ir,1)=0.D0  
      vtop(ir,2)=0.D0  
      !
    ELSE
      !
      ! vbot is => cutoff
      !
      vtop(ir,1) = vtop(ir,1) / vbot(ir)
      vtop(ir,2) = vtop(ir,2) / vbot(ir)
      !
    ENDIF
    !
  ENDDO
  !
  v = REAL(vtop,KIND=DP)
  !
  DEALLOCATE( wfcatomr )
  !
  RETURN
  !
END SUBROUTINE calc_hirshfeld_v

SUBROUTINE EPCDFT_FORCE(force,rho)
  !--------------------------------------------------------------------------
  ! 
  !  Gamma only
  USE kinds,         ONLY : DP
!   USE scf,           ONLY : rho
  USE constants,     ONLY : fpi, eps8, e2, au_debye
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, zv
  USE cell_base,     ONLY : alat, at, omega, bg, tpiba2,tpiba
  USE force_mod,     ONLY : lforce
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_images,     ONLY : intra_image_comm
  USE mp_bands,      ONLY : me_bgrp, intra_bgrp_comm
  USE fft_base,      ONLY : dfftp, dfftp, dffts
  USE mp,            ONLY : mp_bcast, mp_sum
  USE mp_pools,             ONLY : inter_pool_comm, intra_pool_comm, me_pool, &
                                    nproc_pool
  USE control_flags, ONLY : iverbosity
  !
  ! ... Coulomb USE
  ! 
  USE wvfct,                ONLY : npwx, g2kin
  USE klist,                ONLY : xk, nks, igk_k,ngk
  USE gvecw,                ONLY : ecutwfc
  !USE wvfct,                ONLY : npw, npwx, ecutwfc, igk, g2kin
  USE gvect,                ONLY : nl, ngm, g, nlm
  USE uspp_param,           ONLY : upf
  USE gvecs,                ONLY : nls, nlsm
  USE fft_interfaces,       ONLY : invfft,fwfft
  USE wavefunctions_module, ONLY : psic
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_bcast, mp_sum
  USE io_global,     ONLY : stdout,ionode, ionode_id
  USE klist,         ONLY : nelec, xk, nks
  USE control_flags,        ONLY : gamma_only
  USE cell_base,     ONLY : omega
  USE io_files,      ONLY : tmp_dir, prefix

  ! most important control variables are not in epcdft_mod
  USE epcdft, ONLY : do_epcdft,epcdft_locs,epcdft_guess,nconstr_epcdft,epcdft_type
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP) :: dv                              ! volume element
  REAL(DP),INTENT(IN) :: rho(dfftp%nnr,nspin) ! the density whose dipole is computed
  REAL(DP),INTENT(OUT) :: force(3,nat) !
  !
  ! local variables
  !
  INTEGER :: index0, i, j, k, n, ipol, ir, na, ip, ik, npw,ig
  REAL(DP) :: length, vamp, value, sawarg, e_dipole, ion_dipole,gvec
  !
  ! ... Coulomb Vars
  !

  REAL( DP )   :: dist, mindonor, minacceptor
  REAL( DP )   :: r( 3 ), myr(3), s( 3 ), cm(3)
  REAL( DP )   :: inv_nr1, inv_nr2, inv_nr3
  REAL( DP )   :: thresh
  

  INTEGER :: nt, nb, l, m, nwfcm ,nwfc,icon
  COMPLEX(DP), ALLOCATABLE :: wfcatomg (:,:) ! atomic wfcs in g
  COMPLEX(DP), ALLOCATABLE :: wfcatomr (:) ! atomic wfcs in r
  COMPLEX(DP), ALLOCATABLE :: total_atom_rho_r (:) ! atomic wfcs in r - collected
  COMPLEX(DP), ALLOCATABLE :: shifted_atom1 (:,:) ! atomic wfcs in r -shifted and collected
  COMPLEX(DP), ALLOCATABLE :: shifted_atom2 (:,:) ! atomic wfcs in r -shifted and collected
  INTEGER :: orbi ! orbital index for wfcatom
  INTEGER :: nfuncs, lmax ! number of functions as derived from lmax
  REAL(DP) :: orboc ! occupation of orbital
  REAL(DP) :: dx ! finite shift in direction
  REAL(DP) :: force_idir ! force in given direction
  COMPLEX(DP) :: v(dfftp%nnr,2) ! hirshfeld weighting function
  COMPLEX(DP) :: vtop(dfftp%nnr,2) ! top of the hirshfeld potential fraction
  COMPLEX(DP) :: vbot(dfftp%nnr) ! bottom of the hirshfeld potential fraction

  LOGICAL :: write_debug_cubes=.false.

  COMPLEX(DP) :: normfac, ival1,ival2
  COMPLEX(DP) :: cutoff
  COMPLEX(DP) :: vbottot, val
  character(len=1024) :: filename
  LOGICAL :: hirshfeld=.true.
  LOGICAL :: do_analytical_gradient=.false.
  
  !
  force=0._dp
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )

  ! set up to be in bohr via one_atom_shifted_wfc
  dx=1e-6
  
  ALLOCATE( wfcatomr(dfftp%nnr) )
  ALLOCATE( total_atom_rho_r(dfftp%nnr) )
  ALLOCATE( shifted_atom1(3,dfftp%nnr) )
  ALLOCATE( shifted_atom2(3,dfftp%nnr) )
  !
  n=dfftp%nnr
  ik=1
  npw=ngk(ik)
  wfcatomr = 0.D0
  total_atom_rho_r=0.D0
  shifted_atom1=0.d0
  shifted_atom2=0.d0
  vtop = 0.D0
  vbot = 0.D0
  v = 0.D0
  psic = 0.D0
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  cutoff = 1.D-6
  !
  ! load atomic wfcs
  CALL gk_sort (xk (1, 1), ngm, g, ecutwfc / tpiba2, npw, igk_k(1,ik), g2kin)
  !
  ! construct hirshfeld looping over all atomic states
  !
  DO na = 1, nat ! for each atom
    nt = ityp (na) ! get atom type
    nwfc=sum(upf(nt)%oc(:))
    ALLOCATE( wfcatomg(npwx, nwfc) )

    CALL one_atom_wfc (1, wfcatomg, na,nwfc) 
    orbi = 0 
    DO nb = 1, upf(nt)%nwfc ! for each orbital
      if (upf(nt)%oc(nb) > 0.d0) then
        l = upf(nt)%lchi(nb) ! get l 
        DO m = 1, 2*l+1 ! mag num
          ! add all orbitals 
          ! each state weighted by orboc
          orbi = orbi + 1 ! update orbital index
          !  
          orboc = REAL( upf(nt)%oc(nb) ,KIND=DP) / REAL( 2*l+1 ,KIND=DP)
          !
          wfcatomr = 0.D0
          wfcatomr(nl(igk_k(1:npw,1)))  = wfcatomg(1:npw,orbi)
          IF(gamma_only) wfcatomr(nlm(igk_k(1:npw,1))) = CONJG( wfcatomg(1:npw,orbi) )
          ! convert atomic(g) -> |atomic(r)|^2
          CALL invfft ('Dense', wfcatomr(:), dfftp)
          wfcatomr(:) = wfcatomr(:) * CONJG(wfcatomr(:))
          !
          vbot(:) = vbot(:) + orboc * wfcatomr( : ) 
          !
        ENDDO ! m
      ENDIF
    ENDDO ! l 
    DEALLOCATE( wfcatomg )
  ENDDO ! atom

  ! potential normalization
  vbottot = SUM(vbot) 
  CALL mp_sum( vbottot, intra_bgrp_comm )
  normfac=REAL(nelec,KIND=DP)/(REAL(vbottot,KIND=DP)*dv)
  vbot(:) = normfac * vbot(:)

  DO icon=1,nconstr_epcdft
    DO na = 1, nat ! for each atom
      total_atom_rho_r=0.D0
      nt = ityp (na) ! get atom type
      nwfc=sum(upf(nt)%oc(:))
      ALLOCATE( wfcatomg(npwx, nwfc) )
      CALL one_atom_wfc (1, wfcatomg, na,nwfc) 
      orbi = 0 
      DO nb = 1, upf(nt)%nwfc ! for each orbital
      !
        IF (upf(nt)%oc(nb) > 0.d0) THEN ! if occupied
          !
          l = upf(nt)%lchi(nb) ! get l 
          !
          DO m = 1, 2*l+1 ! mag num
            !
            ! add all orbitals, each state weighted by orboc
            orbi = orbi + 1 ! update orbital index
            !  
            orboc = REAL( upf(nt)%oc(nb) ,KIND=DP) / REAL( 2*l+1 ,KIND=DP)
            !
            ! bring to real space
            !
            wfcatomr = 0.D0
            wfcatomr(nl(igk_k(1:npw,1)))  = wfcatomg(1:npw,orbi)
            IF(gamma_only) wfcatomr(nlm(igk_k(1:npw,1))) = CONJG( wfcatomg(1:npw,orbi) )
            !
            ! convert atomic(g) -> |atomic(r)|^2
            !
            CALL invfft ('Dense', wfcatomr(:), dfftp)
            wfcatomr(:) = wfcatomr(:) * CONJG(wfcatomr(:))
            !
            ! add orbital density to total density of atom
            !
            total_atom_rho_r(:) = total_atom_rho_r(:) + orboc * wfcatomr(:)
            !
            if (write_debug_cubes) THEN
	            if (na < 10) then
	              write(filename,"(A6,I1,A1,I1,A1,I1)") "atomic",na,"_",l,"_",m
	            else
	              write(filename,"(A6,I2,A1,I1,A1,I1)") "atomic",na,"_",l,"_",m
	            ENDIF
	            CALL write_cube_r ( 84332, filename,  REAL(wfcatomr,KIND=DP))
	        ENDIF ! write_debug_cube
            !
          ENDDO ! m
        ENDIF ! end if occupied
      ENDDO ! nb 
      SELECT CASE( epcdft_type(icon) )
        CASE('charge')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - total_atom_rho_r( : ) 
          ENDIF
        CASE('spin')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + total_atom_rho_r( : ) 
          ENDIF
        CASE('delta_charge')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - total_atom_rho_r( : ) 
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,1) = vtop(:,1) + total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + total_atom_rho_r( : ) 
            !
          ENDIF   
        CASE('delta_spin')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + total_atom_rho_r( : ) 
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,1) = vtop(:,1) + total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - total_atom_rho_r( : ) 
            !
          ENDIF
        CASE('delta_alpha')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - total_atom_rho_r( : ) 
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,1) = vtop(:,1) + total_atom_rho_r( : ) 
            !
          ENDIF
        CASE('delta_beta')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,2) = vtop(:,2) - total_atom_rho_r( : ) 
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,2) = vtop(:,2) + total_atom_rho_r( : ) 
            !
          ENDIF
      END SELECT
      DEALLOCATE( wfcatomg )
    END DO !atom

    vtop = normfac * vtop 

    v(:,1) = vtop(:,1) / vbot
    v(:,2) = vtop(:,2) / vbot
    DO ir = 1, n
      if (ABS(REAL(vbot(ir))).lt.REAL(cutoff)) THEN
        v(ir,1)=0.D0  
        v(ir,2)=0.D0  
      ENDIF
      ! NAN CHECK
      if (v(ir,1) /= v(ir,1))  v(ir,1)=0.D0
      if (v(ir,2) /= v(ir,2))  v(ir,2)=0.D0
    ENDDO

    if (write_debug_cubes) THEN
      write(filename,*) "hirshfeld_v"
      CALL write_cube_r ( 9519395, filename, REAL(v(:,1)))

      write(filename,*) "rho_up"
      CALL write_cube_r ( 9519395, filename, rho(:,1))

      write(filename,*) "rho_down"
      CALL write_cube_r ( 9519395, filename, rho(:,2))

      write(filename,*) "vbot"
      CALL write_cube_r ( 9519395, filename, REAL(vbot(:)))
    ENDIF


    DO na = 1, nat ! for each atom
      !
      ALLOCATE( wfcatomg(npwx, nwfc) )
      !
      if (do_analytical_gradient) THEN 
        total_atom_rho_r=0.D0
        nt = ityp (na) ! get atom type
        nwfc=sum(upf(nt)%oc(:))
        wfcatomg=0.D0
        CALL one_atom_wfc (1, wfcatomg, na,nwfc) 
        orbi = 0 
        DO nb = 1, upf(nt)%nwfc ! for each orbital
        !
        IF (upf(nt)%oc(nb) > 0.d0) THEN ! if occupied
            !
            l = upf(nt)%lchi(nb) ! get l 
            !
            DO m = 1, 2*l+1 ! mag num
              !
              ! add all orbitals, each state weighted by orboc
              orbi = orbi + 1 ! update orbital index
              !  
              orboc = REAL( upf(nt)%oc(nb) ,KIND=DP) / REAL( 2*l+1 ,KIND=DP)
              !
              ! bring to real space
              !
              wfcatomr = 0.D0
              wfcatomr(nl(igk_k(1:npw,1)))  = wfcatomg(1:npw,orbi)
              IF(gamma_only) wfcatomr(nlm(igk_k(1:npw,1))) = CONJG( wfcatomg(1:npw,orbi) )
              !
              ! convert atomic(g) -> |atomic(r)|^2
              !
              CALL invfft ('Dense', wfcatomr(:), dfftp)
              wfcatomr(:) = wfcatomr(:) * CONJG(wfcatomr(:))
              !
              ! add orbital density to total density of atom
              !
              total_atom_rho_r(:) = total_atom_rho_r(:) + orboc * wfcatomr(:)
              !
            ENDDO ! m
          ENDIF ! end if occupied
        ENDDO ! nb 
        total_atom_rho_r(:)=total_atom_rho_r(:)*normfac

        if (write_debug_cubes) THEN 
	        if (na < 10) then
	          write(filename,"(A8,I1)") "atomic_0",na
	        else
	          write(filename,"(A7,I2)") "atomic_",na
	        ENDIF
	        CALL write_cube_r ( 9519395, filename, REAL(total_atom_rho_r(:) ))
	    ENDIF ! write_debug_cubes

        shifted_atom1=0.D0!       
        CALL gradient(dfftp%nnr,total_atom_rho_r,ngm,g,nl,shifted_atom1)

        DO ipol=1,3
          shifted_atom1(ipol,:)=(shifted_atom1(ipol,:) - total_atom_rho_r(:))/vbot(:)
        END DO

      ELSE  ! DO DERIVATIVE VIA FINITE DIFFERENCE
        DO ipol=1,3        
          ! Negative displacement, stored in total_atom_rho_r
          shifted_atom1(ipol,:)=0.d0
          wfcatomg=0.D0
          CALL one_atom_shifted_wfc (1, wfcatomg, na, nwfc, ipol, -dx) 
          orbi = 0 
          DO nb = 1, upf(nt)%nwfc ! for each orbital
            if (upf(nt)%oc(nb) > 0.d0) then
              l = upf(nt)%lchi(nb) ! get l 
              DO m = 1, 2*l+1 ! mag num

                ! add all orbitals 
                ! each state weighted by orboc
                orbi = orbi + 1 ! update orbital index
                orboc = REAL( upf(nt)%oc(nb) ,KIND=DP) / REAL( 2*l+1 ,KIND=DP)
                wfcatomr = 0.D0
                wfcatomr(nl(igk_k(1:npw,1)))  = wfcatomg(1:npw,orbi)
                IF(gamma_only) wfcatomr(nlm(igk_k(1:npw,1))) = CONJG( wfcatomg(1:npw,orbi) )

                ! convert atomic(g) -> |atomic(r)|^2
                CALL invfft ('Dense', wfcatomr(:), dfftp)
                wfcatomr(:) = wfcatomr(:) * CONJG(wfcatomr(:))  
                shifted_atom1(ipol,:)=shifted_atom1(ipol,:)+ orboc * wfcatomr(:)
              ENDDO ! m
            ENDIF
          ENDDO ! nb
          shifted_atom1(ipol,:)=shifted_atom1(ipol,:)*normfac

          if (write_debug_cubes) THEN
            if (na < 10) then
              write(filename,"(A8,I1,A2,I1)") "atomic_0",na,"_-",ipol
            else
              write(filename,"(A7,I2,A2,I1)") "atomic_",na,"_-",ipol
            ENDIF
            CALL write_cube_r ( 9519395, filename, REAL(shifted_atom1(ipol,:)))
          ENDIF

          ! Positive displacement
          shifted_atom2(ipol,:)=0.d0
          wfcatomg=0.D0
          CALL one_atom_shifted_wfc (1, wfcatomg, na, nwfc, ipol, dx) 
          orbi = 0 
          DO nb = 1, upf(nt)%nwfc ! for each orbital
            if (upf(nt)%oc(nb) > 0.d0) then
              l = upf(nt)%lchi(nb) ! get l 
              DO m = 1, 2*l+1 ! mag num

                ! add all orbitals 
                ! each state weighted by orboc
                orbi = orbi + 1 ! update orbital index
                orboc = REAL( upf(nt)%oc(nb) ,KIND=DP) / REAL( 2*l+1 ,KIND=DP)
                wfcatomr = 0.D0
                wfcatomr(nl(igk_k(1:npw,1)))  = wfcatomg(1:npw,orbi)
                IF(gamma_only) wfcatomr(nlm(igk_k(1:npw,1))) = CONJG( wfcatomg(1:npw,orbi) )

                ! convert atomic(g) -> |atomic(r)|^2
                CALL invfft ('Dense', wfcatomr(:), dfftp)
                wfcatomr(:) = wfcatomr(:) * CONJG(wfcatomr(:))
                shifted_atom2(ipol,:)=shifted_atom2(ipol,:)+ orboc * wfcatomr(:)
              ENDDO ! m
            ENDIF
          ENDDO ! nb 

          shifted_atom2(ipol,:)=shifted_atom2(ipol,:)*normfac

          if (write_debug_cubes) THEN
            if (na < 10) then
              write(filename,"(A8,I1,A1,I1)") "atomic_0",na,"_",ipol
            else
              write(filename,"(A7,I2,A1,I1)") "atomic_",na,"_",ipol
            ENDIF
            CALL write_cube_r ( 9519395, filename, REAL(shifted_atom2(ipol,:)))
          ENDIF

          ! shifted_atom now has drho_i/tot_rho in it, calculated using the two-point center formula
          shifted_atom1(ipol,:)=(shifted_atom2(ipol,:) - shifted_atom1(ipol,:))/(2*dx*vbot(:))

          ! CLEAN IT UP AS ONE CLEANS UP THE POTENTIAL SHAPE
          DO ir = 1, n
            if (ABS(REAL(vbot(ir))).lt.REAL(cutoff)) THEN
              shifted_atom1(ipol,ir)=0.D0  
            ENDIF
            ! NAN CHECK
            if (shifted_atom1(ipol,ir) /= shifted_atom1(ipol,ir))  shifted_atom1(ipol,ir)=0.D0
          ENDDO
        END DO ! ipol
      ENDIF

      DO ipol=1,3 
        SELECT CASE( epcdft_type(icon) )
        ! -V rho w*(-drho_i/tot_rho) 
        
        CASE('charge')
          force_idir=REAL(SUM((v(:,1)*rho(:,1)+v(:,2)*rho(:,2))*shifted_atom1(ipol,:)))
        CASE('spin')
          force_idir=REAL(SUM((-v(:,1)*rho(:,1)+v(:,2)*rho(:,2))*shifted_atom1(ipol,:)))
        CASE('delta_charge')
          force_idir=REAL(SUM((v(:,1)*rho(:,1)+v(:,2)*rho(:,2))*shifted_atom1(ipol,:)))
        CASE('delta_spin')
          force_idir=REAL(SUM((-v(:,1)*rho(:,1)+v(:,2)*rho(:,2))*shifted_atom1(ipol,:)))
        CASE('delta_alpha')
          force_idir=REAL(SUM((-v(:,1)*rho(:,1)*shifted_atom1(ipol,:))))
        CASE('delta_beta')
          force_idir=REAL(SUM((-v(:,2)*rho(:,2)*shifted_atom1(ipol,:))))
        END SELECT

        ! +/- (-V) rho (-drho_i/tot_rho)
        SELECT CASE( epcdft_type(icon) )
        CASE('charge')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            force_idir=force_idir+REAL(SUM((rho(:,1)+rho(:,2))*shifted_atom1(ipol,:)))
          ENDIF
        CASE('spin')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            force_idir=force_idir+REAL(SUM((-rho(:,1)+rho(:,2))*shifted_atom1(ipol,:)))
          ENDIF
        CASE('delta_charge')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            force_idir=force_idir+REAL(SUM((rho(:,1)+rho(:,2))*shifted_atom1(ipol,:)))
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN ! atom in donor
            force_idir=force_idir-REAL(SUM((rho(:,1)+rho(:,2))*shifted_atom1(ipol,:)))
          ENDIF
        CASE('delta_spin')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            force_idir=force_idir+REAL(SUM((-rho(:,1)+rho(:,2))*shifted_atom1(ipol,:)))
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN ! atom in donor
            force_idir=force_idir-REAL(SUM((-rho(:,1)+rho(:,2))*shifted_atom1(ipol,:)))
          ENDIF
        CASE('delta_alpha')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            force_idir=force_idir+REAL(SUM((-rho(:,1))*shifted_atom1(ipol,:)))
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            force_idir=force_idir-REAL(SUM((-rho(:,1))*shifted_atom1(ipol,:)))
          ENDIF
        CASE('delta_beta')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            force_idir=force_idir+REAL(SUM((-rho(:,2))*shifted_atom1(ipol,:)))
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            force_idir=force_idir-REAL(SUM((-rho(:,2))*shifted_atom1(ipol,:)))
          ENDIF
        END SELECT
        CALL mp_sum( force_idir, intra_bgrp_comm )
        force(ipol,na)=force_idir*dv*epcdft_guess(icon)
      END DO ! ipol
      DEALLOCATE(wfcatomg)
    ENDDO ! atom
  ENDDO ! icon

  DEALLOCATE( wfcatomr )
  DEALLOCATE( shifted_atom1 )
  DEALLOCATE( shifted_atom2 )
  DEALLOCATE( total_atom_rho_r )
  !
  RETURN
  !
	

  CONTAINS
	  !
	SUBROUTINE write_cube_r ( iu, fname, wfc_distr )
	  ! -----------------------------------------------------------------
	  !
	  ! For debugging
	  !
	  USE pwcom,                 ONLY : npw,npwx
	  USE kinds,                 ONLY : DP
	  USE cell_base,             ONLY : celldm, at, bg
	  USE ions_base,             ONLY : nat, tau, atm, ityp
	  USE fft_base,              ONLY : dfftp !, grid_gather
	  USE scatter_mod,           ONLY : gather_grid
	  USE mp_global,             ONLY : me_bgrp,root_bgrp
	  !
	  IMPLICIT NONE
	  !
	  ! I/O 
	  !
	  INTEGER,INTENT(IN) :: iu
	  CHARACTER(LEN=256),INTENT(IN) :: fname
	  REAL(DP),INTENT(IN) :: wfc_distr(dfftp%nnr)
	  ! 
	  ! Workspace
	  !
	  REAL(DP)         :: alat
	  INTEGER          :: nr1, nr2, nr3, nr1x, nr2x, nr3x
	  INTEGER          :: i, nt, i1, i2, i3, at_num, ir
	  INTEGER, EXTERNAL  :: atomic_number
	  REAL(DP)    :: at_chrg, tpos(3), inpos(3)
	  REAL(DP) :: wfc_gat(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x)
	  !
	  wfc_gat=0.0_DP
	  CALL gather_grid(dfftp,wfc_distr,wfc_gat)
	  !
	  !      WRITE A FORMATTED 'DENSITY-STYLE' CUBEFILE VERY SIMILAR
	  !      TO THOSE CREATED BY THE GAUSSIAN PROGRAM OR THE CUBEGEN UTILITY.
	  !      THE FORMAT IS AS FOLLOWS (LAST CHECKED AGAINST GAUSSIAN 98):
	  ! 
	  !      LINE   FORMAT      CONTENTS
	  !      ===============================================================
	  !       1     A           TITLE
	  !       2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
	  !       3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
	  !       4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
	  !       #ATOMS LINES OF ATOM COORDINATES:
	  !       ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
	  !       REST: 6E13.5      CUBE DATA
	  ! 
	  !      ALL COORDINATES ARE GIVEN IN ATOMIC UNITS.
	  !
	  alat = celldm(1)
	  nr1 = dfftp%nr1
	  nr2 = dfftp%nr2
	  nr3 = dfftp%nr3
	  nr1x= dfftp%nr1x
	  nr2x= dfftp%nr2x
	  nr3x= dfftp%nr3x
	  !
	!  wfc_gat(:)=wfc_gat(:)/DBLE( nr1*nr2*nr3 )
	  !
	  IF( me_bgrp == root_bgrp) THEN
	     OPEN(UNIT=iu,FILE=TRIM(ADJUSTL(fname))//".cube")
	     !
	     WRITE(iu,*) 'Cubfile created from PWScf calculation'
	     WRITE(iu,*) ' Total SCF Density'
	     !                        origin is forced to (0.0,0.0,0.0)
	     WRITE(iu,'(I5,3F12.6)') nat, 0.0_DP, 0.0_DP, 0.0_DP
	     WRITE(iu,'(I5,3F12.6)') nr1, (alat*at(i,1)/DBLE(nr1),i=1,3)
	     WRITE(iu,'(I5,3F12.6)') nr2, (alat*at(i,2)/DBLE(nr2),i=1,3)
	     WRITE(iu,'(I5,3F12.6)') nr3, (alat*at(i,3)/DBLE(nr3),i=1,3)
	     !
	     DO i=1,nat
	        nt = ityp(i)
	        ! find atomic number for this atom.
	        at_num = atomic_number(TRIM(atm(nt)))
	        at_chrg= DBLE(at_num)
	        ! at_chrg could be alternatively set to valence charge
	        ! positions are in cartesian coordinates and a.u.
	        !
	        ! wrap coordinates back into cell.
	        tpos = MATMUL( TRANSPOSE(bg), tau(:,i) )
	        tpos = tpos - NINT(tpos - 0.5_DP)
	        inpos = alat * MATMUL( at, tpos )
	        WRITE(iu,'(I5,5F15.6)') at_num, at_chrg, inpos
	     ENDDO
	     !
	     CALL actual_write_cube(wfc_gat,nr1,nr2,nr3,iu)
	     !
	     CLOSE(iu)
	  ENDIF
	  !
	END SUBROUTINE
	!
	SUBROUTINE actual_write_cube(func,nr1,nr2,nr3,ounit) 
	  !
	  ! For debugging
	  !
	  USE kinds,  ONLY : DP
	  !
	  IMPLICIT NONE
	  !
	  REAL(DP) :: func(nr1,nr2,nr3)
	  INTEGER :: ounit,nr1,nr2,nr3
	  !
	  INTEGER :: i1,i2,i3
	  !
	  DO i1=1,nr1
	     DO i2=1,nr2
	        WRITE(ounit,'(6E17.5E3)') (func(i1,i2,i3),i3=1,nr3)
	     ENDDO
	  ENDDO
	  !
	END SUBROUTINE
	!
	SUBROUTINE pprint( strg, num, ustrg, a)
	  !
	  ! write strg, num and ustrg
	  ! a is used to determine the format of num
	  !
	  USE kinds, ONLY:DP
	  !
	  IMPLICIT NONE
	  REAL(DP),INTENT(IN)::num
	  CHARACTER(LEN=*),INTENT(IN)::strg, ustrg
	  CHARACTER(LEN=1),INTENT(IN)::a
	  CHARACTER(LEN=32)::trimstrg
	  CHARACTER(LEN=6)::utrimstrg
	  !
	  trimstrg = strg
	  utrimstrg = ustrg
	  !
	  SELECT CASE(a)
	    CASE('f','F')
	      WRITE(*,1)trimstrg,num,utrimstrg
	    CASE('e','E')
	      WRITE(*,2)trimstrg,num,utrimstrg
	  END SELECT
	  !
	  1 FORMAT(5x,A32,':',1x,F23.16,1x,A6)
	  2 FORMAT(5x,A32,':',1x,e23.16,1x,A6)
	  !
	END SUBROUTINE
END SUBROUTINE EPCDFT_FORCE
