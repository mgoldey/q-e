!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... poorly written by N. P. Brawand
! ... and poorly modified by M. B. Goldey
!
!
!--------------------------------------------------------------------------
SUBROUTINE add_epcdft_efield(vpoten,iflag)
  !--------------------------------------------------------------------------
  !
  !   This routine adds constraining potential for cdft to the local potential. 
  !
  !   Voronoi cells
  !
  !   User defined well around single atom
  !
  !   Hirshfeld partitioning following :
  !
  !     J. Chem. Phys. 133, 244105 (2010); http://dx.doi.org/10.1063/1.3507878 
  !     eqs. 6 & 7
  !
  USE kinds,         ONLY : DP
  USE epcdft,        ONLY : do_epcdft
  USE io_global,     ONLY : stdout,ionode
  USE lsda_mod,      ONLY : nspin
  USE fft_base,      ONLY : dfftp !, grid_scatter, grid_gather
  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP),INTENT(INOUT) :: vpoten(dfftp%nnr,nspin)! ef is added to this potential
  LOGICAL,INTENT(IN)     :: iflag            ! set to true to force recalculation of field
  !
  ! local variables
  !
  REAL(DP) :: vconstr (dfftp%nnr,nspin) ! constraining potential shape*weight
  !
  LOGICAL :: first=.TRUE.
  LOGICAL :: hirshfeld=.TRUE.

  SAVE first

  !---------------------
  !  Execution control
  !---------------------

  IF (.NOT. do_epcdft) RETURN
  IF ((.NOT.first) .AND..NOT. iflag) RETURN
  first=.FALSE.

  ! Necessary for restart/pp runs
  CALL init_at_1
  
  ! efield only needs to be added on the first iteration (of each SCF call)
  ! note that for relax calculations it has to be added
  ! again on subsequent relax steps. 
  
  ! calculate hirshfeld potential array
  !
  IF (hirshfeld) THEN 
    vconstr =0.D0
    CALL calc_hirshfeld_v(vconstr,0)
    vpoten= vpoten + vconstr 
    RETURN
  ENDIF
  !
END SUBROUTINE add_epcdft_efield
!
!
!
!--------------------------------------------------------------------------
SUBROUTINE calc_hirshfeld_v( v,iconstraint)
  !--------------------------------------------------------------------------
  ! 
  !  Gamma only
  !
  !  calculate hirshfeld potential and put into v following :
  !
  !     J. Chem. Phys. 133, 244105 (2010); http://dx.doi.org/10.1063/1.3507878 
  !     eqs. 6 & 7
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE wvfct,                ONLY : npw, npwx, ecutwfc, igk, g2kin
  USE klist,                ONLY : xk, nks
  USE gvect,                ONLY : ngm, g
  USE cell_base,            ONLY : tpiba2
  USE uspp_param,           ONLY : upf
  USE lsda_mod,      ONLY : nspin
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE fft_base,             ONLY : dfftp, dffts
  USE gvect,                ONLY : nl, nlm
  USE gvecs,                ONLY : nls, nlsm
  USE fft_interfaces,       ONLY : invfft
  USE wavefunctions_module, ONLY : psic
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_bcast, mp_sum
  USE io_global,     ONLY : stdout,ionode, ionode_id
  USE klist,         ONLY : nelec
  USE cell_base,     ONLY : omega

  ! most important control variables are not in epcdft_mod
  USE epcdft
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: v(dfftp%nnr,nspin) ! the hirshfeld potential
  INTEGER, INTENT(IN) :: iconstraint
  ! 
  INTEGER :: na, nt, nb, l, m,n,ir, nwfcm ,nwfc,icon
  COMPLEX(DP), ALLOCATABLE :: wfcatomg (:,:) ! atomic wfcs in g
  COMPLEX(DP), ALLOCATABLE :: wfcatomr (:) ! atomic wfcs in r
  COMPLEX(DP), ALLOCATABLE :: total_atom_rho_r (:) ! total density of an atom in r
  INTEGER :: orbi ! orbital index for wfcatom
  INTEGER :: nfuncs, lmax ! number of functions as derived from lmax
  REAL(DP) :: orboc ! occupation of orbital
  COMPLEX(DP) :: vtop(dfftp%nnr,nspin) ! top of the hirshfeld potential fraction
  COMPLEX(DP) :: vbot(dfftp%nnr)       ! bottom of the hirshfeld potential fraction
  COMPLEX(DP) :: normfac
  COMPLEX(DP) :: cutoff
  COMPLEX(DP) :: vbottot
  REAL(DP) :: dv
  character(len=1024) :: filename
  
  n=dfftp%nnr ! Just to be safe
  ! get the max number of wfs any atom will have
  nwfcm=0
  lmax=0
  DO na=1,nat
    DO nb = 1, upf(nt)%nwfc ! for each orbital
      if (upf(nt)%oc(nb) > 0.d0) then
        l = upf(nt)%lchi(nb) ! get l 
        if (l .gt. lmax) lmax=l
      END IF
    END DO
    nfuncs=0
    do l=0,lmax
      nfuncs=nfuncs+(2*l+1)
    end do
    if (nfuncs .gt. nwfcm) nwfcm = nfuncs
  END DO
  !
  ALLOCATE( wfcatomr(n) )
  ALLOCATE( total_atom_rho_r(n) )
  !
  wfcatomr = 0.D0
  vtop = 0.D0
  vbot = 0.D0
  v = 0.D0
  psic = 0.D0
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  cutoff = 1.D-6
  total_atom_rho_r = 0.D0
  !
  ! load atomic wfcs
  CALL gk_sort (xk (1, 1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
  !
  ! construct hirshfeld looping over all atomic states
  !
  DO na = 1, nat ! for each atom
    !
    nt = ityp (na) ! get atom type
    nwfc=sum(upf(nt)%oc(:))
    !
    ALLOCATE( wfcatomg(npwx, nwfc) )
    wfcatomg = 0.D0
    ! if (ionode) write(*,*) "Atom ",na
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
          wfcatomr(nl(igk(1:npw)))  = wfcatomg(1:npw,orbi)
          IF(gamma_only) wfcatomr(nlm(igk(1:npw))) = CONJG( wfcatomg(1:npw,orbi) )
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
          ! if (na < 10) then
          !   write(filename,"(A6,I1,A1,I1,A1,I1)") "atomic",na,"_",l,"_",m
          ! else
          !   write(filename,"(A6,I2,A1,I1,A1,I1)") "atomic",na,"_",l,"_",m
          ! ENDIF
          ! CALL write_cube_r ( 84332, filename,  REAL(wfcatomr,KIND=DP))
          !
        ENDDO ! m
      ENDIF ! end if occupied
    ENDDO ! nb 
    !
    ! add this atoms density to the hirshfeld potential, as appropriate
    IF (iconstraint .eq. 0) THEN 
      !write(*,*) iconstraint,nconstr_epcdft
      DO icon=1,nconstr_epcdft
        !write(*,*) na,icon,nconstr_epcdft,epcdft_type(icon),epcdft_guess(icon), &
        !  epcdft_locs(1,icon),epcdft_locs(2,icon),epcdft_locs(3,icon),epcdft_locs(4,icon)
        SELECT CASE( epcdft_type(icon) )
        CASE('charge')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - epcdft_guess(icon)*total_atom_rho_r( : ) 
          ENDIF
        CASE('spin')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + epcdft_guess(icon)*total_atom_rho_r( : ) 
          ENDIF
        CASE('delta_charge')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            !write(*,*) na, epcdft_locs(1,icon),epcdft_locs(2,icon)
            vtop(:,1) = vtop(:,1) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - epcdft_guess(icon)*total_atom_rho_r( : ) 
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            !write(*,*) na, epcdft_locs(3,icon),epcdft_locs(4,icon)
            vtop(:,1) = vtop(:,1) + epcdft_guess(icon)*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + epcdft_guess(icon)*total_atom_rho_r( : ) 
            !
          ENDIF   
        CASE('delta_spin')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) + epcdft_guess(icon)*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - epcdft_guess(icon)*total_atom_rho_r( : ) 
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,1) = vtop(:,1) + epcdft_guess(icon)*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - epcdft_guess(icon)*total_atom_rho_r( : ) 
            !
          ENDIF
        END SELECT
      END DO
    ELSE 
      icon=iconstraint
      !write(*,*) icon,nconstr_epcdft,epcdft_type(icon)
      SELECT CASE( epcdft_type(icon) )
        CASE('charge')
          IF( na >= epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - total_atom_rho_r( : ) 
          ENDIF
        CASE('spin')
          IF( na >= epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + total_atom_rho_r( : ) 
          ENDIF
        CASE('delta_charge')
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
        CASE('delta_spin')
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
  ! M.G.s bar
  !CALL mp_sum( orboc, intra_image_comm )
  !  
  ! call write_cube_r ( 84332, 'vtop.cube',  REAL(vtop,KIND=DP))
  ! call write_cube_r ( 84332, 'vbot.cube',  REAL(vbot,KIND=DP))
  !  
  ! hirshfeld normalization
  !
  vbottot = SUM(vbot)
  CALL mp_sum( vbottot, intra_bgrp_comm )
  normfac=REAL(nelec,KIND=DP)/(REAL(vbottot,KIND=DP)*dv)
  vtop = normfac * vtop
  vbot = normfac * vbot
  
  !
  vtop(:,1) = vtop(:,1) / vbot
  vtop(:,2) = vtop(:,2) / vbot
  DO ir = 1, n
    if (ABS(REAL(vbot(ir))).lt.REAL(cutoff)) THEN
      vtop(ir,1)=0.D0  
      vtop(ir,2)=0.D0  
    ENDIF
    if (vtop(ir,1) /= vtop(ir,1))  vtop(ir,1)=0.D0
    if (vtop(ir,2) /= vtop(ir,2))  vtop(ir,2)=0.D0
  ENDDO
  !write(*,*) real(sum(vtop))
  !
  !
  !call write_cube_r ( 84332, 'v.cube',  REAL(vtop,KIND=DP))
  v = REAL(vtop,KIND=DP)
  !
  !call write_cube_r ( 84332, 'vup.cube',  v(:,1) )
  !call write_cube_r ( 84332, 'vdown.cube',  v(:,2) )
  !
  
  DEALLOCATE( wfcatomr )
  !
  RETURN
  !
END SUBROUTINE calc_hirshfeld_v
!
!

SUBROUTINE EPCDFT_FORCE(force,rho)
  !--------------------------------------------------------------------------
  ! 
  !  Gamma only
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8, e2, au_debye
  USE ions_base,     ONLY : nat, ntyp => nsp, ityp, tau, zv
  USE cell_base,     ONLY : alat, at, omega, bg, saw, tpiba2
  USE force_mod,     ONLY : lforce
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_images,     ONLY : intra_image_comm
  USE mp_bands,      ONLY : me_bgrp, intra_bgrp_comm
  USE fft_base,      ONLY : dfftp, grid_scatter, grid_gather, dfftp, dffts
  USE mp,            ONLY : mp_bcast, mp_sum
  USE control_flags, ONLY : iverbosity
  !
  ! ... Coulomb USE
  ! 
  USE wvfct,                ONLY : npw, npwx, ecutwfc, igk, g2kin
  USE gvect,                ONLY : ngm, g
  USE uspp_param,           ONLY : upf
  USE gvect,                ONLY : nl, nlm
  USE gvecs,                ONLY : nls, nlsm
  USE fft_interfaces,       ONLY : invfft
  USE wavefunctions_module, ONLY : psic
  USE mp_bands,      ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_bcast, mp_sum
  USE io_global,     ONLY : stdout,ionode, ionode_id
  USE klist,         ONLY : nelec, xk, nks
  USE control_flags,        ONLY : gamma_only
  USE cell_base,     ONLY : omega

  ! most important control variables are not in epcdft_mod
  USE epcdft, ONLY : do_epcdft,epcdft_locs,epcdft_guess,nconstr_epcdft,epcdft_type
  !
  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP) :: dv                             ! volume element
  REAL(DP),INTENT(IN) :: rho(dfftp%nnr,nspin) ! the density whose dipole is computed
  REAL(DP),INTENT(OUT) :: force(3,nat) !
  !
  ! local variables
  !
  INTEGER :: index0, i, j, k, n, ipol, ir, na, ip
  REAL(DP) :: length, vamp, value, sawarg, e_dipole, ion_dipole
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
  COMPLEX(DP), ALLOCATABLE :: shifted_atom (:) ! atomic wfcs in r -shifted and collected
  INTEGER :: orbi ! orbital index for wfcatom
  INTEGER :: nfuncs, lmax ! number of functions as derived from lmax
  REAL(DP) :: orboc ! occupation of orbital
  REAL(DP) :: dx ! finite shift in direction
  REAL(DP) :: force_idir ! force in given direction
  COMPLEX(DP) :: v(dfftp%nnr,2) ! hirshfeld weighting function
  COMPLEX(DP) :: v2(dfftp%nnr) ! hirshfeld weighting function * rho
  COMPLEX(DP) :: vtop(dfftp%nnr,2) ! top of the hirshfeld potential fraction
  COMPLEX(DP) :: vbot(dfftp%nnr) ! bottom of the hirshfeld potential fraction

! shifted
  COMPLEX(DP) :: sv(dfftp%nnr,2) ! hirshfeld weighting function for shifted atom
  COMPLEX(DP) :: svtop(dfftp%nnr,2) ! top of the hirshfeld potential fraction
  COMPLEX(DP) :: svbot(dfftp%nnr) ! bottom of the hirshfeld potential fraction

  COMPLEX(DP) :: normfac, ival1,ival2
  COMPLEX(DP) :: cutoff
  COMPLEX(DP) :: vbottot, val
  character(len=1024) :: filename
  LOGICAL :: hirshfeld=.true.
  
  !
  force=0._dp
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )

  ! set up to be in bohr via one_atom_shifted_wfc
  dx=0.001

  ! get the max number of wfs any atom will have
  nwfcm=0
  lmax=0
  DO na=1,nat
    DO nb = 1, upf(nt)%nwfc ! for each orbital
      if (upf(nt)%oc(nb) > 0.d0) then
        l = upf(nt)%lchi(nb) ! get l 
        if (l .gt. lmax) lmax=l
      END IF
    END DO
    nfuncs=0
    do l=0,lmax
      nfuncs=nfuncs+(2*l+1)
    end do
    if (nfuncs .gt. nwfcm) nwfcm = nfuncs
  END DO

  
  ALLOCATE( wfcatomr(dfftp%nnr) )
  ALLOCATE( total_atom_rho_r(dfftp%nnr) )
  ALLOCATE( shifted_atom(dfftp%nnr) )
  !
  n=dfftp%nnr
  wfcatomg = 0.D0
  wfcatomr = 0.D0
  total_atom_rho_r=0.D0
  shifted_atom=0.d0
  vtop = 0.D0
  vbot = 0.D0
  svtop = 0.D0
  svbot = 0.D0
  v = 0.D0
  v2 = 0.D0
  sv = 0.D0
  psic = 0.D0
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  cutoff = 1.D-6
  !
  ! load atomic wfcs
  CALL gk_sort (xk (1, 1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
  !
  ! construct hirshfeld looping over all atomic states
  !
  DO na = 1, nat ! for each atom
    nt = ityp (na) ! get atom type
    nwfc=sum(upf(nt)%oc(:))
    ALLOCATE( wfcatomg(npwx, nwfc) )
!     if (ionode) write(*,*) "Atom ",na
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
          wfcatomr(nl(igk(1:npw)))  = wfcatomg(1:npw,orbi)
          IF(gamma_only) wfcatomr(nlm(igk(1:npw))) = CONJG( wfcatomg(1:npw,orbi) )

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
  vbottot = SUM(vbot)*dv
  CALL mp_sum( vbottot, intra_bgrp_comm )
  normfac=REAL(nelec,KIND=DP)/(REAL(vbottot,KIND=DP)*dv)
  vbot(:) = normfac * vbot(:)

  DO icon=1,nconstr_epcdft
    DO na = 1, nat ! for each atom
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
            wfcatomr(nl(igk(1:npw)))  = wfcatomg(1:npw,orbi)
            IF(gamma_only) wfcatomr(nlm(igk(1:npw))) = CONJG( wfcatomg(1:npw,orbi) )
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
            ! if (na < 10) then
            !   write(filename,"(A6,I1,A1,I1,A1,I1)") "atomic",na,"_",l,"_",m
            ! else
            !   write(filename,"(A6,I2,A1,I1,A1,I1)") "atomic",na,"_",l,"_",m
            ! ENDIF
            ! CALL write_cube_r ( 84332, filename,  REAL(wfcatomr,KIND=DP))
            !
          ENDDO ! m
        ENDIF ! end if occupied
      ENDDO ! nb 
      SELECT CASE( epcdft_type(icon) )
        CASE('charge')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - normfac*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - normfac*total_atom_rho_r( : ) 
          ENDIF
        CASE('spin')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - normfac*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + normfac*total_atom_rho_r( : ) 
          ENDIF
        CASE('delta_charge')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            !write(*,*) na, epcdft_locs(1,icon),epcdft_locs(2,icon)
            vtop(:,1) = vtop(:,1) - normfac*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - normfac*total_atom_rho_r( : ) 
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            !write(*,*) na, epcdft_locs(3,icon),epcdft_locs(4,icon)
            vtop(:,1) = vtop(:,1) + normfac*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + normfac*total_atom_rho_r( : ) 
            !
          ENDIF   
        CASE('delta_spin')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            vtop(:,1) = vtop(:,1) - normfac*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) + normfac*total_atom_rho_r( : ) 
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            vtop(:,1) = vtop(:,1) + normfac*total_atom_rho_r( : ) 
            vtop(:,2) = vtop(:,2) - normfac*total_atom_rho_r( : ) 
            !
          ENDIF
      END SELECT
    END DO !atom

    vtop(:,1) = vtop(:,1) / vbot
    vtop(:,2) = vtop(:,2) / vbot
    DO ir = 1, n
      if (ABS(REAL(vbot(ir))).lt.REAL(cutoff)) THEN
        vtop(ir,1)=0.D0  
        vtop(ir,2)=0.D0  
      ENDIF
      if (vtop(ir,1) /= vtop(ir,1))  vtop(ir,1)=0.D0
      if (vtop(ir,2) /= vtop(ir,2))  vtop(ir,2)=0.D0
    ENDDO
    v2(:)=v(:,1)*rho(:,1)+v(:,2)*rho(:,2)

    DO na = 1, nat ! for each atom
      total_atom_rho_r=0.d0
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
            orboc = REAL( upf(nt)%oc(nb) ,KIND=DP) / REAL( 2*l+1 ,KIND=DP)
            wfcatomr = 0.D0
            wfcatomr(nl(igk(1:npw)))  = wfcatomg(1:npw,orbi)
            IF(gamma_only) wfcatomr(nlm(igk(1:npw))) = CONJG( wfcatomg(1:npw,orbi) )

            ! convert atomic(g) -> |atomic(r)|^2
            CALL invfft ('Dense', wfcatomr(:), dfftp)
            wfcatomr(:) = wfcatomr(:) * CONJG(wfcatomr(:))
            total_atom_rho_r(:)=total_atom_rho_r(:)+wfcatomr(:)
          ENDDO ! m
        ENDIF
      ENDDO ! nb
      total_atom_rho_r(:)=total_atom_rho_r(:)*normfac

      DO ipol=1,3
        ! reset everything
        svtop=vtop
        svbot=vbot

        shifted_atom=0.d0
        CALL one_atom_shifted_wfc (1, wfcatomg, na,nwfc,ipol,dx) 
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
              wfcatomr(nl(igk(1:npw)))  = wfcatomg(1:npw,orbi)
              IF(gamma_only) wfcatomr(nlm(igk(1:npw))) = CONJG( wfcatomg(1:npw,orbi) )

              ! convert atomic(g) -> |atomic(r)|^2
              CALL invfft ('Dense', wfcatomr(:), dfftp)
              wfcatomr(:) = wfcatomr(:) * CONJG(wfcatomr(:))
              shifted_atom(:)=shifted_atom(:)+wfcatomr(:)
            ENDDO ! m
          ENDIF
        ENDDO ! l 
        shifted_atom(:)=shifted_atom(:)*normfac
        SELECT CASE( epcdft_type(icon) )
        CASE('charge')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            svtop(:,1) = svtop(:,1) - (shifted_atom(:) - total_atom_rho_r(:))
            svtop(:,2) = svtop(:,2) - (shifted_atom(:) - total_atom_rho_r(:))
          ENDIF
        CASE('spin')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            svtop(:,1) = svtop(:,1) - (shifted_atom(:) - total_atom_rho_r(:))
            svtop(:,2) = svtop(:,2) + (shifted_atom(:) - total_atom_rho_r(:))
          ENDIF
        CASE('delta_charge')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            !write(*,*) na, epcdft_locs(1,icon),epcdft_locs(2,icon)
            svtop(:,1) = svtop(:,1) - (shifted_atom(:) - total_atom_rho_r(:))
            svtop(:,2) = svtop(:,2) - (shifted_atom(:) - total_atom_rho_r(:))
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            !write(*,*) na, epcdft_locs(3,icon),epcdft_locs(4,icon)
            svtop(:,1) = svtop(:,1) + (shifted_atom(:) - total_atom_rho_r(:))
            svtop(:,2) = svtop(:,2) + (shifted_atom(:) - total_atom_rho_r(:))
            !
          ENDIF   
        CASE('delta_spin')
          IF( na .ge. epcdft_locs(1,icon) .and. na .le. epcdft_locs(2,icon) )THEN ! atom in acceptor
            !
            svtop(:,1) = svtop(:,1) - (shifted_atom(:) - total_atom_rho_r(:))
            svtop(:,2) = svtop(:,2) + (shifted_atom(:) - total_atom_rho_r(:))
          ELSE IF( na .ge. epcdft_locs(3,icon) .and. na .le. epcdft_locs(4,icon) )THEN
            ! 
            svtop(:,1) = svtop(:,1) + (shifted_atom(:) - total_atom_rho_r(:))
            svtop(:,2) = svtop(:,2) - (shifted_atom(:) - total_atom_rho_r(:))
            !
          ENDIF
        END SELECT
        svbot(:) = vbot(:) + (shifted_atom(:)-total_atom_rho_r(:))
            

        sv(:,1) = svtop(:,1) / svbot
        sv(:,2) = svtop(:,2) / svbot
        DO ir = 1, n
          if (ABS(REAL(svbot(ir))).lt.REAL(cutoff)) THEN
            sv(ir,1)=0.D0  
            sv(ir,2)=0.D0  
          ENDIF
          if (sv(ir,1) /= sv(ir,1))  sv(ir,1)=0.D0
          if (sv(ir,2) /= sv(ir,2))  sv(ir,2)=0.D0
        ENDDO
        sv(:,1)=sv(:,1)*rho(:,1)+sv(:,2)*rho(:,2)

        force_idir=-epcdft_guess(icon)*dv*(REAL(sum(sv(:,1)) ,KIND=DP)-real(sum(v2) ,KIND=DP))/dx
        CALL mp_sum(force_idir, intra_bgrp_comm )
        force(ipol,na)=force(ipol,na)+force_idir
      ENDDO ! ipol
      DEALLOCATE( wfcatomg )
    ENDDO ! atom
  ENDDO ! icon
!     CALL write_cube_r ( 84332, "v_hirsh.cub",  REAL(v,KIND=DP))
!     CALL write_cube_r ( 84332, "rhoup.cub",  REAL(rho(:,1),KIND=DP))
!     CALL write_cube_r ( 84332, "rhodown.cub",  REAL(rho(:,2),KIND=DP))
  DEALLOCATE( wfcatomr )
  !
  RETURN
  !
END SUBROUTINE EPCDFT_FORCE
