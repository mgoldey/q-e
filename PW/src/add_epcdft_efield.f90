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
!--------------------------------------------------------------------------
SUBROUTINE add_epcdft_efield(vpoten,iflag)
  !--------------------------------------------------------------------------
  !
  !   This routine adds the constraining potential for cdft. 
  !
  USE kinds,            ONLY : DP
  USE epcdft,           ONLY : do_epcdft, epcdft_surface
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
  REAL(DP),INTENT(INOUT) :: vpoten(dfftp%nnr,nspin)! ef is added to this potential
  LOGICAL,INTENT(IN)     :: iflag            ! set to true to force recalculation of field
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
  ! efield only needs to be added on the first iteration (of each SCF call)
  ! note that for relax calculations it has to be added
  ! again on subsequent relax steps. 
  !
  ! calculate hirshfeld potential array
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
!
!--------------------------------------------------------------------------
SUBROUTINE calc_epcdft_surface_field( vin, x0, qq, dipole )
  !--------------------------------------------------------------------------
  ! 
  ! add the monopole and dipole potential due to an image charge
  !
  ! vin is zeroed and surface energy is returned.
  !
  !
  !                  metal@z=0
  !                    |
  !                    |
  !  (image)-----z-----|-----z-----(charge_center)
  !                    |
  !                    |
  ! 
  !
  USE kinds,     ONLY : DP
  USE scf,       ONLY : rho
  USE lsda_mod,  ONLY : nspin
  USE fft_base,  ONLY : dfftp
  USE cell_base, ONLY : alat, at
  USE constants, ONLY : au_debye
  USE ions_base, ONLY : nat, ityp, tau, zv
  USE io_global, ONLY : stdout,ionode
  USE mp_bands,  ONLY : me_bgrp
  USE io_global, ONLY : ionode
  USE constants, ONLY : e2
  USE epcdft,    ONLY : epcdft_surface_cm_start, epcdft_surface_cm_end
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: vin(dfftp%nnr,nspin) ! the field due to slab
  REAL(DP), INTENT(INOUT) :: x0(3) ! center of charge of system
  REAL(DP), INTENT(INOUT) :: dipole(3)!, quadrupole(3) ! total dips
  REAL(DP), INTENT(INOUT) :: qq ! total charge
  !
  ! local variables
  !
  INTEGER :: index0, i, j, k, ip, ia, ir
  REAL( DP )   :: inv_nr1, inv_nr2, inv_nr3
  REAL(DP) :: zvia
  REAL( DP )   :: r( 3 ), rmag ! position in cell and |r|
  REAL(DP) :: zvtot ! tot ion charge
  REAL(DP) :: e_dipole(0:3), e_quadrupole(3) ! electronic monopole dipole and quadrupole
  REAL(DP) :: quadrupole(3) 
  REAL(DP) :: dipole_ion(3), quadrupole_ion(3) ! ion dips
  REAL(DP), EXTERNAL :: ddot
  !
  ! debug
  !
  LOGICAL :: debug = .false.
  !
  !---------------------
  !  Variable initialization
  !---------------------
  !
  vin = 0.D0
  r = 0.D0
  qq = 0.D0
  zvtot = 0.D0
  x0 = 0.D0
  e_dipole = 0.D0
  e_quadrupole = 0.D0
  dipole_ion = 0.D0
  quadrupole_ion = 0.D0
  dipole = 0.D0
  quadrupole = 0.D0
  zvia = 0.D0
  rmag = 0.D0
  !
  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
  !
  ! get center of charge
  !
  ! if user defined center of charge use atoms from input
  !
  IF( (epcdft_surface_cm_start.gt.0) .and. (epcdft_surface_cm_end.gt.0) )THEN
    !
    zvtot = 0.D0
    x0 = 0.D0
    DO ia = epcdft_surface_cm_start, epcdft_surface_cm_end
       !
       zvtot = zvtot + zv(ityp(ia))
       !
       x0(:) = x0(:) + tau(:,ia)*zv(ityp(ia))
       !
    END DO
    !
    x0(:) = x0(:) / zvtot
    !
  ELSE ! not user defined use center of nucs
    !
    zvtot = 0.D0
    x0 = 0.D0
    DO ia = 1, nat
       !
       zvtot = zvtot + zv(ityp(ia))
       !
       x0(:) = x0(:) + tau(:,ia)*zv(ityp(ia))
       !
    END DO
    !
    x0(:) = x0(:) / zvtot
    !
  ENDIF
  !
  ! get total nuc charge
  !
  zvtot = 0.D0
  DO ia = 1, nat
     !
     zvtot = zvtot + zv(ityp(ia))
     !
  END DO
  !
  ! get electronic dipole
  !
  CALL compute_dipole( dfftp%nnr, nspin, rho%of_r, x0, e_dipole, e_quadrupole )
  !
  ! compute ionic+electronic total charge 
  !
  qq = -e_dipole(0) + zvtot
  !
  ! compute ion dipole moments
  !
  DO ia = 1, nat
     !
     zvia = zv(ityp(ia))
     !
     !zvtot = zvtot + zvia
     !
     DO ip = 1, 3
        !
        dipole_ion(ip) = dipole_ion(ip) + &
                         zvia*( tau(ip,ia) - x0(ip) )*alat
        quadrupole_ion(ip) = quadrupole_ion(ip) + &
                         zvia*( ( tau(ip,ia) - x0(ip) )*alat )**2
        !
     END DO
  END DO
  !
  ! compute ionic+electronic dipole and quadrupole moments
  !
  dipole(:)  = -e_dipole(1:3) + dipole_ion(:)
  !
  quadrupole = -e_quadrupole  + quadrupole_ion
  !
  ! convert center of charge to bohr
  !
  x0(:) = x0(:)*alat
  !
  IF(ionode)THEN
    !
    WRITE( stdout, &
           '(/5X,"reference position (x0):",5X,3F14.8," bohr")' ) x0(:)
    WRITE( stdout, '( 5X,"Total charge",3F9.4," au (Ha)")' )  qq 
    !
    WRITE( stdout, '( 5X,"Total dipole moment",3F9.4," au (Ha),", 3F9.4," Debye")' ) &
        ( dipole(ip),    ip = 1, 3), ( dipole(ip)*au_debye,    ip = 1, 3 )
    WRITE( stdout, * ) 
    !
  ENDIF
  !
  !---------------------
  !  Add potential
  !---------------------
  !
  ! Index for parallel summation
  !
  index0 = 0
#if defined (__MPI)
  !
  DO i = 1, me_bgrp
     index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  END DO
  !
#endif
  !
  ! ... Loop over position grid
  !
  DO ir = 1, dfftp%nnr
    !
    ! ... three dimensional indexes
    !
    i = index0 + ir - 1
    k = i / (dfftp%nr1x*dfftp%nr2x)
    i = i - (dfftp%nr1x*dfftp%nr2x)*k
    j = i / dfftp%nr1x
    i = i - dfftp%nr1x*j
    !
    ! calculate position
    !
    DO ip = 1, 3
      !
      r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
              DBLE( j )*inv_nr2*at(ip,2) + &
              DBLE( k )*inv_nr3*at(ip,3)
      !
    END DO
    !
    r = r * alat
    !
    ! shift image to center of charge
    !
    r = r - x0
    !
    ! shift to other side of xy plane
    !
    r(3) = r(3) + 2.D0*x0(3) 
    !
    rmag = DSQRT( r(1)**2 + r(2)**2 + r(3)**2 )
    !
    ! monopole
    !
    vin(ir,:) = vin(ir,:) + e2*qq / rmag
    !
    ! dipole
    !
    vin(ir,:) = vin(ir,:) + e2*DDOT(3, dipole, 1, r, 1) / rmag**3
    !
  END DO 
  !
END SUBROUTINE calc_epcdft_surface_field
!
!--------------------------------------------------------------------
SUBROUTINE print_epcdft_surface_energy_and_warning ( )
  !--------------------------------------------------------------------
  !
  ! print warning message about using surface method and the starting energy
  !
  USE kinds,     ONLY : DP
  USE io_global, ONLY : stdout, ionode
  USE lsda_mod,  ONLY : nspin
  USE fft_base,  ONLY : dfftp
  !
  IMPLICIT NONE
  !
  REAL(DP) field(dfftp%nnr, nspin)
  REAL(DP) energy
  !
  IF (ionode) THEN
    !
    WRITE( stdout,'(5x,"":)')
    WRITE( stdout,'(5x,"Adding image charge field due to neutral metal slab at origin in XY plane.":)')
    WRITE( stdout,'(5x,"Including monopole and dipole terms.":)')
    WRITE( stdout,'(5x,"System should not overlap with cell edges.":)')
    WRITE( stdout,'(5x,"":)')
    !
  ENDIF
  !
  field = 0.D0
  energy = 0.D0
  !
  CALL epcdft_surface_energy(field, energy)
  !
  IF(ionode)THEN
    WRITE(*,*)
    WRITE(*,'(5x,"First surface image interaction energy:",e10.3,"[Ry]")') energy
    WRITE(*,*)
  ENDIF
  !
  !
END SUBROUTINE print_epcdft_surface_energy_and_warning
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
  USE wvfct,                ONLY : npwx, g2kin
  USE klist,                ONLY : xk, nks, igk_k,ngk
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
  COMPLEX(DP) :: cutoff
  COMPLEX(DP) :: vbottot
  REAL(DP) :: dv
  CHARACTER(len=1024) :: filename
  !
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
  psic = 0.D0
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
  normfac=REAL(nelec,KIND=DP)/(REAL(vbottot,KIND=DP)*dv)
  vtop = normfac * vtop
  vbot = normfac * vbot
  !
  vtop(:,1) = vtop(:,1) / vbot
  vtop(:,2) = vtop(:,2) / vbot
  !
  DO ir = 1, n
    !
    IF (ABS(REAL(vbot(ir))).lt.REAL(cutoff)) THEN
      vtop(ir,1)=0.D0  
      vtop(ir,2)=0.D0  
    ENDIF
    !
    IF (vtop(ir,1) /= vtop(ir,1))  vtop(ir,1)=0.D0
    IF (vtop(ir,2) /= vtop(ir,2))  vtop(ir,2)=0.D0
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
