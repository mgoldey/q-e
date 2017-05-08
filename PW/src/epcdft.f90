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
!----------------------------------------------------------------------------
SUBROUTINE epcdft_controller()
  !----------------------------------------------------------------------------
  !
  !   CDFT
  !
  !   This routine calculates the correction to the total energy
  !   due to the constraining potential from cdft and prints it. 
  !   It also evaluates the charge difference between the donor
  !   and acceptor. If this number is not equal to that of the 
  !   target charge within epcdft_tol, the amplitude of the well
  !   is changed and the scf loop is restarted. 
  !
  !
  USE io_global,     ONLY : stdout, ionode
  USE kinds,         ONLY : DP
  USE plugin_flags
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8, e2, au_debye
  USE ions_base,     ONLY : nat, ityp, zv
  USE cell_base,     ONLY : alat, at, omega, bg, saw
  USE extfield,      ONLY : tefield, dipfield, edir, eamp, emaxpos, &
                            eopreg, forcefield
  USE force_mod,     ONLY : lforce
  USE io_global,     ONLY : stdout,ionode, ionode_id
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_images,     ONLY : intra_image_comm
  USE mp_bands,      ONLY : me_bgrp, intra_bgrp_comm
  USE mp_pools,      ONLY : inter_pool_comm
  USE mp_world,      ONLY : world_comm

  USE fft_base,      ONLY : dfftp
  USE mp,            ONLY : mp_bcast, mp_sum
  USE control_flags, ONLY : iverbosity, conv_elec
  USE scf,           ONLY : rho
  USE io_files,      ONLY : tmp_dir, prefix
  USE epcdft,        ONLY : do_epcdft, conv_epcdft,epcdft_shift, epcdft_locs,&
                            epcdft_guess,nconstr_epcdft,epcdft_tol,epcdft_type,&
                            epcdft_target, epcdft_delta_fld,epcdft_update_intrvl,&
                            reset_field, epcdft_surface_shift, epcdft_surface,&
                            epcdft_surface_tol
  !
  IMPLICIT NONE
  !
  ! local variables
  !
  INTEGER  :: i, is, iatom,iconstraint
  INTEGER  :: ictr=0
  LOGICAL  :: first =.true.
  REAL(DP) :: tmp
  REAL(DP) :: dv
  REAL(DP) :: acharge, dcharge                     ! acceptor/donor charge 
  REAL(DP) :: achargep, dchargep                   ! acceptor/donor charge parallel
  REAL(DP) :: einwell                              ! number of electrons in well
  REAL(DP) :: einwellp                             ! number of electrons in well
  REAL(DP) :: epcdft_shiftp                        ! size of shift
  REAL(DP) :: enumerr                              ! epcdft_charge - einwell  (e number error)
  LOGICAL  :: elocflag                             ! true if charge localization condition is satisfied
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: vpotenp ! ef is added to this potential
  REAL(DP) :: next_epcdft_amp,epcdft_amp           ! guess of amp for next iteration
  CHARACTER(LEN=256) :: filename                   ! cube file for final cdft potential
  REAL(DP) :: surface_echange                      !change in image interaction energy for surface
  !
  ! SAVED VARS
  !
  REAL(DP), DIMENSION(:), allocatable,save :: last_epcdft_amp(:) ! used to determine next guess for potential
  REAL(DP), DIMENSION(:), allocatable,save :: last_einwell(:)    ! used to determine next guess for potential
  SAVE ictr,first
  !
  IF (ictr.eq.0) THEN
    !
    allocate(last_epcdft_amp(nconstr_epcdft))
    allocate(last_einwell(nconstr_epcdft))
    last_epcdft_amp=0.D0
    last_einwell=0.D0
    !
  ENDIF
  !
  reset_field=.false.
  !
  ! Don't update potential too often or CDFT will never converge
  !
  ictr=ictr+1
  !
  !
  !
  IF ((.NOT. conv_elec) .AND. (mod( ictr, epcdft_update_intrvl ) .NE. 0) ) RETURN
  !
  ! restart the counter but dont start at 0 or u will need to reallocate
  !
  ictr = 1
  !
  ! calc is converged lets compute and print the correction
  !
  ! first setup vars
  !
  ALLOCATE(vpotenp(dfftp%nnr, nspin))
  !
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  tmp      = 0.D0
  epcdft_shift = 0.D0
  epcdft_shiftp = 0.D0
  !
  conv_epcdft=.true. ! set to false whenever anything is not satisfied
  !
  ! calculate energy due to image charges
  ! if energy is not converged within epcdft_tol  
  ! set conv_epcdft = false
  !
  IF(epcdft_surface)THEN
    !
    vpotenp = 0.D0
    tmp = 0.D0
    !
    CALL epcdft_surface_energy(vpotenp, tmp)
    !
    surface_echange = tmp - epcdft_surface_shift
    !
    IF( abs(surface_echange) > epcdft_surface_tol)THEN
      epcdft_surface_shift = tmp
      conv_epcdft = .false.
    ENDIF
    !
  ENDIF
  !
  ! moving onto epcdft constraint fields
  !
  DO iconstraint=1,nconstr_epcdft
    !
    acharge = 0.D0
    dcharge = 0.D0
    achargep = 0.D0
    dchargep = 0.D0
    elocflag=.true.
    einwell  = 0.D0
    einwellp  = 0.D0
    epcdft_shiftp=0.D0
    next_epcdft_amp = 0.D0
    vpotenp=0.D0
    !
    epcdft_amp=epcdft_guess(iconstraint)
    !
    ! calculate hirshfeld potential
    !
    CALL calc_hirshfeld_v(vpotenp,iconstraint)
    !
    ! begin calculation of the correction 
    !
    DO i=1, dfftp%nnr
      !
      ! calculate energy correction
      !
      IF (nspin.eq.1) THEN
        epcdft_shiftp = epcdft_shiftp + epcdft_amp * vpotenp(i,1) * rho%of_r(i,1) * dv
      ELSE
        epcdft_shiftp = epcdft_shiftp + epcdft_amp * vpotenp(i,1) * rho%of_r(i,1) * dv &
                                      + epcdft_amp * vpotenp(i,2) * rho%of_r(i,2) * dv
      ENDIF
      !
      ! count number of electrons in well - this must depend on the well
      !
      SELECT CASE( epcdft_type(iconstraint) )
      CASE('charge','delta_charge')
        !
        einwellp = einwellp - vpotenp(i,1) * rho%of_r(i,1) * dv - vpotenp(i,2) * rho%of_r(i,2) * dv
        !
        IF(vpotenp(i,1)<0.D0) THEN
          !
          achargep = achargep - ABS( vpotenp(i,1)  * rho%of_r(i,1) +vpotenp(i,2) * rho%of_r(i,2)) * dv
          !
        ELSE IF(vpotenp(i,1)>0.D0) THEN
          !
          dchargep = dchargep - ABS( vpotenp(i,1)  * rho%of_r(i,1) +vpotenp(i,2) * rho%of_r(i,2)) * dv
          !
        ENDIF
        !
      CASE('spin','delta_spin')
        !
        einwellp = einwellp - vpotenp(i,1) * rho%of_r(i,1) * dv - vpotenp(i,2) * rho%of_r(i,2) * dv
        !
        IF(vpotenp(i,1)<0.D0) THEN
          !
          achargep = achargep + ABS( vpotenp(i,1)  * rho%of_r(i,1) + vpotenp(i,2) * rho%of_r(i,2)) * dv
          !
        ELSE IF(vpotenp(i,1)>0.D0) THEN
          !
          dchargep = dchargep + ABS( vpotenp(i,1)  * rho%of_r(i,1) + vpotenp(i,2) * rho%of_r(i,2)) * dv
          !
        ENDIF
        !
      CASE('delta_alpha')
        !
        einwellp = einwellp - vpotenp(i,1) * rho%of_r(i,1) * dv  
        !
        IF(vpotenp(i,1)<0.D0) THEN
          !
          achargep = achargep + ABS( vpotenp(i,1)  * rho%of_r(i,1)) * dv
          !
        ELSE IF(vpotenp(i,1)>0.D0) THEN
          !
          dchargep = dchargep + ABS( vpotenp(i,1)  * rho%of_r(i,1)) * dv
          !
        ENDIF
        !
      CASE('delta_beta')
        !
        einwellp = einwellp - vpotenp(i,2) * rho%of_r(i,2) * dv  
        !
        IF(vpotenp(i,2)<0.D0) THEN
          !
          achargep = achargep + ABS( vpotenp(i,2)  * rho%of_r(i,2)) * dv
          !
        ELSE IF(vpotenp(i,2)>0.D0) THEN
          !
          dchargep = dchargep + ABS( vpotenp(i,2)  * rho%of_r(i,2)) * dv
          !
        ENDIF
        !
      END SELECT
      !
    END DO ! i over nnr
    !
#ifdef __MPI
    !
    CALL MP_SUM(epcdft_shiftp,intra_image_comm) 
    CALL MP_SUM(einwellp,intra_image_comm) 
    epcdft_shift=epcdft_shiftp
    einwell=einwellp
    !
    CALL MP_SUM(achargep,intra_image_comm) 
    acharge=achargep
    CALL MP_SUM(dchargep,intra_image_comm)
    dcharge=dchargep
    !
#else
    !
    epcdft_shift=epcdft_shiftp
    einwell=einwellp
    !
    acharge=achargep
    dcharge=dchargep
    !
#endif
    !
    ! count nuc charge
    !
    SELECT CASE( epcdft_type(iconstraint) )
      CASE('charge','delta_charge')
        !
        DO iatom=1, nat
          !
          IF ((iatom.ge.epcdft_locs(1,iconstraint)) .AND. (iatom.le.epcdft_locs(2,iconstraint)) ) THEN
            !
            einwell = einwell - zv(ityp(iatom))
            acharge = acharge + zv(ityp(iatom)) 
            !
          ELSE IF ((iatom.ge.epcdft_locs(3,iconstraint)) .AND. (iatom.le.epcdft_locs(4,iconstraint)) ) THEN
            !
            einwell = einwell + zv(ityp(iatom))
            dcharge = dcharge + zv(ityp(iatom))
            !
          ENDIF
          !
        ENDDO ! atoms
        !
    END SELECT
    !
    ! is the localization condition satisfied?
    !
    enumerr = epcdft_target(iconstraint) - einwell 
    !
    !
    IF(ABS(enumerr) .GE. epcdft_tol) THEN
      elocflag=.false.
    ELSE 
      elocflag=.true.
    ENDIF
    !
    IF(.NOT. elocflag .or. .not. conv_elec) THEN
      !
      conv_epcdft=.FALSE.
      reset_field=.true.
      !
      ! update applied field and restart scf
      !
      IF(ionode) THEN
        !
        ! find next guess for epcdft_amp
        !
        IF(first) THEN
           !
           next_epcdft_amp = epcdft_amp + SIGN(MIN(0.001D0,epcdft_delta_fld), enumerr)*SIGN(1.0D0,einwell)
           first=.false.
           !
        ELSE
           !
           CALL secant_method(next_epcdft_amp, epcdft_amp, last_epcdft_amp(iconstraint), &
                              einwell, last_einwell(iconstraint), epcdft_target(iconstraint))
           !
           ! abs of the change in amp must be <= |delta_fld|
           !
           IF ( ABS(next_epcdft_amp - epcdft_amp) .gt. ABS(epcdft_delta_fld) ) THEN
             !
             next_epcdft_amp = epcdft_amp + SIGN(1.D0, next_epcdft_amp - epcdft_amp) * epcdft_delta_fld
             !
           ENDIF
           !
        ENDIF
        !
        ! save this iteration's einwell and amp
        ! for the next iteration
        !
        last_einwell(iconstraint)    = einwell
        last_epcdft_amp(iconstraint) = epcdft_amp
        !
        epcdft_guess(iconstraint) = next_epcdft_amp
        !
      ENDIF !ionode
      !
      CALL mp_bcast( epcdft_guess(iconstraint), ionode_id, intra_image_comm )
      !
    ENDIF !update lambda
    !
    IF(ionode) THEN
      !
      IF (.not. elocflag .or. .not. conv_elec) THEN
      !
      ! CDFT
      !
      WRITE(*,'(5x,a3,a10,a2,a10,a12,a12,a12,a12)') "I","D",'  ','A','Val','Target','Old','New'
      WRITE(*,'(5x,i3,e10.3,a2,e10.3,e12.3,e12.3,e12.3,e17.8)') &
      !
      iconstraint, dcharge,'  ',acharge,einwell,epcdft_target(iconstraint), &
      last_epcdft_amp(iconstraint),epcdft_guess(iconstraint)
      !
      WRITE(*,'(5x,a,2x,e17.8)') 'CDFT charge error:', enumerr
      !
      ! surface
      !
      IF(epcdft_surface)THEN
        !
        WRITE(*,'(5x,"Image interaction energy: ",e15.8," [Ry] Change: ",e10.3," [Ry] Tolerance: ",e10.3," [Ry]")')&
        epcdft_surface_shift, surface_echange, epcdft_surface_tol
        !
      ENDIF
      !
      ELSE
        !
        ! CDFT
        !
        WRITE(*,'(5x,a3,a10,a2,a10,a12,a12,a15)') "I","D",'  ','A','Val','Target','Str'
        WRITE(*,'(5x,i3,e10.3,a2,e10.3,e12.3,e12.3,e15.6)') &
        iconstraint, dcharge,'  ',acharge,einwell,epcdft_target(iconstraint),epcdft_guess(iconstraint)
        WRITE(*,'(5x,a,2x,e17.8)') 'CDFT charge error:', enumerr
        !
        ! surface
        !
        IF(epcdft_surface)THEN
          !
          WRITE(*,'(5x,"Image interaction energy: ",e15.8," [Ry] Change: ",e10.3," [Ry] Tolerance: ",e10.3," [Ry]")')&
          epcdft_surface_shift, surface_echange, epcdft_surface_tol
          !
        ENDIF
        !
      ENDIF
    ENDIF
    !
  END DO ! iconstraint
  !
  DEALLOCATE(vpotenp)
  !
  ! the correction is - of the energy
  !
  epcdft_shift = -1.D0 * epcdft_shift
  !
END SUBROUTINE epcdft_controller
!
SUBROUTINE secant_method(vnext, v, vold, e, eold, egoal)
  !----------------------------------------------------------------------------
  !
  ! Secant method to determine the next 
  ! estimation of the root of an equation.
  !
  !     vnext = v - \left( \frac{v-vold}{e-eold} \right) * (e-egoal)
  !
  USE kinds,            ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: vnext
  REAL(DP), INTENT(IN) :: v, vold, e, eold, egoal
  !
  vnext = v - ( (v-vold) / (e-eold)  ) * (e-egoal)
  !
END SUBROUTINE secant_method
!
!----------------------------------------------------------------------------
SUBROUTINE epcdft_surface_energy(vin, eout)
  !----------------------------------------------------------------------------
  !
  ! Calculate energy due to image charge
  !
  USE kinds,         ONLY : DP
  USE lsda_mod,      ONLY : nspin
  USE scf,           ONLY : rho
  USE fft_base,      ONLY : dfftp
  USE mp_images,     ONLY : intra_image_comm
  USE mp,            ONLY : mp_sum
  USE ions_base,     ONLY : nat, ityp, zv, tau
  USE cell_base,     ONLY : omega, alat
  USE constants,     ONLY : e2
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: vin(dfftp%nnr, nspin) ! surface field
  REAL(DP), INTENT(INOUT) :: eout ! energy due to image charges
  !
  REAL(DP) :: eoutp
  INTEGER :: i, is, ia
  REAL(DP) :: dv
  REAL(DP) :: x0(3) ! center of charge of system
  REAL(DP) :: qq ! total charge
  REAL(DP) :: dipole(3)!, quadrupole(3) ! total dips
  REAL(DP) :: r(3), rmag, mono_e, dip_e, pdotr
  !
  vin = 0.D0
  eout = 0.D0
  eoutp = 0.D0
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  x0 = 0.D0
  qq = 0.D0
  dipole = 0.D0
  !
  ! calc field due to images
  !
  CALL calc_epcdft_surface_field(vin, x0, qq, dipole)
  !
  ! calculate electronic energy correction
  !
  DO is=1, nspin
    DO i=1, dfftp%nnr
      !
      eoutp = eoutp + vin(i,is) * rho%of_r(i,is)
      !
    ENDDO!grid
  ENDDO!spin
  !
  eoutp = eoutp*dv
  !
#ifdef __MPI
    CALL mp_sum(eoutp, intra_image_comm) 
    eout=eoutp
#else
    eout=eoutp
#endif
  !
  !
  ! ionic part of energy
  !
  ! r = x0-R
  !
  ! mono_e = Z*qq / |r|
  !
  ! dip_e = Z*(p.r) / |r|^3
  !
  ! convert to image charge
  !
  x0(3) = x0(3) - 2.D0*x0(3)
  qq = -1.D0*qq
  dipole(:) = -1.D0*dipole(:)
  !
  DO ia=1, nat
    !
    r(:) = tau(:,ia)*alat - x0(:)
    !
    rmag = DSQRT( r(1)**2 + r(2)**2 + r(3)**2 )
    !
    ! mono part
    !
    mono_e = zv(ityp(ia)) * e2*qq / rmag
    !
    ! dip part
    !
    pdotr = dipole(1)*r(1) + dipole(2)*r(2) + dipole(3)*r(3)
    !
    dip_e = zv(ityp(ia)) * e2*pdotr / rmag**3
    !
    eout = eout + mono_e + dip_e
    !
  ENDDO ! atoms
  !
  eout = eout / 2.D0
  !
END SUBROUTINE epcdft_surface_energy
