!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_print_energies()
  !----------------------------------------------------------------------------
  !
  ! This routine is used for printing energy contrib from plugins
  ! DO NOT REMOVE THE TAGS ! ***ADDSON_NAME KIND_OF_PATCH***
  !
  !   This routine calculates the correction to the total energy
  !   and prints it. It also calculates the total number of 
  !   electrons within the applied well. If this number is not
  !   equal to that of eopreg the amplitude of the well is changed
  !   and the scf loop is restarted. 
  !
  !   edir - atom to center potential well around 
  !   emaxpos - radius of potential well in alat
  !   eamp - strength of potential in Ry a.u.
  !   eopreg - number of electrons that should be in well
  !
  !
  !
  USE io_global,        ONLY : stdout, ionode
  USE kinds,            ONLY : DP
  USE io_files,         ONLY : tmp_dir
  !
  USE plugin_flags
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8, e2, au_debye
  USE ions_base,     ONLY : nat, ityp, zv
  USE cell_base,     ONLY : alat, at, omega, bg, saw
  USE extfield,      ONLY : tefield, dipfield, edir, eamp, emaxpos, &
                            eopreg, forcefield, etotefield
  USE epcdft,        ONLY : do_epcdft, fragment_atom1, &
                            fragment_atom2, epcdft_electrons, &
                            epcdft_amp, epcdft_width, epcdft_shift, &
                            epcdft_thr, hirshfeld
  USE force_mod,     ONLY : lforce
  USE io_global,     ONLY : stdout,ionode, ionode_id
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_images,     ONLY : intra_image_comm
  USE mp_bands,      ONLY : me_bgrp, intra_bgrp_comm
  USE fft_base,      ONLY : dfftp, grid_gather
  USE mp,            ONLY : mp_bcast, mp_sum
  USE control_flags, ONLY : iverbosity, conv_elec, conv_ions
  USE scf,           ONLY : rho
  !
  IMPLICIT NONE
  !
  ! local variables
  !
  INTEGER  :: i, is, iatom
  REAL(DP) :: tmp
  REAL(DP) :: dv
  REAL(DP) :: einwell                              ! number of electrons in well
  REAL(DP) :: enumerr                              ! epcdft_electrons - einwell  (e number error)
  LOGICAL  :: elocflag                             ! true if charge localization condition is satisfied
  LOGICAL  :: ZERO
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: vpotens   ! ef is added to this potential serial
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: vpotenp   ! ef is added to this potential parll
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rhosup    ! rho serial
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rhosdown  ! rho serial
  REAL(DP) :: next_epcdft_amp ! guess of amp for next iteration
  !
  ! SAVED VARS
  !
  REAL(DP) :: last_epcdft_amp = 0.D0 ! used to determine next guess for potential
  REAL(DP) :: last_einwell = 0.D0    ! used to determine next guess for potential
  LOGICAL :: first = .TRUE.          ! used to determine next guess for potential
  SAVE last_epcdft_amp
  SAVE last_einwell 
  SAVE first
  !
  ! dont do anything unless the calculation is converged
  !
  IF(.NOT.conv_elec)RETURN
  !
  ! calc is converged lets compute and print the correction
  !
  ! first setup vars
  ALLOCATE(vpotenp(dfftp%nnr, nspin))
  ALLOCATE(vpotens(dfftp%nr1x * dfftp%nr2x * dfftp%nr3x,nspin ))
  ALLOCATE(rhosup( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x ))
  ALLOCATE(rhosdown( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x ))
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  vpotenp  = 0.D0
  vpotens  = 0.D0
  rhosup   = 0.D0
  rhosdown = 0.D0
  einwell  = 0.D0
  tmp      = 0.D0
  elocflag = .TRUE.
  etotefield = 0.D0
  epcdft_shift = 0.D0
  zero =.false.
  next_epcdft_amp = 0.D0
  !
  ! gather the grids for serial calculation

  
  if (epcdft_amp .eq. 0d0) THEN
    zero=.true.
    epcdft_amp=1d0
  ENDIF
! this call only calulates vpoten
  CALL add_efield(vpotenp(:,1), epcdft_shift, rho%of_r, .true. )
  CALL add_efield(vpotenp(:,2), epcdft_shift, rho%of_r, .true. )
  if (zero) epcdft_amp=0d0

#ifdef __MPI
    CALL grid_gather ( vpotenp(:,1), vpotens(:,1) )
    CALL grid_gather ( rho%of_r(:,1), rhosup )
    IF(nspin > 1)THEN
      CALL grid_gather ( rho%of_r(:,2), rhosdown )
      CALL grid_gather ( vpotenp(:,2), vpotens(:,2) )
    ENDIF
#else
    vpotens(:,1)=vpotenp(:,1)
    rhosup(:) = rho%of_r(:,1)
    IF(nspin > 1)THEN
      rhosdown(:) = rho%of_r(:,2)
      vpotens(:,2)=vpotenp(:,2)
    ENDIF
#endif

  !
  ! begin calculation of the correction 
  IF(ionode) THEN
    !
    ! combine up and down parts of rho
    rhosup(:) = rhosup(:) + rhosdown(:) 
    ! vpotens(:,1)=vpotens(:,1)+vpotens(:,2)
    !
    DO i=1, dfftp%nr1x*dfftp%nr2x*dfftp%nr3x
      !
      ! calculate energy correction
      !
      epcdft_shift = epcdft_shift + (vpotens(i,1)) * rhosup(i) * dv
      !
      ! count number of electrons in well for localization condition check
      !
      IF(hirshfeld) THEN
        ! need number of electrons on 
        ! acceptor so negetive of vpotens is used
        IF (vpotens(i,1).lt.0d0) THEN
          einwell = einwell +  rhosup(i) * dv
        ELSE 
          einwell = einwell -  rhosup(i) * dv
        ENDIF
        !
      ELSE
        !
        IF(vpotens(i,1) .NE. 0.D0)THEN
          einwell = einwell + rhosup(i) * dv
        ELSE
          einwell = einwell - rhosup(i) * dv
        ENDIF
        !
      ENDIF
      !
    ENDDO
	DO iatom=1, nat
		IF (fragment_atom2.ne.0) THEN
			IF ((iatom.ge.fragment_atom1) .AND. (iatom.le.fragment_atom2) ) THEN
				einwell = einwell - zv(ityp(iatom))
			ELSE
				einwell = einwell + zv(ityp(iatom))
			ENDIF
		ELSE
			IF (iatom.eq.fragment_atom1) THEN
				einwell = einwell - zv(ityp(iatom))
			ELSE
				einwell = einwell + zv(ityp(iatom))
			ENDIF
		ENDIF
	ENDDO
    !
    ! the correction is - of the energy
    epcdft_shift = -1.D0 * epcdft_shift
    !
    if (.not.zero) WRITE(*,*)"    E field correction : ",epcdft_shift," Ry"
    WRITE(*,*)"    #e's   in well     : ",einwell," electrons"
    !
    ! is there a localization condition?
    !
    IF(epcdft_amp .NE. 0.D0)THEN
      !
      ! is the localization condition satisfied?
      !
      enumerr = epcdft_electrons - einwell 
        ! conv_ions = false will restart scf
      IF( ABS(enumerr) .GE. epcdft_thr .and. .not. zero) THEN
       conv_ions = .FALSE.
       WRITE(*,*) "    Surplus/deficit electrons    :  ",enumerr, "electrons"
       WRITE(*,*) "    epcdft_thr                   :  ",epcdft_thr, "electrons"
      ELSE 
      !  conv_ions = .TRUE.
      ENDIF
      !
    ENDIF
    !
  ENDIF
  ! CALL mp_bcast( enumerr, ionode_id, intra_bgrp_comm )

  if (zero) THEN
      ! IF(ionode) write(*,*) "All except for number of electrons is meaningless - EXITING NOW"
    conv_ions =.true.
  ENDIF
  !
  CALL mp_bcast( conv_ions, ionode_id, intra_image_comm )
  !
  DEALLOCATE(vpotenp)
  DEALLOCATE(vpotens)
  DEALLOCATE(rhosup)
  DEALLOCATE(rhosdown)
  !
  ! if the charge is not localized
  !
  IF(.NOT.conv_ions .AND. epcdft_amp .NE. 0)THEN
    !
    ! update applied field and restart scf
    !
    IF(ionode)THEN
      !
      WRITE(*,*)""
      WRITE(*,*)"    -----------------------------------------------"
      WRITE(*,*)""
      WRITE(*,*)"                     SCF converged but..."
      WRITE(*,*)""
      WRITE(*,*)"    electron localization condition NOT satisfied"
      WRITE(*,*)"    restarting scf with different applied potential"
      WRITE(*,*)""
      !
      ! find next guess for epcdft_amp
      !
      IF(first) THEN
         !
         first = .FALSE.
         !
         next_epcdft_amp = epcdft_amp - SIGN(0.001D0, enumerr) 
         !
      ELSE
         !
         CALL secant_method(next_epcdft_amp, epcdft_amp,   last_epcdft_amp, &
                            einwell,         last_einwell, epcdft_electrons)
         if (abs(next_epcdft_amp) .gt. abs(epcdft_amp)*1.1) THEN
         	next_epcdft_amp=epcdft_amp*1.1
         ENDIF
         !
      ENDIF
      !
      ! save this iteration's einwell and amp
      ! for the next iteration
      last_einwell    = einwell
      last_epcdft_amp = epcdft_amp
      !
      ! The old iteration has passed away;
      ! behold, the new iteration has come 
      epcdft_amp = next_epcdft_amp
      !
    ENDIF
    !
    CALL mp_bcast( epcdft_amp, ionode_id, intra_image_comm ) ! what is best bcast this or enumerr?
    !
    IF(ionode) WRITE(*,*)"    New field Amp      : ",epcdft_amp," Ry"
    !
    ! epcdft_shift = 0.D0 ! this var is added to etot before 
    !                   ! this routine is called during next scf loop
  ENDIF
  !
END SUBROUTINE plugin_print_energies
!----------------------------------------------------------------------------
!
!
SUBROUTINE secant_method(vnext, v, vold, e, eold, egoal)
  !----------------------------------------------------------------------------
  !
  ! This routine uses secant method to determine the next 
  ! estimation of the root of an equation.
  !
  ! Numerical Mathematics and Computing Sixth Edition
  !  Ward Cheney & David Kincaid Pg 112, eqn. 3
  !
  !     vnext = v - \left( \frac{v-vold}{e-eold} \right) * (e-egoal)
  !
  ! In this case the function we want to minmize is:
  !          e*(v) = e(v) - egoal
  ! this is why the formual has egoal in it
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
!----------------------------------------------------------------------------
