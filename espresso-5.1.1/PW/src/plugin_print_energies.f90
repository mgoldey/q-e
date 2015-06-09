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
                            epcdft_amp, epcdft_width, epcdft_shift
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
  INTEGER  :: i, is 
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
  !
  ! gather the grids for serial calculation

  
  if (epcdft_amp .eq. 0d0) THEN
    zero=.true.
    epcdft_amp=1d0
  ENDIF
! this call only calulates vpoten
  CALL add_efield(vpotenp, epcdft_shift, rho%of_r, .true. )
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
      IF(vpotens(i,1) .NE. 0.D0)THEN
        !
        einwell = einwell + rhosup(i) * dv
        !
      ENDIF
      !
    ENDDO
    !
    ! the correction is - of the energy
    epcdft_shift = -1.D0 * epcdft_shift
    !
    WRITE(*,*)"    E field correction : ",epcdft_shift," Ry"
    WRITE(*,*)"    #e's   in well     : ",einwell," electrons"
    !
    ! is there a localization condition?
    !
    IF(epcdft_electrons .NE. 0.D0)THEN
      !
      ! is the localization condition satisfied?
      !
      enumerr = epcdft_electrons - einwell 
        ! conv_ions = false will restart scf
      IF( ABS(enumerr) .GE. 1.D-2 ) THEN
       conv_ions = .FALSE.
       WRITE(*,*) "    Numerical error    :  ",enumerr, "electrons"
      ELSE 
      !  conv_ions = .TRUE.
      ENDIF
      !
    ENDIF
    !
  ENDIF
  ! CALL mp_bcast( enumerr, ionode_id, intra_bgrp_comm )

  if (zero) THEN
    write(*,*) "All except for number of electrons is meaningless - EXITING NOW"
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
  IF(.NOT.conv_ions .AND. epcdft_electrons .NE. 0)THEN
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
    ENDIF
    !
    epcdft_amp = epcdft_amp - enumerr * ABS(epcdft_amp) 
    CALL mp_bcast( epcdft_amp, ionode_id, intra_image_comm )

    IF(ionode) WRITE(*,*)"    New field Amp      : ",epcdft_amp," Ry"
    !
    ! epcdft_shift = 0.D0 ! this var is added to etot before 
    !                   ! this routine is called during next scf loop
  ENDIF
  !
END SUBROUTINE plugin_print_energies
