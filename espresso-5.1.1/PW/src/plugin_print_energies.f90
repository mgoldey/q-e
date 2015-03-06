!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_print_energies(etotefield)
  !----------------------------------------------------------------------------
  !
  ! This routine is used for printing energy contrib from plugins
  ! DO NOT REMOVE THE TAGS ! ***ADDSON_NAME KIND_OF_PATCH***
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
                            eopreg, forcefield
  USE force_mod,     ONLY : lforce
  USE io_global,     ONLY : stdout,ionode
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_images,     ONLY : intra_image_comm
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dfftp, grid_gather
  USE mp,            ONLY : mp_bcast, mp_sum
  USE control_flags, ONLY : iverbosity, conv_elec
  USE scf,           ONLY : rho
  !
  IMPLICIT NONE
  !
  REAL(DP),INTENT(INOUT) :: etotefield       ! contribution to etot due to ef
  !
  ! local variables
  !
  INTEGER :: i, is 
  REAL(DP) :: tmp
  REAL(DP) :: dv
  REAL(DP), DIMENSION(:), ALLOCATABLE :: vpotens ! ef is added to this potential serial
  REAL(DP), DIMENSION(:), ALLOCATABLE :: vpotenp ! ef is added to this potential parll
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
  ALLOCATE(vpotenp(dfftp%nnr))
  ALLOCATE(vpotens( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x ))
  ALLOCATE(rhosup( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x ))
  ALLOCATE(rhosdown( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x ))
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  vpotenp = 0.D0
  vpotens = 0.D0
  rhosup = 0.D0
  rhosdown = 0.D0
  !
  !
  ! this call only calulates vpoten
  CALL add_efield(vpotenp, etotefield, rho%of_r, .true. )
  !
  ! gather the grids for serial calculation
#ifdef __MPI
    CALL grid_gather ( vpotenp, vpotens )
    CALL grid_gather ( rho%of_r(:,1), rhosup )
    IF(nspin > 1)THEN
      CALL grid_gather ( rho%of_r(:,2), rhosdown )
    ENDIF
#else
    vpotens(:)=vpotenp(:)
    rhosup(:) = rho%of_r(:,1)
    IF(nspin > 1)THEN
      rhosdown(:) = rho%of_r(:,2)
    ENDIF
#endif
  !
  ! begin calculation of the correction 
  IF(ionode) THEN
    !
    ! combine up and down parts of rho
    rhosup(:) = rhosup(:) + rhosdown(:) 
    !
    etotefield = 0.D0
    DO i=1, dfftp%nr1x*dfftp%nr2x*dfftp%nr3x
      !
      etotefield = etotefield + vpotens(i) * rhosup(i) * dv
      !
    ENDDO
    !
    ! the correction is - of the energy
    etotefield = -1.D0 * etotefield
    !
    WRITE(*,*)"    E field correction : ",etotefield," Ry"
    !
  ENDIF
  !
  DEALLOCATE(vpotenp)
  DEALLOCATE(vpotens)
  DEALLOCATE(rhosup)
  DEALLOCATE(rhosdown)
  !
END SUBROUTINE plugin_print_energies
