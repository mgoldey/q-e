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
  USE fft_base,      ONLY : dfftp
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
  REAL(DP), DIMENSION(:), ALLOCATABLE :: vpoten ! ef is added to this potential
  !
  ! dont do anything unless the calculation is converged
  !
  IF(.NOT.conv_elec)RETURN
  !
  ! calc is converged lets compute and print the correction
  !
  ALLOCATE(vpoten(dfftp%nnr))
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  vpoten = 0.D0
  !
  !
  ! this call only calulates vpoten
  CALL add_efield(vpoten, etotefield, rho%of_r, .true. )
  !
  etotefield = 0.D0
  DO is=1, nspin
    !
    DO i=1, dfftp%nnr
      !
      etotefield = etotefield + vpoten(i) * rho%of_r(i,is) * dv
      !
    ENDDO
    !
  ENDDO
  !
  WRITE(*,*)"    E field correction : ",etotefield," Ry"
  !
  DEALLOCATE(vpoten)
  !
END SUBROUTINE plugin_print_energies
