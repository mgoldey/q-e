!
! Copyright (C) 2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
SUBROUTINE plugin_print_energies(etotefield,rho)
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
  USE control_flags, ONLY : iverbosity
!
!
! ***Environ MODULES BEGIN***
! ***Environ MODULES END***
!
  IMPLICIT NONE
  !
  REAL(DP),INTENT(INOUT) :: etotefield       ! contribution to etot due to ef
  REAL(DP),INTENT(IN)    :: rho(dfftp%nnr,nspin) ! the density whose dipole is computed
  !
  ! local variables
  !
  INTEGER :: i, is 
  REAL(DP) :: tmp
  REAL(DP), DIMENSION(:), ALLOCATABLE :: vpoten ! ef is added to this potential
  !
  !
  ALLOCATE(vpoten(dfftp%nnr))
  !
  !
  ! this call only calulates vpoten
  CALL add_efield(vpoten, etotefield, rho, .true. )
  !
  etotefield = 0.D0
  DO is=1, nspin
    DO i=1, dfftp%nnr
      etotefield = etotefield + vpoten(i) * rho(i,is)
    ENDDO
  ENDDO
  !
  DEALLOCATE(vpoten)
  !
END SUBROUTINE plugin_print_energies
