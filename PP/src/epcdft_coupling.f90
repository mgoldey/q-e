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
!-----------------------------------------------------------------------
PROGRAM epcdft_coupling 
  !-----------------------------------------------------------------------
  ! 
  USE environment, ONLY : environment_start, environment_end
  USE mp_global, ONLY : mp_startup
  USE epcdft_mod, ONLY : debug
  !
  IMPLICIT NONE
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'epcdft_coupling' )
  !
  CALL epcdft_setup ( ) ! allocate and read vars
  !
  CALL epcdft_get_s ( ) ! create overlap matrix
  !
  CALL epcdft_get_w ( ) ! create weight matrix
  !
  CALL epcdft_get_h ( ) ! create hamiltonian
  !
  CALL epcdft_get_diabatic ( ) ! create orthogonal diabatic hamiltonian
  !
  CALL epcdft_print ( ) ! print results
  !
  CALL environment_end ( 'epcdft_coupling' )
  !
  CALL stop_pp
  !
END PROGRAM epcdft_coupling
