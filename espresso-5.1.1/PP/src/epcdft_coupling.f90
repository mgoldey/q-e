! 
!-----------------------------------------------------------------------
PROGRAM epcdft_coupling 
  !-----------------------------------------------------------------------
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
  IF(debug) CALL epcdft_check_w ( ) ! check the <A/B|w|A/B> = c1/c2
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
!  RETURN
  !
END PROGRAM epcdft_coupling
