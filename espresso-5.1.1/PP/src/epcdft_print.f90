!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_print
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE epcdft_mod, ONLY : free1, free2, wmat, smat, hc, ohc
  USE klist, ONLY : nks
  USE io_global, ONLY : stdout, ionode
  !
  IMPLICIT NONE
  !
  LOGICAL :: epcdft_debug
  CHARACTER(LEN=18) title
  !
  epcdft_debug = .true.
  !
  IF(ionode) THEN
    !
    WRITE(stdout,1)
    WRITE(stdout,*)""
    !
    IF(epcdft_debug) THEN
      !
      WRITE(stdout,4) free1, free2
      WRITE(*,*)""
      !
      title='S up' ; CALL pcmat( title, smat(:,:,1) )
      title='S down' ; CALL pcmat( title, smat(:,:,2) )
      !
      title='W up' ; CALL pcmat( title, wmat(:,:,1) )
      title='W down' ; CALL pcmat( title, wmat(:,:,2) )
      !
      title='H up' ; CALL pcmat( title, hc(:,:,1) )
      title='H down' ; CALL pcmat( title, hc(:,:,2) )
      !
      title='Horth up' ; CALL pcmat( title, ohc(:,:,1) )
      title='Horth down' ; CALL pcmat( title, ohc(:,:,2) )
      !
    ENDIF
    !
    ! Sab = S_12_alpha * S_12_beta 
    WRITE(stdout,6) ABS( smat(1,2,1) * smat(1,2,2) )
    !
    ! Hab = H_12_alpha * S_12_beta + H_12_beta * S_12_alpha
    WRITE(stdout,5) ABS( ohc(1,2,1)*smat(1,2,2) + ohc(1,2,2)*smat(1,2,1) )
    !
    WRITE(stdout,*)""
    WRITE(stdout,1)
    !
  ENDIF
  !
  1 FORMAT('     =======================================================================')
  4 FORMAT(5x,'F1 = ',F12.5,3x,'F2 = ',F12.5)
  5 FORMAT(5x,'|Hab| = ',F14.6)
  6 FORMAT(5x,'|Sab| = ',F14.6)
  !
END SUBROUTINE epcdft_print
!
!-----------------------------------------------------------------------
SUBROUTINE pcmat (title, m)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=18), INTENT(IN) :: title
  COMPLEX(DP), INTENT(IN) :: m(2,2)
  !
  INTEGER :: i, j
  !
  WRITE(*,*)""
  WRITE(stdout,2) title
  WRITE(stdout,3) ( ( m(i,j) , j=1,2 ) , i=1,2)
  WRITE(*,*)""
  !
  2 FORMAT(5x,A18)
  !3 FORMAT(5x,'(',E12.4,',',E12.4,')',2x,'(',E12.4,',',E12.4,')')
  !3 FORMAT(5x,'(',F14.6,',',F14.6,')',2x,'(',F14.6,',',F14.6,')')
  3 FORMAT(5x,'(',e23.16,',',e23.16,')',2x,'(',e23.16,',',e23.16,')')
  !
END SUBROUTINE pcmat
