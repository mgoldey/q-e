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
  CHARACTER(LEN=18) title
  !
  IF(ionode) THEN
    !
    WRITE(stdout,1)
    WRITE(*,*)""
    !
    WRITE(stdout,4) free1, free2
    WRITE(*,*)""
    !
    title='S up' ; CALL pcmat( title, smat(:,:,1) )
    title='S down' ; CALL pcmat( title, smat(:,:,2) )
    !
    title='W1=v1*w1 up' ; CALL pcmat( title, wmat(:,:,1,1) )
    title='W1 down' ; CALL pcmat( title, wmat(:,:,1,2) )
    title='W2 up' ; CALL pcmat( title, wmat(:,:,2,1) )
    title='W2 down' ; CALL pcmat( title, wmat(:,:,2,2) )
    title='W1+W2 up' ; CALL pcmat( title, wmat(:,:,1,1)+wmat(:,:,2,1) )
    title='W1+W2 down' ; CALL pcmat( title, wmat(:,:,1,2)+wmat(:,:,2,2) )
    !
    title='H up' ; CALL pcmat( title, hc(:,:,1) )
    title='H down' ; CALL pcmat( title, hc(:,:,2) )
    !
    title='Horth up' ; CALL pcmat( title, ohc(:,:,1) )
    title='Horth down' ; CALL pcmat( title, ohc(:,:,2) )
    !
    WRITE(*,*)""
    WRITE(stdout,1)
    !
  ENDIF
  !
  1 FORMAT('     =======================================================================')
  4 FORMAT(5x,'F1 = ',F12.5,3x,'F2 = ',F12.5)
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
