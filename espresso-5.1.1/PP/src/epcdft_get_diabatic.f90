!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_diabatic
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE epcdft_mod, ONLY : free1, free2, wmat, smat, hc, ohc
  USE klist, ONLY : nks
  !
  IMPLICIT NONE
  !
  INTEGER i, j, s
  REAL(DP) :: fmat(2)
  !
  hc = 0.D0
  ohc = 0.D0
  fmat(1) = free1
  fmat(2) = free2
  !
  DO s = 1 , nks
    DO i = 1 , 2
      DO j = 1 , 2
        hc(i,j,s) = 0.5D0 * ( ( fmat(i) + fmat(j) ) * smat(i,j,s) - ( wmat(i,j,1,s) + wmat(i,j,2,s) ) )
      ENDDO
    ENDDO
  ENDDO
  !
END SUBROUTINE epcdft_get_diabatic
