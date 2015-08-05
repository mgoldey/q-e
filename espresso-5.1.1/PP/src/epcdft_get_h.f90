!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_h
  !-----------------------------------------------------------------------
  !
  ! hc(1,1,s) = ecor1
  ! hc(1,2,s) = 0.5 ( (F1+F2) Sab - Wab )
  ! hc(2,1,s) = 0.5 ( (F2+F1) Sba - Wba )
  ! hc(2,2,s) = ecor2
  !
  USE kinds, ONLY : DP
  USE epcdft_mod, ONLY : free1, free2, wmat, smat, hc, ohc
  USE klist, ONLY : nks
  !
  IMPLICIT NONE
  !
  INTEGER i, j, s
  REAL(DP) :: e(2,2)  
  REAL(DP) :: ecor1, ecor2
  !
  hc = 0.D0
  !
  WRITE(*,*)"    H done"
  !
END SUBROUTINE epcdft_get_h
