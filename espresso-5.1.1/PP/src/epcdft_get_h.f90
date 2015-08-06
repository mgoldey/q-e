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
  USE epcdft_mod, ONLY : free1, free2, wmat, smat, hc, cor1, cor2
  USE klist, ONLY : nks
  !
  IMPLICIT NONE
  !
  INTEGER i, j, s
  REAL(DP) :: core(2) ! corrected energy
  REAL(DP) :: ftot ! F1+F2, Sup*Sdown, W_up+W_down
  COMPLEX(DP) :: stot(2,2), wtot(2,2) ! F1+F2, Sup*Sdown, W_up+W_down
  !
  hc = 0.D0
  core(1) = free1 + cor1 
  core(2) = free2 + cor2
  ftot = free1 + free2
  stot(:,:) = smat(:,:,1) * smat(:,:,2)
  wtot(:,:) = wmat(:,:,1) + wmat(:,:,2)
  !
!  DO s = 1, nks
    DO i = 1, 2
      DO j = 1, 2
        !
        IF(i==j)THEN
          hc(i,j) = core(i)
        ELSE
          hc(i,j) = 0.5D0 * ( ftot * stot(i,j) - wtot(i,j) )
        ENDIF
        !
      ENDDO
    ENDDO
!  ENDDO
  !
  WRITE(*,*)"    H done Note : assuming Sab=Sab_up*Sab_down and Wab = Wab_up + Wab_down"
  !
END SUBROUTINE epcdft_get_h
