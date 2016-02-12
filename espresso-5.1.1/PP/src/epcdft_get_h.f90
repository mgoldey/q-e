!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_h
  !-----------------------------------------------------------------------
  !
  ! Here we calculate H and take the average of the off diags of H.
  !
  ! hc(1,1,s) = ecor1
  ! hc(1,2,s) = 0.5 ( (F1+F2) Sab - (Wab + Wba) )
  ! hc(2,1,s) = 0.5 ( (F2+F1) Sba - (Wba + Wab) )
  ! hc(2,2,s) = ecor2
  !
  ! Sab = Sab_up*Sab_down
  !
  ! Wab = Wab_up*det(Sab_down) + Wab_down*det(Sab_up)
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode, stdout
  USE epcdft_mod, ONLY : free1, free2, wmat, smat, hc, cor1, cor2
  USE klist,      ONLY : nks
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
  !
!  DO s = 1, nks
    DO i = 1, 2
      DO j = 1, 2
        !
        wtot(i,j) = ( wmat(i,j,1) + wmat(j,i,1) ) * smat(i,j,2) + &
                    ( wmat(i,j,2) + wmat(j,i,2) ) * smat(i,j,1)
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
  IF( ionode ) WRITE( stdout, * )"    H done"
  !
END SUBROUTINE epcdft_get_h
