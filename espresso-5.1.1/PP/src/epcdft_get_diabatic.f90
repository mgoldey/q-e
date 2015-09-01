!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_diabatic
  !-----------------------------------------------------------------------
  !
  USE epcdft_mod, ONLY : eig_of_w
  !
  IMPLICIT NONE
  !
  IF(eig_of_w) THEN
    !
    CALL epcdft_get_diabatic_w()
    !
  ELSE
    !
    CALL epcdft_get_diabatic_lowdin()
    !
  ENDIF
  !
END SUBROUTINE epcdft_get_diabatic
!
!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_diabatic_w()
  !-----------------------------------------------------------------------
  !
  !  use eigen vecs of W following :
  !
  !    The Journal of Chemical Physics 133, 244105 (2010); doi: 10.1063/1.3507878
  !
  USE kinds, ONLY : DP
  USE epcdft_mod, ONLY : wmat, smat, ohc, hc
  !
  IMPLICIT NONE
  !
  INTEGER i, j, s
  COMPLEX(DP) :: stot(2,2), wtot(2,2) ! Sup*Sdown, W_up+W_down
  COMPLEX(DP) :: w_eigvec(2,2) 
  REAL(DP) :: w_eigval(2)
  !
  w_eigval = 0.D0
  ohc = (0.D0,0.D0)
  stot = (0.D0,0.D0)
  wtot = (0.D0,0.D0)
  w_eigvec = (0.D0,0.D0)
  !
  stot(:,:) = smat(:,:,1) * smat(:,:,2)
  !Not 100% sure this is how you combine up and down for W
  wtot(:,:) = wmat(:,:,1) * smat(:,:,2) + wmat(:,:,2) * smat(:,:,1) 
  !
  ! Solve W.V = S.V * eig
  !
  CALL cdiaghg( 2, 2, wtot, stot, 2, w_eigval, w_eigvec )
  !
  ! H_orthg =  v^(dagger).H.v
  !
  ohc = MATMUL( TRANSPOSE(CONJG(w_eigvec)) , hc )
  ohc = MATMUL( ohc , w_eigvec ) 
  !
END SUBROUTINE epcdft_get_diabatic_w
!
!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_diabatic_lowdin()
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE epcdft_mod, ONLY : smat, hc, ohc
  USE klist, ONLY : nks
  !
  IMPLICIT NONE
  !
  INTEGER i
  REAL(DP) :: l(2)            ! sev eigvals of S matrix
  COMPLEX(DP) :: u(2,2)       ! U  S eigenvector matrix sev(i,j) 
  COMPLEX(DP) :: d(2,2)       ! SD^-1/2 is diagonal S^-1/2
  COMPLEX(DP) :: invssqr(2,2) ! S-1/2
  COMPLEX(DP) :: smattot(2,2)   ! Sup*Sdown
  !
  ohc = 0.D0
  l = 0.D0
  u = 0.D0
  d = 0.D0
  invssqr = 0.D0
  smattot(:,:) = smat(:,:,1)*smat(:,:,2)
  !
  !
  ! get U and U^-1 and eigenvals (l) of S
  CALL cdiagh(2,smattot,2,l,u)
  !
  ! get diag S^-1/2
  DO i = 1, 2
    d(i,i) = l( i )**(-0.5D0)
  ENDDO
  !
  !
  ! get S^-1/2 = U^-1 . SD^-1/2 . U
  invssqr(:,:) = MATMUL( u(:,:), d(:,:) )
  invssqr(:,:) = MATMUL( invssqr(:,:), Transpose(u(:,:)) )
  !
  ! CALL print_cmat("S**-.5",invssqr,2)
  !
  ! get oHc = S^-1/2 . Hc . S^-1/2
  ohc(:,:) = MATMUL( invssqr(:,:), hc(:,:) )
  ohc(:,:) = MATMUL( ohc(:,:) , invssqr(:,:) )
  !
  !
END SUBROUTINE epcdft_get_diabatic_lowdin
!
