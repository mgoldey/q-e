!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_diabatic
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
END SUBROUTINE epcdft_get_diabatic
!