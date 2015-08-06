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
  INTEGER i, s
  REAL(DP) :: l(2)            ! sev eigvals of S matrix
  COMPLEX(DP) :: u(2,2)       ! U  S eigenvector matrix sev(i,j) 
  COMPLEX(DP) :: invu(2,2)    !U^-1 inverse of S eigenvector matrix invsev(i,j) 
  COMPLEX(DP) :: d(2,2)       ! SD^-1/2 is diagonal S^-1/2
  COMPLEX(DP) :: invssqr(2,2) ! S-1/2
  COMPLEX(DP) :: smattot(2,2)   ! Sup*Sdown
  !
  ohc = 0.D0
  l = 0.D0
  u = 0.D0
  invu = 0.D0
  d = 0.D0
  invssqr = 0.D0
  smattot(:,:) = smat(:,:,1)*smat(:,:,2)
  !
!  DO s = 1, nks
    !
    ! get U and U^-1 and eigenvals (l) of S
    CALL get_evs(smattot(:,:),2,u(:,:),l)
    CALL get_inv(u(:,:),invu(:,:))
    !
    ! get diag S^-1/2
    DO i = 1, 2
      d(i,i) = l( 2/i )**(-0.5D0) !eigenvals need to be switched thus 2/i
    ENDDO
    !
    !
    ! get S^-1/2 = U^-1 . SD^-1/2 . U
    invssqr(:,:) = MATMUL( invu(:,:), d(:,:) )
    invssqr(:,:) = MATMUL( invssqr(:,:), u(:,:) )
    !
    ! get oHc = S^-1/2 . Hc . S^-1/2
    ohc(:,:) = MATMUL( invssqr(:,:), hc(:,:) )
    ohc(:,:) = MATMUL( ohc(:,:) , invssqr(:,:) )
    !
!  ENDDO
  !
END SUBROUTINE epcdft_get_diabatic
!
!-----------------------------------------------------------------------
SUBROUTINE get_evs(a,n,z,w)
  !---------------------------------------------------------------------
  !
  ! This ugly bit of code is from the lapack software.intel.com
  ! given square complex Hermitian matrix a(n,n) it will give the eigenvectors
  ! in matrix z, eigenvalues are real and in w
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: N
  INTEGER          LDA, LDZ
  INTEGER          LWMAX
  PARAMETER        ( LWMAX = 1000 )
  INTEGER          INFO, LWORK, IL, IU, M
  DOUBLE PRECISION ABSTOL, VL, VU
  INTEGER          IWORK( 5*N ), IFAIL( N )
  DOUBLE PRECISION RWORK( 7*N )
  COMPLEX*16, INTENT(IN)::  A( N,N )
  COMPLEX*16, INTENT(INOUT)::  Z( N,N )
  DOUBLE PRECISION,INTENT(INOUT):: W( N )
  COMPLEX*16 ::  WORK( LWMAX )
  EXTERNAL         ZHEEVX
  EXTERNAL         PRINT_MATRIX, PRINT_RMATRIX
  INTRINSIC        INT, MIN
  !
  LDA = N
  LDZ = N 
  ABSTOL = -1.0
  VL = 0.0
  VU = 100.0
  LWORK = -1
  CALL ZHEEVX( 'Vectors', 'Values', 'Lower', N, A, LDA, VL, VU, IL,&
              IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, &
              IFAIL, INFO )                                        
  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )                            
  CALL ZHEEVX( 'Vectors', 'Values', 'Lower', N, A, LDA, VL, VU, IL, &
              IU, ABSTOL, M, W, Z, LDZ, WORK, LWORK, RWORK, IWORK, &
              IFAIL, INFO )                                        
  IF( INFO.GT.0 ) THEN
     WRITE(*,*)'The algorithm failed to compute eigenvalues.'
     STOP
  END IF
END SUBROUTINE get_evs
!-----------------------------------------------------------------------
SUBROUTINE get_inv(ain,bout)
  !---------------------------------------------------------------------
  ! gives inverse of 2x2 matrix ain
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: ain(2,2)
  COMPLEX(DP), INTENT(INOUT) :: bout(2,2)
  COMPLEX(DP) :: a, b, c, d
  !
  a = ain(1,1)
  b = ain(2,1)
  c = ain(1,2)
  d = ain(2,2)
  !
  bout(1,1) = d/(-b*c + a*d)
  bout(2,1) = -(b/(-b*c + a*d))
  bout(1,2) = -(c/(-b*c + a*d))
  bout(2,2) = a/(-b*c + a*d)
  !
END SUBROUTINE get_inv
