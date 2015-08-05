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
  INTEGER i, j, s
  REAL(DP) :: l(2,2) ! sev eigvals of s matrix (i, spin)
  COMPLEX(DP) :: smato(2,2,2) !(i,j,spin)
  COMPLEX(DP) :: u(2,2,2) ! S eigenvector matrix sev(i,j,spin) 
  COMPLEX(DP) :: invu(2,2,2) ! inverse of S eigenvector matrix invsev(i,j,spin) 
  COMPLEX(DP) :: d(2,2,2) ! diagonal S-1/2
  COMPLEX(DP) :: invssqr(2,2,2) ! S-1/2
  COMPLEX(DP) :: work(2,2,2) 
  !
  ohc = 0.D0
  u = 0.D0
  l = 0.D0
  invu = 0.D0
  d = 0.D0
  invssqr = 0.D0
  !
  DO s = 1, nks
    !
    ! get U and U^-1 and eigenvals of s
    CALL get_evs(smat(:,:,s),2,u(:,:,s),l(:,s))
    CALL get_inv(u(:,:,s),invu(:,:,s))
    !
    ! get diag s^-1/2
    DO i = 1, 2
      d(i,i,s) = SIGN(1.D0,l(i,s))*l(i,s)**(-0.5D0)
    ENDDO
    !
    !
    ! get s^-1/2
    CALL trans(invu(:,:,s),work(:,:,s))
    invssqr(:,:,s) = MATMUL( work(:,:,s), d(:,:,s) )
    invssqr(:,:,s) = MATMUL( invssqr(:,:,s), u(:,:,s) )
    !
    ! get ohc
    CALL trans(invssqr(:,:,s),work(:,:,s))
    ohc(:,:,s) = MATMUL( work(:,:,s), hc(:,:,s) )
    ohc(:,:,s) = MATMUL( ohc(:,:,s) , invssqr(:,:,s) )
    !
  ENDDO
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
!  WRITE(*,'(/A,I2)')' The total number of eigenvalues found:', M
!  CALL PRINT_RMATRIX( 'Selected eigenvalues', 1, M, W, 1 )
!  CALL PRINT_MATRIX( 'Selected eigenvectors (stored columnwise)',&
!                    N, M, Z, LDZ )
END SUBROUTINE get_evs
!
SUBROUTINE PRINT_MATRIX( DESC, M, N, A, LDA )
  CHARACTER(*)    DESC
  INTEGER          M, N, LDA
  COMPLEX*16       A( LDA, * )
  INTEGER          I, J
  WRITE(*,*)
  WRITE(*,*) DESC
  DO I = 1, M
     WRITE(*,9998) ( A( I, J ), J = 1, N )
  END DO
 9998 FORMAT( 11(:,1X,'(',F6.2,',',F6.2,')') )
  RETURN
END SUBROUTINE
!
SUBROUTINE PRINT_RMATRIX( DESC, M, N, A, LDA )
  CHARACTER(*)    DESC
  INTEGER          M, N, LDA
  DOUBLE PRECISION A( LDA, * )
  INTEGER          I, J
  WRITE(*,*)
  WRITE(*,*) DESC
  DO I = 1, M
     WRITE(*,9998) ( A( I, J ), J = 1, N )
  END DO
 9998 FORMAT( 11(:,1X,F6.2) )
  RETURN
END SUBROUTINE
!
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
  b = ain(1,2)
  c = ain(2,1)
  d = ain(2,2)
  !
  bout(1,1) = d/(-b*c + a*d)
  bout(1,2) = -(b/(-b*c + a*d))
  bout(2,1) = -(c/(-b*c + a*d))
  bout(2,2) = a/(-b*c + a*d)
  !
END SUBROUTINE get_inv
!
!-----------------------------------------------------------------------
SUBROUTINE trans(a,b)
  !---------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(IN) :: a(2,2)
  COMPLEX(DP), INTENT(OUT) :: b(2,2)
  INTEGER :: i, j
  !
  DO i = 1, 2
    DO j = 1, 2
      b(i,j)=a(j,i)
    ENDDO
  ENDDO
  !
END SUBROUTINE trans
