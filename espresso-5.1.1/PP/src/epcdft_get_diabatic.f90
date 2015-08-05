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
  REAL(DP) :: fmat(2), sev(2,2) ! sev eigvals of s matrix (i, spin)
  COMPLEX(DP) :: smato(2,2,2) !(i,j,spin)
  COMPLEX(DP) :: sevec(2,2,2) ! S eigenvector matrix sev(i,j,spin) 
  COMPLEX(DP) :: invsevec(2,2,2) ! inverse of S eigenvector matrix invsev(i,j,spin) 
  COMPLEX(DP) :: sd(2,2,2) ! diagonal S 
  COMPLEX(DP) :: sdinvsqr(2,2,2) ! diagonal S-1/2
  COMPLEX(DP) :: sinvsqr(2,2,2) ! S-1/2
  !
  ohc = 0.D0
  smato = 0.D0
  sevec = 0.D0
  sev = 0.D0
  sd = 0.D0
  sdinvsqr = 0.D0
  sinvsqr = 0.D0
  invsevec = 0.D0
  fmat(1) = free1
  fmat(2) = free2
  !
  ! WORK ZONE 
  ! WARNING CRAP BELOW
  CALL get_evs(smat(:,:,1),2,sevec(:,:,1),sev(:,1))
  CALL get_inv(sevec(:,:,1),invsevec(:,:,1))
  !
  sd(:,:,1) = MATMUL(sevec(:,:,1),smat(:,:,1))
  sd(:,:,1) = MATMUL(sd(:,:,1),invsevec(:,:,1))
  !
  DO i = 1, 2
    DO j = 1, 2
      IF(i/=j)THEN
        sdinvsqr(i,j,1) = 0.D0
      ELSE
        sdinvsqr(i,j,1) = sev(i,1)**(-0.5D0)
      ENDIF
    ENDDO
  ENDDO
  !
  sinvsqr(:,:,1) = MATMUL(invsevec(:,:,1), sdinvsqr(:,:,1))
  sinvsqr(:,:,1) = MATMUL(sinvsqr(:,:,1),sevec(:,:,1))
  !
  ohc(:,:,1) = MATMUL(sinvsqr(:,:,1), hc(:,:,1))
  ohc(:,:,1) = MATMUL(ohc(:,:,1),sinvsqr(:,:,1))
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
  !
  USE kinds, ONLY : DP
  ! gives inverse of 2x2 matrix ain
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
