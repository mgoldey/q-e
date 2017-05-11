!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors to this file:
!   Nicholas P. Brawand (nicholasbrawand@gmail.com)
!   Matthew B. Goldey (matthew.goldey@gmail.com)
!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_diabatic
  !-----------------------------------------------------------------------
  !
  IMPLICIT NONE
  !
  CALL epcdft_get_diabatic_lowdin()
  !
END SUBROUTINE epcdft_get_diabatic
!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_diabatic_lowdin()
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE epcdft_mod, ONLY : smat, hc, ohc, debug2
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
  CHARACTER(LEN=256) :: fname
  INTEGER :: filunit=3234873
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
  !
  CALL cdiagh(2,smattot,2,l,u)
  !
  ! get diag S^-1/2
  !
  DO i = 1, 2
    d(i,i) = l( i )**(-0.5D0)
  ENDDO
  !
  !
  ! get S^-1/2 = U . SD^-1/2 . U^-1
  ! get S^-1/2 = U . SD^-1/2 . U^dagger
  !
  invssqr(:,:) = MATMUL( u(:,:), d(:,:) )
  invssqr(:,:) = MATMUL( invssqr(:,:), TRANSPOSE(DCONJG(u(:,:))) )
  !
  IF(debug2)THEN
    fname='sinvoh'
    CALL realpart_dumpmat(fname,filunit,invssqr,2,2)
  ENDIF
  !
  ! get oHc = S^-1/2 . Hc . S^-1/2
  !
  ohc(:,:) = MATMUL( invssqr(:,:), hc(:,:) )
  ohc(:,:) = MATMUL( ohc(:,:) , invssqr(:,:) )
  !
  IF(debug2)THEN
    fname='HorthLow'
    CALL realpart_dumpmat(fname,filunit,ohc,2,2)
  ENDIF
  !
END SUBROUTINE epcdft_get_diabatic_lowdin
!
SUBROUTINE invrts(s,sinv)
  !
  ! compute s^-1
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  COMPLEX(DP),INTENT(IN)::s(2,2)
  COMPLEX(DP),INTENT(INOUT)::sinv(2,2)
  COMPLEX(DP)::det
  !
  sinv = (0.D0,0.D0)
  det = (0.D0,0.D0)
  !
  det = (s(1,1)*s(2,2) - s(1,2)*s(2,1))
  sinv(1,1) = s(2,2)
  sinv(1,2) = -s(1,2)
  sinv(2,1) = -s(2,1)
  sinv(2,2) = s(1,1)
  sinv(:,:) = sinv(:,:)/det
  !
ENDSUBROUTINE invrts
!
SUBROUTINE findnormvec(v)
  !
  ! normalize vector
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  COMPLEX(DP),INTENT(INOUT)::v(2)
  COMPLEX(DP)::norm
  !
  norm = v(1)*CONJG(v(1)) + v(2)*CONJG(v(2))
  v(:) = v(:)/SQRT(norm)
  !
ENDSUBROUTINE findnormvec
