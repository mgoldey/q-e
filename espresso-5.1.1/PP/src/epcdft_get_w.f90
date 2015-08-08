!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_w
  !-----------------------------------------------------------------------
  !
  !  <A|W|B> = <A|VA+VB|B> = N \sum_{i,j}   \langle i | W | j \rangle  (S^{-1} (det(S)I))^T_{i,j} 
  !                                                                     ^-----cofactor----------^
  !
  USE kinds, ONLY : DP
  USE klist, ONLY : nks
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global, ONLY : ionode, stdout
  USE wavefunctions_module, ONLY : evc
  USE wvfct, ONLY : nbnd, npwx, igk, npw , g2kin, ecutwfc
  USE cell_base, ONLY : tpiba2
  USE klist, ONLY : xk
  USE epcdft_mod, ONLY : evc2, iunwfc2, occup1, occdown1, wmat, w, smat
  USE gvect, ONLY : ngm, g
  USE fft_base, ONLY : dfftp
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, i, j
  INTEGER :: occ(2) ! number of occupied states for that spin
  COMPLEX(DP) :: wevc(npwx, nbnd) ! wtot*evc
  COMPLEX(DP) :: wevc2(npwx, nbnd) ! wtot*evc2
  REAL(DP) :: wtot(dfftp%nnr) ! tot W  w(:,1) + w(:,2)
  COMPLEX(DP), EXTERNAL :: dot
  COMPLEX(DP) :: cofc(nbnd,nbnd,2,2) ! cofactor ( i, j, ab ba, spin up down )
  !
  occ(1) = occup1
  occ(2) = occdown1
  wtot(:) = w(:,1) + w(:,2)
  wmat = 0.D0
  !
  ! create S matrix
  !
  DO ik = 1 , nks 
    !
    ! read wfcs for system 1
    CALL gk_sort( xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
    CALL davcio( evc, 2*nwordwfc, iunwfc, ik, -1 )
    !
    ! read wfcs for system 2
    CALL gk_sort( xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
    CALL davcio( evc2, 2*nwordwfc, iunwfc2, ik, -1 )
    !
    ! S^-1  ab 
    !
    CALL get_s_invs(evc, evc2, cofc(:,:,1,ik), nbnd)
    !
    ! S^-1 * det(S) ab 
    !
    cofc(:,:,1,ik) = cofc(:,:,1,ik) * smat(1,2,ik)
    !
    ! C = ( S^-1 * det(S))^T   ab & ba 
    !
    cofc(:,:,2,ik) = CONJG(cofc(:,:,1,ik))
    cofc(:,:,1,ik) = TRANSPOSE(cofc(:,:,1,ik))
    !
    ! aux = w*phiA
    CALL w_psi(evc, wtot, wevc) 
    CALL w_psi(evc2, wtot, wevc2) 
    !
    DO i = 1, occ(ik)
      DO j = 1, occ(ik)
        !
        ! <B|W|A>                                 
        wmat(1,2,ik) = wmat(1,2,ik) + dot(evc2(:,i), wevc(:,j)) * cofc(i,j,1,ik)
        !
        ! <A|W|B>
        wmat(2,1,ik) = wmat(2,1,ik) + dot(evc(:,i), wevc2(:,j)) * cofc(i,j,2,ik)
        !
      ENDDO
    ENDDO
    !
  ENDDO !ik
  !
  WRITE(*,*)"    W done Note"
  !
END SUBROUTINE epcdft_get_w
!
!-----------------------------------------------------------------------
SUBROUTINE w_psi(evc, w, auxg)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : nbnd, npwx, igk, npw , g2kin, ecutwfc
  USE gvect,                ONLY : nl, nlm
  USE gvecs,                ONLY : nls, nlsm
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE fft_base,             ONLY : dfftp, dffts
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), INTENT(INOUT) :: auxg(npwx, nbnd)
  COMPLEX(DP), INTENT(IN) :: evc(npwx, nbnd)
  REAL(DP), INTENT(IN) :: w(dfftp%nnr)
  !
  INTEGER :: ib
  COMPLEX(DP) :: aux(dfftp%nnr, nbnd)
  !
  aux = 0.D0
  auxg = 0.D0
  !
  DO ib = 1, nbnd
    aux( nl(igk(1:npw)), ib ) = evc( 1:npw , ib )
    IF(gamma_only) aux( nlm(igk(1:npw)), ib ) = CONJG( evc( 1:npw, ib ) )
    CALL invfft ('Wave', aux(:,ib), dfftp)
    aux(:,ib) = aux(:,ib) * w(:)
    CALL fwfft  ('Wave', aux(:,ib), dfftp)
    auxg(1:npw, ib) = aux( nl(igk(1:npw)), ib ) 
  ENDDO
  !
END SUBROUTINE w_psi
!
!
!-----------------------------------------------------------------------
FUNCTION dot(a, b) result(c)
  !-----------------------------------------------------------------------
  !
  ! return dot prod of two complex vecs
  !
  USE kinds,      ONLY : DP 
  USE pwcom,      ONLY : npwx
  USE control_flags, ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),INTENT(IN) :: a(npwx), b(npwx)
  COMPLEX(DP), EXTERNAL :: zdotc
  REAL(DP), EXTERNAL :: gdot ! gamma dot prod function
  COMPLEX(DP) :: c
  !
  IF (gamma_only) THEN
     c = CMPLX( gdot(a, b) , 0.D0, KIND=DP )
  ELSE
     c = zdotc( npwx, a, 1, b, 1 )
  ENDIF
  !
END FUNCTION dot
!
!-----------------------------------------------------------------------
FUNCTION gdot(a, b) result(c)
  !-----------------------------------------------------------------------
  !
  ! return real dot prod of two complex vec using gamma point tricks
  !
  USE kinds,      ONLY : DP 
  USE pwcom,      ONLY : npwx
  USE gvect,      ONLY : gstart
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),INTENT(IN) :: a(npwx), b(npwx)
  REAL(DP), EXTERNAL :: DDOT
  REAL(DP) :: c
  !
  c = 2._DP * DDOT(2*npwx, a, 1, b, 1)
  !
  IF(gstart==2) THEN
     c = c - REAL(a(1),KIND=DP)*REAL(b(1),KIND=DP)
  ENDIF
  !
END FUNCTION gdot
!
!-----------------------------------------------------------------------------
SUBROUTINE get_s_invs(evc, evc2, sinvs, occ)
  !--------------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE control_flags, ONLY : gamma_only
  USE wvfct, ONLY : npwx, npw
  USE becmod, ONLY : calbec
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: occ
  COMPLEX(DP), INTENT(IN) :: evc(npwx, occ), evc2(npwx, occ)
  COMPLEX(DP), INTENT(INOUT) :: sinvs(occ, occ)
  REAL(DP) :: r_s_aux(occ, occ)
  COMPLEX(DP) :: c_s_aux(occ, occ)
  COMPLEX(DP) :: filler
  !
  r_s_aux = 0.D0
  c_s_aux = 0.D0
  sinvs = 0.D0
  !
  IF( gamma_only ) THEN
      CALL calbec ( npw, evc, evc2, r_s_aux, occ ) ! get over laps of each state
      c_s_aux = CMPLX(r_s_aux, 0.D0) 
      CALL invmat_complex (occ, c_s_aux, sinvs, filler)
  ELSE
      CALL calbec ( npw, evc, evc2, c_s_aux, occ )
      CALL invmat_complex (occ, c_s_aux, sinvs, filler)
  ENDIF
  !
END SUBROUTINE get_s_invs
