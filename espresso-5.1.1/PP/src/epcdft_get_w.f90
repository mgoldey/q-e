!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_w
  !-----------------------------------------------------------------------
  !
  !  Compute W matrix
  !
  !    <A|W|A> = <A|Wa|A>  (divide by the amplitude)
  !    <A|W|B> = <A|Wb|B> = N \sum_{i,j}   \langle i | W | j \rangle  (S^{-1} (det(S)I))^T_{i,j} 
  !                                                                     ^-----cofactor----------^
  !
  USE kinds,                ONLY : DP
  USE klist,                ONLY : nks
  USE io_global,            ONLY : ionode, stdout
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npwx
  USE epcdft_mod,           ONLY : evc1, evc2, occup1, occdown1, wmat, w, smat
  USE fft_base,             ONLY : dfftp
  USE mp,                   ONLY : mp_sum
  USE mp_images,            ONLY : intra_image_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, i, j
  INTEGER :: occ(2)                  ! number of occupied states for that spin
  COMPLEX(DP) :: wevc(npwx, nbnd)    ! wtot*evc
  COMPLEX(DP) :: wevc2(npwx, nbnd)   ! wtot*evc2
  REAL(DP) :: wtot(dfftp%nnr)        ! tot W  w(:,1) + w(:,2)
  COMPLEX(DP), EXTERNAL :: dot
  COMPLEX(DP) :: cofc(nbnd,nbnd,2,2) ! cofactor ( i, j, ab ba, spin up down )
  !
  occ(1) = occup1
  occ(2) = occdown1
  wmat = ( 0.D0, 0.D0 )
  !
  ! create S matrix
  !
  DO ik = 1 , nks 
    !
    ! S^-1  ab 
    !
    CALL get_s_invs(evc1(:,:,ik), evc2(:,:,ik), cofc(:,:,1,ik), nbnd)
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
    ! wevc = w*evc
    !
    wtot(:) = w(:,1)
    CALL w_psi(evc1(:,:,ik), wtot, wevc) 
    !
    wtot(:) = w(:,2)
    CALL w_psi(evc2(:,:,ik), wtot, wevc2) 
    !
    DO i = 1, occ(ik)
      DO j = 1, occ(ik)
        !
        ! <A|W|A>                                 
        IF(i==j) wmat(1,1,ik) = wmat(1,1,ik) + dot(evc1(:,i,ik), wevc(:,j))
        !
        ! <B|Wa|A>                                 
        wmat(1,2,ik) = wmat(1,2,ik) + dot(evc2(:,j,ik), wevc(:,i)) * cofc(i,j,1,ik)
        !
        ! <A|Wb|B>                                 
        wmat(2,1,ik) = wmat(2,1,ik) + dot(evc1(:,j,ik), wevc2(:,i)) * cofc(i,j,2,ik)
        !
        ! <B|W|B>                                 
        IF(i==j) wmat(2,2,ik) = wmat(2,2,ik) + dot(evc2(:,i,ik), wevc2(:,j))
        !
      ENDDO !j
    ENDDO !i
    !
  ENDDO !ik
  !
  CALL mp_sum(wmat,intra_image_comm)
  !
  IF( ionode ) WRITE( stdout,* )"    W done"
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
  USE pwcom,      ONLY : npwx, npw
  USE gvect,      ONLY : gstart
  !
  IMPLICIT NONE
  !
  COMPLEX(DP),INTENT(IN) :: a(npwx), b(npwx)
  REAL(DP), EXTERNAL :: DDOT
  REAL(DP) :: c
  !
  c = 2._DP * DDOT(2*npw, a, 1, b, 1)
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
      c_s_aux = CMPLX(r_s_aux, 0.D0, KIND=DP) 
      CALL invmat_complex (occ, c_s_aux, sinvs, filler)
  ELSE
      CALL calbec ( npw, evc, evc2, c_s_aux, occ )
      CALL invmat_complex (occ, c_s_aux, sinvs, filler)
  ENDIF
  !
END SUBROUTINE get_s_invs
