!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_w
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE klist, ONLY : nks
  USE io_files, ONLY : nwordwfc, iunwfc
  USE io_global, ONLY : ionode, stdout
  USE wavefunctions_module, ONLY : evc
  USE wvfct, ONLY : nbnd, npwx, igk, npw , g2kin, ecutwfc
  USE cell_base, ONLY : tpiba2
  USE klist, ONLY : xk
  USE epcdft_mod, ONLY : evc2, iunwfc2, occup1, occdown1, wmat, w
  USE gvect, ONLY : ngm, g
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, ib
  INTEGER :: occ(2) ! number of occupied states for that spin
  COMPLEX(DP), ALLOCATABLE :: c_s_aux(:,:) ! complex overlap matrix S_ij = <system 1_i | system 2_j>
  REAL(DP), ALLOCATABLE :: r_s_aux(:,:) ! real overlap matrix S_ij = <system 1_i | system 2_j>
  COMPLEX(DP) :: aux(npwx, nbnd)
  !
  ALLOCATE( c_s_aux(nbnd,nbnd) )
  ALLOCATE( r_s_aux(nbnd,nbnd) )
  !
  ! bounds for spin channels in S matrix
  !
  occ(1) = occup1
  occ(2) = occdown1
  !
  ! create S matrix
  !
  DO ik = 1 , nks 
    !
    ! read wfcs for system 1
    CALL gk_sort( xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
    CALL davcio( evc, 2*nwordwfc, iunwfc, ik, -1 )
    !
    ! read wfcs for system 1
    CALL gk_sort( xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
    CALL davcio( evc2, 2*nwordwfc, iunwfc2, ik, -1 )
    !
    !
    !
    ! aux = w1*evc1
    !
    CALL w_psi(evc, w(:,1), aux) 
    !
    ! <1|w1|1>
    CALL get_det( evc, aux, r_s_aux, c_s_aux, occ(ik), wmat(1,1,1,ik) )
    !
    ! <2|w1|1>
    CALL get_det( evc2, aux, r_s_aux, c_s_aux, occ(ik), wmat(2,1,1,ik) )
    !
    !
    !
    ! aux = w1*evc2
    !
    CALL w_psi(evc2, w(:,1), aux) 
    !
    ! <1|w1|2>
    CALL get_det( evc, aux, r_s_aux, c_s_aux, occ(ik), wmat(1,2,1,ik) )
    !
    ! <2|w1|2>
    CALL get_det( evc2, aux, r_s_aux, c_s_aux, occ(ik), wmat(2,2,1,ik) )
    !
    !
    !
    ! aux = w2*evc2
    CALL w_psi(evc2, w(:,2), aux) 
    ! 
    ! <1|w2|2>
    CALL get_det( evc, aux, r_s_aux, c_s_aux, occ(ik), wmat(1,2,2,ik) )
    !
    ! <2|w2|2>
    CALL get_det( evc2, aux, r_s_aux, c_s_aux, occ(ik), wmat(2,2,2,ik) )
    !
    !
    !
    ! aux = w2*evc
    CALL w_psi(evc, w(:,2), aux) 
    ! 
    ! <1|w2|1>
    CALL get_det( evc, aux, r_s_aux, c_s_aux, occ(ik), wmat(1,1,2,ik) )
    !
    ! <2|w2|1>
    CALL get_det( evc2, aux, r_s_aux, c_s_aux, occ(ik), wmat(2,1,2,ik) )
    !
  ENDDO !ik
  !
  WRITE(*,*)"    W done"
  !
  ! close shop
  !
  DEALLOCATE( r_s_aux )
  DEALLOCATE( c_s_aux )
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
