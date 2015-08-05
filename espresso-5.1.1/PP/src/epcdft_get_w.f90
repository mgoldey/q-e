!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_w
  !-----------------------------------------------------------------------
  !
  !  <A|W|B> = <A|VA+VB|B> = N \sum_{i,j} (-1)^{i+j} <i|W|j> det( S_minor_i,j )
  !
  !  <i|W|j> is caclualted by sum_g conj(i(g)) Wj(g) not the same thing as sum_g conj(i(g)) w(g) j(g)
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
  USE fft_base, ONLY : dfftp
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, i, j
  INTEGER :: occ(2) ! number of occupied states for that spin
  COMPLEX(DP) :: wevc(npwx, nbnd) ! wtot*evc
  COMPLEX(DP) :: wevc2(npwx, nbnd) ! wtot*evc2
  REAL(DP) :: wtot(dfftp%nnr) ! tot W  w(:,1) + w(:,2)
  COMPLEX(DP), EXTERNAL :: dot, det_s_minor ! dot product and det of S minor
  !
  occ(1) = occup1
  occ(2) = occdown1
  wtot(:) = w(:,1) + w(:,2)
  wmat = 0.D0
  !
  WRITE(*,*)"     NOTE : det_s_minor in epcdft_get_w calculates S every call."
  WRITE(*,*)"     much time can be saved by saving S and just deleting the rows for S minor"
  WRITE(*,*)"     we can change this later."
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
    ! aux = w*phiA
    CALL w_psi(evc, wtot, wevc) 
    CALL w_psi(evc2, wtot, wevc2) 
    !
    DO i = 1, occ(ik)
      DO j = 1, occ(ik)
      !
      ! <B|W|A>                                        <phi_ib|w|phi_ja>
      wmat(1,2,ik) = wmat(1,2,ik) + (-1.D0)**(i+j) * dot(evc2(:,i), wevc(:,j)) * det_s_minor(evc2,evc,i,j,occ(ik))
      !
      ! <A|W|B>
      wmat(2,1,ik) = wmat(1,2,ik) + (-1.D0)**(i+j) * dot(evc(:,i), wevc2(:,j)) * det_s_minor(evc,evc2,i,j,occ(ik))
      !
      ENDDO
    ENDDO
    !
    wmat(1,2,ik) = occ(ik) * wmat(1,2,ik)
    wmat(2,1,ik) = occ(ik) * wmat(2,1,ik)
    !
  ENDDO !ik
  !
  WRITE(*,*)"    W done"
  !
  ! close shop
  !
END SUBROUTINE epcdft_get_w
!
!-----------------------------------------------------------------------
FUNCTION det_s_minor(evc1,evc2,i,j,occ) RESULT(c)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE wvfct, ONLY : npwx 
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: i, j, occ
  COMPLEX(DP) :: c
  COMPLEX(DP) :: c_s_aux(occ,occ) ! work array
  REAL(DP) :: r_s_aux(occ,occ)  ! work array
  COMPLEX(DP), INTENT(INOUT) :: evc1(npwx,occ), evc2(npwx,occ) ! these should be unchanged at end
  COMPLEX(DP) :: aux1(npwx), aux2(npwx) ! work arrays
  !
  aux1(:)=evc1(:,i) ! keep the deleted row/colmn in aux
  aux2(:)=evc2(:,j)
  !
  evc1(:,i)=evc1(:,occ) ! copy last row/colmn to deleted row/colmn
  evc2(:,j)=evc2(:,occ)
  !
  CALL get_det( evc1, evc2, r_s_aux, c_s_aux, occ-1, c ) ! here we only get det S without last row/colmn
  !
  evc1(:,i)=aux1(:) ! put the deleted row/colmn back
  evc2(:,j)=aux2(:)
  !
END FUNCTION det_s_minor
!
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
     c = CMPLX( gdot(a, b) , 0.D0 )
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
!
!-----------------------------------------------------------------------
!SUBROUTINE epcdft_get_s
!  !-----------------------------------------------------------------------
!  !
!  USE kinds, ONLY : DP
!  USE klist, ONLY : nks
!  USE io_files, ONLY : nwordwfc, iunwfc
!  USE io_global, ONLY : ionode, stdout
!  USE wavefunctions_module, ONLY : evc
!  USE wvfct, ONLY : nbnd, npwx, igk, npw , g2kin, ecutwfc
!  USE cell_base, ONLY : tpiba2
!  USE klist, ONLY : xk
!  USE epcdft_mod, ONLY : evc2, iunwfc2, occup1, occdown1, smat
!  USE gvect, ONLY : ngm, g
!  !
!  IMPLICIT NONE
!  !
!  INTEGER :: ik
!  INTEGER :: occ(2) ! number of occupied states for that spin
!  COMPLEX(DP), ALLOCATABLE :: c_s_aux(:,:) ! complex overlap matrix S_ij = <system 1_i | system 2_j>
!  REAL(DP), ALLOCATABLE :: r_s_aux(:,:) ! real overlap matrix S_ij = <system 1_i | system 2_j>
!  !
!  ALLOCATE( c_s_aux(nbnd,nbnd) )
!  ALLOCATE( r_s_aux(nbnd,nbnd) )
!  !
!  ! bounds for spin channels in S matrix
!  !
!  occ(1) = occup1
!  occ(2) = occdown1
!  !
!  ! create S matrix
!  !
!  DO ik = 1 , nks 
!    !
!    ! read wfcs for system 1
!    CALL gk_sort( xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
!    CALL davcio( evc, 2*nwordwfc, iunwfc, ik, -1 )
!    !
!    ! read wfcs for system 1
!    CALL gk_sort( xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
!    CALL davcio( evc2, 2*nwordwfc, iunwfc2, ik, -1 )
!    !
!    ! <1|1>
!    CALL get_det( evc, evc, r_s_aux, c_s_aux, occ(ik), smat(1,1,ik) )
!    ! 
!    ! <1|2>
!    CALL get_det( evc, evc2, r_s_aux, c_s_aux, occ(ik), smat(1,2,ik) )
!    !
!    ! <2|1>
!    CALL get_det( evc2, evc, r_s_aux, c_s_aux, occ(ik), smat(2,1,ik) )
!    !
!    ! <2|2>
!    CALL get_det( evc2, evc2, r_s_aux, c_s_aux, occ(ik), smat(2,2,ik) )
!    !
!  ENDDO !ik
!  !
!  WRITE(*,*)"    S done"
!  !
!  ! close shop
!  !
!  DEALLOCATE( r_s_aux )
!  DEALLOCATE( c_s_aux )
!  !
!END SUBROUTINE epcdft_get_s
!!
