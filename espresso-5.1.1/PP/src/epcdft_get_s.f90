!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_s
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
  USE epcdft_mod, ONLY : evc2, iunwfc2, occup1, occdown1, smat
  USE gvect, ONLY : ngm, g
  USE becmod, ONLY : calbec
  USE control_flags, ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  INTEGER :: occ(2) ! number of occupied states for that spin
  COMPLEX(DP), ALLOCATABLE :: c_s_aux(:,:) ! complex overlap matrix S_ij = <system 1_i | system 2_j>
  REAL(DP), ALLOCATABLE :: r_s_aux(:,:) ! real overlap matrix S_ij = <system 1_i | system 2_j>
  !
  ALLOCATE( c_s_aux(nbnd,nbnd) )
  IF( gamma_only ) ALLOCATE( r_s_aux(nbnd,nbnd) )
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
    IF( gamma_only ) THEN 
      ! 
      ! <1|2>
      CALL calbec ( npw, evc, evc2, r_s_aux, occ(ik) ) ! get over laps of each state
      c_s_aux = CMPLX(r_s_aux, 0.D0) ! pass real to complex
      CALL get_det( c_s_aux, occ(ik), smat(1,2,ik) ) ! find det of overlap matrix
      !
      ! <2|1>
      CALL calbec ( npw, evc2, evc, r_s_aux, occ(ik) )
      c_s_aux = CMPLX(r_s_aux, 0.D0)
      CALL get_det( c_s_aux, occ(ik), smat(2,1,ik) )
      !
      ! <1|1>
      CALL calbec ( npw, evc, evc, r_s_aux, occ(ik) )
      c_s_aux = CMPLX(r_s_aux, 0.D0)
      CALL get_det( c_s_aux, occ(ik), smat(1,1,ik) )
      !
      ! <2|2>
      CALL calbec ( npw, evc2, evc2, r_s_aux, occ(ik) )
      c_s_aux = CMPLX(r_s_aux, 0.D0)
      CALL get_det( c_s_aux, occ(ik), smat(2,2,ik) )
      !
    ELSE
      !
      CALL calbec ( npw, evc, evc2, c_s_aux, occ(ik) )
      CALL get_det( c_s_aux, occ(ik), smat(1,2,ik) )
      !
      CALL calbec ( npw, evc2, evc, c_s_aux, occ(ik) )
      CALL get_det( c_s_aux, occ(ik), smat(2,1,ik) )
      !
      CALL calbec ( npw, evc, evc, c_s_aux, occ(ik) )
      CALL get_det( c_s_aux, occ(ik), smat(1,1,ik) )
      !
      CALL calbec ( npw, evc2, evc2, c_s_aux, occ(ik) )
      CALL get_det( c_s_aux, occ(ik), smat(2,2,ik) )
      !
    ENDIF
    !
  ENDDO !ik
  !
  ! close shop
  !
  IF( gamma_only ) DEALLOCATE( r_s_aux )
  DEALLOCATE( c_s_aux )
  !
END SUBROUTINE epcdft_get_s
!
!-----------------------------------------------------------------------------
SUBROUTINE get_det(a, n, outdet)
  !--------------------------------------------------------------------------
  !
  !     this routine will solve for the determinant of a matrix "a"
  !     and store it in outdet
  !
  USE kinds, ONLY : DP
  USE io_global,            ONLY : ionode
  !
  IMPLICIT NONE
  !
  INTEGER,      INTENT(IN)    :: n
  COMPLEX(DP),  INTENT(INOUT) :: outdet
  COMPLEX(DP),  INTENT(IN)    :: a(n,n)
  !
  INTEGER                  :: i
  INTEGER                  :: info
  INTEGER                  :: ivpt(n)
  EXTERNAL                 :: zgefa
  EXTERNAL                 :: zgedi
  COMPLEX(DP)              :: work(n)
  COMPLEX(DP) :: det(2)
  !
  det = 0.D0
  !
  ! factor "a" by gaussian elimination
  ! a will be upper trianglular
  CALL zgefa(a, n, n, ivpt, info)
  ! write(*,*) "Info is ", info
  !
  CALL errore( 'epcdft_coupling', 'error in zgefa', abs(info))
  !
  ! compute det of "a" using factors computed by zgefa
  !
  ! !!! the det = det(1) * 10**det(2) !!!
  !
  CALL zgedi(a, n, n, ivpt, det, work, 10)
  !
  ! compute det from output of zgedi
  !
  outdet = det(1) * 10.d0**( det(2) )
  !
  RETURN
  !
END SUBROUTINE get_det
