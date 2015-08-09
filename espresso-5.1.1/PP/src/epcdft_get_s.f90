!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_s
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE klist,      ONLY : nks
  USE io_global,  ONLY : ionode, stdout
  USE wvfct,      ONLY : nbnd
  USE epcdft_mod, ONLY : evc1, evc2, occup1, occdown1, smat
  !
  IMPLICIT NONE
  !
  INTEGER :: ik
  INTEGER :: occ(2)                        ! number of occupied states for that spin
  COMPLEX(DP), ALLOCATABLE :: c_s_aux(:,:) ! complex overlap matrix S_ij = <system 1_i | system 2_j>
  REAL(DP), ALLOCATABLE :: r_s_aux(:,:)    ! real overlap matrix S_ij = <system 1_i | system 2_j>
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
    ! <1|1>
    CALL get_det( evc1(:,:,ik), evc1(:,:,ik), r_s_aux, c_s_aux, occ(ik), smat(1,1,ik) )
    ! 
    ! <1|2>
    CALL get_det( evc1(:,:,ik), evc2(:,:,ik), r_s_aux, c_s_aux, occ(ik), smat(1,2,ik) )
    !
    ! <2|1>
    CALL get_det( evc2(:,:,ik), evc1(:,:,ik), r_s_aux, c_s_aux, occ(ik), smat(2,1,ik) )
    !
    ! <2|2>
    CALL get_det( evc2(:,:,ik), evc2(:,:,ik), r_s_aux, c_s_aux, occ(ik), smat(2,2,ik) )
    !
  ENDDO !ik
  !
  IF( ionode ) WRITE( stdout,*)"    S done"
  !
  ! close shop
  !
  DEALLOCATE( r_s_aux )
  DEALLOCATE( c_s_aux )
  !
END SUBROUTINE epcdft_get_s
!
!-----------------------------------------------------------------------------
SUBROUTINE get_det(evc, evc2, r_s_aux, c_s_aux, occ, outdet)
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
  REAL(DP), INTENT(INOUT) :: r_s_aux(occ, occ)
  COMPLEX(DP), INTENT(INOUT) :: c_s_aux(occ, occ)
  COMPLEX(DP), INTENT(INOUT) :: outdet
  !
  r_s_aux = 0.D0
  c_s_aux = 0.D0
  outdet = 0.D0
  !
  IF( gamma_only ) THEN 
      CALL calbec ( npw, evc, evc2, r_s_aux, occ ) ! get over laps of each state
      c_s_aux = CMPLX(r_s_aux, 0.D0) ! pass real to complex
      CALL zgedi_wrap( c_s_aux, occ, outdet ) ! find det of overlap matrix
  ELSE
      CALL calbec ( npw, evc, evc2, c_s_aux, occ )
      CALL zgedi_wrap( c_s_aux, occ, outdet )
  ENDIF
  !
  !
END SUBROUTINE get_det
!
!-----------------------------------------------------------------------------
SUBROUTINE zgedi_wrap(a, n, outdet)
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
END SUBROUTINE zgedi_wrap
