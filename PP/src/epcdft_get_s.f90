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
SUBROUTINE epcdft_get_s
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE klist,      ONLY : nks
  USE io_global,  ONLY : ionode, stdout
  USE wvfct,      ONLY : nbnd
  USE epcdft_mod, ONLY : evc1, evc2, occup1, occdown1, smat, debug2
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
  DEALLOCATE( r_s_aux )
  DEALLOCATE( c_s_aux )
  !
  IF(debug2) CALL epcdft_write_s ( smat ) ! create overlap matrix
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
  USE io_global,  ONLY : ionode, stdout
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: occ
  COMPLEX(DP), INTENT(IN) :: evc(npwx, occ), evc2(npwx, occ)
  REAL(DP), INTENT(INOUT) :: r_s_aux(occ, occ)
  COMPLEX(DP), INTENT(INOUT) :: c_s_aux(occ, occ)
  COMPLEX(DP), INTENT(INOUT) :: outdet
  integer :: i,j
  !
  r_s_aux = 0.D0
  c_s_aux = 0.D0
  outdet = 0.D0
  !
  IF( gamma_only ) THEN 
      CALL calbec ( npw, evc, evc2, r_s_aux, occ ) ! get over laps of each state
      c_s_aux = CMPLX(r_s_aux, 0.D0, KIND=DP) ! pass real to complex
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
!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_write_s(smat)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE klist,      ONLY : nks
  USE wvfct,      ONLY : nbnd
  USE epcdft_mod, ONLY : evc1, evc2, occup1, occdown1
  USE becmod, ONLY : calbec
  USE control_flags, ONLY : gamma_only
  USE wvfct, ONLY : npwx, npw
  USE wvfct,      ONLY : nbnd
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, filunit, i, j
  INTEGER :: occ(2)                        ! number of occupied states for that spin
  COMPLEX(DP), ALLOCATABLE :: c_s_aux(:,:) ! complex overlap matrix S_ij = <system 1_i | system 2_j>
  REAL(DP), ALLOCATABLE :: r_s_aux(:,:)    ! real overlap matrix S_ij = <system 1_i | system 2_j>
  COMPLEX(DP), INTENT(IN) :: smat(2,2,nks)
  CHARACTER(LEN=256) :: fname
  !
  ALLOCATE( c_s_aux(nbnd,nbnd) )
  ALLOCATE( r_s_aux(nbnd,nbnd) )
  !
  ! bounds for spin channels in S matrix
  !
  occ(1) = occup1
  occ(2) = occdown1
  filunit=234976
  !fname = 's'
  !
  ! create S matrix
  !
  IF(gamma_only)THEN
    DO ik = 1 , nks 
      !
      ! <1|1>
      r_s_aux = 0.d0
      CALL calbec( npw, evc1(:,:,ik), evc1(:,:,ik), r_s_aux, occ(ik))
      c_s_aux(:,:) = CMPLX(r_s_aux(:,:),0.d0,KIND=DP)
      WRITE(unit=fname,fmt=*) ik
      fname="Saa"//TRIM(ADJUSTL(fname))
      CALL realpart_dumpmat(fname,filunit,c_s_aux,nbnd,occ(ik))
      ! 
      ! <1|2>
      r_s_aux = 0.d0
      CALL calbec( npw, evc1(:,:,ik), evc2(:,:,ik), r_s_aux, occ(ik))
      c_s_aux = CMPLX(r_s_aux,KIND=DP)
      WRITE(unit=fname,fmt=*) ik
      fname="Sab"//TRIM(ADJUSTL(fname))
      CALL realpart_dumpmat(fname,filunit,c_s_aux,nbnd,occ(ik))
      !
      ! <2|1>
      r_s_aux = 0.d0
      CALL calbec( npw, evc2(:,:,ik), evc1(:,:,ik), r_s_aux, occ(ik))
      c_s_aux = CMPLX(r_s_aux,KIND=DP)
      WRITE(unit=fname,fmt=*) ik
      fname="Sba"//TRIM(ADJUSTL(fname))
      CALL realpart_dumpmat(fname,filunit,c_s_aux,nbnd,occ(ik))
      !
      ! <2|2>
      r_s_aux = 0.d0
      CALL calbec( npw, evc2(:,:,ik), evc2(:,:,ik), r_s_aux, occ(ik))
      c_s_aux = CMPLX(r_s_aux,KIND=DP)
      WRITE(unit=fname,fmt=*) ik
      fname="Sbb"//TRIM(ADJUSTL(fname))
      CALL realpart_dumpmat(fname,filunit,c_s_aux,nbnd,occ(ik))
      !
      WRITE(unit=fname,fmt=*) ik
      fname="detS"//TRIM(ADJUSTL(fname))
      CALL realpart_dumpmat(fname,filunit,smat(:,:,ik),2,2)
      !
    ENDDO !ik
  ENDIF
  !
  DEALLOCATE( r_s_aux )
  !
END SUBROUTINE epcdft_write_s
!
!-----------------------------------------------------------------------
SUBROUTINE realpart_dumpmat(fname,filunit,mat,m,n)
  !-----------------------------------------------------------------------
  !
  ! will print real
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j
  INTEGER,INTENT(IN) :: filunit, n,m
  COMPLEX(DP),INTENT(IN) :: mat(m,m)
  CHARACTER(LEN=256),INTENT(IN) :: fname
  !
  !OPEN(UNIT=filunit,FILE=TRIM(ADJUSTL(fname))//".cub")
  IF( ionode ) THEN
    OPEN(UNIT=filunit,FILE=fname)
    DO i = 1, n
      WRITE(filunit, *) ( REAL(mat(i,j),KIND=DP) , j = 1, n )
    ENDDO
  ENDIF
  !-----------------------------------------------------------------------
END SUBROUTINE realpart_dumpmat
