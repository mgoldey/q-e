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
SUBROUTINE epcdft_get_w
  !-----------------------------------------------------------------------
  !
  !  Compute W matrix
  !
  !    <A|W|A> = <A|Wa|A> 
  !    <A|W|B> = <A|Wb|B> = \sum_{i,j}   \langle i | W | j \rangle  (S^{-1} (det(S)I))^T_{i,j} 
  !                                                                     ^-----cofactor----------^
  !
  !  Below I sacrifice speed for assurance by avoiding simplifications
  !
  !
  USE kinds,                ONLY : DP
  USE klist,                ONLY : nks, igk_k
  USE io_global,            ONLY : ionode, stdout
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npwx, npw
  USE epcdft_mod,           ONLY : evc1, evc2, occup1, occdown1, wmat, w, smat, debug, debug2
  USE fft_base,             ONLY : dfftp
  USE mp,                   ONLY : mp_sum
  USE mp_images,            ONLY : intra_image_comm
  USE becmod,               ONLY : calbec
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, i, j, filunit
  INTEGER :: occ(2)                  ! number of occupied states for that spin
  REAL(DP)    :: wtot(dfftp%nnr)     ! tot W  w(:,1) + w(:,2)
  COMPLEX(DP) :: wevc(npwx, nbnd)    ! wtot*evc
  COMPLEX(DP) :: wevc2(npwx, nbnd)   ! wtot*evc2
  COMPLEX(DP) :: cofc(nbnd,nbnd,2,2,2) ! cofactor ( i, j, row a b, colm a b , spin up down )
  CHARACTER(LEN=256) :: fname
  !
  occ(1) = occup1
  occ(2) = occdown1
  wmat = ( 0.D0, 0.D0 )
  !
  ! first build cofactor matrix and apply w to psi then 
  ! calculate wmat for each spin
  !
  DO ik = 1 , nks 
    !
    ! S^-1  
    !
    CALL get_s_invs(evc1(:,1:occ(ik),ik), evc1(:,1:occ(ik),ik), cofc(1:occ(ik),1:occ(ik),1,1,ik), occ(ik) )
    CALL get_s_invs(evc1(:,1:occ(ik),ik), evc2(:,1:occ(ik),ik), cofc(1:occ(ik),1:occ(ik),1,2,ik), occ(ik) )
    CALL get_s_invs(evc2(:,1:occ(ik),ik), evc1(:,1:occ(ik),ik), cofc(1:occ(ik),1:occ(ik),2,1,ik), occ(ik) )
    CALL get_s_invs(evc2(:,1:occ(ik),ik), evc2(:,1:occ(ik),ik), cofc(1:occ(ik),1:occ(ik),2,2,ik), occ(ik) )
    !
    ! C = (S^-1 * det(S))^T
    !
    DO i = 1, 2
      DO j = 1, 2
        !
        cofc(:,:,i,j,ik) = cofc(:,:,i,j,ik) * smat(i,j,ik)
        !
        cofc(:,:,i,j,ik) = TRANSPOSE(cofc(:,:,i,j,ik))
        !
      ENDDO
    ENDDO
    !
    ! write cofac to file
    !
    IF(debug2)THEN
      !
      WRITE(unit=fname,fmt=*) ik
      fname="Caa"//TRIM(ADJUSTL(fname))
      CALL realpart_dumpmat(fname, filunit, cofc(:,:,1,1,ik), nbnd, occ(ik))
      WRITE(unit=fname,fmt=*) ik
      fname="Cab"//TRIM(ADJUSTL(fname))
      CALL realpart_dumpmat(fname, filunit, cofc(:,:,1,2,ik), nbnd, occ(ik))
      WRITE(unit=fname,fmt=*) ik
      fname="Cba"//TRIM(ADJUSTL(fname))
      CALL realpart_dumpmat(fname, filunit, cofc(:,:,2,1,ik), nbnd, occ(ik))
      WRITE(unit=fname,fmt=*) ik
      fname="Cbb"//TRIM(ADJUSTL(fname))
      CALL realpart_dumpmat(fname, filunit, cofc(:,:,2,2,ik), nbnd, occ(ik))
      !
    ENDIF
    !
    ! wevc = FFT w_A,ik(r) * evc_A,ik(r)
    !
    CALL w_psi(evc1(:,:,ik), w(:,1,ik), wevc) 
    CALL w_psi(evc2(:,:,ik), w(:,2,ik), wevc2) 
    !
    ! calculate wmat using cofc and w_psi from above
    !
    IF(gamma_only)THEN
      CALL calc_w_real(ik, occ, wevc, wevc2, cofc)
    ELSE
      CALL calc_w_img(ik, occ, wevc, wevc2, cofc)
    ENDIF
    !
    ! write wmat to file
    !
    IF(debug2)THEN
      !
      WRITE(unit=fname,fmt=*) ik
      fname="W"//TRIM(ADJUSTL(fname))
      CALL realpart_dumpmat(fname, filunit, wmat(:,:,ik), 2, 2)
      !
    ENDIF
    !
  ENDDO ! ik
  !
  IF( ionode ) WRITE( stdout,* )"    W done"
  !
  !IF( debug ) THEN
  !  IF( ionode ) WRITE( stdout,* )"    Check W"
  !  IF( ionode ) WRITE( stdout,* )"      W_AA,up + W_AA,down should be = C1 ", wmat(1,1,1) + wmat(1,1,2)
  !  IF( ionode ) WRITE( stdout,* )"      W_BB,up + W_BB,down should be = C2 ", wmat(2,2,1) + wmat(2,2,2)
  !  IF( ionode ) WRITE( stdout,* )"    END Check W"
  !ENDIF
  !
END SUBROUTINE epcdft_get_w
!
!-----------------------------------------------------------------------
SUBROUTINE w_psi(evc, w, auxg)
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE wvfct,                ONLY : nbnd, npwx, npw , g2kin
  USE klist,                ONLY : igk_k
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
    aux( nl(igk_k(1:npw,1)), ib ) = evc( 1:npw , ib )
    IF(gamma_only) aux( nlm(igk_k(1:npw,1)), ib ) = CONJG( evc( 1:npw, ib ) )
    CALL invfft ('Wave', aux(:,ib), dfftp)
    aux(:,ib) = aux(:,ib) * w(:)
    CALL fwfft  ('Wave', aux(:,ib), dfftp)
    auxg(1:npw, ib) = aux( nl(igk_k(1:npw,1)), ib ) 
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
  USE matrix_inversion
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
      CALL invmat (occ, c_s_aux, sinvs)
  ELSE
      CALL calbec ( npw, evc, evc2, c_s_aux, occ )
      CALL invmat (occ, c_s_aux, sinvs)
  ENDIF
  !
END SUBROUTINE get_s_invs
!
!-----------------------------------------------------------------------------
SUBROUTINE calc_w_real(ik, occ, wevc, wevc2, cofc)
  !-----------------------------------------------------------------------------
  !
  !    <A_i|W|B_j> = \sum_{i,j}   \langle evc_A,i | wevc_B,j \rangle cofc_{i,j} 
  !
  USE kinds,                ONLY : DP
  USE klist,                ONLY : nks
  USE io_global,            ONLY : ionode, stdout
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npwx, npw
  USE epcdft_mod,           ONLY : evc1, evc2, occup1, occdown1, wmat, w, smat, debug2
  USE fft_base,             ONLY : dfftp
  USE mp,                   ONLY : mp_sum
  USE mp_images,            ONLY : intra_image_comm
  USE becmod,               ONLY : calbec
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j, filunit
  INTEGER, INTENT(IN) :: ik
  INTEGER, INTENT(IN) :: occ(2)                  ! number of occupied states for that spin
  COMPLEX(DP), INTENT(IN) :: wevc(npwx, nbnd)    ! wtot*evc
  COMPLEX(DP), INTENT(IN) :: wevc2(npwx, nbnd)   ! wtot*evc2
  COMPLEX(DP), INTENT(IN) :: cofc(nbnd,nbnd,2,2,2) ! cofactor ( i, j, row a b, colm a b , spin up down )
  REAL(DP) :: single_electron_w(nbnd,nbnd,2,2,2) ! <phi_i|w|phi_j> 
  !                                              ! single electron orbitals (i,j,row a b, colm a b, spin up down)
  CHARACTER(LEN=256) :: fname
  single_electron_w = 0.D0
  !
  CALL calbec ( npw, evc1(:,:,ik), wevc,  single_electron_w(:,:,1,1,ik), nbnd )
  CALL calbec ( npw, evc1(:,:,ik), wevc2, single_electron_w(:,:,1,2,ik), nbnd )
  CALL calbec ( npw, evc2(:,:,ik), wevc,  single_electron_w(:,:,2,1,ik), nbnd )
  CALL calbec ( npw, evc2(:,:,ik), wevc2, single_electron_w(:,:,2,2,ik), nbnd )
  !
  ! write w single electron to file
  !
  IF(debug2)THEN
    !
    WRITE(unit=fname,fmt=*) ik
    fname="Waa"//TRIM(ADJUSTL(fname))
    CALL real_dumpmat(fname, filunit, single_electron_w(:,:,1,1,ik), nbnd, occ(ik))
    WRITE(unit=fname,fmt=*) ik
    fname="Wab"//TRIM(ADJUSTL(fname))
    CALL real_dumpmat(fname, filunit, single_electron_w(:,:,1,2,ik), nbnd, occ(ik))
    WRITE(unit=fname,fmt=*) ik
    fname="Wba"//TRIM(ADJUSTL(fname))
    CALL real_dumpmat(fname, filunit, single_electron_w(:,:,2,1,ik), nbnd, occ(ik))
    WRITE(unit=fname,fmt=*) ik
    fname="Wbb"//TRIM(ADJUSTL(fname))
    CALL real_dumpmat(fname, filunit, single_electron_w(:,:,2,2,ik), nbnd, occ(ik))
    !
  ENDIF
  !
  !
  DO i = 1, occ(ik)
    DO j = 1, occ(ik)
      !
      ! <A|W|A>                                 
      wmat(1,1,ik) = wmat(1,1,ik) + single_electron_w(i,j,1,1,ik) * cofc(i,j,1,1,ik)
      !
      ! <B|Wa|A>                                 
      wmat(1,2,ik) = wmat(1,2,ik) + single_electron_w(i,j,1,2,ik) * cofc(i,j,1,2,ik)
      !
      ! <A|Wb|B>                                 
      wmat(2,1,ik) = wmat(2,1,ik) + single_electron_w(i,j,2,1,ik) * cofc(i,j,2,1,ik)
      !
      ! <B|W|B>                     
      wmat(2,2,ik) = wmat(2,2,ik) + single_electron_w(i,j,2,2,ik) * cofc(i,j,2,2,ik)
      !
    ENDDO !j
  ENDDO !i
  !
  !
END SUBROUTINE calc_w_real
!
!
!-----------------------------------------------------------------------------
SUBROUTINE calc_w_img(ik, occ, wevc, wevc2, cofc)
  !-----------------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE klist,                ONLY : nks
  USE io_global,            ONLY : ionode, stdout
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npwx, npw
  USE epcdft_mod,           ONLY : evc1, evc2, occup1, occdown1, wmat, w, smat, debug2
  USE fft_base,             ONLY : dfftp
  USE mp,                   ONLY : mp_sum
  USE mp_images,            ONLY : intra_image_comm
  USE becmod,               ONLY : calbec
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  !
  INTEGER :: i, j, filunit
  INTEGER, INTENT(IN) :: ik
  INTEGER, INTENT(IN) :: occ(2)                  ! number of occupied states for that spin
  COMPLEX(DP), INTENT(IN) :: wevc(npwx, nbnd)    ! wtot*evc
  COMPLEX(DP), INTENT(IN) :: wevc2(npwx, nbnd)   ! wtot*evc2
  COMPLEX(DP), INTENT(IN) :: cofc(nbnd,nbnd,2,2,2) ! cofactor ( i, j, row a b, colm a b , spin up down )
  COMPLEX(DP) :: single_electron_w(nbnd,nbnd,2,2,2)  ! <phi_i|w|phi_j> 
  !                                                  ! single electron orbitals (i,j,row a b, colm a b, spin up down)
  CHARACTER(LEN=256) :: fname
  !
  single_electron_w = 0.D0
  !
  CALL calbec ( npw, evc1(:,:,ik), wevc,  single_electron_w(:,:,1,1,ik), nbnd )
  CALL calbec ( npw, evc1(:,:,ik), wevc2, single_electron_w(:,:,1,2,ik), nbnd )
  CALL calbec ( npw, evc2(:,:,ik), wevc,  single_electron_w(:,:,2,1,ik), nbnd )
  CALL calbec ( npw, evc2(:,:,ik), wevc2, single_electron_w(:,:,2,2,ik), nbnd )
  !
  ! write w single electron to file
  !
  IF(debug2)THEN
    !
    WRITE(unit=fname,fmt=*) ik
    fname="Waa"//TRIM(ADJUSTL(fname))
    CALL realpart_dumpmat(fname, filunit, single_electron_w(:,:,1,1,ik), nbnd, occ(ik))
    WRITE(unit=fname,fmt=*) ik
    fname="Wab"//TRIM(ADJUSTL(fname))
    CALL realpart_dumpmat(fname, filunit, single_electron_w(:,:,1,2,ik), nbnd, occ(ik))
    WRITE(unit=fname,fmt=*) ik
    fname="Wba"//TRIM(ADJUSTL(fname))
    CALL realpart_dumpmat(fname, filunit, single_electron_w(:,:,2,1,ik), nbnd, occ(ik))
    WRITE(unit=fname,fmt=*) ik
    fname="Wbb"//TRIM(ADJUSTL(fname))
    CALL realpart_dumpmat(fname, filunit, single_electron_w(:,:,2,2,ik), nbnd, occ(ik))
    !
  ENDIF
  !
  !
  DO i = 1, occ(ik)
    DO j = 1, occ(ik)
      !
      ! <A|W|A>                                 
      wmat(1,1,ik) = wmat(1,1,ik) + single_electron_w(i,j,1,1,ik) * cofc(i,j,1,1,ik)
      !
      ! <B|Wa|A>                                 
      wmat(1,2,ik) = wmat(1,2,ik) + single_electron_w(i,j,1,2,ik) * cofc(i,j,1,2,ik)
      !
      ! <A|Wb|B>                                 
      wmat(2,1,ik) = wmat(2,1,ik) + single_electron_w(i,j,2,1,ik) * cofc(i,j,2,1,ik)
      !
      ! <B|W|B>                     
      wmat(2,2,ik) = wmat(2,2,ik) + single_electron_w(i,j,2,2,ik) * cofc(i,j,2,2,ik)
      !
    ENDDO !j
  ENDDO !i
  !
  !
END SUBROUTINE calc_w_img
!
!
!-----------------------------------------------------------------------
SUBROUTINE real_dumpmat(fname,filunit,mat,m,n)
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
  REAL(DP),INTENT(IN) :: mat(m,m)
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
END SUBROUTINE real_dumpmat
