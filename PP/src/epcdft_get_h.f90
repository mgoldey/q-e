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
SUBROUTINE epcdft_get_h
  !-----------------------------------------------------------------------
  !
  ! Here we calculate H and take the average of the off diags of H.
  !
  ! hc(1,1,s) = ecor1
  ! hc(1,2,s) = 0.5 ( (F1+F2) Sab - (Wab + Wba) )
  ! hc(2,1,s) = 0.5 ( (F2+F1) Sba - (Wba + Wab) )
  ! hc(2,2,s) = ecor2
  !
  ! Sab = Sab_up*Sab_down
  !
  ! Wab = Wab_up*det(Sab_down) + Wab_down*det(Sab_up)
  !
  USE kinds,      ONLY : DP
  USE io_global,  ONLY : ionode, stdout
  USE epcdft_mod, ONLY : free1, free2, wmat, smat, hc, cor1, cor2, debug2
  USE klist,      ONLY : nks
  !
  IMPLICIT NONE
  !
  INTEGER i, j, s
  REAL(DP) :: core(2) ! <|H_KS|> = Free - V int rho(r) w(r)
  REAL(DP) :: ftot ! F1+F2, Sup*Sdown, W_up+W_down
  COMPLEX(DP) :: stot(2,2), wtot(2,2) ! F1+F2, Sup*Sdown, W_up+W_down
  !
  hc = 0.D0
  core(1) = free1 - cor1 
  core(2) = free2 - cor2
  stot(:,:) = smat(:,:,1) * smat(:,:,2)
  !
  ! Wtot = Wab + Wba
  !
  DO i = 1, 2
    DO j = 1, 2
      wtot(i,j) =  wmat(i,j,1) * smat(i,j,2) + wmat(i,j,2) * smat(i,j,1)
    ENDDO
  ENDDO
  !
  hc(1,1) = core(1)
  hc(1,2) = ( free2 * stot(1,2) - wtot(1,2) ) !<a|h|b>
  hc(2,1) = ( free1 * stot(2,1) - wtot(2,1) ) !<b|h|a>
  hc(2,2) = core(2)
  !
  ! take average of off diags
  !
  hc(1,2) = 0.5D0*( hc(1,2) + hc(2,1) )
  hc(2,1) = hc(1,2) 
  !
  IF( debug2 ) CALL epcdft_print_h(hc)
  IF( ionode ) WRITE( stdout, * )"    H done"
  !
END SUBROUTINE epcdft_get_h
!
!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_print_h(hc)
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE epcdft_mod, ONLY : free1, free2, cor1, cor2
  !
  IMPLICIT NONE
  !
  REAL(DP) :: core(2,2)
  COMPLEX(DP), INTENT(IN) :: hc(2,2)
  CHARACTER(LEN=256) :: fname
  INTEGER :: filunit = 3240903
  !
  core(1,1) = free1
  core(1,2) = cor1 
  core(2,1) = free2 
  core(2,2) = cor2
  !
  fname="FandC"
  CALL real_dumpmat(fname,filunit,core,2,2)
  !
  fname="H"
  CALL realpart_dumpmat(fname,filunit,hc,2,2)
  !
END SUBROUTINE epcdft_print_h
