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
SUBROUTINE epcdft_print
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE epcdft_mod, ONLY : free1, free2, wmat, smat, hc, ohc, cor1, cor2, debug
  USE klist, ONLY : nks
  USE io_global, ONLY : stdout, ionode
  USE constants,  ONLY : rytoev
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) title
  !
  IF(ionode) THEN
    !
    WRITE(stdout,1)
    WRITE(stdout,*)""
    !
    IF(debug) THEN
      !
      WRITE(stdout,4) free1, free2
      WRITE(stdout,7) free1+cor1, free2+cor2
      WRITE(stdout,8) cor1, cor2
      WRITE(*,*)""
      !
      title='S up' ; CALL pcmat( title, smat(:,:,1) )
      title='S down' ; CALL pcmat( title, smat(:,:,2) )
      !
      title='W up' ; CALL pcmat( title, wmat(:,:,1) )
      title='W down' ; CALL pcmat( title, wmat(:,:,2) )
      !
      title='H' ; CALL pcmat( title, hc(:,:) )
      !
      title='H using lowdin orthogonalization'  ; CALL pcmat( title, ohc(:,:) )
      !
    ENDIF
    !
    ! Sab = S_12_alpha * S_12_beta 
    WRITE(stdout,6) ABS( smat(1,2,1) * smat(1,2,2) )
    !
    ! |Hab| 
    WRITE(stdout,5) ABS( ohc(1,2) )
    WRITE(stdout,9) ABS( ohc(1,2) ) * 0.5D0
    WRITE(stdout,10) ABS( ohc(1,2) ) * rytoev
    !
    WRITE(stdout,*)""
    WRITE(stdout,1)
    WRITE(stdout,*)""
    WRITE(stdout,*)"    Thank you for using CDFT!"
    WRITE(stdout,*)""
    WRITE(stdout,*)"    Please cite: "
    WRITE(stdout,*)"        DOI: 10.1021/acs.chemmater.6b04631 and 10.1021/acs.jctc.7b00088"
    WRITE(stdout,*)"    "
    WRITE(stdout,*)""
    !
  ENDIF
  !
  1 FORMAT('     =======================================================================')
  4 FORMAT(5x,'F1 = ',F12.5,3x,'F2 = ',F12.5)
  7 FORMAT(5x,'E1 = ',F12.5,3x,'E2 = ',F12.5)
  8 FORMAT(5x,'C1 = ',F12.5,3x,'C2 = ',F12.5)
  5 FORMAT(5x,'|Hab| = ',E14.6,' Ry')
  9 FORMAT(5x,'|Hab| = ',E14.6,' Ha')
  10 FORMAT(5x,'|Hab| = ',E14.6,' eV')
  6 FORMAT(5x,'|Sab| = ',E14.6)
  !
END SUBROUTINE epcdft_print
!
!-----------------------------------------------------------------------
SUBROUTINE pcmat (title, m)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE io_global, ONLY : stdout
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256), INTENT(IN) :: title
  COMPLEX(DP), INTENT(IN) :: m(2,2)
  !
  INTEGER :: i, j
  !
  WRITE(*,*)""
  WRITE(stdout,2) title
  WRITE(stdout,3) ( ( m(i,j) , j=1,2 ) , i=1,2)
  WRITE(*,*)""
  !
  2 FORMAT(5x,A100)
  !3 FORMAT(5x,'(',E12.4,',',E12.4,')',2x,'(',E12.4,',',E12.4,')')
  !3 FORMAT(5x,'(',F14.6,',',F14.6,')',2x,'(',F14.6,',',F14.6,')')
  3 FORMAT(5x,'(',e23.16,',',e23.16,')',2x,'(',e23.16,',',e23.16,')')
  !
END SUBROUTINE pcmat
