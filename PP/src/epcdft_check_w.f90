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
SUBROUTINE epcdft_check_w
  !-----------------------------------------------------------------------
  !
  !  Check W matrix
  !
  !    <A|W|A> = <A|Wa|A> = N \sum_{i,j}   \langle i | W | j \rangle  (S^{-1} (det(S)I))^T_{i,j} 
  !            = correction energy
  !                                                                     ^-----cofactor----------^
  !
  USE kinds,                ONLY : DP
  USE klist,                ONLY : nks
  USE io_global,            ONLY : ionode, stdout
  USE wavefunctions_module, ONLY : evc
  USE wvfct,                ONLY : nbnd, npwx
  USE epcdft_mod,           ONLY : evc1, evc2, occup1, occdown1, smat, w
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
  COMPLEX(DP) :: wmat, wmat2
  !
  occ(1) = occup1
  occ(2) = occdown1
  wmat = ( 0.D0, 0.D0 )
  wmat2 = ( 0.D0, 0.D0 )
  !
  ! create S matrix
  !
  DO ik = 1 , nks 
    !
    ! S^-1  ab 
    !
    CALL get_s_invs(evc1(:,:,ik), evc1(:,:,ik), cofc(:,:,1,ik), nbnd)
    CALL get_s_invs(evc2(:,:,ik), evc2(:,:,ik), cofc(:,:,2,ik), nbnd)
    !
    ! S^-1 * det(S) ab 
    !
    cofc(:,:,1,ik) = cofc(:,:,1,ik) * smat(1,1,ik)
    cofc(:,:,2,ik) = cofc(:,:,2,ik) * smat(2,2,ik)
    !
    ! C = ( S^-1 * det(S))^T   ab & ba 
    !
    cofc(:,:,2,ik) = TRANSPOSE(cofc(:,:,2,ik))
    cofc(:,:,1,ik) = TRANSPOSE(cofc(:,:,1,ik))
    !
    ! wevc = w*evc
    !
    !wtot(:) = -1.D0*w(:,1)
    wtot(:) = w(:,1,ik)
    CALL w_psi(evc1(:,:,ik), wtot, wevc) 
    !
    !wtot(:) = -1.D0*w(:,2)
    wtot(:) = w(:,2,ik)
    CALL w_psi(evc2(:,:,ik), wtot, wevc2) 
    !
    DO i = 1, occ(ik)
      DO j = 1, occ(ik)
        !
        ! <A|W|A>                                 
        wmat = wmat + dot(evc1(:,i,ik), wevc(:,j)) * cofc(i,j,1,ik)
        !
        ! <B|W|B>                                 
        wmat2 = wmat2 + dot(evc2(:,i,ik), wevc2(:,j)) * cofc(i,j,2,ik)
        !
      ENDDO !j
    ENDDO !i
    !
  ENDDO !ik
  !
  
  CALL mp_sum(wmat,intra_image_comm)
  CALL mp_sum(wmat2,intra_image_comm)
  !
  IF( ionode ) WRITE( stdout,* ) "    Check W1 and W2 (should be C1 and C2)"
  IF( ionode ) WRITE( stdout,* ) "    W1",wmat
  IF( ionode ) WRITE( stdout,* ) "    W2",wmat2
  IF( ionode ) WRITE( stdout,* ) "    END check W"
  !
END SUBROUTINE epcdft_check_w
