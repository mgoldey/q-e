!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_match_phases
  !-----------------------------------------------------------------------
  !
  ! two separate runs using gamma point tricks can differ
  ! by a factor of e^{i*\pi}. 
  !
  USE kinds,      ONLY : DP
  USE klist,      ONLY : nks
  USE io_global,  ONLY : ionode, stdout
  USE wvfct,      ONLY : nbnd
  USE epcdft_mod, ONLY : evc1, evc2, occup1, occdown1
  USE mp_images,  ONLY : intra_image_comm
  USE mp,         ONLY : mp_sum
  !
  IMPLICIT NONE
  !
  INTEGER :: ik, i
  INTEGER :: occ(2)                        ! number of occupied states for that spin
  COMPLEX(DP), EXTERNAL :: dot
  REAL(DP), ALLOCATABLE :: olaps(:,:)
  !
  ALLOCATE( olaps(nbnd, nks) )
  !
  ! bounds for spin channels in S matrix
  !
  occ(1) = occup1
  occ(2) = occdown1
  olaps = 0.D0
  !
  !
  DO ik = 1, nks 
    DO i = 1, occ(ik) 
      !
      olaps(i,ik) = dot( evc1(:,i,ik), evc2(:,i,ik) )
      !
    ENDDO
  ENDDO !ik
  !
  CALL mp_sum(olaps, intra_image_comm)
  !
  DO ik = 1, nks 
    DO i = 1, occ(ik) 
      !
      IF( olaps(i,ik) < 0.D0 ) evc2(:,i,ik) = -1.D0*evc2(:,i,ik)  
      !WRITE(*,*) ik, i, olaps(i,ik)
      !
    ENDDO
  ENDDO !ik
  !
  IF( ionode ) WRITE( stdout,*)"    Phases corrected"
  !
  ! close shop
  !
END SUBROUTINE epcdft_match_phases
