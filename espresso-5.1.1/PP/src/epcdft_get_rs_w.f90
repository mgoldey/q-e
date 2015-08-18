!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_rs_w
  !-----------------------------------------------------------------------
  !
  !  <A|W|B> = <A|VA+VB|B> = N \sum_{i,j}   \langle i | W | j \rangle  (S^{-1} (det(S)I))^T_{i,j} 
  !-----------------------------------------------------------------------
  !
  USE kinds,      ONLY : DP
  USE klist,      ONLY : nks
  USE io_global,  ONLY : ionode, stdout
  USE wvfct,      ONLY : nbnd
  USE lsda_mod,   only: nspin
  USE epcdft_mod, ONLY : evc1, evc2, occup1, occdown1, smat,w
  USE fft_base,             ONLY : dfftp
  USE control_flags, ONLY : gamma_only
  USE wvfct, ONLY : npwx, npw
  USE becmod, ONLY : calbec
  USE mp,                   ONLY : mp_sum
  USE mp_images,            ONLY : intra_image_comm
  !
  IMPLICIT NONE
  !
  INTEGER :: ik,i,j
  INTEGER :: occ(2)                        ! number of occupied states for that spin
  COMPLEX(DP), ALLOCATABLE :: c_s_aux(:,:) ! complex overlap matrix S_ij = <system 1_i | system 2_j>
  REAL(DP), ALLOCATABLE :: r_s_aux(:,:)    ! real overlap matrix S_ij = <system 1_i | system 2_j>
  COMPLEX(DP) :: wevc1(npwx, nbnd) ! wtot*evc1
  COMPLEX(DP) :: wevc2(npwx, nbnd) ! wtot*evc2
  COMPLEX(DP), EXTERNAL :: dot
  REAL(DP) :: wtot(dfftp%nnr) ! tot W  w(:,1) + w(:,2)
  REAL(DP) :: energy_offset
  !
  ALLOCATE( c_s_aux(nbnd,nbnd) )
  ALLOCATE( r_s_aux(nbnd,nbnd) )
  !
  ! bounds for spin channels in S matrix
  !
  occ(1) = occup1
  occ(2) = occdown1
  wtot(:) = w(:,1) + w(:,2)
  !
  r_s_aux = 0.D0
  c_s_aux = 0.D0
  energy_offset =0.D0

  
  !
  DO ik = 1 , nks 
  !
    ! create S matrix
    CALL get_s_invs(evc1(:,:,ik), evc2(:,:,ik), c_s_aux, nbnd)
    CALL w_psi(evc1(:,:,ik), wtot, wevc1) 
    !CALL w_psi(evc2(:,:,ik), wtot, wevc2) 
    DO i = 1, occ(ik)
      DO j = 1, occ(ik)
      !
      ! <B|W|A>                                 
      energy_offset=energy_offset + dot(evc2(:,i,ik), wevc1(:,j)) * c_s_aux(i,j)
      !
      ENDDO
    ENDDO
    !
  ENDDO !ik
  !
  CALL mp_sum(energy_offset,intra_image_comm)
  IF (ionode) WRITE(*,*)"    W is ",energy_offset

  !
  IF (ionode) WRITE(*,*)"    W done in real space "
  !
  ! close shop
  !
END SUBROUTINE epcdft_get_rs_w