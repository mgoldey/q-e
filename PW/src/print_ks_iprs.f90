!
! Copyright (C) 2007-2015 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!  
!
!----------------------------------------------------------------------------
SUBROUTINE print_ks_iprs()
  !----------------------------------------------------------------------------
  !
  ! ... printout of Kohn-Sham inverse participation ratios
  !
  USE kinds,                ONLY : DP
  USE constants,            ONLY : rytoev
  USE io_global,            ONLY : stdout
  USE klist,                ONLY : xk, ngk, nks, nkstot, wk, lgauss, ltetra, &
                                   two_fermi_energies
  USE fixed_occ,            ONLY : one_atom_occupations
  USE ener,                 ONLY : ef, ef_up, ef_dw 
  USE lsda_mod,             ONLY : lsda, nspin
  USE spin_orb,             ONLY : lforcet
  USE wvfct,                ONLY : nbnd, et, wg
  USE control_flags,        ONLY : conv_elec, lbands, iverbosity
  USE mp_bands,             ONLY : root_bgrp, intra_bgrp_comm, inter_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  USE mp_pools,             ONLY : inter_pool_comm 
  USE mp_world,             only : mpime,nproc
  !
  IMPLICIT NONE
  !
  ! ... a few local variables
  !  
  INTEGER, ALLOCATABLE :: &
      ngk_g(:)       ! number of plane waves summed on all nodes
  REAL(DP), DIMENSION(nbnd,nkstot) :: &
      ipr ! IPR
  INTEGER :: &
      i,            &! counter on polarization
      ik,           &! counter on k points
      ibnd           ! counter on bands
  !
     !
     ALLOCATE ( ngk_g (nkstot) ) 
     !
     ngk_g(1:nks) = ngk(:)
     ipr=0d0
!      write(*,*) "Test inside ipr printer", mpime, nproc
!     CALL mp_sum( ngk_g(1:nks), intra_bgrp_comm )
!     CALL ipoolrecover( ngk_g, 1, nkstot, nks )
!     CALL mp_bcast( ngk_g, root_bgrp, intra_bgrp_comm )
!     CALL mp_bcast( ngk_g, root_bgrp, inter_bgrp_comm )
     !
     DO ik = 1, nkstot
      !
      IF ( lsda ) THEN
         !
         IF ( ik == 1 ) WRITE( stdout, 9015)
         IF ( ik == ( 1 + nkstot / 2 ) ) WRITE( stdout, 9016)
         !
      END IF
      !
      IF ( conv_elec ) THEN
         WRITE( stdout, 9021 ) ( xk(i,ik), i = 1, 3 ), ngk_g(ik)
      ELSE
         WRITE( stdout, 9020 ) ( xk(i,ik), i = 1, 3 )
      END IF
      !

      do ibnd = 1, nbnd, 1
        CALL calculate_ipr(ik,ibnd,ipr(ibnd,ik))                
      end do        

      WRITE( stdout, 9030 ) (ipr(ibnd,ik),ibnd=1,nbnd)
      
      !
      IF( iverbosity > 0 .AND. .NOT. lbands ) THEN
         !
         WRITE( stdout, 9032 )
         IF (ABS(wk(ik))>1.d-10) THEN
            WRITE( stdout, 9030 ) ( wg(ibnd,ik)/wk(ik), ibnd = 1, nbnd )
         ELSE
            WRITE( stdout, 9030 ) ( wg(ibnd,ik), ibnd = 1, nbnd )
         ENDIF
         !
      END IF
      !
   END DO
   !
   DEALLOCATE ( ngk_g )
   !
  !
  FLUSH( stdout )
  RETURN
  !
  ! ... formats
  !

9015 FORMAT(/' ------ SPIN UP ------------'/ )
9016 FORMAT(/' ------ SPIN DOWN ----------'/ )
9020 FORMAT(/'          k =',3F7.4,' Inverse participation ratios (1/Angstrom^3):'/ )
9021 FORMAT(/'          k =',3F7.4,' (',I6,' PWs)   Inverse participation ratios (1/Angstrom^3): '/ )
9030 FORMAT( '   ',8E12.4 )
9032 FORMAT(/'     occupation numbers ' )
!
  !
END SUBROUTINE print_ks_iprs
!
!----------------------------------------------------------------------------

SUBROUTINE calculate_ipr(ik,ibnd,ipr)
  USE mp_bands,             ONLY : root_bgrp,root_bgrp_id, intra_bgrp_comm, inter_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  USE mp_pools,             ONLY : inter_pool_comm
  USE mp_images,            ONLY : intra_image_comm
  USE mp_global,            ONLY : inter_image_comm,my_image_id
  USE mp_world,             ONLY : world_comm
  USE fft_base,             ONLY : dfftp
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega ! in cubic bohr
  USE control_flags,        ONLY : gamma_only
  USE wvfct,                ONLY : npwx, nbnd, current_k
  USE klist,                ONLY : nks,xk,ngk, igk_k
  USE gvect,                ONLY : nl,nlm
  USE gvecs,                ONLY : nls,nlsm
  USE fft_interfaces,       ONLY : fwfft,invfft
  USE fft_types,            ONLY : fft_type_descriptor
  USE wavefunctions_module, ONLY : evc, psic,psic_nc
  USE io_global,            ONLY : stdout
  USE io_files,             ONLY : nwordwfc, iunwfc, iunoldwfc, &
                                   iunoldwfc2, diropn
  USE constants,            ONLY : ANGSTROM_AU
  USE realus,               ONLY : real_space, &
            invfft_orbital_gamma, calbec_rs_gamma, &
            invfft_orbital_k, calbec_rs_k
  USE buffers,               ONLY : get_buffer, save_buffer

  INTEGER, INTENT(IN) :: ik,ibnd
  REAL(DP), INTENT(INOUT) :: ipr
  REAL(DP) :: top,bottom,fac
  INTEGER :: ig, npw,npw_k


  ! Clear work array of psic
  psic(:)=0d0
  ipr=0d0
  top=0d0
  bottom=0d0
  npw = ngk(ik)
  npw_k=ngk(ik)
  current_k=ik

  IF ( nks > 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, ik )
  CALL mp_bcast( evc, root_bgrp_id, inter_bgrp_comm )

  IF (gamma_only) THEN
    CALL invfft_orbital_gamma(evc,ibnd,ibnd,.true.) 
  ELSE
    CALL invfft_orbital_k(evc,ibnd,ibnd,ik,.true.) 
  ENDIF

  psic(:)=CONJG(psic)*psic

  fac=omega/(ANGSTROM_AU**3*(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x))
  top=SUM(ABS(psic(:))*ABS(psic(:)))*fac
  bottom=SUM(ABS(psic(:)))*fac

  
  CALL MP_SUM(top,intra_image_comm)
  CALL MP_SUM(bottom,intra_image_comm)

!   CALL MP_SUM(top,world_comm)
!   CALL MP_SUM(bottom,world_comm)

  ipr=top/(bottom*bottom)

  RETURN
END SUBROUTINE  calculate_ipr