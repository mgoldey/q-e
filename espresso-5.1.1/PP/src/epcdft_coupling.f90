!
! Copyright (C) 2013 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! -----------------------------------------------------------------
! This program reads the prefix1.wfc and prefix2.wfc in G-space written by QE and 
! and computes S_1,2. Where S_1,2 is the overlap of the slater det's of system 1 and 2.
! 
!
! Nicholas Brawand nicholasbrawand@gmail.com
!
! Notes:
!        1 = (2-delta_{NumOfSpin,2}) * 1/TotalGridPnts * <phi|phi>
! 
!-----------------------------------------------------------------------
PROGRAM epcdft_coupling 
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_files,             ONLY : prefix, tmp_dir, diropn
  USE mp_global,            ONLY : npool, mp_startup,  intra_image_comm
  USE wvfct,                ONLY : nbnd, npwx, igk, npw , g2kin, ecutwfc
  USE klist,                ONLY : xk     , nks
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast, mp_barrier
  USE mp_world,             ONLY : world_comm
  USE wavefunctions_module, ONLY : evc
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE gvect,                ONLY : ngm, g 
  USE gvecs,                ONLY : nls
  USE noncollin_module,     ONLY : npol, nspin_mag, noncolin
  USE cell_base,            ONLY : tpiba2, omega
  USE environment,          ONLY : environment_start, environment_end
  USE fft_base,             ONLY : dffts, cgather_smooth
  USE fft_interfaces,       ONLY : invfft
  !
  IMPLICIT NONE
  CHARACTER (len=256)          :: outdir
  CHARACTER (len=256)          :: outdir2
  CHARACTER (len=256)          :: tmp_dir2
  CHARACTER (len=256)          :: prefix2
  CHARACTER(LEN=256), external :: trimcheck
  CHARACTER(len=256)           :: filename
  INTEGER                      :: iunitout,ios,ik,i,iuwfcr,lrwfcr,ibnd, ig, is
  LOGICAL                      :: exst
  COMPLEX(DP), ALLOCATABLE     :: evc_r(:,:), dist_evc_r(:,:)
  REAL(DP)                     :: dtmp  ! temp var
  REAL(DP)                     :: dv    ! for ints over realspace 
  !
  NAMELIST / inputpp / outdir, prefix, prefix2, outdir2
  !
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'epcdft_coupling' )
  !
  ! init input variables
  !
  prefix = 'pwscf'
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  !
  IF ( npool > 1 ) CALL errore('bands','pools not implemented',npool)
  !
  ! Read input file
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     ! 
     READ (5, inputpp, err = 200, iostat = ios)
200  CALL errore ('epcdft_coupling', 'reading inputpp namelist', ABS (ios) )
     !
     tmp_dir = trimcheck (outdir)
     tmp_dir2 = trimcheck (outdir2)
     ! 
  END IF
  !
  ! ... Broadcast variables
  !
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( tmp_dir2, ionode_id, world_comm )
  CALL mp_bcast( prefix2, ionode_id, world_comm )
!WRITE(*,*)"dir1 prefix1 dir2 prefix2",tmp_dir," ",prefix," ",tmp_dir2," ",prefix2," "
  !
  !   Now allocate space for pwscf variables, read and check them.
  !
  CALL read_file
  CALL openfil_pp
  dtmp = 0.0
  dv = omega / DBLE( dffts%nr1*dffts%nr2*dffts%nr3 )
  exst=.false.
  filename='wfc_r'
  write(6,*) 'filename=',filename
  iuwfcr=877
  lrwfcr = 2 * dffts%nr1x*dffts%nr2x*dffts%nr3x * npol
  ! lrwfc = 2 * nbnd * npwx * npol
  write(6,*) dffts%nnr, npwx
  write(6,*) 'length of wfc in real space/per band', nks*lrwfcr*8
  write(6,*) 'length of wfc in k space', 2*nbnd*npwx*nks*8
  CALL init_us_1
  !
  !define lrwfcr
  !
  IF (ionode) CALL diropn (iuwfcr, filename, lrwfcr, exst)
  ALLOCATE ( evc_r(dffts%nnr,npol) )
  ALLOCATE ( dist_evc_r(dffts%nr1x*dffts%nr2x*dffts%nr3x,nspin_mag) )
  !
  ! 
  !
  DO ik = 1,nks
     !
     !    prepare the indices of this k point
     !
     CALL gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, &
          igk, g2kin)
     !
     CALL davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
     DO ibnd=1,nbnd 
        !
        ! Fourier transform to realspace
        !
        evc_r = cmplx(0.d0, 0.d0)     
        do ig = 1, npw
           evc_r( nls( igk(ig) ), 1 ) = evc(ig,ibnd)
        enddo
        CALL invfft ('Wave', evc_r(:,1), dffts)
        IF (noncolin) THEN
           DO ig = 1, npw
              evc_r (nls(igk(ig)),2) = evc (ig+npwx, ibnd)
           ENDDO
           CALL invfft ('Wave', evc_r(:,2), dffts)
        ENDIF
        !
        dist_evc_r=CMPLX(0.d0,0.d0)
        !
#if defined (__MPI)
        DO is = 1, nspin_mag
           !
           CALL cgather_smooth( evc_r(:,is), dist_evc_r(:,is) )
           !
        END DO
#else
        dist_evc_r(1:dffts%nnr,:)=evc_r(1:dffts%nnr,:)
#endif
        !
        ! overlap check
        !
        IF(ionode) THEN
           !  
           dtmp = 0.0
           DO i = 1, dffts%nr1x*dffts%nr2x*dffts%nr3x
              dtmp = dtmp + ( dist_evc_r(i,1) * CONJG(dist_evc_r(i,1)) )
           ENDDO
           dtmp = DOT_PRODUCT( dist_evc_r(:,1) , dist_evc_r(:,1) )
           WRITE(*,*)"overlap of k = ",ik," and ibnd =",ibnd," is = ",2*dtmp*dv/omega
           !
        ENDIF
        !
     ENDDO
     !
  ENDDO
  !
  if(ionode) close(iuwfcr)
  DEALLOCATE (evc_r)
  !
  CALL environment_end ( 'epcdft_coupling' )
  !
  CALL stop_pp
  !
  STOP
  !
END PROGRAM epcdft_coupling
