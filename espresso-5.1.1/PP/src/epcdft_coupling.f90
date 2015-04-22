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
  INTEGER                      :: ios,ik,i,ibnd, ig, is
  REAL(DP)                     :: dtmp  ! temp variable
  COMPLEX(DP)                  :: ztmp  ! temp variable
  COMPLEX(DP), EXTERNAL :: zdotc
  !
  NAMELIST / inputpp / outdir, prefix, prefix2, outdir2
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'epcdft_coupling' )
  !
  ! setup vars and consistency checks
  !
  CALL get_env( 'ESPRESSO_TMPDIR', outdir ) ! set env variable
  IF ( TRIM( outdir ) == ' ' ) outdir = './' ! remove spaces from end
  IF ( TRIM( outdir2 ) == ' ' ) outdir2 = './'
  !
  IF ( npool > 1 ) CALL errore('bands','pools not implemented',npool)
  !
  ! Read input file
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file() ! check command line for input file and attach to 5 
     ! 
     READ (5, inputpp, err = 200, iostat = ios)
200  CALL errore ('epcdft_coupling', 'reading inputpp namelist', ABS (ios) )
     !
     tmp_dir  = trimcheck (outdir) ! make sure dir ends with '/' and remove end space
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
  ! read & construct everything from outfile geometry, wfc, eigenvals, potential...
  CALL read_file() 
  !
  CALL openfil_pp() ! open all files for scf run set filenames/units
  !
  CALL init_us_1 ! compute pseduo pot stuff
  !
  dtmp = 0.0
  ztmp = 0.0
  !
  DO ik = 1,nks
     !
     ! prepare the indices of this k point
     CALL gk_sort( xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
     !
     ! read phi_ik
     CALL davcio( evc, 2*nwordwfc, iunwfc, ik, -1 ) 
     !
     DO ibnd=1,nbnd 
        !
        dtmp = REAL(zdotc(npw, evc(:,ibnd), 1, evc(:,ibnd), 1 ), DP)
        WRITE(*,*)"Overlap of k = ",ik," and ibnd =",ibnd," is = ",2.0*dtmp
        !
     ENDDO
     !
  ENDDO
  !
  CALL environment_end ( 'epcdft_coupling' )
  !
  CALL stop_pp
  !
  STOP
  !
END PROGRAM epcdft_coupling
