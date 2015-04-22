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
!        S( j , i ) matrix = < evc1(i) | evc2(j) > 
!                                 evc1 -->
!                                e
!                                v
!                                c
!                                2
!                                |
!                                v
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
  LOGICAL                      :: exst2
  CHARACTER (len=20)          :: FMT
  CHARACTER (len=256)          :: outdir
  CHARACTER (len=256)          :: outdir2
  CHARACTER (len=256)          :: tmp_dir2
  CHARACTER (len=256)          :: prefix2
  CHARACTER(LEN=256), external :: trimcheck
  INTEGER                      :: ios,ik,i,j,ibnd, ig, is
  INTEGER                      :: ik1, ik2, ibnd1, ibnd2
  INTEGER                      :: iunwfc2 = 3636 ! unit for 2nd set of wfcs
  REAL(DP)                     :: dtmp  ! temp variable
  COMPLEX(DP)                  :: ztmp  ! temp variable
  COMPLEX(DP), EXTERNAL        :: zdotc
  COMPLEX(DP), ALLOCATABLE     :: evc2(:,:)
  COMPLEX(DP),    ALLOCATABLE  :: smat(:,:)  ! S_ij matrix <wfc_j|wfc2_i>
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
  CALL read_file() ! read , wfc, eigenvals, potential
  !
  CALL openfil_pp() ! open all files for scf run set filenames/units
  !
  CALL diropn_gw ( iunwfc2, tmp_dir2, trim( prefix2 )//'.wfc', 2*nwordwfc, exst2, 1 )! diropn with more control
  !WRITE(*,*)"Found file with 2nd set of wfcs: ",exst2
  !
  CALL init_us_1 ! compute pseduo pot stuff
  !
  i = 0
  j = 0
  dtmp = 0.0
  ztmp = 0.0
  ALLOCATE( evc2( npwx, nbnd ) )
  ALLOCATE( smat( nks*nbnd, nks*nbnd ) )
  smat = 0.0
  !WRITE(*,*)"Size of ecv dim1 ",SIZE(evc,1)," size of evc2 dim1 ",SIZE(evc2,1)
  !WRITE(*,*)"Size of ecv dim2 ",SIZE(evc,2)," size of evc2 dim2 ",SIZE(evc2,2)
  !
  DO ik1 = 1, nks
     !
     ! prepare the indices & read evc1
     CALL gk_sort( xk(1,ik1), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
     CALL davcio( evc, 2*nwordwfc, iunwfc, ik1, -1 ) 
     !
     DO ibnd1 = 1, nbnd
        !
        i = i + 1 ! S matrix counter
        !
        DO ik2 = 1, nks
           !
           ! prepare the indices & read evc2
           CALL gk_sort( xk(1,ik2), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
           CALL davcio( evc2, 2*nwordwfc, iunwfc2, ik2, -1 ) 
           !
           DO ibnd2= 1, nbnd 
              !
              j = j + 1 ! S matrix counter
              !
              ! evc is being conjg
              smat(i, j) = zdotc(npw, evc(:,ibnd1), 1, evc2(:,ibnd2), 1 )
              !
           ENDDO ! end ibnd2
           !
        ENDDO ! end ik2
        !
        j = 0 ! S matrix inner counter reset
        !
     ENDDO ! end ibnd1
     !
  ENDDO ! end ik1
  !
  ! for spin = 1
  IF(.NOT.noncolin) smat = 2.0 * smat
  !
  ! Print Results
  !
  WRITE(*,*)""
  WRITE(*,*)"  REAL S_ij"
  DO j = 1, nbnd*nks
    WRITE(*,"(4F8.3)")REAL(smat(j,:))
  ENDDO
  !
  WRITE(*,*)""
  WRITE(*,*)"  IMG S_ij"
  DO j = 1, nbnd*nks
    WRITE(*,"(4F8.3)")AIMAG(smat(j,:))
  ENDDO
  !
  CALL environment_end ( 'epcdft_coupling' )
  !
  CALL stop_pp
  !
  STOP
  !
  !
END PROGRAM epcdft_coupling
!
!
!
!-----------------------------------------------------------------------
SUBROUTINE diropn_gw (unit, tmp_dir2, filename, recl, exst, ndnmbr )
  !-----------------------------------------------------------------------
  !
  !     this routine (taken from PP/src/pw2gw.f90 and modded) opens a file in tmp_dir2 for direct I/O access
  !     If appropriate, the node number is added to the file name
  !
  USE kinds
  USE io_files
  IMPLICIT NONE

  !
  !    first the input variables
  !
  CHARACTER(len=*) :: filename
  ! input: name of the file to open
  INTEGER :: unit, recl
  ! input: unit of the file to open
  ! input: length of the records
  LOGICAL :: exst
  ! output: if true the file exists
  !
  INTEGER :: ndnmbr
  !
  !    local variables
  !
  CHARACTER(len=256) :: tempfile
  ! complete file name
  CHARACTER(len=80) :: assstr
  INTEGER :: ios, unf_recl, ierr
  ! used to check I/O operations
  ! length of the record
  ! error code
  LOGICAL :: opnd
  ! if true the file is already opened
  CHARACTER (len=256), INTENT(IN)          :: tmp_dir2
  CHARACTER(len=256) :: strnum
  !
  !WRITE( strnum, '(i10)' ) ndnmbr
  IF(ndnmbr>-1 .and. ndnmbr<10)   WRITE( strnum, '(I1)' ) ndnmbr
  IF(ndnmbr>9 .and. ndnmbr<100)   WRITE( strnum, '(I2)' ) ndnmbr
  IF(ndnmbr>99 .and. ndnmbr<1000) WRITE( strnum, '(I3)' ) ndnmbr
  !
  IF (unit < 0) CALL errore ('diropn', 'wrong unit', 1)
  !
  !    we first check that the file is not already openend
  !
  ios = 0
  INQUIRE (unit = unit, opened = opnd)
  IF (opnd) CALL errore ('diropn', "can't open a connected unit", abs(unit))
  !
  !      then we check the filename
  !

  IF (filename == ' ') CALL errore ('diropn', 'filename not given', 2)
  tempfile = trim(tmp_dir2) // trim(filename) // trim(strnum)

  INQUIRE (file = tempfile, exist = exst)
  !
  !      the unit for record length is unfortunately machine-dependent
  !
#define DIRECT_IO_FACTOR 8
  unf_recl = DIRECT_IO_FACTOR * recl
  IF (unf_recl <= 0) CALL errore ('diropn', 'wrong record length', 3)
  !
  OPEN ( unit, file = trim(tempfile), iostat = ios, form = 'unformatted', &
       status = 'unknown', access = 'direct', recl = unf_recl )

  IF (ios /= 0) CALL errore ('diropn', 'error opening '//filename, unit)
  RETURN
END SUBROUTINE diropn_gw
