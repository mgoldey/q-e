!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_setup
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_files,             ONLY : prefix, tmp_dir, diropn
  USE mp_global,            ONLY : npool  
  USE wvfct,                ONLY : nbnd, npwx
  USE klist,                ONLY : nks
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE fft_base,             ONLY : dfftp,  grid_gather
  USE scf,                  ONLY : rho
  USE epcdft_mod  
  USE epcdft
  !
  IMPLICIT NONE
  !
  LOGICAL :: exst2 
  CHARACTER (len=256) :: tmp_dir2 ! temp variable to store system 2's dir info 
  CHARACTER (len=256) :: tmp_dir_pass   ! used to store tmp_dir during pass for reading two systems
  INTEGER  :: iunwfc_pass    ! same as above but diff var
  CHARACTER (len=256) :: prefix_pass
  CHARACTER(LEN=256), external :: trimcheck
  INTEGER :: ios
  REAL(DP) :: dtmp ! temp variable
  REAL(DP), ALLOCATABLE :: aux(:)
  !
  NAMELIST / inputpp / outdir, prefix, prefix2, outdir2, occup1, occup2, occdown1, occdown2, &
                       debug,  s_spin, det_by_zgedi, do_epcdft, fragment1_atom1, fragment1_atom2,&
                       fragment2_atom1, fragment2_atom2, fragment1_amp, fragment2_amp
  !
  ! setup vars and consistency checks
  !
  CALL get_env( 'ESPRESSO_TMPDIR', outdir )  ! set env variable
  IF ( TRIM( outdir ) == ' ' ) outdir = './' ! remove spaces from end
  IF ( TRIM( outdir2 ) == ' ' ) outdir2 = './'
  !
  IF ( npool > 1 ) CALL errore('bands','pools not implemented',npool)
  !
  ! Read input file
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file()        ! check command line for input file and attach to 5 
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
  !
  ! first read system 2 and store in system 1's variables 
  ! then we will restore vars to right place and read sys 1's stuff
  !
  tmp_dir_pass = tmp_dir ! store system 1's dir for later reading
  iunwfc_pass = iunwfc    
  prefix_pass = prefix
  !
  tmp_dir = tmp_dir2     ! the great exchange 
  iunwfc = iunwfc2
  prefix = prefix2
  !
  WRITE(*,*)"    ======================================================================= "
  WRITE(*,*) "    SYSTEM 2 INFO"
  do_epcdft=.false.
  CALL read_file()  ! for system 2
  !
  ! deallocate to avoid reallocation of sys 1 vars
  CALL clean_pw( .TRUE. )
  !
  ! re open the wfc file for system 2
  CALL diropn_gw ( iunwfc2, tmp_dir2, trim( prefix2 )//'.wfc', 2*nwordwfc, exst2, 1 ) ! below 
  !
  ! restore sys 1's vars and read sys 1's data
   WRITE(*,*)"    ======================================================================= "
  WRITE(*,*) "    SYSTEM 1 INFO"
  tmp_dir = tmp_dir_pass
  iunwfc = iunwfc_pass
  prefix = prefix_pass
  CALL read_file()  
   WRITE(*,*)"    ======================================================================= "
  !
  CALL openfil_pp() ! open all files for scf run set filenames/units
  !
  CALL init_us_1    ! compute pseduo pot stuff
  !
  ALLOCATE( evc2 ( npwx, nbnd ) )
  ALLOCATE( smat ( 2 , 2 , nks) )
  ALLOCATE( wmat ( 2 , 2 , nks) )
  ALLOCATE( w ( dfftp%nr1x*dfftp%nr2x*dfftp%nr3x , nks ) )
  ALLOCATE( aux( dfftp%nnr ) )
  !
  evc2 = 0.d0
  aux = 0.d0
  w = 0.d0
  smat = 0.d0
  w = 0.D0
  !
  CALL print_checks_warns(prefix, tmp_dir, prefix2, tmp_dir2, nks, nbnd, &
                          occup1, occdown1, occup2, occdown2, debug,  s_spin, det_by_zgedi )
  !
  ! setup weight function for system 1
  !
  do_epcdft=.true.
  ! 
  aux = 0.D0
  fragment_atom1=fragment1_atom1
  fragment_atom2=fragment1_atom2
  epcdft_amp=fragment1_amp
  CALL add_epcdft_efield( aux, dtmp, rho%of_r, .true. )
  !
#ifdef __MPI
    CALL grid_gather ( aux, w(:,1))
#else
    w(:,1)= aux(:)
#endif
  !
  ! setup weight function for system 1
  !
  aux = 0.D0
  fragment_atom1=fragment2_atom1
  fragment_atom2=fragment2_atom2
  epcdft_amp=fragment2_amp
  CALL add_epcdft_efield( aux, dtmp, rho%of_r, .true. )
  !
#ifdef __MPI
    CALL grid_gather ( aux, w(:,2))
#else
    w(:,2)= aux(:)
#endif
  !
END SUBROUTINE epcdft_setup
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
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------------
SUBROUTINE print_checks_warns(prefix, tmp_dir, prefix2, tmp_dir2, nks, nbnd, occup1, &
                              occdown1, occup2, occdown2, debug,  s_spin, det_by_zgedi )
  !--------------------------------------------------------------------------
  !
  !     this routine prints warnings and some data from the input file 
  !     
  !
  USE io_global,            ONLY : ionode
  !
  IMPLICIT NONE
  !
  INTEGER,      INTENT(IN)    :: nks
  INTEGER,      INTENT(IN)    :: nbnd
  CHARACTER(*), INTENT(IN)    :: prefix
  CHARACTER(*), INTENT(IN)    :: prefix2
  CHARACTER(*), INTENT(IN)    :: tmp_dir
  CHARACTER(*), INTENT(IN)    :: tmp_dir2
  INTEGER,      INTENT(IN)    :: occup1, occdown1, occup2, occdown2
  LOGICAL,      INTENT(IN)    :: debug, s_spin, det_by_zgedi
  !
  IF(ionode)THEN
     !
     WRITE(*,*)" "
     WRITE(*,*)"    ======================================================================= "
     WRITE(*,*)" "
     WRITE(*,*)"      EPCDFT_Coupling Code warnings:"
     WRITE(*,*)"      1) Only works with norm conserving pseudos."
     WRITE(*,*)"      2) (which implies) that smooth grid = dense grid."
     WRITE(*,*)"      3) Run must be in serial."
     WRITE(*,*)"      4) No K-points."
     WRITE(*,*)"      5) Make sure your grids/cutoffs... are the same for both systems."
     WRITE(*,*)"      6) Previous PW runs must NOT use parallelization over k points."
     WRITE(*,*)"      7) occup1+occdown1 == occup2+occdown2."
     WRITE(*,*)"      8) if s_spin = .true. then  occup1 must = occup2."
     WRITE(*,*)"      9) if s_spin = .true. then occdown1 must = occdown2."
     WRITE(*,*)" "
     WRITE(*,*)"    ======================================================================= "
     WRITE(*,*)" "
     WRITE(*,*)" "
     WRITE(*,*)"    Data from input file :"
     WRITE(*,*)" "
     WRITE(*,*)"    prefix1      :", prefix
     WRITE(*,*)"    outdir1      :", tmp_dir
     WRITE(*,*)"    prefix2      :", prefix2
     WRITE(*,*)"    outdir2      :", tmp_dir2
     WRITE(*,*)"    debug        :", debug
     WRITE(*,*)"    s_spin       :", s_spin
     WRITE(*,*)"    det_by_zgedi :", det_by_zgedi
     WRITE(*,*)" "
     WRITE(*,*)"    # of spins   :", nks
     WRITE(*,*)"    # of bands   :", nbnd
     WRITE(*,*)"    # occ up states in sys 1   :", occup1 
     WRITE(*,*)"    # occ down states in sys 1 :", occdown1 
     WRITE(*,*)"    # occ up states in sys 2   :", occup2 
     WRITE(*,*)"    # occ down states in sys 2 :", occdown2 
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE print_checks_warns 
