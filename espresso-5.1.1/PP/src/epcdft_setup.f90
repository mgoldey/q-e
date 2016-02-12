!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_setup
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_files,             ONLY : prefix, tmp_dir, diropn
  USE mp_global,            ONLY : npool  
  USE wvfct,                ONLY : nbnd, npwx, igk, npw, g2kin, ecutwfc
  USE klist,                ONLY : nks
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE fft_base,             ONLY : dfftp, grid_gather
  USE scf,                  ONLY : rho
  USE klist,                ONLY : nks
  USE klist,                ONLY : xk
  USE gvect,                ONLY : ngm, g
  USE wavefunctions_module, ONLY : evc
  USE cell_base,            ONLY : tpiba2
  USE scf,                  ONLY : rho
  USE epcdft_mod  
  USE epcdft
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: fil ! variable name to read cdft potential from output dir of run
  LOGICAL :: exst2 
  CHARACTER (len=256) :: tmp_dir2 ! temp variable to store system 2's dir info 
  CHARACTER (len=256) :: tmp_dir_pass   ! used to store tmp_dir during pass for reading two systems
  INTEGER  :: iunwfc_pass, ik   ! same as above but diff var
  CHARACTER (len=256) :: prefix_pass
  CHARACTER(LEN=256), external :: trimcheck
  INTEGER :: ios
  REAL(DP) :: dtmp ! temp variable
  !
  NAMELIST / inputpp / outdir, prefix, prefix2, outdir2, occup1, occup2, occdown1, occdown2, &
                       debug,  s_spin, free1, free2,&
                       cor1, cor2, eig_of_w
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
  CALL mp_bcast (ios, ionode_id, world_comm )
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( tmp_dir2, ionode_id, world_comm )
  CALL mp_bcast( outdir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( prefix2, ionode_id, world_comm )
  CALL mp_bcast( outdir2, ionode_id, world_comm )
  CALL mp_bcast( occup1, ionode_id, world_comm )
  CALL mp_bcast( occup2, ionode_id, world_comm )
  CALL mp_bcast( occdown1, ionode_id, world_comm )
  CALL mp_bcast( occdown2, ionode_id, world_comm )
  CALL mp_bcast( debug, ionode_id, world_comm )
  CALL mp_bcast( s_spin, ionode_id, world_comm )
  CALL mp_bcast( free1, ionode_id, world_comm )
  CALL mp_bcast( free2, ionode_id, world_comm )
  CALL mp_bcast( cor1, ionode_id, world_comm )
  CALL mp_bcast( cor2, ionode_id, world_comm )
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
  IF( ionode ) WRITE( stdout,*)"    ======================================================================= "
  IF( ionode ) WRITE( stdout,*) "    SYSTEM 2 INFO"
  do_epcdft=.false.
  CALL read_file()  ! for system 2
  CALL openfil_pp() ! open all files for scf run set filenames/units
  ALLOCATE( evc2 ( npwx, nbnd, nks ) )
  evc2 = ( 0.D0, 0.D0 )
  !
  DO ik = 1, nks
    CALL gk_sort( xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
    CALL davcio( evc2(:,:,ik), 2*nwordwfc, iunwfc, ik, -1 )
  ENDDO
  !
  ALLOCATE( w ( dfftp%nnr , 2 ) )
  w = 0.D0
  !
  ! setup weight function for system 2
  CALL add_epcdft_efield(w(:,2),.TRUE.)
  !
  ! get fld str
  !
  wamp2 = 1.D0
  DO ik = 1, nconstr_epcdft
    wamp2 = wamp2 * epcdft_guess(ik)
  ENDDO
  !fil =  TRIM( tmp_dir ) // TRIM( prefix ) // 'v_cdft.cub'
  !CALL read_cube(239841274, fil, w(:,2) )
  !
  ! deallocate to avoid reallocation of sys 1 vars
  CALL clean_pw( .TRUE. )
  !
  ! restore sys 1's vars and read sys 1's data
  IF( ionode ) WRITE( stdout,*)"    ======================================================================= "
  IF( ionode ) WRITE( stdout,*) "    SYSTEM 1 INFO"
  do_epcdft=.false.
  tmp_dir = tmp_dir_pass
  iunwfc = iunwfc_pass
  prefix = prefix_pass
  CALL read_file()  
  CALL openfil_pp() 
  !
  ALLOCATE( evc1 ( npwx, nbnd, nks ) )
  evc1 = ( 0.D0, 0.D0 )
  !
  DO ik = 1, nks
    CALL gk_sort( xk(1,ik), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
    CALL davcio( evc1(:,:,ik), 2*nwordwfc, iunwfc, ik, -1 )
  ENDDO
  !
  DEALLOCATE( evc )
  !
  ! setup weight function for system 1
  CALL add_epcdft_efield(w(:,1),.TRUE.)
  !
  ! get fld str
  !
  wamp1 = 1.D0
  DO ik = 1, nconstr_epcdft
    wamp1 = wamp1 * epcdft_guess(ik)
  ENDDO
  !fil =  TRIM( tmp_dir ) // TRIM( prefix ) // 'v_cdft.cub'
  !CALL read_cube(239841275, fil, w(:,1) )
  !
  IF( ionode ) WRITE( stdout,*)"    ======================================================================= "
  !
  !CALL init_us_1    ! compute pseduo pot stuff
  !
  ALLOCATE( smat ( 2 , 2 , nks) )
  ALLOCATE( wmat ( 2 , 2, nks) )
  !
  smat = ( 0.D0, 0.D0 )
  wmat = ( 0.D0, 0.D0 )
  !
  CALL print_checks_warns(prefix, tmp_dir, prefix2, tmp_dir2, nks, nbnd, &
                          occup1, occdown1, occup2, occdown2, debug,  s_spin )
  !
  IF( ionode ) WRITE( stdout,*)" "
  IF( ionode ) WRITE( stdout,*)"    ======================================================================= "
  IF( ionode ) WRITE( stdout,*)"    Progress : "
  IF( ionode ) WRITE( stdout,*)"    Setup done"
  !
END SUBROUTINE epcdft_setup
!
!-----------------------------------------------------------------------------
SUBROUTINE print_checks_warns(prefix, tmp_dir, prefix2, tmp_dir2, nks, nbnd, occup1, &
                              occdown1, occup2, occdown2, debug,  s_spin )
  !--------------------------------------------------------------------------
  !
  !     this routine prints warnings and some data from the input file 
  !     
  !
  USE io_global,            ONLY : ionode
  USE io_global,            ONLY : ionode, stdout
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
  LOGICAL,      INTENT(IN)    :: debug, s_spin
  !
  IF(ionode)THEN
     !
     IF( ionode ) WRITE( stdout,*)" "
     IF( ionode ) WRITE( stdout,*)"    ======================================================================= "
     IF( ionode ) WRITE( stdout,*)" "
     IF( ionode ) WRITE( stdout,*)"      EPCDFT_Coupling Code warnings:"
     IF( ionode ) WRITE( stdout,*)"      1) Only works with norm conserving pseudos."
     IF( ionode ) WRITE( stdout,*)"      2) (which implies) that smooth grid = dense grid."
     IF( ionode ) WRITE( stdout,*)"      3) Run must be in serial."
     IF( ionode ) WRITE( stdout,*)"      4) No K-points."
     IF( ionode ) WRITE( stdout,*)"      5) Make sure your grids/cutoffs... are the same for both systems."
     IF( ionode ) WRITE( stdout,*)"      6) Previous PW runs must NOT use parallelization over k points."
     IF( ionode ) WRITE( stdout,*)"      7) occup1+occdown1 == occup2+occdown2."
     IF( ionode ) WRITE( stdout,*)"      8) if s_spin = .true. then  occup1 must = occup2."
     IF( ionode ) WRITE( stdout,*)"      9) if s_spin = .true. then occdown1 must = occdown2."
     IF( ionode ) WRITE( stdout,*)" "
     IF( ionode ) WRITE( stdout,*)"    ======================================================================= "
     IF( ionode ) WRITE( stdout,*)" "
     IF( ionode ) WRITE( stdout,*)" "
     IF( ionode ) WRITE( stdout,*)"    Data from input file :"
     IF( ionode ) WRITE( stdout,*)" "
     IF( ionode ) WRITE( stdout,*)"    prefix1      :", prefix
     IF( ionode ) WRITE( stdout,*)"    outdir1      :", tmp_dir
     IF( ionode ) WRITE( stdout,*)"    prefix2      :", prefix2
     IF( ionode ) WRITE( stdout,*)"    outdir2      :", tmp_dir2
     IF( ionode ) WRITE( stdout,*)"    debug        :", debug
     IF( ionode ) WRITE( stdout,*)"    s_spin       :", s_spin
     IF( ionode ) WRITE( stdout,*)" "
     IF( ionode ) WRITE( stdout,*)"    # of spins   :", nks
     IF( ionode ) WRITE( stdout,*)"    # of bands   :", nbnd
     IF( ionode ) WRITE( stdout,*)"    # occ up states in sys 1   :", occup1 
     IF( ionode ) WRITE( stdout,*)"    # occ down states in sys 1 :", occdown1 
     IF( ionode ) WRITE( stdout,*)"    # occ up states in sys 2   :", occup2 
     IF( ionode ) WRITE( stdout,*)"    # occ down states in sys 2 :", occdown2 
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE print_checks_warns 
!
!
!-----------------------------------------------------------------------
SUBROUTINE read_cube(iunps, appfile, cubedatout)
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE fft_base,             ONLY : dfftp
  USE mp,            ONLY : mp_bcast
  USE mp_images,     ONLY : intra_image_comm
  USE io_global,     ONLY : stdout, ionode, ionode_id
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: iunps
  CHARACTER(LEN=256), INTENT(IN) :: appfile
  REAL(DP), INTENT(INOUT) :: cubedatout(dfftp%nnr)
  !
  INTEGER :: cnr1, cnr2, cnr3, cnat, i, j, ierr, k
  REAL(DP) :: cubedat(dfftp%nr1x,dfftp%nr2x,dfftp%nr3x)
  !
  ! check for input file and handle errors
  OPEN  (unit = iunps, file = appfile, status = 'old', &
         form = 'formatted', action='read', iostat = ierr)
  !
  IF( ierr /= 0) THEN
    !
    ! if file failed stop QE
    CALL mp_bcast( ierr, ionode_id, intra_image_comm )
    CALL errore( 'add_efield ', &
                 'can not find input cubefile for applied field', ierr )
    !
  ENDIF 
  !
  ! Found file continue to read
  ! 
  ! skip header
  READ(iunps,*)
  READ(iunps,*)
  !
  ! read cube data
  READ(iunps, *) cnat
  READ(iunps, *) cnr1
  READ(iunps, *) cnr2
  READ(iunps, *) cnr3
  !
  ! cube data must match scf data
  IF(cnr1/=dfftp%nr1x .OR. cnr2/=dfftp%nr2x .OR. cnr3/=dfftp%nr3x) THEN
    CALL mp_bcast( ierr, ionode_id, intra_image_comm )
    CALL errore( 'read_cube', &
                 'Number of grid points in cube != dense FFT grid', 0 )
  ENDIF 
  !
  ! skip atoms
  DO i=1, cnat
    READ(iunps,*) 
  ENDDO
  !
  ! init array and store cube data 
  DO i=1, cnr1
    !
    DO j=1, cnr2
      !
      READ(iunps,'(6E13.5)') (cubedat(i,j,k),k=1,cnr3)
      !
    ENDDO
    !
  ENDDO
  !
  CALL scat_wrap(cubedat, dfftp%nr1x*dfftp%nr2x*dfftp%nr3x, cubedatout, dfftp%nnr)
  !
END SUBROUTINE read_cube
!
SUBROUTINE scat_wrap(srl,ns,par,np)
  !
  USE kinds, ONLY : DP
  USE fft_base,             ONLY : grid_scatter
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: ns, np
  REAL(DP), INTENT(IN) :: srl(ns)
  REAL(DP), INTENT(INOUT) :: par(np)
  !
#ifdef __MPI
  CALL grid_scatter(srl, par)
#else
  par = srl
#endif
  !
END SUBROUTINE 
