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
!        1) 1 = (2-delta_{NumOfSpin,2}) * 1/TotalGridPnts * <phi|phi>
!
!           
!        2) S(i,j) matrix is smat_ij = < evc1(i) | evc2(j) > 
!                                 evc2(j) -->
!                                e
!                                v
!                                c
!                                1
!                               (i)
!
!                                |
!                                v
!
!        3) V is applied to psi in realspace. CODE ONLY WORKS FOR Dense = Smooth
!           Vx1(r)_dense_parallel -> Vx1(r)_dense_serial ->
!           Vx1(r)_smooth_serial with (psi(G)->psi(r)) =  vpsi1(r)=|Vx1(r)*psi>
!           vpsi1(r) -> vpsi1(G)
!
!        4) Only Vx1 is working right now.
!
!        5) vex1_smat(i,j) = < Vex1*evc1(i) | evc2(j) > 
!                                 evc2(j) -->
!                                V
!                                e
!                                x
!                                e
!                                *
!                                v
!                                c
!                                1
!                               (i)
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
  USE wavefunctions_module, ONLY : evc, psic
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE gvect,                ONLY : ngm, g, gstart
  USE gvecs,                ONLY : nls, nlsm
  USE noncollin_module,     ONLY : npol, nspin_mag, noncolin
  USE cell_base,            ONLY : tpiba2, omega
  USE environment,          ONLY : environment_start, environment_end
  USE fft_base,             ONLY : dfftp, dffts, cgather_smooth, grid_gather
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE scf,                  ONLY : rho, v
  USE control_flags,        ONLY : gamma_only
  USE fft_base            
  USE gvect           
  !
  IMPLICIT NONE
  LOGICAL                      :: exst2
  LOGICAL                      :: debug
  CHARACTER (len=20)          :: FMT
  CHARACTER (len=256)          :: outdir
  CHARACTER (len=256)          :: outdir2
  CHARACTER (len=256)          :: tmp_dir2
  CHARACTER (len=256)          :: prefix2
  CHARACTER(LEN=256), external :: trimcheck
  INTEGER                      :: ios,ik,i,j, ig, is, itmp
  INTEGER                      :: ik1, ik2, ibnd1, ibnd2
  INTEGER                      :: iunwfc2 = 3636 ! unit for 2nd set of wfcs
  INTEGER,    ALLOCATABLE      :: ivpt(:) ! pivot indices for zgefa 
  INTEGER                      :: info ! for zgefa to stop zgedi 
  REAL(DP)                     :: dtmp  ! temp variable
  COMPLEX(DP)                  :: ztmp  ! temp variable
  COMPLEX(DP)                  :: smatdet ! determinant of smat
  COMPLEX(DP)                  :: vex1_smatdet ! determinant of vex1_smat
  COMPLEX(DP), EXTERNAL        :: zdotc
  COMPLEX(DP), ALLOCATABLE     :: evc2(:,:)
  COMPLEX(DP), ALLOCATABLE     :: vex1_evc1(:) ! |Vex1_evc1>
  COMPLEX(DP), ALLOCATABLE     :: work(:)
  COMPLEX(DP), ALLOCATABLE  :: smat(:,:)     ! S_ij matrix <wfc1_i|wfc2_j>
  COMPLEX(DP),    ALLOCATABLE  :: vex1_smat(:,:)! vex1*s_ij matrix <vex1*wfc1_i|wfc2_j>
  REAL(DP), DIMENSION(:), ALLOCATABLE :: vxs1   ! Vx1 is added to this potential (serial)
  REAL(DP), DIMENSION(:), ALLOCATABLE :: vxp1   ! Vx1 is added to this potential (parallel)
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
  ALLOCATE(vxp1(dfftp%nnr))
  ALLOCATE(vxs1( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x ))
  ALLOCATE( evc2( npwx, nbnd ) )
  ALLOCATE( vex1_evc1( npwx ) )
  ALLOCATE( smat( nks*nbnd, nks*nbnd ) )
  ALLOCATE( vex1_smat( nks*nbnd, nks*nbnd ) )
  ALLOCATE( ivpt( nks*nbnd ) )
  ALLOCATE( work( nks*nbnd ) )
  evc2 = 0.0
  vex1_evc1 = 0.0
  vxp1 = 0.0
  vxs1 = 0.0
  smat = 0.0
  vex1_smat = 0.0
  ivpt = 0
  smatdet = 0.0
  vex1_smatdet = 0.0
  psic = 0.0
  debug = .true.
  !
  ! Print checks and errors
  !
  IF(ionode)THEN
     !
     WRITE(*,*)" "
     WRITE(*,*)"    ======================================================================= "
     WRITE(*,*)" "
     WRITE(*,*)"      EPCDFT_Coupling Code only works:"
     WRITE(*,*)"      1) with norm conserving pseudos."
     WRITE(*,*)"      2) (which implies) when smooth grid = dense grid."
     WRITE(*,*)"      3) in serial. (Lucky Charms form the best medium with your espresso)"
     WRITE(*,*)"      4) Need to check rows and colmns make sure in order (use mathematica)"
     WRITE(*,*)"      5) Need to check Vx1 is applied correctly."
     WRITE(*,*)" "
     WRITE(*,*)"    ======================================================================= "
     WRITE(*,*)" "
     !
  ENDIF
  !
  ! Print Debug Info
  !
  IF(debug)THEN
     !
     WRITE(*,*)" "
     WRITE(*,*)" "
     WRITE(*,*)"      EPCDFT_Coupling Code Stats:"
     WRITE(*,*)" "
     WRITE(*,*)"      Dense Grid Stuff"
     WRITE(*,*)"      size(nl)       : ", size(nl)
     WRITE(*,*)"      size(nlm)      : ", size(nlm)
     WRITE(*,*)"      dfftp%nnr      : ", dfftp%nnr
     WRITE(*,*)"      dfftp%nr1      : ", dfftp%nr1
     WRITE(*,*)"      dfftp%nr2      : ", dfftp%nr2
     WRITE(*,*)"      dfftp%nr3      : ", dfftp%nr3
     WRITE(*,*)"      dfftp%nr1x     : ", dfftp%nr1x
     WRITE(*,*)"      dfftp%nr2x     : ", dfftp%nr2x
     WRITE(*,*)"      dfftp%nr3x     : ", dfftp%nr3x
     WRITE(*,*)"      local g ngm    : ", ngm
     WRITE(*,*)"      all gs ngm_g   : ", ngm_g
     WRITE(*,*)"      # shell ngl    : ", ngl
     WRITE(*,*)" "
     WRITE(*,*)"      Smooth Grid Stuff"
     WRITE(*,*)"      size(nls)      : ", size(nls)
     WRITE(*,*)"      size(nlsm)     : ", size(nlsm)
     WRITE(*,*)"      dffts%nnr      : ", dffts%nnr
     WRITE(*,*)"      dffts%nr1      : ", dffts%nr1
     WRITE(*,*)"      dffts%nr2      : ", dffts%nr2
     WRITE(*,*)"      dffts%nr3      : ", dffts%nr3
     WRITE(*,*)"      dffts%nr1x     : ", dffts%nr1x
     WRITE(*,*)"      dffts%nr2x     : ", dffts%nr2x
     WRITE(*,*)"      dffts%nr3x     : ", dffts%nr3x
     WRITE(*,*)" "
     WRITE(*,*)"      Plane Wave Stuff"
     WRITE(*,*)"      npw            : ", npw 
     WRITE(*,*)"      npwx           : ", npwx 
     WRITE(*,*)"      size(igk)      : ", size(igk)
     WRITE(*,*)" "
     WRITE(*,*)"      Phys Objects"
     WRITE(*,*)"      nbnd           : ", nbnd 
     WRITE(*,*)"      size(psic)     : ", size(psic)
     WRITE(*,*)"      size(evc)      : ", size(evc)
     WRITE(*,*)"      size(evc,1)    : ", size(evc,1)
     WRITE(*,*)"      size(evc,2)    : ", size(evc,2)
     WRITE(*,*)"      size(rho%of_r) : ", size(rho%of_r)
     WRITE(*,*)"      size(v%of_r)   : ", size(v%of_r)
     WRITE(*,*)" "
     WRITE(*,*)" "
     !
  ENDIF
  !
  ! this call only calulates vpoten
  CALL add_efield( vxp1, dtmp, rho%of_r, .true. )
  !
  ! gather the potentials
#ifdef __MPI
    CALL grid_gather ( vxp1, vxs1 )
#else
    vxs1(:)=vxp1(:)
#endif
  !
  !WRITE(*,*)"Size of ecv dim1 ",SIZE(evc,1)," size of evc2 dim1 ",SIZE(evc2,1)
  !WRITE(*,*)"Size of ecv dim2 ",SIZE(evc,2)," size of evc2 dim2 ",SIZE(evc2,2)
  !
  i = 0
  DO ik1 = 1, nks
     !
     ! prepare the indices & read evc1
     CALL gk_sort( xk(1,ik1), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
     CALL davcio( evc, 2*nwordwfc, iunwfc, ik1, -1 ) 
     !
     DO ibnd1 = 1, nbnd
        !
        i = i + 1 ! S matrix counter goes with evc1
        !
        ! create | vxs1 * psic > go to realspace, multiply then back to G space
        !
        psic = 0.0
        !
        ! setup psic = evc_ik_iband
        IF(gamma_only)THEN
           !
           psic( nls (igk(1:npw)) ) =        evc( 1:npw, ibnd1 )
           psic( nlsm(igk(1:npw)) ) = CONJG( evc( 1:npw, ibnd1 ) )
           !
        ELSE
           !
           psic( nls(igk(1:npw)) ) = evc(1:npw,ibnd1)
           !
        ENDIF
        !
        ! go to real space
        CALL invfft ('Wave', psic, dffts)
        !
        ! apply Vex1 to psi in realspace
        DO itmp = 1, dffts%nnr
           !
           psic(itmp) = psic(itmp) * vxs1(itmp)
           !
        ENDDO
        !
        ! go back to g space
        CALL fwfft ('Wave', psic, dffts)
        !
        ! going back to the plane wave sphere
        IF(gamma_only)THEN
           !
           vex1_evc1( 1:npw )  = psic( nls (igk(1:npw)) ) 
           vex1_evc1( 1:npw ) = CONJG( psic( nlsm(igk(1:npw)) ) ) 
           !
        ELSE
           !
           vex1_evc1(1:npw) = psic( nls(igk(1:npw)) ) 
           !
        ENDIF
        !
        ! calculate <psi1|psi2> and <vxs1*psi1|psi2>
        !
        DO ik2 = 1, nks
           !
           ! prepare the indices & read evc2
           CALL gk_sort( xk(1,ik2), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
           CALL davcio( evc2, 2*nwordwfc, iunwfc2, ik2, -1 ) 
           !
           DO ibnd2= 1, nbnd 
              !
              j = j + 1 ! S matrix counter goes with evc2
              !
              IF(gamma_only) THEN
                 !
                 CALL gamma_dot( gstart, npwx, evc(:,ibnd1), evc2(:,ibnd2), smat(i, j) )
                 !
                 CALL gamma_dot( gstart, npwx, vex1_evc1(:), evc2(:,ibnd2), vex1_smat(i, j) )
                 !
              ELSE ! if not gamma_only
                 !
                 smat(i, j)      = zdotc( npw, evc(:,ibnd1), 1, evc2(:,ibnd2), 1 )
                 !
                 vex1_smat(i, j) = zdotc( npw, vex1_evc1(:), 1, evc2(:,ibnd2), 1 )
                 !
              ENDIF ! end if gamma_only
              !
           ENDDO ! end ibnd2
           !
        ENDDO ! end ik2
        !
        j = 0 ! S matrix inner counter reset for evc2
        !
     ENDDO ! end ibnd1
     !
  ENDDO ! end ik1
  !
  ! for spin = 1
  !IF(.NOT.noncolin) smat = 2.0 * smat
  !IF(.NOT.noncolin) vex1_smat = 2.0 * vex1_smat
  !
  ! Print smat
  !
  WRITE(*,*)""
  WRITE(*,*)"  REAL S_row,col <psi(row)|psi(col)>"
  WRITE(*,*)"-------------------------------------"
  DO j = 1, nbnd*nks
    WRITE(*,1)REAL(smat(j,:))
  ENDDO
  !
  WRITE(*,*)""
  WRITE(*,*)"  IMG S_row,col <psi(row)|psi(col)>"
  WRITE(*,*)"-------------------------------------"
  DO j = 1, nbnd*nks
    WRITE(*,1)AIMAG(smat(j,:))
  ENDDO
  !
  ! Print vex1_smat
  !
  WRITE(*,*)""
  WRITE(*,*)"  REAL <Vex1*psi(row)|psi(col)>"
  WRITE(*,*)"-------------------------------------"
  DO j = 1, nbnd*nks
    WRITE(*,1)REAL(vex1_smat(j,:))
  ENDDO
  !
  WRITE(*,*)""
  WRITE(*,*)"  IMG <Vex1*psi(row)|psi(col)>"
  WRITE(*,*)"-------------------------------------"
  DO j = 1, nbnd*nks
    WRITE(*,1)AIMAG(vex1_smat(j,:))
  ENDDO
  !
  ! calculate determinant of smat
  !
  CALL zgefa(smat,nbnd*nks,nbnd*nks,ivpt,info) ! prep matrix for det
  CALL errore('epcdft_coupling','error in zgefa',abs(info))
  CALL zgedi(smat,nbnd*nks,nbnd*nks,ivpt,smatdet,work,10) ! get det of Smat from zgefa ivpt's
  !
  ! calculate determinant of vex1_smat
  !
  ivpt = 0
  CALL zgefa(vex1_smat,nbnd*nks,nbnd*nks,ivpt,info) ! prep matrix for det
  CALL errore('epcdft_coupling','error in zgefa',abs(info))
  CALL zgedi(vex1_smat,nbnd*nks,nbnd*nks,ivpt,vex1_smatdet,work,10) ! get det of Smat from zgefa ivpt's
  !
  ! Print determinant smat
  !
  WRITE(*,*)""
  WRITE(*,*)"  Det( S_ij )"
  WRITE(*,*)"----------------"
  WRITE(*,1)smatdet
  !
  ! Print determinant vex1_smat
  !
  WRITE(*,*)""
  WRITE(*,*)"  Det( <Vex1*psi(row)|psi(col)> )"
  WRITE(*,*)"----------------"
  WRITE(*,1)vex1_smatdet
  !
  CALL environment_end ( 'epcdft_coupling' )
  !
  CALL stop_pp
  !
  STOP
  !
  1 FORMAT(4E12.3)
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
!
!-----------------------------------------------------------------------
SUBROUTINE gamma_dot (gstart, n, a, b, c)
  !-----------------------------------------------------------------------
  !
  ! Calculate <a|b> using gamma point tricks.
  ! Return the result in c. 
  ! n is the length of a and b.   
  !     
  !
  USE kinds, ONLY : DP
  !
  IMPLICIT NONE
  !
  COMPLEX(DP), EXTERNAL      :: zdotc
  INTEGER,     INTENT(IN)    :: gstart
  INTEGER,     INTENT(IN)    :: n
  COMPLEX(DP), INTENT(IN)    :: a(n)
  COMPLEX(DP), INTENT(IN)    :: b(n)
  COMPLEX(DP), INTENT(INOUT) :: c
  !
  ! sum over first half of G vectors
  c = zdotc( n, a, 1, b, 1 )
  !
  ! sum over 2nd half of the G vectors a* -> a and b -> b*
  c = c + zdotc( n, b, 1, a, 1 )
  !
  ! remove the double count at G=0 from above, 
  ! they are both real at G = 0 so Conjg doesn't matter
  IF(gstart == 2) c = c - a(1) * b(1)
  !
  !
END SUBROUTINE gamma_dot
