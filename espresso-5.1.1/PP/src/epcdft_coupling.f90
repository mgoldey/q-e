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
! It also computes det(<Vex1*Phi1|Phi2>) where Phi is a slater det and Vex1 is the 
! applied field to system 1.
! 
!
! Nicholas Brawand nicholasbrawand@gmail.com
!
! Notes:
!
!        1) V is applied to psi in realspace. CODE ONLY WORKS FOR Dense = Smooth
!           Vx1(r)_dense_parallel -> Vx1(r)_dense_serial ->
!           Vx1(r)_smooth_serial with (psi(G)->psi(r)) =  vpsi1(r)=|Vx1(r)*psi>
!           vpsi1(r) -> vpsi1(G)
!
!        2) Only Vx1 is working right now.
!
!        3) vex1_smat(i,j) = < Vex1*evc1(i) | evc2(j) > 
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
  USE lsda_mod,      ONLY : nspin
  USE epcdft, only : do_epcdft, fragment_atom1, fragment_atom2, epcdft_electrons,&
          epcdft_amp, epcdft_shift, epcdft_width
  USE cell_base,            ONLY : tpiba2, omega
  USE environment,          ONLY : environment_start, environment_end
  USE fft_base,             ONLY : dfftp, dffts, cgather_smooth, grid_gather
  USE fft_interfaces,       ONLY : invfft, fwfft
  USE scf,                  ONLY : rho, v
  USE control_flags,        ONLY : gamma_only
  !
  IMPLICIT NONE
  LOGICAL                      :: exst2
  LOGICAL                      :: debug
  LOGICAL                      :: s_spin         ! calculate S matrix for each spin separately
  LOGICAL                      :: det_by_zgedi   ! will use zgedi to get determinants if true
  CHARACTER (len=256)          :: outdir
  CHARACTER (len=256)          :: outdir2
  CHARACTER (len=256)          :: tmp_dir2
  CHARACTER (len=256)          :: tmp_dir_pass   ! used to store tmp_dir during pass for reading two systems
  INTEGER                      :: iunwfc_pass    ! same as above but diff var
  CHARACTER (len=256)          :: prefix_pass
  CHARACTER (len=256)          :: prefix2
  CHARACTER(LEN=256), external :: trimcheck
  INTEGER                      :: ios,ik,i,j, ig, is, itmp
  INTEGER                      :: ik1, ik2, ibnd1, ibnd2
  INTEGER                      :: iunwfc2 = 3636 ! unit for 2nd set of wfcs
  INTEGER                      :: info           ! for zgefa to stop zgedi 
  INTEGER                      :: occup1         ! occupied up states for system 1 
  INTEGER                      :: occup2         ! occupied up states for system 2 
  INTEGER                      :: occdown1       
  INTEGER                      :: occdown2       
  INTEGER                      :: occ(2,2)       ! array that holds all the occupations above (system,spin)
  INTEGER                      :: occtmp1        ! temp var will hold num occ states for system 1 in main loop 
  INTEGER                      :: occtmp2        ! temp var will hold num occ states for system 2 in main loop 
  REAL(DP)                     :: dtmp           ! temp variable
  REAL(DP),    EXTERNAL        :: ddot
  COMPLEX(DP)                  :: ztmp           ! temp variable
  COMPLEX(DP)                  :: vex1_test      ! sum_i <i|vex1|i> for testing 
  COMPLEX(DP)                  :: smatdet(2)     ! determinant of smat
  COMPLEX(DP)                  :: smatdet_spinup(2)! determinant of smat_spin
  COMPLEX(DP)                  :: smatdet_spindown(2)! determinant of smat_spin
  COMPLEX(DP)                  :: vex1_smatdet(2)! determinant of vex1_smat
  COMPLEX(DP), EXTERNAL        :: zdotc
  COMPLEX(DP), ALLOCATABLE     :: evc2(:,:)      ! will store 2nd vecs for dot prods
  COMPLEX(DP), ALLOCATABLE     :: vex1_evc1(:)   ! |Vex1_evc1>
  COMPLEX(DP), ALLOCATABLE     :: smat(:,:)      ! S_ij matrix <wfc1_i|wfc2_j>
  COMPLEX(DP), ALLOCATABLE     :: smat_spinup(:,:)! S_ij matrix <wfc1_i|wfc2_j> for each spin
  COMPLEX(DP), ALLOCATABLE     :: smat_spindown(:,:)! S_ij matrix <wfc1_i|wfc2_j> for each spin
  COMPLEX(DP), ALLOCATABLE     :: vex1_smat(:,:)   ! vex1*s_ij matrix <vex1*wfc1_i|wfc2_j>
  REAL(DP), DIMENSION(:),   ALLOCATABLE :: vxs1    ! Vx1 is added to this potential (serial)
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: vxp1    ! Vx1 is added to this potential (parallel)
  INTEGER  :: fragment1_atom1, fragment1_atom2, fragment2_atom1, fragment2_atom2
  REAL(DP) :: fragment1_amp, fragment2_amp
  !
  NAMELIST / inputpp / outdir, prefix, prefix2, outdir2, occup1, occup2, occdown1, occdown2, &
                       debug,  s_spin, det_by_zgedi, do_epcdft, fragment1_atom1, fragment1_atom2,&
                       fragment2_atom1, fragment2_atom2, fragment1_amp, fragment2_amp
  !
#ifdef __MPI
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'epcdft_coupling' )
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
  CALL read_file()  ! read , wfc, eigenvals, potential for system 2
  !
  ! now put system 2's vecs into system 2's vars
  ALLOCATE( evc2( npwx, nbnd ) )
  evc2 = evc
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
  occ(1,1) = occup1
  occ(1,2) = occdown1
  occ(2,1) = occup2
  occ(2,2) = occdown2
  !
  ALLOCATE( vxp1( dfftp%nnr, nspin  ) )
  ALLOCATE( vex1_evc1( npwx  ) )
  ALLOCATE( smat( SUM(occ(1,:)), SUM(occ(2,:)) ) )
  ALLOCATE( vex1_smat( SUM(occ(1,:)), SUM(occ(2,:)) ) )
  IF( s_spin ) ALLOCATE( smat_spinup( occ(1,1), occ(2,1) ) )
  IF( s_spin ) ALLOCATE( smat_spindown( occ(1,2), occ(2,2) ) )
  ALLOCATE( vxs1( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x ))
  !
  i         = 0
  j         = 0
  dtmp      = 0.d0
  ztmp      = 0.d0
  evc2      = 0.d0
  vex1_evc1 = 0.d0
  vxp1      = 0.d0
  vxs1      = 0.d0
  smat      = 0.d0
  IF( s_spin ) smat_spinup = 0.d0
  IF( s_spin ) smat_spindown = 0.d0
  vex1_smat = 0.d0
  smatdet   = 0.d0
  smatdet_spinup   = 0.d0
  smatdet_spindown = 0.d0
  psic      = 0.d0
  vex1_test = 0.d0
  vex1_smatdet = 0.d0
  !
  CALL print_checks_warns(prefix, tmp_dir, prefix2, tmp_dir2, nks, nbnd, &
                          occ(1,1), occ(1,2), occ(2,1), occ(2,2), debug,  s_spin, det_by_zgedi )
  !
  ! this call only calulates vpoten
  ! write(*,*) "CALLING add_efield"
  do_epcdft=.true.
  fragment_atom1=fragment1_atom1
  fragment_atom2=fragment1_atom2
  write(*,*) fragment_atom1, fragment_atom2
  epcdft_amp=fragment1_amp
  CALL add_efield( vxp1, dtmp, rho%of_r, .true. )
  ! write(*,*) "VXP1", vxp1(1:5)
  !
  ! gather the potentials
  !

#ifdef __MPI
    CALL grid_gather ( vxp1(:,1), vxs1(:))
#else
    vxs1(:)=vxp1(:,1)
#endif
  !
  ! Start calculating evc, vex1_evc and evc2 overlaps over kpnts and bands
  !
  i = 0
  !
  DO ik1 = 1, nks
     !
     ! prepare the indices & read evc1
     CALL gk_sort( xk(1,ik1), ngm, g, ecutwfc/tpiba2, npw, igk, g2kin )
     CALL davcio( evc, 2*nwordwfc, iunwfc, ik1, -1 ) 
     !
     ! set number of states to include in S matrix
     IF(ik1 == 1) occtmp1 = occ(1,1) 
     IF(ik1 == 2) occtmp1 = occ(1,2) 
     !
     DO ibnd1 = 1, occtmp1
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
           ! set number of states to include in S matrix
           IF(ik2 == 1) occtmp2 = occ(2,1) 
           IF(ik2 == 2) occtmp2 = occ(2,2) 
           !
           DO ibnd2= 1, occtmp2
              !
              j = j + 1 ! S matrix counter goes with evc2
              !
              IF(ik1 == ik2) THEN ! if spins match
                 !
                 IF(gamma_only) THEN
                    !
                    ! take the dot products and remove the double count on G=0 point
                    !
                    smat(i,j) =  2.d0 * ddot( 2*npw, evc(:,ibnd1), 1, evc2(:,ibnd2), 1 )
                    IF(gstart == 2) smat(i,j) = smat(i,j) - DBLE(evc(1,ibnd1)) * DBLE(evc2(1,ibnd2))
                    !
                    !
                    vex1_smat(i,j) = 2.d0 * ddot( 2*npw, vex1_evc1(:), 1, evc2(:,ibnd2), 1 )
                    IF(gstart == 2) vex1_smat(i,j) = vex1_smat(i,j) - DBLE(vex1_evc1(1)) * DBLE(evc2(1,ibnd2))
                    !
                 ELSE ! if not gamma_only
                    !
                    smat(i, j)      = zdotc( npwx, evc(:,ibnd1), 1, evc2(:,ibnd2), 1 )
                    !
                    vex1_smat(i, j) = zdotc( npwx, vex1_evc1(:), 1, evc2(:,ibnd2), 1 )
                    !
                 ENDIF ! end if gamma_only
                 !
              ENDIF ! end if spins match
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
  ! build S for each spin
  !
  IF(s_spin) THEN
     !
     ! spin up
     !
     DO i = 1, occ(2,1)
        !
        DO j = 1, occ(1,1)
           !
           smat_spinup(j,i) = smat( j, i )
           !
        ENDDO
        !
     ENDDO
     !
     ! spin down
     !
     DO i = 1, occ(2,2)
        !
        DO j = 1, occ(1,2)
           !
           smat_spindown(j,i) = smat( j+occ(1,1) , i+occ(2,1) )
           !
        ENDDO
        !
     ENDDO
     !
  ENDIF
  !
  ! Processing matrix results
  !
  IF(.true.) THEN
     !
     ! print matrices
     if (debug) THEN
      CALL print_cmat ( " S_row,col <psi1(row)|psi2(col)>", smat, SUM(occ(1,:)) )
      CALL print_cmat ( " <Vex1*psi1(row)|psi2(col)>", vex1_smat, SUM(occ(1,:)) )
      IF(s_spin) CALL print_cmat (" S_ij_up", smat_spinup, occ(1,1) )
      IF(s_spin) CALL print_cmat (" S_ij_down", smat_spindown, occ(1,2) )
     ENDIF
     !
     ! print trace of coupling matrix
     vex1_test = 0.d0
     DO i = 1, SUM(occ(1,:))
        vex1_test = vex1_test + vex1_smat(i,i)
     ENDDO
     CALL print_cnum (" Sum_i <i|Vx1|i>", vex1_test)
     !
     ! print sum of coupling matrix
     vex1_test = 0.d0
     DO i = 1, SUM(occ(2,:)) 
        DO j = 1, SUM(occ(1,:))
           vex1_test = vex1_test + vex1_smat(j,i)
        ENDDO
     ENDDO
     CALL print_cnum (" Sum_i,j <j|Vx1|i>", vex1_test)
     !
  ENDIF
  !
  ! calculate determinant of smat and vex1_smat
  !
  IF(det_by_zgedi) THEN
     !
     CALL get_det_from_zgedi( smat, SUM(occ(1,:)), smatdet )
     !
     CALL get_det_from_zgedi( vex1_smat, SUM(occ(1,:)), vex1_smatdet )
     !
     IF(s_spin) THEN
        !
        CALL get_det_from_zgedi( smat_spinup, occ(1,1), smatdet_spinup )
        !
        CALL get_det_from_zgedi( smat_spindown, occ(1,2), smatdet_spindown )
        !
     ENDIF
     !
  ELSE
     !
     CALL get_det_from_zgeev( SUM(occ(1,:)), smat, smatdet )
     !
     CALL get_det_from_zgeev( SUM(occ(2,:)), vex1_smat, vex1_smatdet )
     !
     IF(s_spin) THEN
        !
        CALL get_det_from_zgeev( occ(1,1), smat_spinup, smatdet_spinup )
        !
        CALL get_det_from_zgeev( occ(1,2), smat_spindown, smatdet_spindown )
        !
     ENDIF
     !
  ENDIF
  !
  ! Print determinants
  !
  CALL print_cnum( " Det( S_ij )",  smatdet(1) )
  !
  CALL print_cnum( " Det( < psi_i | V_x1 | psi_j > )",  vex1_smatdet(1) )
  !
  IF(s_spin) THEN
     !
     CALL print_cnum( " Det( S_up_ij )",  smatdet_spinup(1) )
     !
     CALL print_cnum( " Det( S_down_ij )",  smatdet_spindown(1) )
     !
  ENDIF
  !
  !
  ! closing shop
  !
  WRITE(*,*)" "
  WRITE(*,*)"  ========================================================================= "
  WRITE(*,*)" "
  !
  DEALLOCATE( vxs1 )
  !DEALLOCATE( vxp1 )
  DEALLOCATE( evc2 )
  DEALLOCATE( smat )
  DEALLOCATE( vex1_evc1 )
  DEALLOCATE( vex1_smat )
  IF( s_spin ) DEALLOCATE( smat_spinup )
  IF( s_spin ) DEALLOCATE( smat_spindown )
  !
  CALL environment_end ( 'epcdft_coupling' )
  !
  CALL stop_pp
  !
  STOP
  !
  1 FORMAT(8E12.3)
  2 FORMAT(2E16.5)
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
!
!
!
!-----------------------------------------------------------------------
SUBROUTINE print_cmat (aname, a, lda)
  !-----------------------------------------------------------------------
  !
  !     this routine prints a complex matrix a 
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN)            :: lda
  CHARACTER(*), INTENT(IN)       :: aname
  COMPLEX(DP), INTENT(IN)        :: a(lda,lda)      
  !
  INTEGER j
  !
  WRITE(*,*)""
  WRITE(*,*)" The REAL part of : "
  WRITE(*,*)aname
  WRITE(*,*)"-------------------------------------"
  !
  DO j = 1, lda
    WRITE(*,1)REAL(a(j,:))
  ENDDO
  !
!  WRITE(*,*)""
!  WRITE(*,*)" The IMG part of : "
!  WRITE(*,*)aname
!  WRITE(*,*)"-------------------------------------"
  !
!  DO j = 1, lda
!    WRITE(*,1)AIMAG(a(j,:))
!  ENDDO
  !
  1 FORMAT(8E12.3)
  !
  RETURN
  !
END SUBROUTINE print_cmat
!
!
!
!
!-----------------------------------------------------------------------
SUBROUTINE get_det_from_zgeev(n, a, det)
  !-----------------------------------------------------------------------
  !
  !     this routine will solve for the determinant of a matrix "a"
  !     and store it in det(1)
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  INTEGER,      INTENT(IN)    :: n
  COMPLEX(DP),  INTENT(INOUT) :: det(n)
  COMPLEX(DP),  INTENT(IN)    :: a(n,n)      
  !
  INTEGER                  :: i
  INTEGER                  :: info 
  EXTERNAL                 :: zgeev
  COMPLEX(DP)              :: VL(1), VR(1) 
  REAL(DP)                 :: zgeevrwork(2*n)
  COMPLEX(DP)              :: zgeevwork(2*n), zgeevw(n)
  !
  !
  !
  ! find and store eigenvalues of a in sgeevw
  !
  CALL ZGEEV( 'N', 'N', n, a, n, zgeevw, VL, 1, VR, 1, &
              zgeevwork, 2*n, zgeevrwork, info )
  !
  CALL errore('epcdft_coupling', 'error in zgeev', abs(info)) 
  !
  ! det(s) is given by prod of eigenvalues
  !
  det(1) = zgeevw(1)
  DO i = 2, n 
     det(1) = det(1) * zgeevw(i)
  ENDDO
  !
  RETURN
  !
END SUBROUTINE get_det_from_zgeev
!
!
!
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
!
!
!
!
!-----------------------------------------------------------------------------
SUBROUTINE get_det_from_zgedi(a, n, det)
  !--------------------------------------------------------------------------
  !
  !     this routine will solve for the determinant of a matrix "a"
  !     and store it in det(1)
  !
  USE kinds
  USE io_global,            ONLY : ionode
  !
  IMPLICIT NONE
  !
  INTEGER,      INTENT(IN)    :: n
  COMPLEX(DP),  INTENT(INOUT) :: det(2)
  COMPLEX(DP),  INTENT(IN)    :: a(n,n)      
  !
  INTEGER                  :: i
  INTEGER                  :: info 
  INTEGER                  :: ivpt(n)
  EXTERNAL                 :: zgefa
  EXTERNAL                 :: zgedi
  COMPLEX(DP)              :: work(n)
  !
  ! factor "a" by gaussian elimination
  ! a will be upper trianglular
  CALL zgefa(a, n, n, ivpt, info)  
  ! write(*,*) "Info is ", info
  !
  CALL errore( 'epcdft_coupling', 'error in zgefa', abs(info))
  !
  ! compute det of "a" using factors computed by zgefa
  !
  ! !!! the det = det(1) * 10**det(2) !!!
  !
  CALL zgedi(a, n, n, ivpt, det, work, 10)
  !
  ! compute det from output of zgedi
  !
  det(1) = det(1) * 10.d0**( det(2) )
  !
  RETURN
  !
END SUBROUTINE get_det_from_zgedi
!
!
!-----------------------------------------------------------------------
SUBROUTINE print_cnum (aname, a)
  !-----------------------------------------------------------------------
  !
  !     this routine prints a complex number "a"
  !
  USE kinds
  !
  IMPLICIT NONE
  !
  CHARACTER(*), INTENT(IN)       :: aname
  COMPLEX(DP),  INTENT(IN)       :: a
  !
  WRITE(*,*)""
  WRITE(*,*)"   ",aname
  WRITE(*,*)"   -------------------------------------"
  WRITE(*,*)"   ",a
  !
  RETURN
  !
END SUBROUTINE print_cnum
!
