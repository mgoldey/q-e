!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Contributors to this file:
!   Nicholas P. Brawand (nicholasbrawand@gmail.com)
!   Matthew B. Goldey (matthew.goldey@gmail.com)
!
!-----------------------------------------------------------------------
!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_setup
  !-----------------------------------------------------------------------
  !
  USE kinds,                ONLY : DP
  USE io_files,             ONLY : prefix, tmp_dir, diropn
  USE mp_global,            ONLY : npool  
  USE wvfct,                ONLY : nbnd, npwx, npw, g2kin
  USE klist,                ONLY : nks, igk_k
  USE io_global,            ONLY : ionode, ionode_id, stdout
  USE mp,                   ONLY : mp_bcast
  USE mp_world,             ONLY : world_comm
  USE io_files,             ONLY : nwordwfc, iunwfc
  USE fft_base,             ONLY : dfftp
  USE gvecw,                ONLY : gcutw
  USE scf,                  ONLY : rho
  USE klist,                ONLY : nks
  USE klist,                ONLY : xk
  USE gvect,                ONLY : ngm, g
  USE wavefunctions_module, ONLY : evc
  USE cell_base,            ONLY : tpiba2
  USE scf,                  ONLY : rho
  USE lsda_mod,             ONLY : nspin
  USE esm,                  ONLY : mill_2d, imill_2d
  USE epcdft_mod  
  USE epcdft
  !
  IMPLICIT NONE
  !
  CHARACTER(LEN=256) :: fil ! variable name to read cdft potential from output dir of run
  LOGICAL :: exst2 
  CHARACTER (len=256) :: tmp_dir2 ! temp variable to store system 2's dir info 
  CHARACTER (len=256) :: tmp_dir_pass   ! used to store tmp_dir during pass for reading two systems
  INTEGER  :: iunwfc_pass, ik, nctmp   ! same as above but diff var
  CHARACTER (len=256) :: prefix_pass
  CHARACTER(LEN=256), external :: trimcheck
  INTEGER :: ios
  REAL(DP) :: dtmp
  !
  NAMELIST / inputpp / outdir, prefix, prefix2, outdir2, &
                       debug,  s_spin, eig_of_w, debug2
  !
  ! setup vars and consistency checks
  !
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( TRIM( outdir ) == ' ' ) outdir = './'
  IF ( TRIM( outdir2 ) == ' ' ) outdir2 = './'
  !
  IF ( npool > 1 ) CALL errore('bands','pools not implemented',npool)
  !
  ! Read input file
  !
  IF ( ionode )  THEN
     !
     CALL input_from_file()
     ! 
     READ (5, inputpp, err = 200, iostat = ios)
200  CALL errore ('epcdft_coupling', 'reading inputpp namelist', ABS (ios) )
     !
     tmp_dir  = trimcheck (outdir)
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
!  CALL mp_bcast( occup1, ionode_id, world_comm )
!  CALL mp_bcast( occup2, ionode_id, world_comm )
!  CALL mp_bcast( occdown1, ionode_id, world_comm )
!  CALL mp_bcast( occdown2, ionode_id, world_comm )
  CALL mp_bcast( debug, ionode_id, world_comm )
  CALL mp_bcast( debug2, ionode_id, world_comm )
  CALL mp_bcast( s_spin, ionode_id, world_comm )
!  CALL mp_bcast( free1, ionode_id, world_comm )
!  CALL mp_bcast( free2, ionode_id, world_comm )
!  CALL mp_bcast( cor1, ionode_id, world_comm )
!  CALL mp_bcast( cor2, ionode_id, world_comm )
  CALL mp_bcast( eig_of_w, ionode_id, world_comm )
  !
  ! first read system 2 and store in system 1's variables 
  ! then read sys 1
  !
  tmp_dir_pass = tmp_dir ! store system 1's dir for later reading
  iunwfc_pass = iunwfc    
  prefix_pass = prefix
  !
  tmp_dir = tmp_dir2 
  iunwfc = iunwfc2
  prefix = prefix2
  !
  IF( ionode ) WRITE( stdout,*)"    ======================================================================= "
  IF( ionode ) WRITE( stdout,*) "    SYSTEM 2 INFO"
  do_epcdft=.false.
  CALL read_file()  ! for system 2
  CALL openfil_pp() 
  ALLOCATE( evc2 ( npwx, nbnd, nks ) )
  evc2 = ( 0.D0, 0.D0 )
  !
  DO ik = 1, nks
    CALL gk_sort( xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin )
    CALL davcio( evc2(:,:,ik), 2*nwordwfc, iunwfc, ik, -1 )
  ENDDO
  !
  ALLOCATE( w ( dfftp%nnr , 2, nspin ) )
  w = 0.D0
  !
  ! setup weight function for system 2
  CALL add_epcdft_efield(w(:,2,:),.TRUE.)
  nctmp = nconstr_epcdft
  ALLOCATE( lm( nconstr_epcdft, 2 ) )
  lm(:,2) = epcdft_guess(:)
  !
  !
  ! get correction = V*int*dr*w(r)*rho(r)
  !
  cor2 = epcdft_shift
  !
  ! get free energy and occupations
  !
  CALL epcdft_init_vars(free2, occup2, occdown2)
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
  !
  ! avoid double allocate esm arrays
  !
  IF( ALLOCATED( mill_2d ) ) DEALLOCATE( mill_2d )
  IF( ALLOCATED( imill_2d ) ) DEALLOCATE( imill_2d )
  !
  CALL read_file()  
  CALL openfil_pp() 
  !
  ALLOCATE( evc1 ( npwx, nbnd, nks ) )
  evc1 = ( 0.D0, 0.D0 )
  !
  DO ik = 1, nks
    CALL gk_sort( xk(1,ik), ngm, g, gcutw, npw, igk_k(1,ik), g2kin )
    CALL davcio( evc1(:,:,ik), 2*nwordwfc, iunwfc, ik, -1 )
  ENDDO
  !
  DEALLOCATE( evc )
  !
  ! setup weight function for system 1
  CALL add_epcdft_efield(w(:,1,:),.TRUE.)
  IF(nctmp .ne. nconstr_epcdft .and. ionode ) WRITE( stdout,*)&
  "    !!!! number of constraints for A and B not the same !!!!"
  !IF( nconstr_epcdft .ne. 1 .and. eig_of_w .and. ionode ) WRITE( stdout,*)&
  !"    !!!! number of constraints need to be 1 for eig_of_w = .true. !!!!"
  lm(:,1) = epcdft_guess(:)
  !
  ! get correction = V*int*dr*w(r)*rho(r)
  !
  cor1 = epcdft_shift
  !
  ! get free energy and occupations
  !
  CALL epcdft_init_vars(free1, occup1, occdown1)
  !
  IF( ionode ) WRITE( stdout,*)"    ======================================================================= "
  !
  ALLOCATE( smat ( 2 , 2 , nks) )
  ALLOCATE( wmat ( 2 , 2, nks) )
  !
  smat = ( 0.D0, 0.D0 )
  wmat = ( 0.D0, 0.D0 )
  !
  CALL print_checks_warns(prefix, tmp_dir, prefix2, tmp_dir2, nks, nbnd, &
                          occup1, occdown1, occup2, occdown2, debug,  s_spin, debug2 )
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
                              occdown1, occup2, occdown2, debug,  s_spin, debug2 )
  !--------------------------------------------------------------------------
  !
  !     this routine prints warnings and some data from the input file 
  !     
  !
  USE io_global,        ONLY : ionode
  USE io_global,        ONLY : ionode, stdout
  USE uspp,             ONLY : okvan
  USE paw_variables,    ONLY : okpaw
  USE control_flags,    ONLY : gamma_only
  USE noncollin_module, ONLY : noncolin
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
  LOGICAL,      INTENT(IN)    :: debug2, debug, s_spin
  INTEGER :: ier=1
  !
  IF(.NOT. gamma_only) CALL errore('add_epcdft_efield', 'CDFT requires gamma_only.', ier)  
  IF(okvan) CALL errore('add_epcdft_efield', 'CDFT: ultrasoft not implemented.', ier)  
  IF(okpaw) CALL errore('add_epcdft_efield', 'CDFT: PAW not implemented.', ier)  
  IF(noncolin) CALL errore('add_epcdft_efield', 'CDFT noncolin not implemented.', ier)  
  !
  IF(ionode)THEN
     !
     IF( ionode ) WRITE( stdout,*)" "
     IF( ionode ) WRITE( stdout,*)"    ======================================================================= "
     IF( ionode ) WRITE( stdout,*)" "
     IF( ionode ) WRITE( stdout,*)"      EPCDFT_Coupling:"
     IF( ionode ) WRITE( stdout,*)"      1) Make sure your grids/cutoffs... are the same for both systems."
     IF( ionode ) WRITE( stdout,*)"      2) Previous PW runs must NOT use parallelization over k points."
     IF( ionode ) WRITE( stdout,*)"      3) occup1+occdown1 == occup2+occdown2."
     IF( ionode ) WRITE( stdout,*)"      4) if s_spin = .true. then  occup1 must = occup2."
     IF( ionode ) WRITE( stdout,*)"      5) if s_spin = .true. then occdown1 must = occdown2."
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
     IF( ionode ) WRITE( stdout,*)"    debug2       :", debug2
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
!-----------------------------------------------------------------------------
SUBROUTINE epcdft_init_vars(inFree, inOccUp, inOccDown)
  !--------------------------------------------------------------------------
  !
  ! Init variables: free energy, occup, occdown
  !
  USE kinds, ONLY : DP
  USE pw_restart_new, ONLY : pw_readschema_file
  USE qes_types_module, ONLY : output_type, parallel_info_type, general_info_type
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: inFree
  INTEGER, INTENT(INOUT) :: inOccUp, inOccDown
  INTEGER :: ierr
  TYPE (output_type) :: output_obj
  TYPE (parallel_info_type) :: parinfo_obj
  TYPE (general_info_type) :: geninfo_obj
  !
  ! read xml and add output to obj to get etot
  !
  CALL pw_readschema_file ( ierr, output_obj, parinfo_obj, geninfo_obj)
  !
  CALL epcdft_build_free_energy(inFree, output_obj)
  !
  CALL epcdft_build_occ(inOccUp, inOccDown, output_obj)
  !
END SUBROUTINE  epcdft_init_vars
!
!
!-----------------------------------------------------------------------------
SUBROUTINE epcdft_build_occ(inOccUp, inOccDown, output_obj)
  !--------------------------------------------------------------------------
  !
  ! calculate occupations for up and down
  !
  USE kinds, ONLY : DP
  USE qes_types_module, ONLY : output_type
  !
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(INOUT)    :: inOccUp, inOccDown
  TYPE (output_type), INTENT(IN) :: output_obj
  !
  INTEGER :: i
  !
  ! count up
  !
  !
  IF(output_obj%band_structure%nbnd_up_ispresent)THEN
    inOccUp = 0
    DO i= 1, output_obj%band_structure%nbnd_up
      inOccUp = inOccUp + output_obj%band_structure%ks_energies(1)%occupations(i)
    ENDDO
  ENDIF
  !
  ! count down
  !
  IF(output_obj%band_structure%nbnd_dw_ispresent)THEN
    inOccDown = 0
    DO i = 1, output_obj%band_structure%nbnd_dw
      inOccDown = inOccDown + output_obj%band_structure%ks_energies(1)%occupations(output_obj%band_structure%nbnd_up+i)
    ENDDO
  ENDIF
  !
END SUBROUTINE  epcdft_build_occ
!
!
!-----------------------------------------------------------------------------
SUBROUTINE epcdft_build_free_energy(inFree, output_obj)
  !--------------------------------------------------------------------------
  !
  ! calculate free energy and put it in inFree
  !
  USE kinds, ONLY : DP
  USE epcdft, ONLY : epcdft_guess, epcdft_target, nconstr_epcdft, epcdft_shift
  USE qes_types_module, ONLY : output_type
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT)    :: inFree
  TYPE (output_type), INTENT(IN) :: output_obj
  !
  REAL(DP) :: etot
  INTEGER :: ik
  !
  ! etot = E_ks + V*int*dr*w(r)*rho(r) - V*N_cdft
  !
  ! free = E_ks + V*int*dr*w(r)*rho(r)
  !      = etot + V*N_cdft 
  !
  !
  ! etot is in Hartree it needs to be in Ry
  !
  etot = 2.D0*output_obj%total_energy%etot
  !
  ! build V*N_cdft 
  !
  inFree = 0.D0
  DO ik = 1, nconstr_epcdft
    inFree = inFree + epcdft_guess(ik)*epcdft_target(ik)
  ENDDO
  !
  inFree = etot + inFree
  !
END SUBROUTINE  epcdft_build_free_energy
!
