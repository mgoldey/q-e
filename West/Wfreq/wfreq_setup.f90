!
! Copyright (C) 2015-2016 M. Govoni 
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! This file is part of WEST.
!
! Contributors to this file: 
! Marco Govoni
!
!-----------------------------------------------------------------------
SUBROUTINE wfreq_setup
  !-----------------------------------------------------------------------
  !
  USE westcom,                ONLY : alphapv_dfpt,npwq0,sqvc,west_prefix,wfreq_dirname,&
                                   & n_pdep_eigen_to_use,n_imfreq,nbnd_occ,l_macropol,macropol_calculation,&
                                   & n_refreq,isz,qp_bandrange
  USE mp,                     ONLY : mp_max
  USE mp_global,              ONLY : intra_bgrp_comm
  USE pwcom,                  ONLY : npw,nbnd
  USE kinds,                  ONLY : DP
  USE gvect,                  ONLY : gstart,g,ig_l2g
  USE io_files,               ONLY : tmp_dir
  USE distribution_center,    ONLY : pert,macropert,ifr,rfr,aband
  USE class_idistribute,      ONLY : idistribute
  USE wavefunctions_module,   ONLY : evc
  USE mod_mpiio,              ONLY : set_io_comm
  !
  IMPLICIT NONE
  !
  REAL(DP) :: q(3)
  REAL(DP) :: qq
  COMPLEX(DP),EXTERNAL :: get_alpha_pv
  INTEGER :: ig
  !
  CALL do_setup ( ) 
  !
  ! Calculate ALPHA_PV
  !
  alphapv_dfpt = get_alpha_pv()
  !
  CALL set_npwq0()
  !
  ALLOCATE(sqvc(npwq0))
  !
  CALL store_sqvc(sqvc,npwq0,1,isz)
  !
  IF(qp_bandrange(1)>nbnd) CALL errore('wfreq_setup','Err: qp_bandrange(1)>nbnd', 1) 
  IF(qp_bandrange(2)>nbnd) CALL errore('wfreq_setup','Err: qp_bandrange(2)>nbnd', 1) 
  !
  CALL set_nbndocc()
  !
  wfreq_dirname = TRIM( tmp_dir ) // TRIM( west_prefix ) // '.wfreq.save'
  CALL my_mkdir( wfreq_dirname )
  !
  pert = idistribute()
  CALL pert%init(n_pdep_eigen_to_use,'i','npdep',.TRUE.)
  macropert = idistribute()
  CALL macropert%init(n_pdep_eigen_to_use+3,'i','npdep+macro',.TRUE.)
  ifr = idistribute()
  CALL ifr%init(n_imfreq,'z','n_imfreq',.TRUE.)
  rfr = idistribute()
  CALL rfr%init(n_refreq,'z','n_refreq',.TRUE.)
  aband = idistribute()
  CALL aband%init(nbnd,'i','nbnd',.TRUE.)
  !
  CALL set_freqlists( )
  !
  SELECT CASE(macropol_calculation)
  CASE('c','C')
     l_macropol = .TRUE.
  END SELECT
  !
  CALL set_io_comm( ) ! this defines communicator between heads of each image (me_bgrp==0) 
  !
END SUBROUTINE 
