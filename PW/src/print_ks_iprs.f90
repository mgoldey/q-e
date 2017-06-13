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
  REAL(DP), DIMENSION(nbnd,nkstot,3) :: &
      wf_locs ! location of wavefunction - expression is meaningless for delocalized wavefunctions 
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
        CALL calculate_wf_loc(ik,ibnd,wf_locs(ibnd,ik,:))    
      end do        

      WRITE( stdout, 9030 ) (ipr(ibnd,ik),ibnd=1,nbnd)

      write(stdout,9022)
      do ibnd=1,nbnd,1
        write(stdout, 9031) ibnd, wf_locs(ibnd,ik,1),wf_locs(ibnd,ik,2),wf_locs(ibnd,ik,3)
      end do
      
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
9022 FORMAT(/'          Wavefunction centers (Angstrom):'/ )
9030 FORMAT( '   ',8E12.4 )
9031 FORMAT(' Band: ',I4,'   ',8E12.4, '   ',8E12.4, '   ',8E12.4,' '/)
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

  ipr=top/(bottom*bottom)

  RETURN
END SUBROUTINE  calculate_ipr

SUBROUTINE calculate_wf_loc(ik,ibnd,locs)
  USE mp_bands,             ONLY : root_bgrp,root_bgrp_id, intra_bgrp_comm, inter_bgrp_comm
  USE mp,                   ONLY : mp_sum, mp_bcast
  USE mp_pools,             ONLY : inter_pool_comm, me_pool
  USE mp_images,            ONLY : intra_image_comm
  USE mp_global,            ONLY : inter_image_comm,my_image_id
  USE mp_world,             ONLY : world_comm
  USE fft_base,             ONLY : dfftp
  USE kinds,                ONLY : DP
  USE cell_base,            ONLY : omega,at,alat ! in cubic bohr
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
  USE constants,            ONLY : ANGSTROM_AU, pi
  USE realus,               ONLY : real_space, &
            invfft_orbital_gamma, calbec_rs_gamma, &
            invfft_orbital_k, calbec_rs_k
  USE buffers,               ONLY : get_buffer, save_buffer

  INTEGER, INTENT(IN) :: ik,ibnd
  REAL(DP), DIMENSION(3), INTENT(INOUT) :: locs
  REAL(DP)   :: x,y,z,a,b,c,posi(3)
  COMPLEX(DP):: x0,y0,z0
  INTEGER    :: npw,npw_k,idx0,ir,ipol,idx,i,j,k
  REAL(DP)   :: inv_nr1, inv_nr2, inv_nr3
  CHARACTER(len=1024) :: filename
  LOGICAL :: debug_write_cubes=.false.

  ! Clear work array of psic
  psic(:)=0d0
  top=0d0
  bottom=0d0
  npw = ngk(ik)
  npw_k=ngk(ik)
  current_k=ik
  inv_nr1 = 1._DP / DBLE(  dfftp%nr1 )
  inv_nr2 = 1._DP / DBLE(  dfftp%nr2 )
  inv_nr3 = 1._DP / DBLE(  dfftp%nr3 )

  IF ( nks > 1 ) CALL get_buffer( evc, nwordwfc, iunwfc, ik )
  CALL mp_bcast( evc, root_bgrp_id, inter_bgrp_comm )

  IF (gamma_only) THEN
    CALL invfft_orbital_gamma(evc,ibnd,ibnd,.true.) 
  ELSE
    CALL invfft_orbital_k(evc,ibnd,ibnd,ik,.true.) 
  ENDIF

  psic(:)=CONJG(psic)*psic
  psic(:)=psic(:)/(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x)

  if (debug_write_cubes .eqv. .true.) THEN
    if (ibnd < 10) then
      write(filename,"(A5,I1)") "band0",ibnd
    else
      write(filename,"(A4,I2)") "band",ibnd
    ENDIF
    CALL write_cube_r ( 84332, filename,  REAL(psic,KIND=DP))
  ENDIF

  ! Loop in the charge array
  ! idx0 = starting index of real-space FFT arrays for this processor
  !
  idx0 = dfftp%nr1x * dfftp%nr2x * dfftp%ipp(me_pool+1)
  !
  a=alat*at(1,1)/ANGSTROM_AU
  b=alat*at(2,2)/ANGSTROM_AU
  c=alat*at(3,3)/ANGSTROM_AU  
  x0=0.D0
  y0=0.D0
  z0=0.D0
  DO ir = 1, dfftp%nr1x*dfftp%nr2x*dfftp%npl
    !
    ! ... three dimensional indexes
    !
    idx = idx0 + ir - 1
    k   = idx / (dfftp%nr1x*dfftp%nr2x)
    idx = idx - (dfftp%nr1x*dfftp%nr2x)*k
    j   = idx / dfftp%nr1x
    idx = idx - dfftp%nr1x*j
    i   = idx

    DO ipol = 1, 3
      posi(ipol) = DBLE( i )*inv_nr1*alat*at(ipol,1)/ANGSTROM_AU + &
                   DBLE( j )*inv_nr2*alat*at(ipol,2)/ANGSTROM_AU + &
                   DBLE( k )*inv_nr3*alat*at(ipol,3)/ANGSTROM_AU
    ENDDO


    ! ... do not include points outside the physical range

    IF ( i >= dfftp%nr1 .OR. j >= dfftp%nr2 .OR. k >= dfftp%nr3 ) CYCLE

    x0=x0+EXP(((2.0*pi)/a)*posi(1)*CMPLX(0.d0,1.d0))*psic(ir)
    y0=y0+EXP(((2.0*pi)/b)*posi(2)*CMPLX(0.d0,1.d0))*psic(ir)
    z0=z0+EXP(((2.0*pi)/c)*posi(3)*CMPLX(0.d0,1.d0))*psic(ir)
  END DO

  CALL MP_SUM(x0,intra_image_comm)
  CALL MP_SUM(y0,intra_image_comm)
  CALL MP_SUM(z0,intra_image_comm)
  locs(1)=REAL(AIMAG(LOG(x0)))/(2.0*pi/a)
  if (locs(1).lt.0)     locs(1)=locs(1)+a
  locs(2)=REAL(AIMAG(LOG(y0)))/(2.0*pi/b)
  if (locs(2).lt.0)     locs(2)=locs(2)+b
  locs(3)=REAL(AIMAG(LOG(z0)))/(2.0*pi/c)
  if (locs(3).lt.0)     locs(3)=locs(3)+c
  RETURN

  CONTAINS
    !
  SUBROUTINE write_cube_r ( iu, fname, wfc_distr )
    ! -----------------------------------------------------------------
    !
    ! For debugging
    !
    USE pwcom,                 ONLY : npw,npwx
    USE kinds,                 ONLY : DP
    USE cell_base,             ONLY : celldm, at, bg
    USE ions_base,             ONLY : nat, tau, atm, ityp
    USE fft_base,              ONLY : dfftp !, grid_gather
    USE scatter_mod,           ONLY : gather_grid
    USE mp_global,             ONLY : me_bgrp,root_bgrp
    !
    IMPLICIT NONE
    !
    ! I/O 
    !
    INTEGER,INTENT(IN) :: iu
    CHARACTER(LEN=256),INTENT(IN) :: fname
    REAL(DP),INTENT(IN) :: wfc_distr(dfftp%nnr)
    ! 
    ! Workspace
    !
    REAL(DP)         :: alat
    INTEGER          :: nr1, nr2, nr3, nr1x, nr2x, nr3x
    INTEGER          :: i, nt, i1, i2, i3, at_num, ir
    INTEGER, EXTERNAL  :: atomic_number
    REAL(DP)    :: at_chrg, tpos(3), inpos(3)
    REAL(DP) :: wfc_gat(dfftp%nr1x*dfftp%nr2x*dfftp%nr3x)
    !
    wfc_gat=0.0_DP
    CALL gather_grid(dfftp,wfc_distr,wfc_gat)
    !
    !      WRITE A FORMATTED 'DENSITY-STYLE' CUBEFILE VERY SIMILAR
    !      TO THOSE CREATED BY THE GAUSSIAN PROGRAM OR THE CUBEGEN UTILITY.
    !      THE FORMAT IS AS FOLLOWS (LAST CHECKED AGAINST GAUSSIAN 98):
    ! 
    !      LINE   FORMAT      CONTENTS
    !      ===============================================================
    !       1     A           TITLE
    !       2     A           DESCRIPTION OF PROPERTY STORED IN CUBEFILE
    !       3     I5,3F12.6   #ATOMS, X-,Y-,Z-COORDINATES OF ORIGIN
    !       4-6   I5,3F12.6   #GRIDPOINTS, INCREMENT VECTOR
    !       #ATOMS LINES OF ATOM COORDINATES:
    !       ...   I5,4F12.6   ATOM NUMBER, CHARGE, X-,Y-,Z-COORDINATE
    !       REST: 6E13.5      CUBE DATA
    ! 
    !      ALL COORDINATES ARE GIVEN IN ATOMIC UNITS.
    !
    alat = celldm(1)
    nr1 = dfftp%nr1
    nr2 = dfftp%nr2
    nr3 = dfftp%nr3
    nr1x= dfftp%nr1x
    nr2x= dfftp%nr2x
    nr3x= dfftp%nr3x
    !
  !  wfc_gat(:)=wfc_gat(:)/DBLE( nr1*nr2*nr3 )
    !
    IF( me_bgrp == root_bgrp) THEN
       OPEN(UNIT=iu,FILE=TRIM(ADJUSTL(fname))//".cube")
       !
       WRITE(iu,*) 'Cubfile created from PWScf calculation'
       WRITE(iu,*) ' Total SCF Density'
       !                        origin is forced to (0.0,0.0,0.0)
       WRITE(iu,'(I5,3F12.6)') nat, 0.0_DP, 0.0_DP, 0.0_DP
       WRITE(iu,'(I5,3F12.6)') nr1, (alat*at(i,1)/DBLE(nr1),i=1,3)
       WRITE(iu,'(I5,3F12.6)') nr2, (alat*at(i,2)/DBLE(nr2),i=1,3)
       WRITE(iu,'(I5,3F12.6)') nr3, (alat*at(i,3)/DBLE(nr3),i=1,3)
       !
       DO i=1,nat
          nt = ityp(i)
          ! find atomic number for this atom.
          at_num = atomic_number(TRIM(atm(nt)))
          at_chrg= DBLE(at_num)
          ! at_chrg could be alternatively set to valence charge
          ! positions are in cartesian coordinates and a.u.
          !
          ! wrap coordinates back into cell.
          tpos = MATMUL( TRANSPOSE(bg), tau(:,i) )
          tpos = tpos - NINT(tpos - 0.5_DP)
          inpos = alat * MATMUL( at, tpos )
          WRITE(iu,'(I5,5F15.6)') at_num, at_chrg, inpos
       ENDDO
       !
       CALL actual_write_cube(wfc_gat,nr1,nr2,nr3,iu)
       !
       CLOSE(iu)
    ENDIF
    !
  END SUBROUTINE
  !
  SUBROUTINE actual_write_cube(func,nr1,nr2,nr3,ounit) 
    !
    ! For debugging
    !
    USE kinds,  ONLY : DP
    !
    IMPLICIT NONE
    !
    REAL(DP) :: func(nr1,nr2,nr3)
    INTEGER :: ounit,nr1,nr2,nr3
    !
    INTEGER :: i1,i2,i3
    !
    DO i1=1,nr1
       DO i2=1,nr2
          WRITE(ounit,'(6E17.5E3)') (func(i1,i2,i3),i3=1,nr3)
       ENDDO
    ENDDO
    !
  END SUBROUTINE


END SUBROUTINE  calculate_wf_loc