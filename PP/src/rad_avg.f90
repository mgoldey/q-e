!
! Copyright (C) 2001-2007 Quantum-Espresso group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! Written by T. Anh Pham
! UC Davis & LLNL: atupham@ucdavis.edu, pham16@llnl.gov
! 
! Rewritten for spherical average by M. Voros
! ANL: vormar@gmail.com
!
! Touched up for QE 6.1 by M. Goldey
! U of Chicago : matthew.goldey@gmail.com 
!
! Original purpose:
! The program is designed to calculate the local density of states (LDOS): D(z,e)
! and local band edges along z-direction of any interface models. 
! Current version is for one-k point, that is for supercell calculation. 
!
! New purpose:
! The program is designed to calculate the local density of states (LDOS): D(r,e)
! and local band edges along a radial direction of any (almost) spherical model:
! eg. nanoparticle embedded in a matrix. 
! Current version is for one-k point, that is for supercell calculation. 
!
!-----------------------------------------------------------------------
PROGRAM do_plan_avg
!-----------------------------------------------------------------------
  USE kinds,     ONLY : DP
  USE run_info,  ONLY: title
  USE cell_base, ONLY : ibrav, celldm, at, alat,tpiba, tpiba2, omega
  USE gvect
  USE wvfct,     ONLY : npw, npwx, nbnd, wg, g2kin,et
  USE klist,     ONLY : xk, nks, nkstot, ngk, igk_k,wk
  USE ions_base, ONLY : nat, ntyp=>nsp, ityp, tau, atm, zv
  USE io_files,  ONLY : tmp_dir, prefix, nd_nmbr, iunwfc, nwordwfc,iunoldwfc, diropn
  USE io_global, ONLY : ionode, ionode_id, stdout
  USE mp,        ONLY : mp_bcast, mp_sum 
  USE mp_global, ONLY : mp_startup
  USE mp_global, ONLY : intra_pool_comm,world_comm
  USE control_flags, ONLY : gamma_only
  USE gvecs,     ONLY : nls, nlsm, doublegrid,dual
  USE constants, ONLY: rytoev
  USE scf,       ONLY: rho, rho_core, rhog_core, vnew
  USE lsda_mod,  ONLY : nspin, isk
  USE wavefunctions_module,  ONLY: evc
  USE becmod,    ONLY : bec_type, becp, calbec, allocate_bec_type, deallocate_bec_type
  USE gvecw,                ONLY : ecutwfc
  USE uspp, ONLY: vkb, nkb
  USE constants, ONLY: tpi !, e2  
  USE fft_base,        ONLY : dffts, dfftp
  USE fft_interfaces,  ONLY : fwfft, invfft
  USE basic_algebra_routines, ONLY : norm
  USE environment,   ONLY : environment_start, environment_end
  !
  IMPLICIT NONE
  !
  CHARACTER(len=256) :: filplot, outdir
  REAL(DP), ALLOCATABLE    :: ravg (:,:,:)
  CHARACTER(LEN=256), EXTERNAL :: trimcheck
  REAL(DP), ALLOCATABLE    :: DOSofE(:)
  !
  ! array for wave functions
  COMPLEX(DP), ALLOCATABLE    :: phi (:)
  ! array for wave functions squared
  REAL(DP), ALLOCATABLE    :: phi1 (:) 
  ! 
  !
  REAL(DP) :: w1, w2, intnorm, rmax,  dr, r, rg, absg, phase 
  REAL(DP) :: E,  Elw, Eup, DeltaE, Emin, Emax, degauss1,degauss
  REAL(DP) :: w0gauss,mp_shift
  external w0gauss
  !
  REAL(DP), PARAMETER      :: x(3) = (/ 0.5D0, 0.0D0, 0.0D0 /), &
                              y(3) = (/ 0.0D0, 0.5D0, 0.0D0 /), &
                              z(3) = (/ 0.0D0, 0.0D0, 0.5D0 /)
  INTEGER :: npts = 100
  
  INTEGER :: iunplot = 4, ios, ibnd, ik, ir, nt, na, i,j,k,ig, ndos, ngauss,ia
  INTEGER :: k1,k2, bndmin, bndmax
  INTEGER  :: ed1, ed2
  REAL(DP) :: Eed1, Eed2
  REAL(DP) :: Ef, delta 
  REAL(DP) :: Ef1, Ef2
  COMPLEX(DP) :: cavg 
  CHARACTER(len=10) :: filename
  LOGICAL :: exst
  REAL(DP) :: r0(3), r01, r02, r03, zvtot, sinxx 
  !
  NAMELIST / inputpp / bndmin, bndmax,outdir, prefix, Emin,Emax,&
      & filplot, Ef, Ef1, Ef2, delta, degauss, DeltaE, r01, r02, r03, npts, mp_shift
  !
#if defined(__MPI)
  CALL mp_startup ( )
#endif
  CALL environment_start ( 'rad_avg' )
  !
  !   set default values for variables in namelist
  !
  prefix = 'pwscf'
  CALL get_environment_variable( 'ESPRESSO_TMPDIR', outdir )
  IF ( trim( outdir ) == ' ' ) outdir = './'
  filplot = 'tmp.pp'
  ios = 0
  mp_shift=0d0
  ! set default bndmin, bndmax
  !
  bndmin=1
  bndmax=nbnd
  !
  ! set default r0 as center of charge
  ! good for molecules in vacuum (molecule can be centered anywhere in the cell)
  !
  zvtot = 0.D0
  r0(:) = 0.D0
  DO ia = 1, nat
     zvtot = zvtot + zv(ityp(ia))
     r0(:) = r0(:) + tau(:,ia)*zv(ityp(ia))
  END DO
  r0(:) = r0(:) / zvtot
  !
  Emin   =-1000000.d0
  Emax   = 1000000.d0

  ! read input file and broadcast
  IF ( ionode )  THEN
     !
     CALL input_from_file ( )
     READ (5, inputpp, iostat = ios)
     tmp_dir = trimcheck (outdir)
     !
  END IF
  !
  CALL mp_bcast( ios, ionode_id ,world_comm)
  IF ( ios /= 0 ) CALL errore ('rad_avg', 'reading inputpp namelist',ABS(ios))
  !
  CALL mp_bcast( tmp_dir, ionode_id, world_comm )
  CALL mp_bcast( prefix, ionode_id, world_comm )
  CALL mp_bcast( filplot, ionode_id, world_comm )
  CALL mp_bcast( Ef, ionode_id, world_comm )
  CALL mp_bcast( Ef1, ionode_id, world_comm )
  CALL mp_bcast( Ef2, ionode_id, world_comm )
  CALL mp_bcast( delta, ionode_id, world_comm )
  CALL mp_bcast( degauss, ionode_id, world_comm )
  CALL mp_bcast( DeltaE, ionode_id, world_comm )
  CALL mp_bcast( r01, ionode_id, world_comm )
  CALL mp_bcast( r02, ionode_id, world_comm )
  CALL mp_bcast( r03, ionode_id, world_comm )
  CALL mp_bcast( npts, ionode_id, world_comm )
  CALL mp_bcast( mp_shift, ionode_id, world_comm )
  CALL mp_bcast( bndmin, ionode_id, world_comm )
  CALL mp_bcast( bndmax, ionode_id, world_comm )
  CALL mp_bcast( Emin, ionode_id, world_comm )
  CALL mp_bcast( Emax, ionode_id, world_comm )
  !
  r0(1)=r01
  r0(2)=r02
  r0(3)=r03

  CALL read_file ( )
  CALL openfil_pp ( )
  !
  IF ( ionode ) THEN
     !
     OPEN (UNIT = iunplot, FILE = filplot, FORM = 'formatted', &
          STATUS = 'unknown', err = 100, IOSTAT = ios)
100  CALL errore ('plan_avg', 'opening file '//TRIM(filplot), abs (ios))
     !
  ENDIF
  !
  ! find sphere with maximum radius that can fit the simulation cell
  !
  rmax = norm( MATMUL( at(:,:), x(:) ) ) 
  rmax = MIN( rmax, norm( MATMUL( at(:,:), y(:) ) ) )
  rmax = MIN( rmax, norm( MATMUL( at(:,:), z(:) ) ) )
  !
  ! bring it to atomic units (Bohr)
  !
  rmax = rmax*alat
  if (IONODE) write(*,'(A,F8.5,A)') "Maximum radial distance is ",rmax," Bohr"
  if (bndmax.eq.0) bndmax=nbnd

  !
  ! step along radial direction
  !
  dr = rmax / DBLE( npts )
  !
  ! allocate array for spherically averaged "DOS"
  !
  ALLOCATE (ravg(npts, nbnd, nks))
  !
  !
  CALL allocate_bec_type ( nkb, nbnd, becp)
  allocate(phi(dfftp%nnr))
  allocate(phi1(dfftp%nnr))
  phi=0.D0
  phi1=0.D0
  ravg(:,:,:)=0.d0

  if (IONODE) write(*,'(A)',advance="yes") "Performing radial average for band "
  do ik = 1, nks
     call gk_sort (xk (1, ik), ngm, g, ecutwfc / tpiba2, npw, igk_k(1,ik),g2kin)
     !call davcio (evc, 2*nwordwfc, iunwfc, ik, - 1)
     call init_us_2 (npw, igk_k(1,ik), xk (1, ik), vkb)
     call calbec ( npw, vkb, evc, becp)
     do ibnd = bndmin,bndmax
         IF (MOD(ibnd-bndmin+1,8).eq.0 .and. (ibnd-bndmin) .gt. 1) THEN
             if (IONODE) write(*,'(I6)',advance="yes") ibnd
         ELSE
             if (IONODE) write(*,'(I6,A)',advance="no") ibnd, " "
         ENDIF
         if (IONODE .and. ibnd.eq.bndmax) WRITE(*,*) ""
        !
        phi = (0.d0,0.d0)
        !
        if ( gamma_only ) then
          do ig = 1, npw
         phi (nls (igk_k(ig,ik) ) ) = evc (ig, ibnd)
         phi (nlsm (igk_k(ig,ik) ) ) = conjg(evc (ig, ibnd))
          enddo
        else
           do ig = 1, npw
              phi (nls (igk_k(ig,ik) ) ) = evc (ig, ibnd)
           enddo
        endif
        call invfft('Wave',phi, dfftp)
        !
        ! calculate 1/V*|psi(r)|^2
        !
        w1 = wk (ik) / omega
        do ir = 1, dfftp%nnr
           phi1(ir) = w1*(REAL(phi(ir))**2.0 + AIMAG(phi(ir))**2.0)
        enddo
        !
        !write (stdout, *),'nrxx',nrxx,'omega',omega,'ngm',ngm,'npw',npw,'nr1s',nr1s,'nrx1s',nrx1s,'nr1',nr1,'nrx1',nrx1
        ! use phi as work array
        phi(:) =  CMPLX ( phi1(:), 0.d0)
        ! 
        CALL fwfft ('Dense', phi, dfftp)

        DO ir = 1, npts
          !
          ! ... r is in atomic units
          !
          r = dr*ir
          !
          !
          DO ig = 1, ngm
             !
             !WRITE( stdout, *), 'ibnd', ibnd,'ir', ir,'ig', ig
             ! ... g vectors are in units of 2pi / alat :
             ! ... to go to atomic units g must be multiplied by 2pi / alat
             !
             absg=tpiba*SQRT( gg(ig) )
             !
             rg = r*absg
             !
             IF ( r == 0.D0 .AND. absg /= 0 ) THEN
                sinxx = 1.D0
             ELSE IF ( absg == 0) THEN  
                sinxx = 0.5D0
             ELSE
                sinxx = SIN( rg ) / rg
             END IF
             !
             ! ... add the phase factor corresponding to the translation of the
             ! ... origin by x0 (notice that x0 is in alat units)
             !
             phase = tpi*( g(1,ig)*r0(1) + g(2,ig)*r0(2) + g(3,ig)*r0(3) )
             cavg = phi(nl(ig))*CMPLX( COS( phase ), SIN( phase ) ,kind=DP)
             !
             if (gamma_only) then
                ravg(ir,ibnd,1) = ravg(ir,ibnd,1) + 2.D0*REAL( cavg )*sinxx
             else
               ravg(ir,ibnd,1) = ravg(ir,ibnd,1) + REAL( cavg )*sinxx
             endif
             !
             !
          enddo
          CALL mp_sum( ravg(ir,ibnd,1), world_comm )
       enddo
     enddo
     !
  enddo
  !
  !WRITE( stdout, * ),'Finished writing wfc squared'
  !
  CALL deallocate_bec_type ( becp)
  deallocate(phi)
  deallocate(phi1)
  !
  !   WRITE( stdout, *), 'finished doing radial average for bands'
  !
  !
  !  WRITE (stdout, '(7i10)') nr1, nr2, nr3, nrx1, nrx2, nrx3, nrxx
  !
  ! DOS
  !

  ! IF UNDEFINED BY USER, SET TO SOMETHING LOGICAL
  ngauss = 0
  if (Emin.eq.-1000000.d0 .and. Emax .eq. 1000000.d0) THEN 
      !
      Elw = et (bndmin, 1) 
      Eup = et (bndmax, 1)  
      DO ik = 2, nks
         Elw = MIN (Elw, et (bndmin, ik) )
         Eup = MAX (Eup, et (bndmax, ik) )
      ENDDO
      IF (degauss.NE.0.d0) THEN
         Eup = Eup + 3.d0 * degauss
         Elw = Elw - 3.d0 * degauss
      ENDIF
      !
      Emin=MAX(Emin/rytoev,Elw) + mp_shift/rytoev 
      Emax=MIN(Emax/rytoev,Eup) + mp_shift/rytoev 
  ELSE 
      Emin=Emin/rytoev
      Emax=Emax/rytoev
  ENDIF
  DeltaE = DeltaE / rytoev
  ndos = NINT ( (Emax - Emin) / DeltaE+0.500001d0)
  ALLOCATE(DOSofE(ndos))
  DOSofE=0.d0
  !
  DO ir = 1,npts
   !
   ! first calculate smeared dos at each point along radial direction
   !
   WRITE(filename,'(i3.3)') ir    
   ! 
   IF ( ionode ) OPEN(200,file=TRIM(filplot)//'_rdos.'//filename,status='unknown')
   !
   DOSofE = 0.0d0
   DO i = 1, ndos
      E = Emin + (i - 1) * DeltaE
      !
      do ik = 1, nks
         do ibnd = bndmin, bndmax
            DOSofE(i) = DOSofE(i) + wk (ik) * ravg(ir,ibnd,ik) * w0gauss ((E-(et (ibnd, ik) + mp_shift/rytoev) )& 
                    / degauss, ngauss)
         enddo
      enddo
      !
      DOSofE(i) = DOSofE(i) / degauss
      !
      IF ( ionode ) WRITE (200, '(f7.3,3e12.4)') E * rytoev, DOSofE(i)/rytoev
      !
    ENDDO
    !
    ! INT\DOS(z,E)dE
    !
    intnorm = 0.0d0
    DO i = 1, ndos
       E = (Emin + (i - 1)*DeltaE)*rytoev
       IF (E.lt.Ef) THEN
       intnorm = intnorm + DOSofE(i)*DeltaE
       ENDIF
       IF (E.lt.Ef1) THEN
       k1 = i
       ENDIF
       IF (E.lt.Ef2) THEN
       k2 = i
       ENDIF
    ENDDO
    intnorm = intnorm*delta
    !
    ! VBM
    !
    w2 = 0.0d0
    DO i = k1,1,-1
       IF (w2.lt.intnorm) THEN
       w2  = w2 + DOSofE(i)*DeltaE
       ed1  = i
       Eed1 = (Emin + (i - 1)*DeltaE)*rytoev
       ENDIF
    ENDDO
    !
    ! CBM
    !
    w2 = 0.0d0
    DO i = k2, ndos, 1
       IF (w2.lt.intnorm) THEN
          w2  = w2 + DOSofE(i)*DeltaE
          ed2 = i
          Eed2 = (Emin + (i - 1)*DeltaE)*rytoev
       ENDIF
    ENDDO
    IF ( ionode ) THEN
       WRITE(iunplot,*) DBLE(ir)/DBLE(npts)*rmax, Eed1, Eed2
       CLOSE(200)
    ENDIF
    !
  ENDDO
  !
  ! free at least some of the allocated arrays...
  ! this is still not totally enough I guess (PWscf has its routines to free memory)
  !
  deallocate(ravg)
  deallocate(DOSofE)
  !
  CALL environment_end ( 'rad_avg' )
  !
  CALL stop_pp
  STOP
  !
END PROGRAM do_plan_avg

