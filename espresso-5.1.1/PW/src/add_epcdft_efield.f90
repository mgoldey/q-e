!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... poorly written by N. P. Brawand
! ... and poorly modified by M. B. Goldey
!
!
!--------------------------------------------------------------------------
SUBROUTINE add_epcdft_efield(vpoten,etotefield,rho,iflag)
  !--------------------------------------------------------------------------
  !
  !   This routine adds constraining potential for cdft to the local potential. 
  !
  !   Voronoi cells:
  !
  !     acceptor_start - first atom in fragment to center potential well around 
  !     acceptor_end - last atom  in fragment to center potential well around 
  !     donor_start - first atom in fragment in donor
  !     donor_end - last atom  in fragment in donor 
  !     epcdft_amp - strength of potential in Ry a.u.
  !     epcdft_charge - number of electrons that should be in well
  !
  !   User defined well around single atom:
  !
  !     acceptor_start - only atom to center potential well around 
  !     acceptor_end - = 0
  !     donor_start - do not use
  !     donor_end - do not use
  !     epccdft_width  - radius of potential well in alat
  !     epcdft_amp - strength of potential in Ry a.u.
  !     epcdft_charge - number of electrons that should be in well
  !
  !   Hirshfeld partitioning following :
  !
  !     J. Chem. Phys. 133, 244105 (2010); http://dx.doi.org/10.1063/1.3507878 
  !     eqs. 6 & 7
  !
  !     acceptor_start - first atom in fragment in acceptor
  !     acceptor_end - last atom  in fragment in acceptor (all others are donors)
  !     donor_start - first atom in fragment in donor
  !     donor_end - last atom  in fragment in donor 
  !     epcdft_amp - strength of potential in Ry a.u. 
  !     epcdft_charge - ?? not sure yet
  !
  !
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8, e2, au_debye
  USE ions_base,     ONLY : nat, ityp, zv
  USE cell_base,     ONLY : alat, at, omega, bg, saw
  USE epcdft,        ONLY : do_epcdft, donor_start, donor_end, &
                            acceptor_start, acceptor_end, &
                            epcdft_amp, epcdft_width, epcdft_shift, &
                            epcdft_charge, hirshfeld
  USE force_mod,     ONLY : lforce
  USE io_global,     ONLY : stdout,ionode
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_images,     ONLY : intra_image_comm
  USE mp_bands,      ONLY : me_bgrp, intra_bgrp_comm
  USE fft_base,      ONLY : dfftp, grid_scatter, grid_gather
  USE mp,            ONLY : mp_bcast, mp_sum
  USE control_flags, ONLY : iverbosity
  !
  ! ... Coulomb USE
  !
  USE ions_base,     ONLY : tau
  
  IMPLICIT NONE
  !
  ! I/O variables
  !
  REAL(DP),INTENT(INOUT) :: vpoten(dfftp%nnr)! ef is added to this potential
  REAL(DP) :: dv                             ! volume element
  REAL(DP) :: einwellp                       ! number of electrons in well
  REAL(DP) :: einwells                       ! number of electrons in well
  REAL(DP) :: enumerr                        ! epcdft_charge - einwell  (e number error)
  REAL(DP) :: oldamp                         ! stores old external potential strength
  SAVE oldamp
  LOGICAL  :: elocflag                       ! true if charge localization condition is satisfied
  REAL(DP),INTENT(INOUT) :: etotefield       ! contribution to etot due to ef
  REAL(DP),INTENT(IN):: rho(dfftp%nnr,nspin) ! the density whose dipole is computed
  LOGICAL,INTENT(IN)     :: iflag            ! set to true to force recalculation of field
  !
  ! local variables
  !
  INTEGER :: index0, i, j, k
  INTEGER :: ir, na, ipol, iatom
  REAL(DP) :: length, vamp, value, sawarg, e_dipole, ion_dipole
  REAL(DP) :: tot_dipole, bmod
  INTEGER :: icfd(-1:1), ir_end, in 
  INTEGER :: ix(-1:1),iy(-1:1),iz(-1:1)
  REAL( DP ), DIMENSION( :, : ), ALLOCATABLE :: gradtmp
  REAL( DP ), DIMENSION( :, : ), ALLOCATABLE :: grad
  REAL(DP) :: hirshv(dfftp%nnr) ! constraining potential from hirshfeld partitioning
  !
  LOGICAL :: first=.TRUE.
  SAVE first

  !
  ! ... Coulomb Vars
  !
  INTEGER      :: ip
  REAL( DP )   :: dist, mindonor, minacceptor
  REAL( DP )   :: r( 3 ), myr(3), s( 3 ), cm(3)
  REAL( DP )   :: inv_nr1, inv_nr2, inv_nr3
  REAL( DP )   :: thresh
  
  LOGICAL :: on_donor    = .TRUE.
  LOGICAL :: on_acceptor = .TRUE.

  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  einwellp = 0.D0
  einwells = 0.D0
  enumerr  = 0.D0 
  oldamp = epcdft_amp

  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )

  !---------------------
  !  Execution control
  !---------------------
  IF (.NOT. do_epcdft) RETURN
  IF (.NOT. iflag) RETURN  !TURN OFF SELF-CONSISTENCY INSIDE SCF CYCLE OR ELSE!
  if (iflag) first=.true.
   
  if (.not. first) RETURN

  ! Necessary for restart/pp runs
  CALL init_at_1
  
  ! efield only needs to be added on the first iteration (of each SCF call)
  ! note that for relax calculations it has to be added
  ! again on subsequent relax steps. 
  !
  !---------------------
  !  Variable initialization
  !---------------------
  !
  cm(:) = 0.D0
  !
  ! print potential type used for cdft
  !
  IF (ionode) THEN
    !
    WRITE( stdout,*)
    !
    IF (hirshfeld) THEN
      !
      WRITE( stdout,'(5x,"Using Hirshfeld")')
      WRITE( stdout,'(8x,"Amplitude [Ry a.u.] : ", es11.4)') epcdft_amp 
      WRITE( stdout,'(8x,"Acceptor Fragment start : ", I11.1)') acceptor_start
      WRITE( stdout,'(8x,"Acceptor Fragment end   : ", I11.1)') acceptor_end
      WRITE( stdout,'(8x,"Donor Fragment start : ", I11.1)') donor_start
      WRITE( stdout,'(8x,"Donor Fragment end   : ", I11.1)') donor_end
      WRITE( stdout,*)     
      !
    ELSE ! not hirshfeld then voronoi or user well
      !
      if (acceptor_end .ne. 0) then
        !
        write( stdout,'(5x,"Using Voronoi cells")')
        WRITE( stdout,'(5x,"Adding the potential well":)')
        WRITE( stdout,'(8x,"Amplitude [Ry a.u.] : ", es11.4)') epcdft_amp 
        WRITE( stdout,'(8x,"Acceptor Fragment start : ", I11.1)') acceptor_start
        WRITE( stdout,'(8x,"Acceptor Fragment end   : ", I11.1)') acceptor_end
        WRITE( stdout,'(8x,"Donor Fragment start : ", I11.1)') donor_start
        WRITE( stdout,'(8x,"Donor Fragment end   : ", I11.1)') donor_end
        WRITE( stdout,*)     
        !
      else
        !
        write( stdout,'(5x,"Using a distance-based cell")')
        WRITE( stdout,'(5x,"Adding the potential well":)')
        WRITE( stdout,'(8x,"Amplitude [Ry a.u.] : ", es11.4)') epcdft_amp 
        WRITE( stdout,'(8x,"Atom selected : ", I11.1)') acceptor_start
        WRITE( stdout,'(8x,"Well size     : ", es11.4)') epcdft_width
        WRITE( stdout,*)     
        !
      endif ! end if voronoi or user well
      !
    ENDIF !end if hirshfeld
    !
  ENDIF 
  !
  ! end print
  !
  ! calculate hirshfeld potential array
  !
  IF (hirshfeld) THEN 
    hirshv=0.D0
    !CALL calc_hirshfeld_v_pointlists(hirshv, dfftp%nnr)
    CALL calc_hirshfeld_v(hirshv, dfftp%nnr)
    vpoten = vpoten + epcdft_amp * hirshv
    !write(*,*) "vpoten is ",sum(vpoten(:))
    RETURN
  ENDIF
  !
  !
  !
  index0 = 0
#if defined (__MPI)
  DO i = 1, me_bgrp
     index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  END DO
  ir_end = MIN(dfftp%nnr,dfftp%nr1x*dfftp%nr2x*dfftp%npl)
#else
  ir_end = dfftp%nnr
#endif
  
  if (sum(rho).lt.1e-3) THEN
  !  write(*,*) "Density is really small. I forget why this matters."
  ENDIF

  ! voronoi
  if (acceptor_end .ne. 0) then
    DO ir = 1, dfftp%nnr
      i = index0 + ir - 1
      k = i / (dfftp%nr1x*dfftp%nr2x)
      i = i - (dfftp%nr1x*dfftp%nr2x)*k
      j = i / dfftp%nr1x
      i = i - dfftp%nr1x*j
      DO ip = 1, 3
        r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                DBLE( j )*inv_nr2*at(ip,2) + &
                DBLE( k )*inv_nr3*at(ip,3)
      END DO
      
      mindonor=5d6
      minacceptor=5d6
      DO iatom= acceptor_start, acceptor_end
        cm(:) = tau(:,iatom)
        myr(:) = r(:) - cm(:)
        s(:) = MATMUL( myr(:), bg(:,:) )
        s(:) = s(:) - ANINT(s(:))
        myr(:) = MATMUL( at(:,:), s(:) )
        dist = SQRT( SUM( myr * myr ) )*alat
        IF(dist .le. minacceptor) THEN
          minacceptor=dist
        END IF
      END DO
      DO iatom= donor_start, donor_end
        cm(:) = tau(:,iatom)
        myr(:) = r(:) - cm(:)
        s(:) = MATMUL( myr(:), bg(:,:) )
        s(:) = s(:) - ANINT(s(:))
        myr(:) = MATMUL( at(:,:), s(:) )
        dist = SQRT( SUM( myr * myr ) )*alat
        IF(dist .le. mindonor) THEN
          mindonor=dist
        END IF
      END DO
      on_donor = .TRUE.
      on_acceptor = .TRUE.
      !write(*,*)" mindist are ",mindonor," ",minacceptor
      DO iatom=1, nat
        IF ((.not. on_donor).and. (.not. on_acceptor)) CYCLE
        IF ((iatom.ge.acceptor_start) .AND. (iatom.le.acceptor_end) ) CYCLE
        IF ((iatom.ge.donor_start) .AND. (iatom.le.donor_end) ) CYCLE
        cm(:) = tau(:,iatom)
        myr(:) = r(:) - cm(:)
        s(:) = MATMUL( myr(:), bg(:,:) )
        s(:) = s(:) - ANINT(s(:))
        myr(:) = MATMUL( at(:,:), s(:) )
        dist = SQRT( SUM( myr * myr ) )*alat
        IF(dist .lt. mindonor) on_donor = .FALSE.
        IF(dist .lt. minacceptor) on_acceptor=.FALSE.        
      END DO
      if (mindonor.lt.minacceptor) on_acceptor=.FALSE.
      if (minacceptor.lt.mindonor) on_donor=.FALSE.
      !write(*,*)" bools are ",on_donor," ",on_acceptor
      IF(on_donor) vpoten(ir) = vpoten(ir) - epcdft_amp
      IF(on_acceptor) vpoten(ir) = vpoten(ir) + epcdft_amp
    END DO 
  else ! APPLY POTENTIAL WITHIN WELL OF SIZE EPCDFT_WIDTH
    DO ir = 1, dfftp%nnr
      i = index0 + ir - 1
      k = i / (dfftp%nr1x*dfftp%nr2x)
      i = i - (dfftp%nr1x*dfftp%nr2x)*k
      j = i / dfftp%nr1x
      i = i - dfftp%nr1x*j
      DO ip = 1, 3
        r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
                DBLE( j )*inv_nr2*at(ip,2) + &
                DBLE( k )*inv_nr3*at(ip,3)
      END DO
      iatom= acceptor_start
      cm(:) = tau(:,iatom)
      myr(:) = r(:) - cm(:)
      s(:) = MATMUL( myr(:), bg(:,:) )
      s(:) = s(:) - ANINT(s(:))
      myr(:) = MATMUL( at(:,:), s(:) )
      dist = SQRT( SUM( myr * myr ) )*alat
      IF(dist .le. epcdft_width) THEN
        vpoten(ir) = vpoten(ir) - epcdft_amp
      ELSE
        vpoten(ir) = vpoten(ir) + epcdft_amp
      ENDIF
    END DO 
  ENDIF
  RETURN
  !
END SUBROUTINE add_epcdft_efield
!
!
!
!--------------------------------------------------------------------------
SUBROUTINE calc_hirshfeld_v( v, n )
  !--------------------------------------------------------------------------
  ! 
  !  Gamma only
  !
  !  calculate hirshfeld potential and put into v following :
  !
  !     J. Chem. Phys. 133, 244105 (2010); http://dx.doi.org/10.1063/1.3507878 
  !     eqs. 6 & 7
  !
  USE kinds,                ONLY : DP
  USE control_flags,        ONLY : gamma_only
  USE basis,                ONLY : natomwfc
  USE wvfct,                ONLY : npw, npwx, ecutwfc, igk, g2kin
  USE klist,                ONLY : xk, nks
  USE gvect,                ONLY : ngm, g
  USE cell_base,            ONLY : tpiba2
  USE uspp_param,           ONLY : upf
  USE lsda_mod,      ONLY : nspin
  USE ions_base,            ONLY : nat, ityp
  USE epcdft,               ONLY : donor_start,donor_end,&
                                   acceptor_start,acceptor_end
  USE fft_base,             ONLY : dfftp, dffts
  USE gvect,                ONLY : nl, nlm
  USE gvecs,                ONLY : nls, nlsm
  USE fft_interfaces,       ONLY : invfft
  USE wavefunctions_module, ONLY : psic
  USE mp_bands,     ONLY : intra_bgrp_comm
  USE mp,            ONLY : mp_bcast, mp_sum
  USE io_global,     ONLY : stdout,ionode, ionode_id
  USE klist,         ONLY : nelec
  USE cell_base,     ONLY : omega
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: v(n) ! the hirshfeld potential
  INTEGER, INTENT(IN) :: n
  ! 
  INTEGER :: s, na, nt, nb, l, m,ir
  COMPLEX(DP), ALLOCATABLE :: wfcatomg (:,:) ! atomic wfcs in g
  COMPLEX(DP), ALLOCATABLE :: wfcatomr (:,:) ! atomic wfcs in r
  INTEGER :: orbi ! orbital index for wfcatom
  REAL(DP) :: orboc ! occupation of orbital
  COMPLEX(DP) :: vtop(n) ! top of the hirshfeld potential fraction
  COMPLEX(DP) :: vbot(n) ! bottom of the hirshfeld potential fraction
  COMPLEX(DP) :: normfac
  COMPLEX(DP) :: cutoff
  COMPLEX(DP) :: vbottot
  REAL(DP) :: dv
  !
  ALLOCATE( wfcatomg(npwx, natomwfc) )
  ALLOCATE( wfcatomr(n, natomwfc) )
  !
  wfcatomg = 0.D0
  wfcatomr = 0.D0
  vtop = 0.D0
  vbot = 0.D0
  v = 0.D0
  psic = 0.D0
  orbi = 0 
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  !
  ! load atomic wfcs
  CALL gk_sort (xk (1, 1), ngm, g, ecutwfc / tpiba2, npw, igk, g2kin)
  CALL atomic_wfc (1, wfcatomg) 
  !
  ! convert atomic(g) -> |atomic(r)|^2
  DO s = 1, natomwfc
    !
    wfcatomr(nl(igk(1:npw)),s)  = wfcatomg(1:npw,s)
    IF(gamma_only) wfcatomr(nlm(igk(1:npw)),s) = CONJG( wfcatomg(1:npw,s) )
    !
    CALL invfft ('Dense', wfcatomr(:,s), dfftp)
    wfcatomr(:,s) = wfcatomr(:,s) * CONJG(wfcatomr(:,s))

    !normfac=sum(wfcatomr(:,s))
    !wfcatomr(:,s)=wfcatomr(:,s)/sqrt(normfac)
    !
    ! Prune low values away
    !
    !DO ir=1, n
    !  if (abs(wfcatomr(ir,s)).lt.1d-6) wfcatomr(ir,s)=0.D0
    !ENDDO  
    !
  ENDDO
  !
  ! construct hirshfeld looping over all atomic states
  !
  !
  DO na = 1, nat ! for each atom
    nt = ityp (na) ! get atom type
    DO nb = 1, upf(nt)%nwfc ! for each orbital
      l = upf(nt)%lchi(nb) ! get l 
      DO m = 1, 2*l+1 ! mag num
        !  
        ! add all orbitals 
        ! each state weighted by orboc
        !  
        orbi = orbi + 1 ! update orbital index
        !  
        orboc = REAL( upf(nt)%oc(nb) ,KIND=DP) / REAL( 2*l+1 ,KIND=DP)
        !
        IF( na >= acceptor_start .and. na <= acceptor_end )THEN ! atom in acceptor
          !
          vtop(:) = vtop(:) - orboc * wfcatomr( : , orbi ) 
          !
        ELSE IF ( na >= donor_start .and. na <= donor_end )THEN ! atom in donor
          ! 
          vtop(:) = vtop(:) + orboc * wfcatomr( : , orbi ) 
          !
        ENDIF 
        !
        vbot(:) = vbot(:) + orboc * wfcatomr( : , orbi ) 
        !
        !
        !
      ENDDO ! m
    ENDDO ! l 
  ENDDO ! atom
  !
  ! M.G.s bar
  !CALL mp_sum( orboc, intra_image_comm )
  !
  !call write_cube_r ( 84332, 'vtop.cube',  REAL(vtop,KIND=DP))
  !call write_cube_r ( 84332, 'vbot.cube',  REAL(vbot,KIND=DP))
  !
  ! force normalization
  !
  vbottot = SUM(vbot)
  CALL mp_sum( vbottot, intra_bgrp_comm )
  normfac=REAL(nelec,KIND=DP)/(REAL(vbottot,KIND=DP)*dv)
  vtop = normfac * vtop
  vbot = normfac * vbot
  !
  cutoff = 1.D-9
  vtop = vtop / vbot
  DO ir = 1, n
    if (ABS(REAL(vbot(ir))).lt.REAL(cutoff)) vtop(ir)=0.D0  
    if (vtop(ir) /= vtop(ir)) vtop(ir)=0.D0
  ENDDO
  !
  !
  !call write_cube_r ( 84332, 'v.cube',  REAL(vtop,KIND=DP))
  v(:) = REAL(vtop(:),KIND=DP)
  !
  !call write_cube_r ( 84332, 'v.cube',  v )
  !
  DEALLOCATE( wfcatomg )
  DEALLOCATE( wfcatomr )
  !
  RETURN
  !
END SUBROUTINE calc_hirshfeld_v
!
!
!--------------------------------------------------------------------------
SUBROUTINE calc_hirshfeld_v_pointlists( v, n )
  !--------------------------------------------------------------------------
  ! 
  !  Gamma only
  !
  !  calculate hirshfeld potential and put into v following :
  !
  !     J. Chem. Phys. 133, 244105 (2010); http://dx.doi.org/10.1063/1.3507878 
  !     eqs. 6 & 7
  !
  USE kinds,            ONLY : DP
  USE mp_bands,         ONLY : intra_bgrp_comm
  USE noncollin_module, ONLY : pointlist, factlist
  USE fft_base,         ONLY : dfftp
  USE cell_base,        ONLY : omega
  USE lsda_mod,         ONLY : nspin
  USE mp,               ONLY : mp_sum
  USE klist,            ONLY : nelec
  USE epcdft,           ONLY : donor_start,donor_end,&
                               acceptor_start,acceptor_end
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: v(n) ! the hirshfeld potential
  INTEGER, INTENT(IN) :: n        !dfftp%nnr
  ! 
  INTEGER :: s, na, nt, ir
  REAL(DP) :: vtop(n)             ! top of the hirshfeld potential fraction
  REAL(DP) :: vbot(n)             ! bottom of the hirshfeld potential fraction
  REAL(DP) :: dv                  ! volume element
  REAL(DP) :: vbot_tot
  REAL(DP) :: rho(n,nspin)        ! atomic rho
  REAL(DP) :: cutoff
  ! 
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  ! 
  ! calculate atomic charge
  ! 
  print *, " is atomic rho saved somewhere? recalcing it now"
  CALL atomic_rho (rho, nspin)
  !
  ! renormalize charge
  !
  vbot_tot = SUM(rho(:,1:nspin)) * dv
  CALL mp_sum( vbot_tot, intra_bgrp_comm )
  rho = rho * REAL(nelec,KIND=DP)/vbot_tot
  DO s = 1, nspin
    rho(:,s) = rho(:,s) * factlist(:) * dv
  ENDDO
  ! 
  ! calc hirsh 
  ! 
  vbot = 0.D0
  vtop = 0.D0
  DO s = 1, nspin
    DO ir = 1, n
      ! 
      ! top part of hirsh 
      ! 
      IF( pointlist(ir) >= acceptor_start .and. pointlist(ir) <= acceptor_end )THEN ! acceptor
        !
        vtop(ir) = vtop(ir) - rho(ir,s) 
        !
      ELSE IF ( pointlist(ir) >= donor_start .and. pointlist(ir) <= donor_end )THEN ! donor
        ! 
        vtop(ir) = vtop(ir) + rho(ir,s)
        !
      ENDIF 
      ! 
      ! top part of hirsh 
      ! 
      vbot(ir) = vbot(ir) + rho(ir,s)
      ! 
    ENDDO
  ENDDO
  !
  ! set to zero places where atomic density is low
  !
  cutoff = 1.D-9
  vtop = vtop / vbot
  DO ir = 1, n
    IF ( vbot(ir) .lt. cutoff ) vtop(ir) = 0.D0  
    IF ( vtop(ir) /= vtop(ir) ) vtop(ir) = 0.D0
  ENDDO
  !
  v = vtop
  !
END SUBROUTINE calc_hirshfeld_v_pointlists
!
