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
SUBROUTINE add_efield(vpoten,etotefield,rho,iflag)
  !--------------------------------------------------------------------------
  !
  !   This routine adds stuff to the local potential. 
  !
  !   fragment_atom1 - first atom in fragment to center potential well around 
  !   fragment_atom2 - last atom  in fragment to center potential well around 
  !   epccdft_width  - radius of potential well in alat
  !   epcdft_amp - strength of potential in Ry a.u.
  !   epcdft_electrons - number of electrons that should be in well
  !
  !
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8, e2, au_debye
  USE ions_base,     ONLY : nat, ityp, zv
  USE cell_base,     ONLY : alat, at, omega, bg, saw
  USE extfield,      ONLY : tefield, dipfield, edir, eamp, emaxpos, &
                            eopreg, forcefield
  USE epcdft,        ONLY : do_epcdft, fragment_atom1, fragment_atom2, &
                            epcdft_amp, epcdft_width, epcdft_shift, &
                            epcdft_electrons
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
  REAL(DP),INTENT(INOUT) :: vpoten(dfftp%nnr,nspin)! ef is added to this potential
  REAL(DP) :: dv                             ! volume element
  REAL(DP) :: einwellp                       ! number of electrons in well
  REAL(DP) :: einwells                       ! number of electrons in well
  REAL(DP) :: enumerr                        ! epcdft_electrons - einwell  (e number error)
  REAL(DP) :: oldamp                         ! stores old external potential strength
  SAVE oldamp
  LOGICAL  :: elocflag                       ! true if charge localization condition is satisfied
  REAL(DP),DIMENSION(:), ALLOCATABLE :: tmpv ! temporary copy of v to check number of electrons
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

  LOGICAL :: first=.TRUE.
  SAVE first
  !
  ! ... Coulomb Vars
  !
  INTEGER      :: ip
  REAL( DP )   :: dist, mindist
  REAL( DP )   :: r( 3 ), myr(3), s( 3 ), cm(3)
  REAL( DP )   :: inv_nr1, inv_nr2, inv_nr3
  REAL( DP )   :: thresh
  
  LOGICAL :: on_frag = .TRUE.

  ALLOCATE(tmpv(dfftp%nnr))
  tmpv=0.D0
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
  !write(*,*) "do_epcdft is ",do_epcdft
  !write(*,*) "tefield is ",tefield

  !WRITE( stdout,'(8x,"Amplitude [Ry a.u.] : ", es11.4)') epcdft_amp 
  !WRITE( stdout,'(8x,"Fragment start : ", I11.1)') fragment_atom1
  !WRITE( stdout,'(8x,"Fragment end   : ", I11.1)') fragment_atom2
  !WRITE( stdout,'(8x,"Well size     : ", es11.4)') epcdft_width

  ! write(*,*) iflag, first, sum(rho(:,1)), sum(vpoten(:,1))
  IF (.NOT. do_epcdft) RETURN
  IF (.NOT. iflag) RETURN  !TURN OFF SELF-CONSISTENCY INSIDE SCF CYCLE OR ELSE!
  if (iflag) first=.true.
   
  if (.not. first) RETURN
  
  ! efield only needs to be added on the first iteration (of each SCF call)
  ! note that for relax calculations it has to be added
  ! again on subsequent relax steps. 

  !---------------------
  !  Variable initialization
  !---------------------
  cm(:) = 0.D0
  IF (ionode) THEN
    !
    WRITE( stdout,*)
    if (fragment_atom2 .ne. 0) then
      write( stdout,'(5x,"Using Voronoi cells")')
      WRITE( stdout,'(5x,"Adding the potential well":)')
      WRITE( stdout,'(8x,"Amplitude [Ry a.u.] : ", es11.4)') epcdft_amp 
      WRITE( stdout,'(8x,"Fragment start : ", I11.1)') fragment_atom1
      WRITE( stdout,'(8x,"Fragment end   : ", I11.1)') fragment_atom2
      WRITE( stdout,*)     
    else
      write( stdout,'(5x,"Using a distance-based cell")')
      WRITE( stdout,'(5x,"Adding the potential well":)')
      WRITE( stdout,'(8x,"Amplitude [Ry a.u.] : ", es11.4)') epcdft_amp 
      WRITE( stdout,'(8x,"Atom selected : ", I11.1)') fragment_atom1
      WRITE( stdout,'(8x,"Well size     : ", es11.4)') epcdft_width
      WRITE( stdout,*)     
    endif
  ENDIF
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
  if (fragment_atom2 .ne. 0) then
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
      mindist=5.D6
      DO iatom= fragment_atom1, fragment_atom2
        cm(:) = tau(:,iatom)
        myr(:) = r(:) - cm(:)
        s(:) = MATMUL( myr(:), bg(:,:) )
        s(:) = s(:) - ANINT(s(:))
        myr(:) = MATMUL( at(:,:), s(:) )
        dist = SQRT( SUM( myr * myr ) )*alat
        IF(dist .le. mindist) THEN
          mindist=dist
        END IF
      END DO
      on_frag = .TRUE.
      DO iatom=1, nat
        IF ((iatom.ge.fragment_atom1) .AND. (iatom.le.fragment_atom2) ) THEN
          CYCLE
        END IF
        cm(:) = tau(:,iatom)
        myr(:) = r(:) - cm(:)
        s(:) = MATMUL( myr(:), bg(:,:) )
        s(:) = s(:) - ANINT(s(:))
        myr(:) = MATMUL( at(:,:), s(:) )
        dist = SQRT( SUM( myr * myr ) )*alat
        IF(dist .lt. mindist) THEN
          mindist = dist
          on_frag = .FALSE.
          EXIT
        ENDIF
      END DO
      IF(on_frag) THEN
          vpoten(ir,1) = vpoten(ir,1) + epcdft_amp
          if (nspin.ne.1) vpoten(ir,2) = vpoten(ir,2) + epcdft_amp
      ENDIF
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
      iatom= fragment_atom1
      cm(:) = tau(:,iatom)
      myr(:) = r(:) - cm(:)
      s(:) = MATMUL( myr(:), bg(:,:) )
      s(:) = s(:) - ANINT(s(:))
      myr(:) = MATMUL( at(:,:), s(:) )
      dist = SQRT( SUM( myr * myr ) )*alat
      IF(dist .le. epcdft_width) THEN
        vpoten(ir,1) = vpoten(ir,1) + epcdft_amp
          if (nspin.ne.1) vpoten(ir,2) = vpoten(ir,2) + epcdft_amp        
      ENDIF
    END DO 
  ENDIF

#ifdef __MPI
  CALL MP_SUM(einwellp,intra_image_comm) ! FIX ME?
  einwells=einwellp
#else
  einwells=einwellp
#endif
  
  !dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )

  RETURN

  oldamp = epcdft_amp
  IF (ionode) WRITE(*,*)"    #e's   in well     : ",einwells," electrons"
  enumerr = epcdft_electrons - einwells
  epcdft_amp = epcdft_amp - enumerr * ABS(epcdft_amp)*.1
  IF(ionode) WRITE(*,*)"    New field Amp      : ",epcdft_amp," Ry" 
  !
  RETURN
  !
END SUBROUTINE add_efield
