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
  !   edir - atom to center potential well around 
  !   emaxpos - radius of potential well in alat
  !   eamp - strength of potential in Ry a.u.
  !   eopreg - number of electrons that should be in well
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
                            epcdft_amp, epcdft_shift
  USE force_mod,     ONLY : lforce
  USE io_global,     ONLY : stdout,ionode
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_images,     ONLY : intra_image_comm
  USE mp_bands,      ONLY : me_bgrp
  USE fft_base,      ONLY : dfftp
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
  REAL(DP),INTENT(INOUT) :: etotefield       ! contribution to etot due to ef
  REAL(DP),INTENT(IN)    :: rho(dfftp%nnr,nspin) ! the density whose dipole is computed
  LOGICAL,INTENT(IN)     :: iflag ! set to true to force recalculation of field
  !
  ! local variables
  !
  INTEGER :: index0, i, j, k
  INTEGER :: ir, na, ipol, iatom
  REAL(DP) :: length, vamp, value, sawarg, e_dipole, ion_dipole
  REAL(DP) :: tot_dipole, bmod
  
  LOGICAL :: first=.TRUE.
  SAVE first
  !
  ! ... Coulomb Vars
  !
  INTEGER      :: ip
  REAL( DP )   :: dist, mindist
  REAL( DP )   :: r( 3 ), myr(3), s( 3 ), cm(3)
  REAL( DP )   :: inv_nr1, inv_nr2, inv_nr3
  
  LOGICAL :: on_frag = .TRUE.

  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
  
  !---------------------
  !  Execution control
  !---------------------

  IF (.NOT. do_epcdft) RETURN
  IF ((.NOT.first) .AND..NOT. iflag) RETURN
  ! efield only needs to be added on the first iteration, but
  ! this should be changed to be self-consistent soon -MBG
  ! note that for relax calculations it has to be added
  ! again on subsequent relax steps. 
  first=.FALSE.

  !---------------------
  !  Variable initialization
  !---------------------
  cm(:) = 0.D0
  !
  ! write the input variables
  !
  IF (ionode) THEN
    !
    WRITE( stdout,*)
    WRITE( stdout,'(5x,"Adding potential well":)')
    WRITE( stdout,'(8x,"Amplitude [Ry a.u.] : ", es11.4)') epcdft_amp 
    WRITE( stdout,'(8x,"Fragment start : ", I11.1)') fragment_atom1
    WRITE( stdout,'(8x,"Fragment end : ", I11.1)') fragment_atom2
    WRITE( stdout,*)     
    !
  ENDIF
  !
  !---------------------
  !  Add potential
  !---------------------
  !
  ! Index for parallel summation
  !
  index0 = 0
#if defined (__MPI)
  !
  DO i = 1, me_bgrp
     index0 = index0 + dfftp%nr1x*dfftp%nr2x*dfftp%npp(i)
  END DO
  !
#endif

  write(*,*) "atom data ",nat, fragment_atom1, fragment_atom2
  !
  ! ... Loop over position grid
  !
  
  DO ir = 1, dfftp%nnr
    !
    ! ... three dimensional indexes
    !
    i = index0 + ir - 1
    k = i / (dfftp%nr1x*dfftp%nr2x)
    i = i - (dfftp%nr1x*dfftp%nr2x)*k
    j = i / dfftp%nr1x
    i = i - dfftp%nr1x*j
    !
    ! calculate position
    DO ip = 1, 3
      !
      r(ip) = DBLE( i )*inv_nr1*at(ip,1) + &
              DBLE( j )*inv_nr2*at(ip,2) + &
              DBLE( k )*inv_nr3*at(ip,3)
      !
    END DO
    !
    mindist=5.D0
    DO iatom= fragment_atom1, fragment_atom2
      ! write(*,*) "atom ",iatom, fragment_atom1, fragment_atom2
      !
      ! GRAB iatom xyzs
      !
      cm(:) = tau(:,iatom)
      ! get 3-space (r) vector between position and atom center
      !
      myr(:) = r(:) - cm(:)
      !
      ! ... minimum image convention for periodic BC
      !
      s(:) = MATMUL( myr(:), bg(:,:) )
      s(:) = s(:) - ANINT(s(:))
      myr(:) = MATMUL( at(:,:), s(:) )
      !
      dist = SQRT( SUM( myr * myr ) ) 
      !
      IF(dist .le. mindist) THEN
        !
        mindist=dist
        !
      END IF
    END DO

    on_frag = .TRUE.

    DO iatom=1, nat
      !
      ! CYCLE if wrong set of atoms
      !
      IF ((iatom.ge.fragment_atom1) .AND. (iatom.le.fragment_atom2) ) THEN
        ! write(*,*) "Skipping this atom"
        CYCLE
      END IF
      !
      ! GRAB iatom xyzs
      !
      cm(:) = tau(:,iatom)
      !
      ! Get 3-space (r) vector between position and atom center
      !
      myr(:) = r(:) - cm(:)
      !
      ! ... minimum image convention for periodic BC
      !
      s(:) = MATMUL( myr(:), bg(:,:) )
      s(:) = s(:) - ANINT(s(:))
      myr(:) = MATMUL( at(:,:), s(:) )
      !
      dist = SQRT( SUM( myr * myr ) ) 
      !
      IF(dist .lt. mindist) THEN
        !
        mindist = dist
        on_frag = .FALSE.
        ! write(*,*) "Here"
        EXIT
        !
      ENDIF
      !
    END DO
    ! if within emaxpos add potential 
    !
    IF(on_frag) THEN
      !
      ! write(*,*) "Applying the potential",mindist
      vpoten(ir) = vpoten(ir) + epcdft_amp
      !
    ELSE
      ! write(*,*) "Not applying",mindist
    ENDIF
    !
  END DO 
  !
  RETURN
  !
END SUBROUTINE add_efield
