!
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
! ... poorly written by N. P. Brawand
!
! Changes 01/02/2015 (ADC) : 
!               Stuff
!
!
!--------------------------------------------------------------------------
SUBROUTINE add_efield(vpoten,etotefield,rho,iflag)
  !--------------------------------------------------------------------------
  !
  !   This routine adds stuff to the local potential. 
  !
  !   edir - atom to center potential well around 
  !   emaxpos - radius of potential well
  !   eamp - strength of potential in Ry a.u.
  !
  !
  !
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8, e2, au_debye
  USE ions_base,     ONLY : nat, ityp, zv
  USE cell_base,     ONLY : alat, at, omega, bg, saw
  USE extfield,      ONLY : tefield, dipfield, edir, eamp, emaxpos, &
                            eopreg, forcefield
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
  INTEGER :: ir, na, ipol
  REAL(DP) :: length, vamp, value, sawarg, e_dipole, ion_dipole
  REAL(DP) :: tot_dipole, bmod
  
  LOGICAL :: first=.TRUE.
  SAVE first
  !
  ! ... Coulomb Vars
  !
  INTEGER      :: ip
  REAL( DP )   :: dist
  REAL( DP )   :: r( 3 ), s( 3 ), cm( 3 )
  REAL( DP )   :: inv_nr1, inv_nr2, inv_nr3
  inv_nr1 = 1.D0 / DBLE( dfftp%nr1 )
  inv_nr2 = 1.D0 / DBLE( dfftp%nr2 )
  inv_nr3 = 1.D0 / DBLE( dfftp%nr3 )
  
  !---------------------
  !  Execution control
  !---------------------

  IF (.NOT.tefield) RETURN
  ! efield only needs to be added on the first iteration, if dipfield
  ! is not used. note that for relax calculations it has to be added
  ! again on subsequent relax steps.
  IF ((.NOT.dipfield).AND.(.NOT.first) .AND..NOT. iflag) RETURN
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
    WRITE( stdout,'(8x,"Amplitude [Ry a.u.] : ", es11.4)') eamp 
    WRITE( stdout,'(8x,"Postion on atom # : ", I11.1)') edir
    WRITE( stdout,'(8x,"Well radius [bohr] : ", es11.4)') emaxpos
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
  !
  ! ... set center of well on edir atom
  !
  cm(:) = cm(:) + tau(:,edir)
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
    ! translate potential well to atom center
    !
    r(:) = r(:) - cm(:)
    !
    ! ... minimum image convention for periodic BC
    !
    s(:) = MATMUL( r(:), bg(:,:) )
    s(:) = s(:) - ANINT(s(:))
    r(:) = MATMUL( at(:,:), s(:) )
    !
    dist = SQRT( SUM( r * r ) ) 
    !
    ! if within emaxpos add potential 
    !
    IF(dist <= emaxpos) THEN
      !
      vpoten(ir) = vpoten(ir) + eamp
      !
    ENDIF
    !
  END DO 
  !
  RETURN
  !
END SUBROUTINE add_efield
