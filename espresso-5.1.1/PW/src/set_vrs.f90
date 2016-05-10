!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!--------------------------------------------------------------------
subroutine set_vrs (vrs, vltot, vr, kedtau, kedtaur,nrxx, nspin, doublegrid)
  !--------------------------------------------------------------------
  ! set the total local potential vrs on the smooth mesh to be used in 
  ! h_psi, adding the (spin dependent) scf (H+xc) part and the sum of 
  ! all the local pseudopotential contributions.
  !
  USE kinds
  USE funct, only : dft_is_meta
  USE fft_base, only : dffts 
  implicit none

  integer :: nspin, nrxx
  ! input: number of spin components: 1 if lda, 2 if lsd, 4 if noncolinear
  ! input: the fft grid dimension
  real(DP) :: vrs (nrxx, nspin), vltot (nrxx), vr (nrxx, nspin), &
              kedtau(dffts%nnr,nspin), kedtaur(nrxx,nspin)
  ! output: total local potential on the smooth grid
  !         vrs=vltot+vr
  ! input: the total local pseudopotential
  ! input: the scf(H+xc) part of the local potential
  logical :: doublegrid
  ! input: true if a doublegrid is used
  !
  CALL sum_vrs( nrxx, nspin, vltot, vr, vrs )
  !
  CALL interpolate_vrs( nrxx, nspin, doublegrid, kedtau, kedtaur, vrs )
  ! 
  return

end subroutine set_vrs
!
!--------------------------------------------------------------------
subroutine sum_vrs ( nrxx, nspin, vltot, vr, vrs )
  !--------------------------------------------------------------------
  ! accumulates local potential contributions in to vrs 
  !
  USE kinds
  USE io_global,     ONLY : ionode
  USE io_global,     ONLY : stdout, ionode
  USE epcdft,    ONLY : do_epcdft, reset_field, epcdft_field, epcdft_surface_shift,&
                        epcdft_surface
  !
  implicit none

  logical :: first=.true.
  integer :: nspin, nrxx
  ! input: number of spin components: 1 if lda, 2 if lsd, 4 if noncolinear
  ! input: the fft grid dimension
  real(DP) :: vrs(nrxx,nspin), vltot(nrxx), vr(nrxx,nspin)
  ! output: total local potential on the smooth grid
  !         vrs=vltot+vr
  ! input: the total local pseudopotential
  ! input: the scf(H+xc) part of the local potential
  !
  REAL(DP) :: x0(3) ! center of charge of system
  REAL(DP) :: qq ! total charge
  REAL(DP) :: dipole(3)!, quadrupole(3) ! total dips
  integer:: is
  REAL(DP) :: tmp
  SAVE first
  !write(*,*) "sum_vrs",sum(vr(:,1))
  x0 = 0.D0
  qq = 0.D0
  dipole = 0.D0
  !
  IF (first .or. .not. allocated(epcdft_field)) THEN
    !
    allocate(epcdft_field(nrxx,nspin))
    reset_field=.true.
    !
  END IF
  !
  IF(first)THEN
    !
    ! set image interaction energy
    !
    IF(epcdft_surface)THEN
      !
      !
      IF (ionode) THEN
        !
        WRITE( stdout,'(5x,"":)')
        WRITE( stdout,'(5x,"Adding image charge field due to neutral metal slab at origin in XY plane.":)')
        WRITE( stdout,'(5x,"Including monopole and dipole terms.":)')
        WRITE( stdout,'(5x,"Quadrupole not implemented yet.":)')
        WRITE( stdout,'(5x,"System should not overlap with cell edges.":)')
        WRITE( stdout,'(5x,"":)')
        !WRITE( stdout,'(8x,"Amplitude [Ry a.u.] : ", es11.4)') eamp 
        !WRITE( stdout,'(8x,"Postion on atom # : ", I11.1)') edir
        !WRITE( stdout,'(8x,"Well radius [bohr] : ", es11.4)') emaxpos * alat
        !
      ENDIF
      !
      epcdft_field = 0.D0
      tmp = 0.D0
      !
      CALL epcdft_surface_energy(epcdft_field, tmp)
      !
      epcdft_surface_shift = tmp
      !
      IF(ionode)THEN
        WRITE(*,*)
        WRITE(*,'(5x,"First surface image interaction energy:",e10.3,"[Ry]")') epcdft_surface_shift
        WRITE(*,*)
      ENDIF
      !
    ENDIF
      !
  ENDIF
  !
  first=.false.
  !
  IF ( do_epcdft .and. reset_field) THEN
    !
    epcdft_field=0.0D0
    !
    ! epcdft surface also added in call below only if epcdft_surface=T
    !
    CALL add_epcdft_efield(epcdft_field,.TRUE.)
    IF(epcdft_surface) CALL calc_epcdft_surface_field( epcdft_field, x0, qq, dipole )
    reset_field=.false.
    !
  ENDIF
  !
  vr=vr+epcdft_field
  !write(*,*) "sum_vrs",sum(vr(:,1))
  ! confirmed the potential was changing here

  do is = 1, nspin
     !
     ! define the total local potential (external + scf) for each spin ...
     !
     if (is > 1 .and. nspin == 4) then
        !
        ! noncolinear case: only the first component contains vltot
        !
        vrs (:, is) = vr (:, is)
     else
        vrs (:, is) = vltot (:) + vr (:, is)
     end if
     !
  enddo
  return

end subroutine sum_vrs
!
!--------------------------------------------------------------------
subroutine interpolate_vrs ( nrxx, nspin, doublegrid, kedtau, kedtaur, vrs )
  !--------------------------------------------------------------------
  ! set the total local potential vrs on the smooth mesh to be used in 
  ! h_psi, adding the (spin dependent) scf (H+xc) part and the sum of 
  ! all the local pseudopotential contributions.
  !
  USE kinds
  USE funct, only : dft_is_meta
  USE fft_base, only : dffts 
  implicit none

  integer :: nspin, nrxx
  ! input: number of spin components: 1 if lda, 2 if lsd, 4 if noncolinear
  ! input: the fft grid dimension
  real(DP) :: vrs (nrxx, nspin), &
              kedtau(dffts%nnr,nspin), kedtaur(nrxx,nspin)
  ! output: total local potential interpolated on the smooth grid
  ! input: the scf(H+xc) part of the local potential
  logical :: doublegrid
  ! input: true if a doublegrid is used

  integer:: is

  do is = 1, nspin
     !
     ! ... and interpolate it on the smooth mesh if necessary
     !
     if (doublegrid) call interpolate (vrs (1, is), vrs (1, is), - 1)
     if (dft_is_meta()) call interpolate(kedtaur(1,is),kedtau(1,is),-1)
  enddo
  return

end subroutine interpolate_vrs
