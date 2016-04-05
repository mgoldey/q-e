!
!    EPCDFT
!
MODULE epcdft
!
! ... The quantities needed in calculations with epcdft
!
USE kinds, ONLY : DP
!
SAVE
!
  PUBLIC :: do_epcdft, conv_epcdft, epcdft_fields, nconstr_epcdft,&
    epcdft_delta_fld, epcdft_tol, epcdft_shift, epcdft_type, &
    epcdft_locs, epcdft_guess, epcdft_target, epcdft_forces, &
    update_intrvl
  !
  LOGICAL :: &
       do_epcdft,       &  ! if .TRUE. a finite electric field is added to the
                           ! local potential
       conv_epcdft,     &  ! localization condition flag
       reset_field=.false. ! reset field (so we don't do this every time)
  !
  INTEGER :: epcdft_fields = 4   ! max number of fields that is allowed to
                                 ! define a constraint
  !
  INTEGER :: epcdft_update_intrvl = 40  ! update the potential every this many steps and 
                                        ! after electrons reach scf
  !
  INTEGER  :: nconstr_epcdft    = 0
  REAL(DP) :: epcdft_delta_fld =1.E-1_DP
  REAL(DP) :: epcdft_tol = 1.E-4_DP
  !
  REAL(DP) :: epcdft_shift    ! energy shift from all fields
  !
  CHARACTER(len=20), ALLOCATABLE :: epcdft_type(:)    ! type
  INTEGER,           ALLOCATABLE :: epcdft_locs(:,:)  ! atoms start end (start end) ()-if delta
  REAL(DP),          ALLOCATABLE :: epcdft_guess(:)   !  guess values of constraints
  REAL(DP),          ALLOCATABLE :: epcdft_target(:)  ! target values of constraints
  REAL(DP),          ALLOCATABLE :: epcdft_field(:,:) ! constraint field
  !
  REAL(DP), ALLOCATABLE :: epcdft_forces(:,:)
  !
CONTAINS
  !
  !-----------------------------------------------------------------------------
  SUBROUTINE allocate_input_epcdft()
  !-----------------------------------------------------------------------------
  !
  IF ( allocated( epcdft_type ) )   DEALLOCATE( epcdft_type )
  IF ( allocated( epcdft_locs ) )   DEALLOCATE( epcdft_locs )
  IF ( allocated( epcdft_guess ) )  DEALLOCATE( epcdft_guess )
  IF ( allocated( epcdft_target ) ) DEALLOCATE( epcdft_target )
  !
  ALLOCATE( epcdft_type(   nconstr_epcdft ) )
  ALLOCATE( epcdft_guess(  nconstr_epcdft ) )
  ALLOCATE( epcdft_target( nconstr_epcdft ) )
  !
  ALLOCATE( epcdft_locs( epcdft_fields, nconstr_epcdft ) )
  !
  epcdft_type   = ' '
  epcdft_locs   = 0
  epcdft_target = 0.0_DP
  epcdft_guess  = 0.0_DP
  !
  RETURN
  !
  END SUBROUTINE allocate_input_epcdft
END MODULE epcdft

