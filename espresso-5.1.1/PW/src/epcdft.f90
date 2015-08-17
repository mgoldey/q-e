!
! ==============================================================
! Developer info: 
! Matthew Goldey & Nicholas Brawand
! The University of Chicago, Institute for Molecular Engineering
! matthew.goldey@gmail.com nicholasbrawand@gmail.com
! ==============================================================
!
!----------------------------------------------------------------------------
SUBROUTINE epcdft_controller()
  !----------------------------------------------------------------------------
  !
  !   CDFT
  !
  !   This routine calculates the correction to the total energy
  !   due to the constraining potential from cdft and prints it. 
  !   It also calculates the total number of 
  !   electrons within the applied well. If this number is not
  !   equal to that of epcdft_charge within epcdft_thr, 
  !   the amplitude of the well is changed and the scf loop is restarted. 
  !
  !   do_epcdft - flag to do cdft
  !
  !   acceptor_start - first atom in voronoi cell or in acceptor (if hirshfeld)
  !
  !   acceptor_end   - last atom in voronoi cell or in acceptor (if hirshfeld)
  !                    if zero user defined well is used (not working right now)
  !
  !   epcdft_charge - number of electrons in voronoi cell or on acceptor
  !
  !   epcdft_amp - amplitude of voronoi cell or lagrange multiplier for hirshfeld
  !
  !   epcdft_width - width of the user defined well (not working right now
  !
  !   epcdft_shift -  energy correction to etot due to constraining potential
  !
  !   epcdft_thr - threshold on number of electrons to match epcdft_charge
  !
  !   hirshfeld - if .true. hirshfeld is used rather than voronoi cells
  !
  !
  !
  USE io_global,     ONLY : stdout, ionode
  USE kinds,         ONLY : DP
  USE io_files,      ONLY : tmp_dir
  USE plugin_flags
  USE kinds,         ONLY : DP
  USE constants,     ONLY : fpi, eps8, e2, au_debye
  USE ions_base,     ONLY : nat, ityp, zv
  USE cell_base,     ONLY : alat, at, omega, bg, saw
  USE extfield,      ONLY : tefield, dipfield, edir, eamp, emaxpos, &
                            eopreg, forcefield, etotefield
  USE epcdft,        ONLY : do_epcdft, donor_start, &
                            acceptor_start, acceptor_end, &
                            donor_end, epcdft_charge, &
                            epcdft_amp, epcdft_width, epcdft_shift, &
                            epcdft_thr, hirshfeld, conv_epcdft, epcdft_delta_fld
  USE force_mod,     ONLY : lforce
  USE io_global,     ONLY : stdout,ionode, ionode_id
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_images,     ONLY : intra_image_comm
  USE mp_bands,      ONLY : me_bgrp, intra_bgrp_comm
  USE fft_base,      ONLY : dfftp, grid_gather
  USE mp,            ONLY : mp_bcast, mp_sum
  USE control_flags, ONLY : iverbosity, conv_elec
  USE scf,           ONLY : rho
  USE io_files,  ONLY : tmp_dir, prefix
  !
  IMPLICIT NONE
  !
  ! local variables
  !
  INTEGER  :: i, is, iatom
  REAL(DP) :: tmp
  REAL(DP) :: dv
  REAL(DP) :: einwell                              ! number of electrons in well
  REAL(DP) :: enumerr                              ! epcdft_charge - einwell  (e number error)
  LOGICAL  :: elocflag                             ! true if charge localization condition is satisfied
  LOGICAL  :: zero                                 ! used to run cdft with zero constraining potential
  REAL(DP) :: safe_epcdft_amp                      ! used for zero constraining potential run 
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: vpotens ! ef is added to this potential serial
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: vpotenp ! ef is added to this potential parll
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rhosup    ! rho serial
  REAL(DP), DIMENSION(:), ALLOCATABLE :: rhosdown  ! rho serial
  REAL(DP) :: next_epcdft_amp                      ! guess of amp for next iteration
  CHARACTER(LEN=256) :: filename                   ! cube file for final cdft potential
  !
  ! SAVED VARS
  !
  REAL(DP) :: last_epcdft_amp = 0.D0 ! used to determine next guess for potential
  REAL(DP) :: last_einwell = 0.D0    ! used to determine next guess for potential
  LOGICAL :: first = .TRUE.          ! used to determine next guess for potential
  SAVE last_epcdft_amp
  SAVE last_einwell 
  SAVE first
  !
  ! dont do anything unless the calculation is converged
  !
  IF(.NOT.conv_elec)RETURN
  !
  ! calc is converged lets compute and print the correction
  !
  ! first setup vars
  ALLOCATE(vpotenp(dfftp%nnr, nspin))
  ALLOCATE(vpotens(dfftp%nr1x * dfftp%nr2x * dfftp%nr3x,nspin ))
  ALLOCATE(rhosup( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x ))
  ALLOCATE(rhosdown( dfftp%nr1x * dfftp%nr2x * dfftp%nr3x ))
  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  vpotenp  = 0.D0
  vpotens  = 0.D0
  rhosup   = 0.D0
  rhosdown = 0.D0
  einwell  = 0.D0
  tmp      = 0.D0
  elocflag = .TRUE.
  etotefield = 0.D0
  epcdft_shift = 0.D0
  zero =.false.
  next_epcdft_amp = 0.D0
  !
  ! gather the grids for serial calculation
  !
  IF (epcdft_amp .eq. 0.D0) THEN
    zero=.true.
    epcdft_amp=1.D0 
  ENDIF
  !
  safe_epcdft_amp = epcdft_amp
  !
  ! this call only calulates vpoten
  !
  CALL add_epcdft_efield(vpotenp(:,1), epcdft_shift, rho%of_r, .true. )
  CALL add_epcdft_efield(vpotenp(:,2), epcdft_shift, rho%of_r, .true. )
  !
  IF (zero) epcdft_amp=0.D0
  !
  ! gather the grids for serial calculation
  !
#ifdef __MPI
    CALL grid_gather ( vpotenp(:,1), vpotens(:,1) )
    CALL grid_gather ( rho%of_r(:,1), rhosup )
    IF(nspin > 1)THEN
      CALL grid_gather ( rho%of_r(:,2), rhosdown )
      CALL grid_gather ( vpotenp(:,2), vpotens(:,2) )
    ENDIF
#else
    vpotens(:,1)=vpotenp(:,1)
    rhosup(:) = rho%of_r(:,1)
    IF(nspin > 1)THEN
      rhosdown(:) = rho%of_r(:,2)
      vpotens(:,2)=vpotenp(:,2)
    ENDIF
#endif
  !
  ! begin calculation of the correction 
  !
  IF(ionode) THEN
    !
    ! combine up and down parts of rho
    rhosup(:) = rhosup(:) + rhosdown(:) 
    ! vpotens(:,1)=vpotens(:,1)+vpotens(:,2)
    !
    DO i=1, dfftp%nr1x*dfftp%nr2x*dfftp%nr3x
      !
      ! calculate energy correction
      !
      epcdft_shift = epcdft_shift + (vpotens(i,1)) * rhosup(i) * dv
      !
      ! count number of electrons in well 
      !
      IF(hirshfeld) THEN
        ! need number of electrons on 
        ! acceptor so negetive of vpotens is used
        einwell = einwell - ( vpotens(i,1) / safe_epcdft_amp ) * rhosup(i) * dv
        !
      ELSE
        !
        IF(vpotens(i,1) .lt. 0.D0)THEN
          einwell = einwell + rhosup(i) * dv
        ELSE IF (vpotens(i,1).gt. 0.D0) THEN
          einwell = einwell - rhosup(i) * dv
        ENDIF
        !
      ENDIF
      !
    ENDDO
    !
    ! count nuc charge
    !
    DO iatom=1, nat
      !
      ! Voronoi and Hirshfeld below
      IF (acceptor_end.ne.0) THEN
        !
        IF ((iatom.ge.acceptor_start) .AND. (iatom.le.acceptor_end) ) THEN
          einwell = einwell - zv(ityp(iatom))
        ELSE IF ((iatom.ge.donor_start) .AND. (iatom.le.donor_end) ) THEN
          einwell = einwell + zv(ityp(iatom))
        ENDIF
        !
      ELSE ! just a well around one atom
        !
        IF (iatom.eq.acceptor_start) THEN
          einwell = einwell - zv(ityp(iatom))
        ELSE
          einwell = einwell + zv(ityp(iatom))
        ENDIF
        !
      ENDIF
      !
    ENDDO ! atoms
    !
    ! the correction is - of the energy
    epcdft_shift = -1.D0 * epcdft_shift
    !
    IF (.not.zero) WRITE(*,*)"    E field correction              :  ",epcdft_shift,"Ry"
    WRITE(*,*)"    Difference of charges           : ",einwell," electrons"
    !
    ! is there a localization condition?
    !
    IF(epcdft_amp .ne. 0.D0)THEN
      !
      ! is the localization condition satisfied?
      !
      enumerr = epcdft_charge - einwell 
      !
      ! conv_epcdft = false will restart scf
      !
      IF( ABS(enumerr) .GE. epcdft_thr .and. .not. zero) THEN
        !
        conv_epcdft = .FALSE.
        !
        WRITE(*,*) "    Surplus(+)/deficit(-) charge    :  ", -enumerr,    "electrons"
        WRITE(*,*) "    epcdft_thr                      :  ", epcdft_thr, "electrons"
        !
      ELSE
        conv_epcdft =.true.
        ! WRITE OUT potential
        !call write_wfc_cube_r ( 84332, 'v',  v )
      ENDIF
      !
    ENDIF ! localization condition
    !
  ENDIF ! io node
  !
  IF (zero) THEN
      ! IF(ionode) write(*,*) "All except for number of electrons is meaningless - EXITING NOW"
    conv_epcdft =.true.
  ENDIF
  !
  CALL mp_bcast( conv_epcdft, ionode_id, intra_image_comm )
  !
  IF(conv_epcdft)THEN  ! cdft done write external pot to cube
    !
    filename =  TRIM( tmp_dir ) // TRIM( prefix ) // 'v_cdft'
    CALL write_cube_r ( 9519395, filename, vpotenp(:,1) )
    !
  ENDIF
  !
  DEALLOCATE(vpotenp)
  DEALLOCATE(vpotens)
  DEALLOCATE(rhosup)
  DEALLOCATE(rhosdown)
  !
  ! if the charge is not localized
  !
  IF(.NOT.conv_epcdft .AND. epcdft_amp .NE. 0)THEN
    !
    ! update applied field and restart scf
    !
    IF(ionode)THEN
      !
      WRITE(*,*)""
      WRITE(*,*)"    -----------------------------------------------"
      WRITE(*,*)""
      WRITE(*,*)"                     SCF converged but..."
      WRITE(*,*)""
      WRITE(*,*)"    electron localization condition NOT satisfied"
      WRITE(*,*)"    restarting scf with different applied potential"
      WRITE(*,*)""
      !
      ! find next guess for epcdft_amp
      !
      IF(first) THEN
         !
         first = .FALSE.
         !
         next_epcdft_amp = epcdft_amp + SIGN(0.001D0, enumerr) 
         !
      ELSE
         !
         CALL secant_method(next_epcdft_amp, epcdft_amp,   last_epcdft_amp, &
                            einwell,         last_einwell, epcdft_charge)
         !
         ! abs of the change in amp must be <= |delta_fld|
         !
         IF ( ABS(next_epcdft_amp - epcdft_amp) .gt. ABS(epcdft_delta_fld) ) THEN
           !
           next_epcdft_amp = epcdft_amp + SIGN(1.D0, next_epcdft_amp - epcdft_amp) * epcdft_delta_fld
           !
         ENDIF
         !
      ENDIF
      !
      ! save this iteration's einwell and amp
      ! for the next iteration
      last_einwell    = einwell
      last_epcdft_amp = epcdft_amp
      !
      ! The old iteration has passed away;
      ! behold, the new iteration has come 
      epcdft_amp = next_epcdft_amp
      !
    ENDIF
    !
    CALL mp_bcast( epcdft_amp, ionode_id, intra_image_comm ) ! what is best bcast this or enumerr?
    !
    IF(ionode) WRITE(*,*)"    New field Amp      : ",epcdft_amp," Ry"
    !
    ! epcdft_shift = 0.D0 ! this var is added to etot before 
    !                   ! this routine is called during next scf loop
  ENDIF
  !
END SUBROUTINE epcdft_controller
!----------------------------------------------------------------------------
!
!
SUBROUTINE secant_method(vnext, v, vold, e, eold, egoal)
  !----------------------------------------------------------------------------
  !
  ! This routine uses secant method to determine the next 
  ! estimation of the root of an equation.
  !
  ! Numerical Mathematics and Computing Sixth Edition
  !  Ward Cheney & David Kincaid Pg 112, eqn. 3
  !
  !     vnext = v - \left( \frac{v-vold}{e-eold} \right) * (e-egoal)
  !
  ! In this case the function we want to minmize is:
  !          e*(v) = e(v) - egoal
  ! this is why the formual has egoal in it
  !
  USE kinds,            ONLY : DP
  !
  IMPLICIT NONE
  !
  REAL(DP), INTENT(INOUT) :: vnext
  REAL(DP), INTENT(IN) :: v, vold, e, eold, egoal
  !
  vnext = v - ( (v-vold) / (e-eold)  ) * (e-egoal)
  !
END SUBROUTINE secant_method
!----------------------------------------------------------------------------
!
! -------------------------------------------------------------------
SUBROUTINE write_cube_r ( iu, fname, wfc_distr )
  ! -----------------------------------------------------------------
  !
  USE pwcom,                 ONLY : npw,npwx
  USE kinds,                 ONLY : DP
  USE cell_base,             ONLY : celldm, at, bg
  USE ions_base,             ONLY : nat, tau, atm, ityp
  USE fft_base,              ONLY : dfftp, grid_gather
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
  CALL grid_gather(wfc_distr,wfc_gat)
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
     OPEN(UNIT=iu,FILE=TRIM(ADJUSTL(fname))//".cub")
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
        WRITE(iu,'(I5,5F12.6)') at_num, at_chrg, inpos
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
        WRITE(ounit,'(6E13.5)') (func(i1,i2,i3),i3=1,nr3)
     ENDDO
  ENDDO
  !
END SUBROUTINE
