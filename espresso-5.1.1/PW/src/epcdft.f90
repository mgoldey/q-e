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
  !   equal to that of the target charge within epcdft_tol, 
  !   the amplitude of the well is changed and the scf loop is restarted. 
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
  USE epcdft,        ONLY : do_epcdft, conv_epcdft,epcdft_shift
  USE force_mod,     ONLY : lforce
  USE io_global,     ONLY : stdout,ionode, ionode_id
  USE control_flags, ONLY : mixing_beta
  USE lsda_mod,      ONLY : nspin
  USE mp_images,     ONLY : intra_image_comm
  USE mp_bands,      ONLY : me_bgrp, intra_bgrp_comm
  USE mp_pools,      ONLY : inter_pool_comm
  USE mp_world,      ONLY : world_comm

  USE fft_base,      ONLY : dfftp, grid_gather
  USE mp,            ONLY : mp_bcast, mp_sum
  USE control_flags, ONLY : iverbosity, conv_elec
  USE scf,           ONLY : rho
  USE io_files,  ONLY : tmp_dir, prefix

  ! here's the meat
  USE input_parameters, ONLY : epcdft_locs,epcdft_guess,nconstr_epcdft, &
                    & epcdft_tol, epcdft_type,epcdft_target, epcdft_delta_fld

  !
  IMPLICIT NONE
  !
  ! local variables
  !
  INTEGER  :: i, is, iatom,iconstraint
  INTEGER  :: ictr=0
  REAL(DP) :: tmp
  REAL(DP) :: dv
  REAL(DP) :: acharge, dcharge                     ! acceptor/donor charge 
  REAL(DP) :: achargep, dchargep                   ! acceptor/donor charge parallel
  REAL(DP) :: einwell                              ! number of electrons in well
  REAL(DP) :: einwellp                             ! number of electrons in well
  REAL(DP) :: epcdft_shiftp                        ! size of shift
  REAL(DP) :: enumerr                              ! epcdft_charge - einwell  (e number error)
  LOGICAL  :: elocflag                             ! true if charge localization condition is satisfied
  LOGICAL  :: zero                                 ! used to run cdft with zero constraining potential
  REAL(DP) :: safe_epcdft_amp                      ! used for zero constraining potential run 
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: vpotenp ! ef is added to this potential
  REAL(DP) :: next_epcdft_amp,epcdft_amp           ! guess of amp for next iteration
  CHARACTER(LEN=256) :: filename                   ! cube file for final cdft potential
  !
  ! SAVED VARS
  !
  REAL(DP), DIMENSION(:), allocatable,save :: last_epcdft_amp(:) ! used to determine next guess for potential
  REAL(DP), DIMENSION(:), allocatable,save :: last_einwell(:) ! used to determine next guess for potential
  LOGICAL :: first = .TRUE.          ! used to determine next guess for potential
  SAVE ictr
  SAVE first
  if (ictr.eq.0) THEN
    allocate(last_epcdft_amp(nconstr_epcdft))
    allocate(last_einwell(nconstr_epcdft))
    last_epcdft_amp=0.D0
    last_einwell=0.D0
  ENDIF
  !
  ! Don't try to update too often or else the number of electrons will sometimes go the wrong way and mess up the solver
  ictr=ictr+1
  if (conv_elec) ictr=8
  if (mod(ictr,10) .ne. 8) RETURN
  !
  ! calc is converged lets compute and print the correction
  !
  ! first setup vars
  ALLOCATE(vpotenp(dfftp%nnr, nspin))

  dv = omega / DBLE( dfftp%nr1 * dfftp%nr2 * dfftp%nr3 )
  tmp      = 0.D0
  elocflag = .TRUE.
  etotefield = 0.D0
  epcdft_shift = 0.D0
  epcdft_shiftp = 0.D0
  zero =.false.
  next_epcdft_amp = 0.D0
  acharge = 0.D0
  dcharge = 0.D0
  achargep = 0.D0
  dchargep = 0.D0
  !
  ! gather the grids for serial calculation
  !
  DO iconstraint=1,nconstr_epcdft
    zero=.false.
    einwell  = 0.D0
    einwellp  = 0.D0
    epcdft_shiftp=0.D0
    epcdft_amp=epcdft_guess(iconstraint)
    IF (epcdft_amp .eq. 0.D0) THEN
      zero=.true.
      epcdft_amp=1.D0 
    ENDIF
    safe_epcdft_amp = epcdft_amp
    !
    ! this call only calulates vpoten - hardcoded to hirshfeld right now
    !
    vpotenp=0.D0
    
    if (.true.) CALL calc_hirshfeld_v(vpotenp,iconstraint)
    !write(*,*) iconstraint,sum(vpotenp(:,1)*rho%of_r(:,1))*dv,sum(vpotenp(:,2)*rho%of_r(:,2))*dv
    !
    IF (zero) epcdft_amp=0.D0
    !
    ! begin calculation of the correction 
    !
    DO i=1, dfftp%nnr
      !
      ! calculate energy correction
      !
      IF (nspin.eq.1) THEN
        epcdft_shiftp = epcdft_shiftp + epcdft_amp * vpotenp(i,1) * rho%of_r(i,1) * dv
      ELSE
        epcdft_shiftp = epcdft_shiftp + epcdft_amp * vpotenp(i,1) * rho%of_r(i,1) * dv &
                                      + epcdft_amp * vpotenp(i,2) * rho%of_r(i,2) * dv
      ENDIF
      !
      ! count number of electrons in well - this must depend on the well
      !
      SELECT CASE( epcdft_type(iconstraint) )
      CASE('charge','delta_charge')
        einwellp = einwellp - vpotenp(i,1) * rho%of_r(i,1) * dv - vpotenp(i,2) * rho%of_r(i,2) * dv
        !
        IF(vpotenp(i,1)<0.D0) THEN
          achargep = achargep - ABS( vpotenp(i,1)  * rho%of_r(i,1) +vpotenp(i,2) * rho%of_r(i,2)) * dv
        ELSE IF(vpotenp(i,1)>0.D0) THEN
          dchargep = dchargep - ABS( vpotenp(i,1)  * rho%of_r(i,1) +vpotenp(i,2) * rho%of_r(i,2)) * dv
        ENDIF
      CASE('spin','delta_spin')
        einwellp = einwellp - vpotenp(i,1) * rho%of_r(i,1) * dv + vpotenp(i,2) * rho%of_r(i,2) * dv
        !
        IF(vpotenp(i,1)<0.D0) THEN
          achargep = achargep - ABS( vpotenp(i,1)  * rho%of_r(i,1) + vpotenp(i,2) * rho%of_r(i,2)) * dv
        ELSE IF(vpotenp(i,1)>0.D0) THEN
          dchargep = dchargep - ABS( vpotenp(i,1)  * rho%of_r(i,1) + vpotenp(i,2) * rho%of_r(i,2)) * dv
        ENDIF
      END SELECT
    END DO
#ifdef __MPI
  CALL MP_SUM(epcdft_shiftp,intra_image_comm) 
  CALL MP_SUM(einwellp,intra_image_comm) 
  epcdft_shift=epcdft_shiftp
  einwell=einwellp
!
  CALL MP_SUM(achargep,intra_image_comm) 
  acharge=achargep
  CALL MP_SUM(dchargep,intra_image_comm)
  dcharge=dchargep
#else
  epcdft_shift=epcdft_shiftp
  einwell=einwellp
  !
  acharge=achargep
  dcharge=dchargep
#endif
    !
    ! count nuc charge
    !
    SELECT CASE( epcdft_type(iconstraint) )
      CASE('charge','delta_charge')
        DO iatom=1, nat
          !
          IF ((iatom.ge.epcdft_locs(1,iconstraint)) .AND. (iatom.le.epcdft_locs(2,iconstraint)) ) THEN
            !write(*,*) iatom, epcdft_locs(1,iconstraint), epcdft_locs(2,iconstraint)
            einwell = einwell - zv(ityp(iatom))
            acharge = acharge + zv(ityp(iatom)) 
          ELSE IF ((iatom.ge.epcdft_locs(3,iconstraint)) .AND. (iatom.le.epcdft_locs(4,iconstraint)) ) THEN
            !write(*,*) iatom, epcdft_locs(3,iconstraint), epcdft_locs(4,iconstraint)
            einwell = einwell + zv(ityp(iatom))
            dcharge = dcharge + zv(ityp(iatom))
          ENDIF
          !
        ENDDO ! atoms
      END SELECT
    !
    !
    ! the correction is - of the energy
    epcdft_shift = -1.D0 * epcdft_shift
    !   


    !
    ! is there a localization condition?
    !
    IF(epcdft_amp .ne. 0.D0) THEN
      !
      ! is the localization condition satisfied?
      !
      enumerr = epcdft_target(iconstraint) - einwell 
      !
      ! conv_epcdft = false will restart scf
      !
      IF( ABS(enumerr) .GE. epcdft_tol .and. .not. zero .and. ionode) THEN
        !
        conv_epcdft = .FALSE.
        !
!         WRITE(*,*)""
!         CALL pprint("epcdft_tol",epcdft_tol,'e','e')
!         CALL pprint("Surplus(+)/deficit(-) charge",-enumerr,'e','e')
        !
      ELSE
        conv_epcdft =.true.
      ENDIF
      !
    ENDIF ! localization condition
    !
    IF (zero) THEN
      conv_epcdft =.true.
    ENDIF
    !
    CALL mp_bcast( conv_epcdft, ionode_id, intra_image_comm )
    !
    IF(conv_epcdft .and. conv_elec )THEN  ! cdft done write external pot to cube
      !
      filename =  TRIM( tmp_dir ) // TRIM( prefix ) // 'rho_up'
      CALL write_cube_r ( 9519395, filename, rho%of_r(:,1) )
      filename =  TRIM( tmp_dir ) // TRIM( prefix ) // 'rho_down'
      CALL write_cube_r ( 9519395, filename, rho%of_r(:,2) )
      filename =  TRIM( tmp_dir ) // TRIM( prefix ) // 'v_cdft_up'
      CALL write_cube_r ( 9519395, filename, vpotenp(:,1) )
      filename =  TRIM( tmp_dir ) // TRIM( prefix ) // 'v_cdft_down'
      CALL write_cube_r ( 9519395, filename, vpotenp(:,2) )
      !
    ENDIF
    !
    ! if the charge is not localized
    !
    IF(.NOT.conv_epcdft .AND. epcdft_amp .NE. 0)THEN
      !
      ! update applied field and restart scf
      !
      IF(ionode) THEN
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
           CALL secant_method(next_epcdft_amp, epcdft_amp,   last_epcdft_amp(iconstraint), &
                              einwell,         last_einwell(iconstraint), epcdft_target(iconstraint))
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
        last_einwell(iconstraint)    = einwell
        last_epcdft_amp(iconstraint) = epcdft_amp
        !
        ! The old iteration has passed away;
        ! behold, the new iteration has come 
        epcdft_guess(iconstraint) = next_epcdft_amp
        !
      ENDIF
      !
      CALL mp_bcast( epcdft_guess(iconstraint), ionode_id, intra_image_comm ) ! what is best bcast this or enumerr?
      !
    ENDIF
    IF(ionode) THEN
      write(*,'(5x,a6,a2,a6,a9,a9,a9)') "D",'  ','A','Diff','Old','New'
      write(*,'(5x,f6.3,a2,f6.3,f9.3,f9.3,f9.3)') dcharge,'  ',acharge,einwell,last_epcdft_amp,next_epcdft_amp
    ENDIF
    if (conv_epcdft .and. conv_elec) &
      write(*,*) "Converged potential ",iconstraint," with strength  ", epcdft_guess(iconstraint)
  END DO
  DEALLOCATE(vpotenp)

  !
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
!
SUBROUTINE pprint( strg, num, ustrg, a)
  !
  ! write strg, num and ustrg
  ! a is used to determine the format of num
  !
  USE kinds, ONLY:DP
  !
  IMPLICIT NONE
  REAL(DP),INTENT(IN)::num
  CHARACTER(LEN=*),INTENT(IN)::strg, ustrg
  CHARACTER(LEN=1),INTENT(IN)::a
  CHARACTER(LEN=32)::trimstrg
  CHARACTER(LEN=6)::utrimstrg
  !
  trimstrg = strg
  utrimstrg = ustrg
  !
  SELECT CASE(a)
    CASE('f','F')
      WRITE(*,1)trimstrg,num,utrimstrg
    CASE('e','E')
      WRITE(*,2)trimstrg,num,utrimstrg
  END SELECT
  !
  1 FORMAT(5x,A32,':',1x,F23.16,1x,A6)
  2 FORMAT(5x,A32,':',1x,e23.16,1x,A6)
  !
END SUBROUTINE
