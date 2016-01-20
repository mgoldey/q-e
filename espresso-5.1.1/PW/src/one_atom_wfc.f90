!
! Copyright (C) 2001-2007 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
SUBROUTINE one_atom_wfc (ik, wfcatom,iatom,nfuncs)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the superposition of atomic wavefunctions
  ! for k-point "ik" - output in "wfcatom"
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi, fpi, pi
  USE cell_base,  ONLY : omega, tpiba, bg
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  !USE basis,      ONLY : natomwfc
  USE gvect,      ONLY : mill, g ! eigts1, eigts2, eigts3
  USE klist,      ONLY : xk
  USE fft_base,      ONLY : dfftp !, grid_scatter, grid_gather
  USE wvfct,      ONLY : npwx, npw, nbnd, igk
  USE us,         ONLY : tab_at, dq
  USE uspp_param, ONLY : upf
  USE io_global,     ONLY : stdout,ionode, ionode_id

  !
  implicit none
  !
  integer, intent(in) :: ik
  integer, intent(in) :: iatom
  integer, intent(in) :: nfuncs
  complex(DP), intent(out) :: wfcatom (npwx, nfuncs)
  !
  integer :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3,  ipol, n1, n2, n3, natomwfc
  real(DP), allocatable :: qg(:), ylm (:,:), chiq (:,:), gk (:,:)
  complex(DP), allocatable :: sk (:), aux(:)
  complex(DP) :: kphase, lphase

  complex(DP) ::      eigts1 ( -dfftp%nr1:dfftp%nr1), &
                      eigts2 ( -dfftp%nr2:dfftp%nr2), &
                      eigts3 ( -dfftp%nr3:dfftp%nr3)

  real(DP) :: arg, px, ux, vx, wx
  real(DP) :: bgtau (3)

  call start_clock ('atomic_wfc')

  ! calculate max angular momentum required in wavefunctions
  nt=ityp(iatom)
  ! natomwfc=SUM(upf(nt)%oc(:))/2
  ! natomwfc=upf(nt)%nwfc
  lmax_wfc = MAXVAL (upf(nt)%lchi(1:upf(nt)%nwfc) )
  !
  na=iatom
  !
  allocate ( ylm (npw,(lmax_wfc+1)**2), chiq(npw,upf(nt)%nwfc), &
             sk(npw), gk(3,npw), qg(npw) )
  !
  do ig = 1, npw
     gk (1,ig) = xk(1, ik) + g(1, igk(ig) )
     gk (2,ig) = xk(2, ik) + g(2, igk(ig) )
     gk (3,ig) = xk(3, ik) + g(3, igk(ig) )
     qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  !  ylm = spherical harmonics
  !
  call ylmr2 ((lmax_wfc+1)**2, npw, gk, qg, ylm)
  !
  ! set now q=|k+G| in atomic units
  !
  do ig = 1, npw
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  
!   if (ionode) write(*,*) "Atom ",iatom
!   if (ionode) write(*,*) "itype ",nt
!   if (ionode) write(*,*) "nwfc ",nfuncs
!    if (ionode) write(*,*) "npw ",npw
!   if (ionode) write(*,*) "x ",tau(1,na)
!   if (ionode) write(*,*) "y ",tau(2,na)
!   if (ionode) write(*,*) "z ",tau(3,na)

  do nb = 1, upf(nt)%nwfc
    if ( upf(nt)%oc (nb) >= 0.d0) then
        do ig = 1, npw
            px = qg (ig) / dq - int (qg (ig) / dq)
            ux = 1.d0 - px
            vx = 2.d0 - px
            wx = 3.d0 - px
            i0 = INT( qg (ig) / dq ) + 1
            i1 = i0 + 1
            i2 = i0 + 2
            i3 = i0 + 3
            chiq (ig, nb) = &
                   tab_at (i0, nb, nt) * ux * vx * wx / 6.d0 + &
                   tab_at (i1, nb, nt) * px * vx * wx / 2.d0 - &
                   tab_at (i2, nb, nt) * px * ux * wx / 2.d0 + &
                   tab_at (i3, nb, nt) * px * ux * vx / 6.d0
        enddo
    endif
  enddo
   !if (ionode) write(*,*) "made it here" ! works

  deallocate (qg, gk)
  !if (ionode) write(*,*) "made it here2"
  allocate ( aux(npw) )
  !
  wfcatom(:,:) = (0.0_dp, 0.0_dp)
  !if (ionode) write(*,*) "made it here" !fail
  !
  arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
  kphase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
  !
  !     sk is the structure factor
  !
  ! from struct_fact
  do ipol = 1, 3
    bgtau (ipol) = bg (1, ipol) * tau (1, na) + &
                   bg (2, ipol) * tau (2, na) + &
                   bg (3, ipol) * tau (3, na)
  enddo
  !if (ionode) write(*,*) "made it here" !fail
  do n1 = - dfftp%nr1, dfftp%nr1
    arg = tpi * n1 * bgtau (1)
    eigts1 (n1) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
  enddo
  do n2 = - dfftp%nr2, dfftp%nr2
    arg = tpi * n2 * bgtau (2)
    eigts2 (n2) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
  enddo
  do n3 = - dfftp%nr3, dfftp%nr3
    arg = tpi * n3 * bgtau (3)
    eigts3 (n3) = CMPLX(cos (arg), - sin (arg) ,kind=DP)
  enddo
!  if (ionode) write(*,*) "made it here" !fail
  do ig = 1, npw
    iig = igk (ig)
    sk (ig) = kphase * eigts1 (mill (1,iig)) * &
                       eigts2 (mill (2,iig)) * &
                       eigts3 (mill (3,iig))
  enddo
   !
  do nb = 1, upf(nt)%nwfc
    if (upf(nt)%oc(nb) > 0.d0) then
       l = upf(nt)%lchi(nb)
!       if (ionode) write(*,*) "l ",l
       lphase = (0.d0,1.d0)**l
       !
       !  the factor i^l MUST BE PRESENT in order to produce
       !  wavefunctions for k=0 that are real in real space
       !
       call atomic_wfc___ ( )
       !
    END IF
    !
  END DO
   !

  !if (n_starting_wfc /= natomwfc) THEN
  !  write(*,*) "nwfc ", n_starting_wfc
  !  write(*,*) "nwfc expected ", natomwfc

    !call errore ('atomic_wfc', &
    !   'internal error: some wfcs were lost ', 1)
  !ENDIF

  deallocate(aux, sk, chiq, ylm)

  call stop_clock ('atomic_wfc')

!  if (ionode) write(*,*) "Atom " ,iatom, " done"
  return

CONTAINS
   SUBROUTINE atomic_wfc___( )
   !
   ! ... LSDA or nonmagnetic case
   !
   DO m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1
!      if (ionode) write(*,*) "wf id ",n_starting_wfc
!      if (ionode) write(*,*) "m ",m


!       if (n_starting_wfc > natomwfc) call errore &
!          ('atomic_wfc___', 'internal error: too many wfcs', 1)
      !
      DO ig = 1, npw
         wfcatom (ig, n_starting_wfc) = lphase * &
            sk (ig) * ylm (ig, lm) * chiq (ig, nb)
      ENDDO
      !
   END DO
   !
   END SUBROUTINE atomic_wfc___
END SUBROUTINE one_atom_wfc
