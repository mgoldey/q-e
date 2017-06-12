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
  USE cell_base,  ONLY : omega, tpiba, bg,alat
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  !USE basis,      ONLY : natomwfc
  USE gvect,      ONLY : mill, g ! eigts1, eigts2, eigts3
  USE klist,      ONLY : xk, igk_k, ngk
  USE fft_base,      ONLY : dfftp !, grid_scatter, grid_gather
  USE wvfct,      ONLY : npwx
  USE us,         ONLY : tab_at, dq
  USE uspp_param, ONLY : upf
  USE io_global,     ONLY : stdout,ionode, ionode_id
  USE mp_bands,   ONLY : inter_bgrp_comm, set_bgrp_indices
  USE mp,         ONLY : mp_sum

  !
  implicit none
  !
  integer, intent(in) :: ik
  integer, intent(in) :: iatom
  integer, intent(in) :: nfuncs
  complex(DP), intent(out) :: wfcatom (npwx, nfuncs)
  !
  integer :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3,  ipol, n1, n2, n3, natomwfc,npw
  integer :: ig_start, ig_end
  real(DP), allocatable :: qg(:), ylm (:,:), chiq (:,:), gk (:,:)
  complex(DP), allocatable :: sk (:)
  complex(DP) :: kphase, lphase

  complex(DP) ::      eigts1 ( -dfftp%nr1:dfftp%nr1), &
                      eigts2 ( -dfftp%nr2:dfftp%nr2), &
                      eigts3 ( -dfftp%nr3:dfftp%nr3)

  real(DP) :: arg, px, ux, vx, wx
  real(DP) :: bgtau (3)

  call start_clock ('atomic_wfc')

  ! calculate max angular momentum required in wavefunctions
  nt=ityp(iatom)
  lmax_wfc = 0
  ! natomwfc=SUM(upf(nt)%oc(:))/2
  ! natomwfc=upf(nt)%nwfc
  lmax_wfc = MAXVAL (upf(nt)%lchi(1:upf(nt)%nwfc) )
  npw = ngk(ik)
  !
  na=iatom
  !
  allocate ( ylm (npw,(lmax_wfc+1)**2), chiq(npw,upf(nt)%nwfc), &
             sk(npw), gk(3,npw), qg(npw) )
  !

  
  do ig = 1, npw
     iig = igk_k (ig,ik)
     gk (1,ig) = xk(1, ik) + g(1,iig)
     gk (2,ig) = xk(2, ik) + g(2,iig)
     gk (3,ig) = xk(3, ik) + g(3,iig)
     qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  !  ylm = spherical harmonics
  !
  call ylmr2 ((lmax_wfc+1)**2, npw, gk, qg, ylm)

  ! from now to the end of the routine the ig loops are distributed across bgrp
  call set_bgrp_indices(npw,ig_start,ig_end)
  !
  ! set now q=|k+G| in atomic units
  !
  do ig = ig_start, ig_end
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  
  do nb = 1, upf(nt)%nwfc
    if ( upf(nt)%oc (nb) >= 0.d0) then
        do ig = ig_start, ig_end
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

  deallocate (qg, gk)

  !
  wfcatom(:,:) = (0.0_dp, 0.0_dp)
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

  do ig = ig_start, ig_end
	iig = igk_k (ig,ik)
    sk (ig) = kphase * eigts1 (mill (1,iig)) * &
                       eigts2 (mill (2,iig)) * &
                       eigts3 (mill (3,iig))
  enddo
   !
  do nb = 1, upf(nt)%nwfc
    if (upf(nt)%oc(nb) > 0.d0) then
       l = upf(nt)%lchi(nb)

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

  deallocate(sk, chiq, ylm)
  ! collect results across bgrp
  call mp_sum(wfcatom, inter_bgrp_comm)

  call stop_clock ('atomic_wfc')


  return

CONTAINS
   SUBROUTINE atomic_wfc___( )
   !
   !
   DO m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1

      DO ig = ig_start, ig_end
         wfcatom (ig, n_starting_wfc) = lphase * &
            sk (ig) * ylm (ig, lm) * chiq (ig, nb)
      ENDDO
      !
   END DO
   !
   END SUBROUTINE atomic_wfc___
END SUBROUTINE one_atom_wfc


SUBROUTINE one_atom_shifted_wfc (ik, wfcatom,iatom,nfuncs,idir,dx)
  !-----------------------------------------------------------------------
  !
  ! This routine computes the superposition of atomic wavefunctions
  ! for k-point "ik" - output in "wfcatom"
  !
  USE kinds,      ONLY : DP
  USE constants,  ONLY : tpi, fpi, pi
  USE cell_base,  ONLY : omega, tpiba, bg, alat
  USE ions_base,  ONLY : nat, ntyp => nsp, ityp, tau
  USE gvect,      ONLY : mill, g ! eigts1, eigts2, eigts3
  USE klist,      ONLY : xk, igk_k, ngk
  USE fft_base,   ONLY : dfftp !, grid_scatter, grid_gather
  USE wvfct,      ONLY : npwx
  USE us,         ONLY : tab_at, dq
  USE uspp_param, ONLY : upf
  USE io_global,  ONLY : stdout,ionode, ionode_id
  USE mp_bands,   ONLY : inter_bgrp_comm, set_bgrp_indices
  USE mp,         ONLY : mp_sum

  !
  implicit none
  !
  integer, intent(in) :: ik,idir
  integer, intent(in) :: iatom
  integer, intent(in) :: nfuncs
  complex(DP), intent(out) :: wfcatom (npwx, nfuncs)
  real(DP), intent(in) :: dx
  !
  integer :: n_starting_wfc, lmax_wfc, nt, l, nb, na, m, lm, ig, iig, &
             i0, i1, i2, i3,  ipol, n1, n2, n3, natomwfc, npw
  integer :: ig_start, ig_end             
  real(DP), allocatable :: qg(:), ylm (:,:), chiq (:,:), gk (:,:)
  complex(DP), allocatable :: sk (:)
  complex(DP) :: kphase, lphase

  complex(DP) ::      eigts1 ( -dfftp%nr1:dfftp%nr1), &
                      eigts2 ( -dfftp%nr2:dfftp%nr2), &
                      eigts3 ( -dfftp%nr3:dfftp%nr3)

  real(DP) :: arg, px, ux, vx, wx, x, y, z
  real(DP) :: bgtau (3)

  call start_clock ('one_atomic_shifted_wfc')

  ! calculate max angular momentum required in wavefunctions
  nt=ityp(iatom)
  lmax_wfc = 0
  lmax_wfc = MAXVAL (upf(nt)%lchi(1:upf(nt)%nwfc) )
  npw = ngk(ik)
  !
  na=iatom
  !
  allocate ( ylm (npw,(lmax_wfc+1)**2), chiq(npw,upf(nt)%nwfc), &
             sk(npw), gk(3,npw), qg(npw) )
  !
  do ig = 1, npw
	 iig = igk_k (ig,ik)
	 gk (1,ig) = xk(1, ik) + g(1,iig)
	 gk (2,ig) = xk(2, ik) + g(2,iig)
	 gk (3,ig) = xk(3, ik) + g(3,iig)
	 qg(ig) = gk(1, ig)**2 +  gk(2, ig)**2 + gk(3, ig)**2
  enddo
  !
  !  ylm = spherical harmonics
  !
  call ylmr2 ((lmax_wfc+1)**2, npw, gk, qg, ylm)


  ! from now to the end of the routine the ig loops are distributed across bgrp
  call set_bgrp_indices(npw,ig_start,ig_end)
  !
  ! set now q=|k+G| in atomic units
  !
  do ig = ig_start, ig_end
     qg(ig) = sqrt(qg(ig))*tpiba
  enddo
  !
  n_starting_wfc = 0
  !
  ! chiq = radial fourier transform of atomic orbitals chi
  !
  
  x=tau(1,na)
  y=tau(2,na)
  z=tau(3,na)

  ! alat in bohr, so 

  ! convert dx from bohr to % of cell
  if (idir.eq.1) x=x+dx/(alat)
  if (idir.eq.2) y=y+dx/(alat)
  if (idir.eq.3) z=z+dx/(alat)

  do nb = 1, upf(nt)%nwfc
    if ( upf(nt)%oc (nb) >= 0.d0) then
        do ig = ig_start, ig_end
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

  deallocate(qg, gk)

  !
  wfcatom(:,:) = (0.0_dp, 0.0_dp)

  !
  arg = (xk(1,ik)*tau(1,na) + xk(2,ik)*tau(2,na) + xk(3,ik)*tau(3,na)) * tpi
  kphase = CMPLX(cos (arg), - sin (arg) ,kind=DP)
  !
  !     sk is the structure factor
  !
  ! from struct_fact
  do ipol = 1, 3
    bgtau (ipol) = bg (1, ipol) * x + &
                   bg (2, ipol) * y + &
                   bg (3, ipol) * z
  enddo

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

  do ig = ig_start, ig_end
    iig = igk_k (ig,ik)
    sk (ig) = kphase * eigts1 (mill (1,iig)) * &
                       eigts2 (mill (2,iig)) * &
                       eigts3 (mill (3,iig))
  enddo
   !
  do nb = 1, upf(nt)%nwfc
    if (upf(nt)%oc(nb) > 0.d0) then
       l = upf(nt)%lchi(nb)
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

  deallocate(sk, chiq, ylm)

  ! collect results across bgrp
  call mp_sum(wfcatom, inter_bgrp_comm)

  call stop_clock ('one_atomic_shifted_wfc')

  return

CONTAINS
   SUBROUTINE atomic_wfc___( )
   DO m = 1, 2 * l + 1
      lm = l**2 + m
      n_starting_wfc = n_starting_wfc + 1

      DO ig = ig_start, ig_end
         wfcatom (ig, n_starting_wfc) = lphase * &
            sk (ig) * ylm (ig, lm) * chiq (ig, nb)
      ENDDO
      !
   END DO
   !
   END SUBROUTINE atomic_wfc___
END SUBROUTINE one_atom_shifted_wfc

