!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_diabatic
  !-----------------------------------------------------------------------
  !
  USE epcdft_mod, ONLY : eig_of_w
  !
  IMPLICIT NONE
  !
  IF(eig_of_w) THEN
    !
    CALL epcdft_get_diabatic_w()
    !
  ELSE
    !
    CALL epcdft_get_diabatic_lowdin()
    !
  ENDIF
  !
END SUBROUTINE epcdft_get_diabatic
!
!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_diabatic_w()
  !-----------------------------------------------------------------------
  !
  !  use eigen vecs of W following :
  !
  !    The Journal of Chemical Physics 133, 244105 (2010); doi: 10.1063/1.3507878
  !
  USE kinds, ONLY : DP
  USE epcdft_mod, ONLY : wmat, smat, ohc, hc, debug2, lm
  !
  IMPLICIT NONE
  !
  COMPLEX(DP) :: stot(2,2), wtot(2,2), stotinv(2,2) ! Sup*Sdown, W_up+W_down
  COMPLEX(DP) :: w_eigvec(2,2), vl(2,2), work(2)
  !REAL(DP) :: w_eigval(2)
  COMPLEX(DP) :: w_eigval(2) 
  !COMPLEX(DP) :: hp ! S^-1*hc
  CHARACTER(LEN=256) :: fname
  INTEGER :: filunit=3234873, info
  REAL(DP)::rwork(4)
  !
  w_eigval = 0.D0
  ohc = (0.D0,0.D0)
  !hp = (0.D0,0.D0)
  stot = (0.D0,0.D0)
  stotinv = (0.D0,0.D0)
  wtot = (0.D0,0.D0)
  w_eigvec = (0.D0,0.D0)
  !
  stot(:,:) = smat(:,:,1) * smat(:,:,2)
 !
 ! S^-1
 !
  CALL invrts(stot,stotinv)
  !
  ! S^-1*hc
  !
  !hp = MATMUL(stotinv,hc)
  !
  ! Solve Wab = wabup*sabdown + wabdown*sabup
  !
  wtot(:,:) = wmat(:,:,1) * smat(:,:,2) + wmat(:,:,2) * smat(:,:,1)
  wtot(:,1) = wtot(:,1)/lm(1,1)
  wtot(:,2) = wtot(:,2)/lm(1,2)
  !
  !
  IF(debug2)THEN
    fname='W'
    CALL realpart_dumpmat(fname,filunit,wtot,2,2)
  ENDIF
  !
  ! Solve W.V = S.V * eig
  !
  ohc = MATMUL(stotinv, wtot) ! use ohc as work array here
  CALL ZGEEV( 'N', 'V', 2, ohc,   2, w_eigval, vl, 2, w_eigvec, 2, work,&
                        4, rwork, info )
  !print *,'info',info
  !CALL findnormvec(w_eigvec(:,1))
  !CALL findnormvec(w_eigvec(:,2))
  !
  IF(debug2)THEN
    fname='Weigvec'
    CALL realpart_dumpmat(fname,filunit,w_eigvec,2,2)
  ENDIF
  !
  ! H_orthg =  v^(dagger).Hc.v
  !
  ohc = MATMUL( TRANSPOSE(DCONJG(w_eigvec)) , hc )
  ohc = MATMUL( ohc , w_eigvec ) 
  !
  IF(debug2)THEN
    fname='Horth'
    CALL realpart_dumpmat(fname,filunit,ohc,2,2)
  ENDIF
  !
END SUBROUTINE epcdft_get_diabatic_w
!
!
!-----------------------------------------------------------------------
SUBROUTINE epcdft_get_diabatic_lowdin()
  !-----------------------------------------------------------------------
  !
  USE kinds, ONLY : DP
  USE epcdft_mod, ONLY : smat, hc, ohc
  USE klist, ONLY : nks
  !
  IMPLICIT NONE
  !
  INTEGER i
  REAL(DP) :: l(2)            ! sev eigvals of S matrix
  COMPLEX(DP) :: u(2,2)       ! U  S eigenvector matrix sev(i,j) 
  COMPLEX(DP) :: d(2,2)       ! SD^-1/2 is diagonal S^-1/2
  COMPLEX(DP) :: invssqr(2,2) ! S-1/2
  COMPLEX(DP) :: smattot(2,2)   ! Sup*Sdown
  !
  ohc = 0.D0
  l = 0.D0
  u = 0.D0
  d = 0.D0
  invssqr = 0.D0
  smattot(:,:) = smat(:,:,1)*smat(:,:,2)
  !
  !
  ! get U and U^-1 and eigenvals (l) of S
  CALL cdiagh(2,smattot,2,l,u)
  !
  ! get diag S^-1/2
  DO i = 1, 2
    d(i,i) = l( i )**(-0.5D0)
  ENDDO
  !
  !
  ! get S^-1/2 = U^-1 . SD^-1/2 . U
  invssqr(:,:) = MATMUL( u(:,:), d(:,:) )
  invssqr(:,:) = MATMUL( invssqr(:,:), Transpose(u(:,:)) )
  !
  ! CALL print_cmat("S**-.5",invssqr,2)
  !
  ! get oHc = S^-1/2 . Hc . S^-1/2
  ohc(:,:) = MATMUL( invssqr(:,:), hc(:,:) )
  ohc(:,:) = MATMUL( ohc(:,:) , invssqr(:,:) )
  !
  !
END SUBROUTINE epcdft_get_diabatic_lowdin
!
SUBROUTINE invrts(s,sinv)
  !
  ! compute s^-1
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  COMPLEX(DP),INTENT(IN)::s(2,2)
  COMPLEX(DP),INTENT(INOUT)::sinv(2,2)
  COMPLEX(DP)::det
  !
  sinv = (0.D0,0.D0)
  det = (0.D0,0.D0)
  !
  det = (s(1,1)*s(2,2) - s(1,2)*s(2,1))
  sinv(1,1) = s(2,2)
  sinv(1,2) = -s(1,2)
  sinv(2,1) = -s(2,1)
  sinv(2,2) = s(1,1)
  sinv(:,:) = sinv(:,:)/det
  !
ENDSUBROUTINE invrts
!
SUBROUTINE findnormvec(v)
  !
  ! normalize vector
  !
  USE kinds, ONLY : DP
  IMPLICIT NONE
  COMPLEX(DP),INTENT(INOUT)::v(2)
  COMPLEX(DP)::norm
  !
  norm = v(1)*CONJG(v(1)) + v(2)*CONJG(v(2))
  v(:) = v(:)/SQRT(norm)
  !
ENDSUBROUTINE findnormvec
