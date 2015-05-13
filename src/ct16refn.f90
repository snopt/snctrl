!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File:  ct16refn.f90
!   Subroutines for adaptive refinement.
!
! 24 Dec 2008: Updated for v4 of snctrl.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module ct16refn
  use precision, only : ip, rp
  use ct02data,  only : ctGrid
  use ct15usr,   only : odeconS, odeconA, odeconD
  implicit none

  private
  public  :: c1refn
  private :: c1grid

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c1refn( iIntf, refTol, nPhs, nY, nU, nP, neJ, n, &
                     x, hs, grid, usrode, nAdded, &
                     cu, lencu, iu, leniu, ru, lenru, &
                     cw, lencw, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: iIntf, nPhs, nY, nU, nP, neJ(nPhs), n, &
                                  lencu, leniu, lenru, lencw, leniw, lenrw
    real(rp),    intent(in)    :: refTol

    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)

    integer(ip), intent(inout), pointer :: hs(:)
    real(rp),    intent(inout), pointer :: x(:)

    integer(ip), intent(out)   :: nAdded

    external             :: usrode
    type(ctGrid)         :: grid

    !===========================================================================
    ! Adaptive Refinement
    ! 1. Interpolate midpoints
    ! 2. Calculate HS estimates
    ! 3. Compare HS with TR
    ! 4. Determine whether a point is added.
    !
    ! On entry:
    !    x, hs are solutions from a previous run.
    ! On exit:
    !    grid contains the new refined grid.
    !    x, hs are new refined estimates of the solution.
    !
    ! 24 Dec 2008: First version of c1refn.
    ! 09 Feb 2010: v5.
    !===========================================================================

    nAdded = 0

    if ( iw(iIntf) == 0 ) then        ! S version
       call c1grid ( refTol, nPhs, nY, nU, nP, neJ, n, x, hs, &
                     grid, nAdded, odeconS, usrode, &
                     cu, lencu, iu, leniu, ru, lenru, &
                     cw, lencw, iw, leniw, rw, lenrw )

    else if ( iw(iIntf) == 1 ) then   ! D version
       call c1grid ( refTol, nPhs, nY, nU, nP, neJ, n, x, hs, &
                     grid, nAdded, odeconD, usrode, &
                     cu, lencu, iu, leniu, ru, lenru, &
                     cw, lencw, iw, leniw, rw, lenrw )

    else if ( iw(iIntf) == 2 ) then   ! A version
       call c1grid ( refTol, nPhs, nY, nU, nP, neJ, n, x, hs, &
                     grid, nAdded, odeconA, usrode, &
                     cu, lencu, iu, leniu, ru, lenru, &
                     cw, lencw, iw, leniw, rw, lenrw )
    end if


  end subroutine c1refn

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c1grid ( refTol, nPhs, nY, nU, nP, neJ, n, x, hs, &
                      grid, nAdded, odecon, usrode, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )

    integer(ip), intent(in) :: nPhs, nY, nU, nP, neJ(nPhs), n, &
                               lencu, leniu, lenru, lencw, leniw, lenrw
    real(rp),    intent(in) :: refTol

    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)

    integer(ip), intent(inout), pointer :: hs(:)
    real(rp),    intent(inout), pointer :: x(:)

    integer(ip), intent(out)   :: nAdded

    external             :: usrode, odecon
    type(ctGrid)         :: grid

    !===========================================================================
    ! Given data from a previous grid, a refined grid is created.
    ! On entry:
    !   nInt, ndPtr, and step are for the old grid.
    ! On exit:
    !   grid is the refined grid.
    !   x, hs are the new estimates.
    !
    ! 24 Dec 2008: First version of c1grid for v4 of snctrl.
    ! 02 Feb 2009: Adapted for all versions of snctrl (usrode).
    ! 09 Feb 2010: v5.
    !=============================================================================
    integer(ip) :: Status, need1, need2, nVar, pind, nrefNd, indx, indrx, &
                   indS, indE, ndS, ndE, nNd, nJ, p, j, i, nn
    real(rp)    :: hstep, val

    integer(ip) :: dflag(nY), refhs((nY+nU)*(2*grid%nInt+1))
    real(rp)    :: midpt(nY+nU), F2(nY,1), hsx(nY), &
                   refstep(2*grid%nInt), refx((nY+nU)*(2*grid%nInt+1))

    integer(ip), pointer :: prevhs(:)
    real(rp),    pointer :: prevx(:)

    integer(ip), allocatable :: rndPtr(:), Jrow(:), Jcol(:)
    real(rp),    allocatable :: F(:,:), Jval(:,:)


    Status = 0
    need1  = 1
    need2  = 0
    nVar   = nY+nU
    pind   = nVar*grid%nNodes + 1

    allocate ( rndPtr(nPhs+1) )

    refhs  = 0
    nrefNd = 0
    indx   = 0
    indrx  = 0
    indS   = 1

    allocate ( Jcol(nY+nU+nP+1) )

    do p = 1, nPhs
       rndPtr(p) = nrefNd + 1

       nJ   = neJ(p)
       ndS  = grid%ndPtr(p)
       ndE  = grid%ndPtr(p+1)
       nNd  = ndE - ndS + 1
       indE = indS-1 + nNd*(nY+nU)

       ! Get user-defined state equations
       allocate( F(nY,ndS:ndE), Jrow(nJ), Jval(nJ,ndS:ndE) )
       call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNd, dFlag, &
                    F, Jrow, Jval, Jcol, nJ, x(indS:indE), x(pind), &
                    1_ip, 0_ip, &
                    cu, lencu, iu, leniu, ru, lenru, &
                    cw, lencw, iw, leniw, rw, lenrw )

       indS = indE + 1

       do j = grid%ndPtr(p), grid%ndPtr(p+1)-1
          hstep = grid%step(j)

          ! Midpoint
          do i = 1, nY
             midpt(i) = 0.5*( x(i+indx) + x(i+indx+nVar) ) + &
                  (hstep/8.0)*( F(i,j) - F(i,j+1) )
          end do
          do i = 1+nY, nU+nY
             midpt(i) = 0.5*( x(i+indx) + x(i+indx+nVar) )
          end do

          ! Get user-defined state equations at midpoints
          call odecon( Status, usrode, p, nPhs, nY, nU, nP, 1_ip, dFlag, &
                       F2, Jrow, Jval, Jcol, nJ, midpt, x(pind), &
                       1_ip, 0_ip, &
                       cu, lencu, iu, leniu, ru, lenru, &
                       cw, lencw, iw, leniw, rw, lenrw )

          !-------------------------------------------------------------------
          ! HS values
          !-------------------------------------------------------------------
          do i = 1, nY
             hsx(i) = x(i+indx)+(hstep/6.0)*( F(i,j)+4.0*F2(i,1)+ F(i,j+1) )
          end do

          val = maxval (abs (hsx - x(1+indx+nVar:nY+indx+nVar)) )

          if ( val > refTol ) then     ! Add a point
             nAdded  = nAdded + 1
             ! Set up new step sizes.
             nrefNd = nrefNd + 1
             refstep(nrefNd) = hstep/2.0

             nrefNd = nrefNd + 1
             refstep(nrefNd) = hstep/2.0

             ! Copy old x, hs values at point.
             refx(1+indrx:nVar+indrx)  = x(1+indx:nVar+indx)
             refhs(1+indrx:nVar+indrx) = hs(1+indx:nVar+indx)
             indrx = indrx + nVar

             ! Copy midpoint values.
             refx(1+indrx:nVar+indrx) = midpt
             refhs(1+indrx:nY+indrx) = 3
             refhs(1+nY+indrx:nU+nY+indrx) = 0
             indrx = indrx + nVar

          else
             nrefNd = nrefNd + 1
             refstep(nrefNd) = hstep

             refx(1+indrx:nVar+indrx)  = x(1+indx:nVar+indx)
             refhs(1+indrx:nVar+indrx) = hs(1+indx:nVar+indx)
             indrx = indrx + nVar

          end if

          indx = indx + nVar
       end do

       refx(1+indrx:nVar+indrx)  =  x(1+indx:nVar+indx)
       refhs(1+indrx:nVar+indrx) = hs(1+indx:nVar+indx)
       indx = indx + nVar   ! for phase continuity vars

       deallocate ( F, Jrow, Jval )
    end do

    deallocate ( Jcol )

    rndPtr(nPhs+1) = nrefNd + 1

    if ( nAdded == 0 ) &
         return

    !-----------------------------------
    ! Set refined grid.
    !   nNodes = nrefNd + 1
    !   nInt   = nrefNd
    !   step   = refstep

    grid%nInt   = nrefNd
    grid%nNodes = nrefNd + 1

    if( allocated(grid%intPtr) ) deallocate( grid%intPtr )
    if( allocated(grid%ndPtr) )  deallocate( grid%ndPtr )
    if( allocated(grid%step) )   deallocate( grid%step )
    allocate ( grid%step(nrefNd), grid%intPtr(nPhs+1), grid%ndPtr(nPhs+1) )

    grid%step   = refstep(1:nrefNd)
    grid%intPtr = rndPtr
    grid%ndPtr  = rndPtr

    deallocate ( rndPtr )

    ! Set up new x, hs values.
    prevx  => x(1:n)
    prevhs => hs(1:n)

    nn = (nY+nU)*grid%nNodes + nP
    nullify  ( x, hs )
    allocate ( x(nn), hs(nn) )

    x(1:nn-nP)  = refx(1:nn-nP)
    hs(1:nn-nP) = refhs(1:nn-nP)

    ! Parameters
    if ( nP > 0 ) then
       x(nn-nP+1:nn)  = prevx(n-nP+1:n)
       hs(nn-nP+1:nn) = prevhs(n-nP+1:n)
    end if

    deallocate ( prevx, prevhs )

  end subroutine c1grid

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module ct16refn
