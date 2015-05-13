!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File:   ct28derv.f90
! Subroutines to set up the derivative structures for the user-defined state
! equations and algebraic constraints.  All routines in this file assume that
! the user-defined Jacobians are in sparse structures and are defined at every
! node.
!
! 30 Dec 2008: First version for v4 of snctrl.
! 18 Jan 2009: Sparse structures for user Jacobians implemented.
! 06 Feb 2010: v5.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module ct28derv
  use precision, only : ip, rp
  implicit none

  private
  public :: odeSetupTR, algSetupTR, odeSetupHS, algSetupHS

  interface

     subroutine ictrode( Status, usrode, curPhs, nPhs, nY, nU, nP, nNodes, &
                         dFlag, F, Jrow, Jval, Jcol, lenJ, dvar, pvar, &
                         needF, needJ, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
       use precision, only : ip, rp
       implicit none

       integer(ip),  intent(in)    :: Status, curPhs, nPhs, nY, nU, nP, nNodes, &
                                      lenJ, needF, needJ, lencu, leniu, lenru, &
                                      lencw, leniw, lenrw
       real(rp),     intent(in)    :: dvar(nY+nU,nNodes), pvar(nP)

       integer(ip),  intent(inout) :: iu(leniu), iw(leniw)
       character(8), intent(inout) :: cu(lencu), cw(lencw)
       real(rp),     intent(inout) :: ru(lenru), rw(lenrw)

       integer(ip),  intent(out)   :: Jrow(lenJ), Jcol(nY+nU+nP+1), dFlag(nY)
       real(rp),     intent(out)   :: F(nY,nNodes), Jval(lenJ,nNodes)
       external                    :: usrode
     end subroutine ictrode


     subroutine ictralg( Status, usralg, curPhs, nPhs, nC, nY, nU, nP, nNodes, &
                         C, Grow, Gval, Gcol, lenG, dvar, pvar, &
                         needC, needG, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
       use precision, only : ip, rp
       implicit none

       integer(ip),  intent(in)    :: Status, curPhs, nPhs, nY, nU, nP, nC, nNodes, &
                                      lenG, needC, needG, &
                                      lencu, leniu, lenru, lencw, leniw, lenrw
       real(rp),     intent(in)    :: dvar(nY+nU,nNodes), pvar(nP)

       integer(ip),  intent(inout) :: iu(leniu), iw(leniw)
       character(8), intent(inout) :: cu(lencu), cw(lencw)
       real(rp),     intent(inout) :: ru(lenru), rw(lenrw)

       integer(ip),  intent(out)   :: Grow(lenG), Gcol(nY+nU+nP+1)
       real(rp),     intent(out)   :: C(nC,nNodes), Gval(lenG,nNodes)
       external                    :: usralg
     end subroutine ictralg
  end interface

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine odeSetupTR( nPhs, nY, nU, nP, objL, ytype, nnY, nlY, &
                         nnC, nlC, neJ, nnCon, ndPtr, step, s, n, m, &
                         lenJ, nindJ, nlocJ, lJcol, lindJ, llocJ, &
                         odecon, usrode, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
    integer(ip),  intent(in)    ::  nPhs, nY, nU, nP, objL(nPhs), ytype(nY,nPhs), &
                                   nnY(nPhs), nlY(nPhs), nnC(nPhs), nlC(nPhs), &
                                   neJ(nPhs), nnCon, ndPtr(nPhs+1), s, lenJ, &
                                   n, m, lencu, leniu, lenru, lencw, leniw, lenrw
    real(rp),     intent(in)    :: step(s)

    integer(ip),  intent(inout) :: iu(leniu), iw(leniw)
    character(8), intent(inout) :: cu(lencu), cw(lencw)
    real(rp),     intent(inout) :: ru(lenru), rw(lenrw)

    integer(ip),  intent(out)   :: nindJ(lenJ), nlocJ(n+1), lindJ(lenJ), llocJ(n+1)
    real(rp),     intent(out)   :: lJcol(lenJ)
    procedure(ictrode)          :: odecon
    external                    :: usrode

    !===========================================================================
    ! Trapezoid discretization.
    ! Sets up the SNOPT structures for the ODE constraint derivatives.
    ! Differentiates linear and nonlinear constraints.
    !
    !           y,u(t0) | y,u(t1) | .... | y,u(t) | y,u(tf) | p
    ! C(t0)
    ! F(t0,t1)     x        x                                 x
    ! C(t1)
    ! F(t1,t2)              x        x                        x
    !
    !   ...
    ! C(t)
    ! F(t, tf)                               x        x       x
    ! C(tf)
    !
    ! 28 Dec 2008: First version of odeSetupTR for snctrl v4.
    ! 18 Jan 2009: Implemented sparse data structures for user
    !              Jacobian.
    ! 22 Jan 2009: Redid phases.
    !===========================================================================
    integer(ip) :: Status, nVar, colind, ncount, lcount, &
                   nrwind, lrwind, nNodes, ndS, ndE, nJ, &
                   yrows(nY,nPhs), dFlag(nY), &
                   p, i, j, k, q, in, il
    real(rp)    :: phstep, hstep, dvar(n-nP), pvar(nP)

    integer(ip), allocatable :: usrrow(:), usrcol(:)
    real(rp),    allocatable :: F(:,:), usrval(:,:)


    Status = 0
    dvar   = 1.0
    pvar   = 1.0
    nVar   = nY+nU

    ! Set up row indices for nonlinear/linear state eqns.
    do p = 1, nPhs
       in = 0
       il = 0
       do i = 1, nY
          if ( ytype(i,p) == 0 ) then
             in = in + 1
             yrows(i,p) = in
          else
             il = il + 1
             yrows(i,p) = il
          end if
       end do
    end do


    ! Initialize variables
    nindJ  = 0
    nlocJ  = 0

    lJcol  = 0
    lindJ  = 0
    llocJ  = 0

    ncount  = 1
    lcount  = 1

    colind  = 0
    nrwind  = 0
    lrwind  = nnCon

    ! And let's go...
    allocate ( usrcol(nY+nU+nP+1) )

    do p = 1, nPhs

       nJ  = neJ(p)
       ndS = ndPtr(p)
       ndE = ndPtr(p+1)
       nNodes = ndE - ndS + 1

       ! Get user-defined Jacobian elements
       allocate ( F(nY,ndS:ndE) )
       allocate ( usrrow(nJ), usrval(nJ,ndS:ndE) )
       call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                    F, usrrow, usrval, usrcol, nJ, dvar, pvar, 0_ip, 1_ip, &
                    cu, lencu, iu, leniu, ru, lenru, &
                    cw, lencw, iw, leniw, rw, lenrw )

       ! At first node of the phase:
       j = ndS
       hstep = step(j)

       nrwind = nrwind + nnC(p)
       lrwind = lrwind + nlC(p)
       if ( p /= 1 ) &
            lrwind = lrwind + nY

       do k = 1, nY+nU
          nlocJ(k+colind) = ncount
          llocJ(k+colind) = lcount

          ! Phase continuity
          if ( p > 1 .and. k <= nY ) then
             lindJ(lcount) = k + lrwind - nY
             lJcol(lcount) = 1d+0
             lcount = lcount + 1
          end if

          ! ODE Constraints
          do q = usrcol(k), usrcol(k+1)-1
             i = usrrow(q)

             if ( ytype(i,p) == 0 ) then
                nindJ(ncount) = yrows(i,p) + nrwind
                ncount = ncount + 1
             else
                lindJ(lcount) = yrows(i,p) + lrwind
                lJcol(lcount) = (hstep/2d+0)*usrval(q,j)
                if ( i == k ) &
                     lJcol(lcount) = lJcol(lcount) + 1d+0
                lcount = lcount + 1
             end if
          end do

          if ( k <= nY ) then
             if ( dFlag(k) == 0 ) then
                if ( ytype(k,p) == 0 ) then
                   nindJ(ncount) = yrows(k,p) + nrwind
                   ncount = ncount + 1
                else
                   lindJ(lcount) = yrows(k,p) + lrwind
                   lJcol(lcount) = 1d+0
                   lcount = lcount + 1
                end if
             end if
          end if
       end do


       ! At inner nodes:
       do j = ndS+1, ndE-1
          colind = colind + nVar
          phstep = hstep
          hstep  = step(j)

          do k = 1, nY+nU
             nlocJ(k+colind) = ncount
             llocJ(k+colind) = lcount

             ! ODE Constraints
             do q = usrcol(k), usrcol(k+1)-1
                i = usrrow(q)

                if (ytype(i,p) == 0) then
                   nindJ(ncount) = yrows(i,p) + nrwind
                   ncount = ncount + 1

                   nindJ(ncount) = yrows(i,p) + nrwind + nnC(p) + nnY(p)
                   ncount = ncount + 1
                else
                   lindJ(lcount) = yrows(i,p) + lrwind
                   lJcol(lcount) = (phstep/2d+0)*usrval(q,j)
                   lcount = lcount + 1

                   lindJ(lcount) = yrows(i,p) + lrwind + nlC(p) + nlY(p)
                   lJcol(lcount) = (hstep/2d+0)*usrval(q,j)
                   lcount = lcount + 1

                   if ( i == k ) then
                      lJcol(lcount-2) = lJcol(lcount-2) - 1.0
                      lJcol(lcount-1) = lJcol(lcount-1) + 1.0
                   end if
                end if
             end do

             if ( k <= nY ) then
                if ( dFlag(k) == 0 ) then
                   if ( ytype(k,p) == 0 ) then
                      nindJ(ncount) = yrows(k,p) + nrwind
                      ncount = ncount + 1

                      nindJ(ncount) = yrows(k,p) + nrwind + nnC(p) + nnY(p)
                      ncount = ncount + 1
                   else
                      lindJ(lcount) = yrows(k,p) + lrwind
                      lJcol(lcount) = -1.0
                      lcount = lcount + 1

                      lindJ(lcount) = yrows(k,p) + lrwind + nlC(p) + nlY(p)
                      lJcol(lcount) = 1.0
                      lcount = lcount + 1
                   end if
                end if
             end if

          end do

          nrwind = nrwind + nnY(p) + nnC(p)
          lrwind = lrwind + nlY(p) + nlC(p)
       end do

       ! At final node of phase:
       j = ndE
       phstep = hstep
       colind = colind + nVar

       do k = 1, nY+nU
          nlocJ(k+colind) = ncount
          llocJ(k+colind) = lcount

          ! ODE Constraints
          do q = usrcol(k), usrcol(k+1)-1
             i = usrrow(q)

             if ( ytype(i,p) == 0 ) then
                nindJ(ncount) = yrows(i,p) + nrwind
                ncount = ncount + 1
             else
                lindJ(lcount) = yrows(i,p) + lrwind
                lJcol(lcount) = (phstep/2.0)*usrval(q,j)
                if ( i == k ) &
                     lJcol(lcount) = lJcol(lcount) - 1d+0
                lcount = lcount + 1
             end if
          end do

          if ( k <= nY ) then
             ! Diagonals
             if ( dFlag(k) == 0 ) then
                if ( ytype(k,p) == 0 ) then
                   nindJ(ncount) = yrows(k,p) + nrwind
                   ncount = ncount + 1
                else
                   lindJ(lcount) = yrows(k,p) + lrwind
                   lJcol(lcount) = -1.0
                   lcount = lcount + 1
                end if
             end if

             ! Equality constraints for continuity at phases
             if ( p < nPhs ) then
                lindJ(lcount) = k + lrwind + nlY(p)
                lJcol(lcount) = -1.0
                lcount = lcount + 1
             end if
          end if

          ! Objective
          if ( k == objL(p) .or. k+nY == objL(p) ) then
             lindJ(lcount) = m
             lJcol(lcount) = 1.0
             lcount = lcount + 1
          end if
       end do

       colind = colind + nVar
       nrwind = nrwind + nnC(p) + nnY(p)
       lrwind = lrwind + nlC(p) + nlY(p)

       deallocate ( F, usrrow, usrval )
    end do


    ! Parameters
    nrwind = 0
    lrwind = nnCon

    do k = 1, nP
       nlocJ(k+colind) = ncount
       llocJ(k+colind) = lcount

       do p = 1, nPhs

          nJ  = neJ(p)
          ndS = ndPtr(p)
          ndE = ndPtr(p+1)
          nNodes = ndE - ndS + 1

          ! Get user-defined state equations
          allocate ( F(nY,ndS:ndE) )
          allocate ( usrrow(nJ), usrval(nJ,ndS:ndE) )
          call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                       F, usrrow, usrval, usrcol, nJ, dvar, pvar, 0_ip, 1_ip, &
                       cu, lencu, iu, leniu, ru, lenru, &
                       cw, lencw, iw, leniw, rw, lenrw )

          nrwind = nrwind + nnC(p)
          lrwind = lrwind + nlC(p)

          do j = ndS, ndE-1
             hstep = step(j)

             do q = usrcol(k+nVar), usrcol(k+nVar+1)-1
                i = usrrow(q)

                if ( ytype(i,p) == 0 ) then
                   nindJ(ncount) = yrows(i,p) + nrwind
                   ncount = ncount + 1
                else
                   lindJ(lcount) = yrows(i,p) + lrwind
                   lJcol(lcount) = hstep*usrval(q,j)
                   lcount = lcount + 1
                end if
             end do
             nrwind = nrwind + nnY(p) + nnC(p)
             lrwind = lrwind + nlY(p) + nlC(p)
          end do

          lrwind = lrwind + nY
          deallocate ( F, usrrow, usrval )

          ! Objective
          if ( k+nY+nU == objL(p) ) then
             lindJ(lcount) = m
             lJcol(lcount) = 1.0
             lcount = lcount + 1
          end if
       end do
    end do

    deallocate ( usrcol )

    nlocJ(n+1) = ncount
    llocJ(n+1) = lcount

  end subroutine odeSetupTR

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine algSetupTR ( nPhs, nY, nU, nP, nC, ctype, nnY, nlY, &
                          nnC, nlC, neG, nnCon, ndPtr, n, lenJ, &
                          nindJ, nlocJ, lJcol, lindJ, llocJ, &
                          algcon, usralg, &
                          cu, lencu, iu, leniu, ru, lenru, &
                          cw, lencw, iw, leniw, rw, lenrw )
    integer(ip),  intent(in)    :: nPhs, nY, nU, nP, nC, ctype(nC,nPhs), &
                                   nnY(nPhs), nlY(nPhs), nnC(nPhs), nlC(nPhs), &
                                   neG(nPhs), nnCon, ndPtr(nPhs+1), n, lenJ, &
                                   lencu, leniu, lenru, lencw, leniw, lenrw

    integer(ip),  intent(inout) :: iu(leniu), iw(leniw)
    character(8), intent(inout) :: cu(lencu), cw(lencw)
    real(rp),     intent(inout) :: ru(lenru), rw(lenrw)

    integer(ip),  intent(out)   :: nindJ(lenJ), nlocJ(n+1), &
                                   lindJ(lenJ), llocJ(n+1)
    real(rp),     intent(out)   :: lJcol(lenJ)
    procedure(ictralg)          :: algcon
    external                    :: usralg

    !===========================================================================
    ! Sets up the SNOPT structures for the algebraic constraint derivatives.
    ! Differentiates linear and nonlinear constraints.
    !
    !         y,u(t0) | y,u(t1) | .... | y,u(tf) | p
    ! C(t0)      x                                 x
    ! C(t1)                x                       x
    !   ...
    ! C(tf)                                x       x
    !
    !
    ! 26 Dec 2008: First version of algSetup for snctrl v4.
    ! 18 Jan 2009: Sparse structures implemented.
    !===========================================================================
    integer(ip) :: Status, nVar, colind, ncount, lcount, nrwind, lrwind, &
                   crows(nC,nPhs), nNodes, ndS, ndE, nG, p, i, j, k, q, in, il
    real(rp)    :: dvar(n-nP), pvar(nP)
    integer(ip), allocatable :: usrrow(:), usrcol(:)
    real(rp),    allocatable :: C(:,:), usrval(:,:)


    Status = 0
    dvar   = 0.0
    pvar   = 0.0
    nVar   = nY+nU

    ! Set up row indices for nonlinear/linear state eqns.
    do p = 1, nPhs
       in = 0
       il = 0
       do i = 1, nC
          if ( ctype(i,p) == 0 ) then
             in = in + 1
             crows(i,p) = in
          else
             il = il + 1
             crows(i,p) = il
          end if
       end do
    end do

    ! Initialize variables
    nindJ = 0
    nlocJ = 0

    lJcol = 0
    lindJ = 0
    llocJ = 0

    nVar   = nY+nU
    ncount = 1
    lcount = 1

    colind = 0
    nrwind = 0
    lrwind = nnCon

    ! Start...
    allocate ( usrcol(nY+nU+nP+1) )
    do p = 1, nPhs

       nG  = neG(p)
       ndS = ndPtr(p)
       ndE = ndPtr(p+1)
       nNodes = ndE - ndS + 1

       ! Get user-defined Jacobian elements
       allocate( C(nC,ndS:ndE) )
       allocate( usrrow(nG), usrval(nG,ndS:ndE) )
       call algcon( Status, usralg, p, nPhs, nC, nY, nU, nP, nNodes, &
                    C, usrrow, usrval, usrcol, nG, dvar, pvar, 0_ip, 1_ip, &
                    cu, lencu, iu, leniu, ru, lenru, &
                    cw, lencw, iw, leniw, rw, lenrw )


       do j = ndS, ndE
          do k = 1, nY+nU
             nlocJ(k+colind) = ncount
             llocJ(k+colind) = lcount

             do q = usrcol(k), usrcol(k+1)-1
                i = usrrow(q)

                if ( ctype(i,p) == 0 ) then
                   nindJ(ncount) = crows(i,p) + nrwind
                   ncount = ncount + 1
                else
                   lindJ(lcount) = crows(i,p) + lrwind
                   lJcol(lcount) = usrval(q,j)
                   lcount = lcount + 1
                end if
             end do
          end do
          colind = colind + nVar

          nrwind = nrwind + nnC(p) + nnY(p)
          lrwind = lrwind + nlC(p) + nlY(p)
       end do

       deallocate ( C, usrrow, usrval )
    end do

    ! Parameters
    nrwind = 0
    lrwind = nnCon

    do k = 1, nP
       nlocJ(k+colind) = ncount
       llocJ(k+colind) = lcount

       do p = 1, nPhs

          nG  = neG(p)
          ndS = ndPtr(p)
          ndE = ndPtr(p+1)
          nNodes = ndE - ndS + 1

          ! Get user-defined Jacobian elements
          allocate( C(nC,ndS:ndE) )
          allocate( usrrow(nG), usrval(nG,ndS:ndE) )
          call algcon( Status, usralg, p, nPhs, nC, nY, nU, nP, nNodes, &
                       C, usrrow, usrval, usrcol, nG, dvar, pvar, 0_ip, 1_ip, &
                       cu, lencu, iu, leniu, ru, lenru, &
                       cw, lencw, iw, leniw, rw, lenrw )

          do j = ndS, ndE
             do q = usrcol(k+nVar), usrcol(k+nVar+1)-1
                i = usrrow(q)

                if ( ctype(i,p) == 0 ) then
                   nindJ(ncount) = crows(i,p) + nrwind
                   ncount = ncount + 1
                else
                   lindJ(lcount) = crows(i,p) + lrwind
                   lJcol(lcount) = usrval(q,j)
                   lcount = lcount + 1
                end if
             end do
             nrwind = nrwind + nnC(p) + nnY(p)
             lrwind = lrwind + nlC(p) + nlY(p)
          end do

          deallocate ( C, usrrow, usrval )
       end do
    end do

    deallocate ( usrcol )

    nlocJ(n+1) = ncount
    llocJ(n+1) = lcount

  end subroutine algSetupTR

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine odeSetupHS( nPhs, nY, nU, nP, objL, ytype, nnY, nlY, &
                         nnC, nlC, neJ, nnCon, ndPtr, step, s, n, m, &
                         lenJ, nindJ, nlocJ, lJcol, lindJ, llocJ, &
                         odecon, usrode, &
                         cu, lencu, iu, leniu, ru, lenru, &
                         cw, lencw, iw, leniw, rw, lenrw )
    integer(ip),  intent(in)    :: nPhs, nY, nU, nP, objL(nPhs), ytype(nY,nPhs), &
                                   nnY(nPhs), nlY(nPhs), nnC(nPhs), nlC(nPhs), &
                                   neJ(nPhs), nnCon, ndPtr(nPhs+1), s, lenJ, &
                                   n, m, &
                                   lencu, leniu, lenru, lencw, leniw, lenrw
    real(rp),     intent(in)    :: step(s)

    integer(ip),  intent(inout) :: iu(leniu), iw(leniw)
    character(8), intent(inout) :: cu(lencu), cw(lencw)
    real(rp),     intent(inout) :: ru(lenru), rw(lenrw)

    integer(ip),  intent(out)   :: nindJ(lenJ), nlocJ(n+1), &
                                   lindJ(lenJ), llocJ(n+1)
    real(rp),     intent(out)   :: lJcol(lenJ)
    procedure(ictrode)          :: odecon
    external                    :: usrode

    !===========================================================================
    ! Hermite-Simpson discretization.
    ! Sets up the SNOPT structures for the ODE constraint derivatives.
    ! Differentiates linear and nonlinear constraints.
    !
    ! 26 Dec 2008: First version of odeSetupHS for v4 of snctrl.
    ! 27 Jan 2009: Sparse structures and phases implemented.
    !===========================================================================
    integer(ip) :: Status, nVar, colind, ncount, lcount, nrwind, lrwind, &
                   yrows(nY,nPhs), nNodes, ndS, ndE, nJ, dFlag(nY), &
                   p, i, j, k, q, il, in
    real(rp)    :: phstep, hstep, dvar(n-nP), pvar(nP)
    integer(ip), allocatable :: usrrow(:), usrcol(:)
    real(rp),    allocatable :: F(:,:), usrval(:,:)


    Status = 0
    dvar   = 0.0
    pvar   = 0.0
    nVar   = nY+nU

    ! Set up row indices for nonlinear/linear state eqns.
    do p = 1, nPhs
       in = 0
       il = 0
       do i = 1, nY
          if ( ytype(i,p) == 0 ) then
             in = in + 1
             yrows(i,p) = in
          else
             il = il + 1
             yrows(i,p) = il
          end if
       end do
    end do

    ! Initialize variables
    nindJ = 0
    nlocJ = 0

    lJcol = 0
    lindJ = 0
    llocJ = 0

    ncount = 1
    lcount = 1

    colind = 0
    nrwind = 0
    lrwind = nnCon

    allocate ( usrcol(nY+nU+nP+1) )
    do p = 1, nPhs

       nJ  = neJ(p)
       ndS = ndPtr(p)
       ndE = ndPtr(p+1)
       nNodes = ndE - ndS + 1

       ! Get user-defined Jacobian elements
       allocate( F(nY,ndS:ndE) )
       allocate( usrrow(nJ), usrval(nJ,ndS:ndE) )
       call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                    F, usrrow, usrval, usrcol, nJ, dvar, pvar, 0_ip, 1_ip, &
                    cu, lencu, iu, leniu, ru, lenru, &
                    cw, lencw, iw, leniw, rw, lenrw )

       ! At first node of the phase:
       j = ndS
       hstep = step( (j+1)/2 )

       nrwind = nrwind + nnC(p)
       lrwind = lrwind + nlC(p)
       if ( p /= 1 ) &
            lrwind = lrwind + nY

       do k = 1, nY+nU
          nlocJ(k+colind) = ncount
          llocJ(k+colind) = lcount

          ! Phase continuity
          if ( p > 1 .and. k <= nY ) then
             lindJ(lcount) = k + lrwind - nY
             lJcol(lcount) = 1d+0
             lcount = lcount + 1
          end if

          ! ODE Constraints 1 and 2
          do q = usrcol(k), usrcol(k+1)-1
             i = usrrow(q)

             if ( ytype(i,p) == 0 ) then
                nindJ(ncount) = yrows(i,p) + nrwind
                ncount = ncount + 1

                nindJ(ncount) = yrows(i,p) + nrwind + nnY(p)
                ncount = ncount + 1
             else
                lindJ(lcount) = yrows(i,p) + lrwind
                lJcol(lcount) = (hstep/6.0d+0)*usrval(q,j)
                if (k == i) &
                     lJcol(lcount) = lJcol(lcount) + 1d+0
                lcount = lcount + 1

                lindJ(lcount) = yrows(i,p) + lrwind + nlY(p)
                lJcol(lcount) = (hstep/8.0d+0)*usrval(q,j)
                if (k == i) &
                     lJcol(lcount) = lJcol(lcount) + 0.5d+0
                lcount = lcount + 1
             end if
          end do

          if ( k <= nY ) then
             if ( dFlag(k) == 0 ) then
                if ( ytype(k,p) == 0 ) then
                   nindJ(ncount) = yrows(k,p) + nrwind
                   ncount = ncount + 1

                   nindJ(ncount) = yrows(k,p) + nrwind + nnY(p)
                   ncount = ncount + 1
                else
                   lindJ(lcount) = yrows(k,p) + lrwind
                   lJcol(lcount) = 1d+0
                   lcount = lcount + 1

                   lindJ(lcount) = yrows(k,p) + lrwind + nlY(p)
                   lJcol(lcount) = 0.5d+0
                   lcount = lcount + 1
                end if
             end if
          else ! Linear Control
             lindJ(lcount) = (k-nY) + lrwind + 2*nlY(p)
             lJcol(lcount) = 0.5d+0
             lcount = lcount + 1
          end if
       end do

       ! Midpoint
       colind = colind + nVar
       do k = 1, nY+nU
          nlocJ(k+colind) = ncount
          llocJ(k+colind) = lcount

          ! ODE Constraint 1 and 2
          do q = usrcol(k), usrcol(k+1)-1
             i = usrrow(q)

             if ( ytype(i,p) == 0 ) then
                nindJ(ncount) = yrows(i,p) + nrwind
                ncount = ncount + 1

                nindJ(ncount) = yrows(i,p) + nrwind + nnY(p)
                ncount = ncount + 1
             else
                lindJ(lcount) = yrows(i,p) + lrwind
                lJcol(lcount) = (2*hstep/3.0d+0)*usrval(q,j+1)
                lcount = lcount + 1

                lindJ(lcount) = yrows(i,p) + lrwind + nlY(p)
                lJcol(lcount) = 0.0
                if (k == i) &
                     lJcol(lcount) = -1d+0
                lcount = lcount + 1
             end if
          end do

          if ( k <= nY ) then
             if ( dFlag(k) == 0 ) then
                if ( ytype(k,p) == 0 ) then
                   nindJ(ncount) = yrows(k,p) + nrwind + nnY(p)
                   ncount = ncount + 1
                else
                   lindJ(lcount) = yrows(k,p) + lrwind + nlY(p)
                   lJcol(lcount) = -1d+0
                   lcount = lcount + 1
                end if
             end if
          else ! Linear Control
             lindJ(lcount) = (k-nY) + lrwind + 2*nlY(p)
             lJcol(lcount) = -1d+0
             lcount = lcount + 1
          end if
       end do


       ! At inner nodes:
       do j = ndS+2, ndE-2, 2
          colind = colind + nVar
          phstep = hstep
          hstep  = step( (j+1)/2 )

          do k = 1, nY+nU
             nlocJ(k+colind) = ncount
             llocJ(k+colind) = lcount

             ! ODE Constraint 1 and 2
             do q = usrcol(k), usrcol(k+1)-1
                i = usrrow(q)

                if ( ytype(i,p) == 0 ) then
                   nindJ(ncount) = yrows(i,p) + nrwind
                   ncount = ncount + 1

                   nindJ(ncount) = yrows(i,p) + nrwind + nnY(p)
                   ncount = ncount + 1


                   nindJ(ncount) = yrows(i,p) + nrwind + 2*nnY(p) + nnC(p)
                   ncount = ncount + 1

                   nindJ(ncount) = yrows(i,p) + nrwind + 3*nnY(p) + nnC(p)
                   ncount = ncount + 1
                else
                   lindJ(lcount) = yrows(i,p) + lrwind
                   lJcol(lcount) = (phstep/6.0d+0)*usrval(q,j)
                   if (k == i) &
                        lJcol(lcount) = lJcol(lcount) - 1d+0
                   lcount = lcount + 1

                   lindJ(lcount) = yrows(i,p) + lrwind + nlY(p)
                   lJcol(lcount) = -(phstep/8.0d+0)*usrval(q,j)
                   if (k == i) &
                        lJcol(lcount) = lJcol(lcount) + 0.5d+0
                   lcount = lcount + 1


                   lindJ(lcount) = yrows(i,p) + lrwind + nU + 2*nlY(p) + nlC(p)
                   lJcol(lcount) = (hstep/6.0d+0)*usrval(q,j)
                   if (k == i) &
                        lJcol(lcount) = lJcol(lcount) + 1d+0
                   lcount = lcount + 1

                   lindJ(lcount) = yrows(i,p) + lrwind + nU + 3*nlY(p) + nlC(p)
                   lJcol(lcount) = (hstep/8.0d+0)*usrval(q,j)
                   if (k == i) &
                        lJcol(lcount) = lJcol(lcount) + 0.5d+0
                   lcount = lcount + 1
                end if
             end do

             if ( k <= nY ) then
                if ( dFlag(k) == 0) then
                   if ( ytype(k,p) == 0 ) then
                      nindJ(ncount) = yrows(k,p) + nrwind
                      ncount = ncount + 1

                      nindJ(ncount) = yrows(k,p) + nrwind + nnY(p)
                      ncount = ncount + 1


                      nindJ(ncount) = yrows(k,p) + nrwind + 2*nnY(p) + nnC(p)
                      ncount = ncount + 1

                      nindJ(ncount) = yrows(k,p) + nrwind + 3*nnY(p) + nnC(p)
                      ncount = ncount + 1
                   else
                      lindJ(lcount) = yrows(k,p) + lrwind
                      lJcol(lcount) = -1d+0
                      lcount = lcount + 1

                      lindJ(lcount) = yrows(k,p) + lrwind + nlY(p)
                      lJcol(lcount) = 0.5d+0
                      lcount = lcount + 1


                      lindJ(lcount) = yrows(k,p) + lrwind + nU + 2*nlY(p) + nlC(p)
                      lJcol(lcount) = 1d+0
                      lcount = lcount + 1

                      lindJ(lcount) = yrows(i,p) + lrwind + nU + 3*nlY(p) + nlC(p)
                      lJcol(lcount) = 0.5d+0
                      lcount = lcount + 1
                   end if
                end if
             else ! Linear Control
                lindJ(lcount) = (k-nY) + lrwind + 2*nlY(p)
                lJcol(lcount) = 0.5d+0
                lcount = lcount + 1

                lindJ(lcount) = (k-nY) + lrwind + nU + 4*nlY(p) + nlC(p)
                lJcol(lcount) = 0.5d+0
                lcount = lcount + 1
             end if
          end do

          nrwind = nrwind + 2*nnY(p) + nnC(p)
          lrwind = lrwind + nU + 2*nlY(p) + nlC(p)

          ! Midpoint
          colind = colind + nVar
          do k = 1, nY+nU
             nlocJ(k+colind) = ncount
             llocJ(k+colind) = lcount

             ! ODE Constraint 1 and 2
             do q = usrcol(k), usrcol(k+1)-1
                i = usrrow(q)

                if ( ytype(i,p) == 0 ) then
                   nindJ(ncount) = yrows(i,p) + nrwind
                   ncount = ncount + 1

                   nindJ(ncount) = yrows(i,p) + nrwind + nnY(p)
                   ncount = ncount + 1
                else
                   lindJ(lcount) = yrows(i,p) + lrwind
                   lJcol(lcount) = (2*hstep/3.0d+0)*usrval(q,j+1)
                   lcount = lcount + 1

                   lindJ(lcount) = yrows(i,p) + lrwind + nlY(p)
                   lJcol(lcount) = 0.0
                   if (k == i) &
                        lJcol(lcount) = -1d+0
                   lcount = lcount + 1
                end if
             end do

             if ( k <= nY ) then
                if ( dFlag(k) == 0 ) then
                   if ( ytype(k,p) == 0 ) then
                      nindJ(ncount) = yrows(k,p) + nrwind + nnY(p)
                      ncount = ncount + 1
                   else
                      lindJ(lcount) = yrows(k,p) + lrwind + nlY(p)
                      lJcol(lcount) = -1d+0
                      lcount = lcount + 1
                   end if
                end if
             else ! Linear Control
                lindJ(lcount) = (k-nY) + lrwind + 2*nlY(p)
                lJcol(lcount) = -1d+0
                lcount = lcount + 1
             end if
          end do
       end do

       ! At final node of phase:
       j = ndE
       phstep = hstep
       colind = colind + nVar
       do k = 1, nY+nU
          nlocJ(k+colind) = ncount
          llocJ(k+colind) = lcount

          ! ODE Constraints 1 and 2
          do q = usrcol(k), usrcol(k+1)-1
             i = usrrow(q)

             if ( ytype(i,p) == 0 ) then
                nindJ(ncount) = yrows(i,p) + nrwind
                ncount = ncount + 1

                nindJ(ncount) = yrows(i,p) + nrwind + nnY(p)
                ncount = ncount + 1
             else
                lindJ(lcount) = yrows(i,p) + lrwind
                lJcol(lcount) = (phstep/6d+0)*usrval(q,j)
                if (i == k) &
                     lJcol(lcount) = lJcol(lcount) - 1d+0
                lcount = lcount + 1

                lindJ(lcount) = yrows(i,p) + lrwind + nlY(p)
                lJcol(lcount) = -(phstep/8d+0)*usrval(q,j)
                if (i == k) &
                     lJcol(lcount) = lJcol(lcount) + 0.5d+0
                lcount = lcount + 1
             end if
          end do

          if ( k <= nY ) then
             if ( dFlag(k) == 0 ) then
                if ( ytype(k,p) == 0 ) then
                   nindJ(ncount) = yrows(k,p) + nrwind
                   ncount = ncount + 1

                   nindJ(ncount) = yrows(k,p) + nrwind + nnY(p)
                   ncount = ncount + 1
                else
                   lindJ(lcount) = yrows(k,p) + lrwind
                   lJcol(lcount) = -1d+0
                   lcount = lcount + 1

                   lindJ(lcount) = yrows(k,p) + lrwind + nlY(p)
                   lJcol(lcount) = 0.5d+0
                   lcount = lcount + 1
                end if
             end if

             ! Equality constraints for continuity at phases
             if ( p < nPhs ) then
                lindJ(lcount) = k + lrwind + 2*nlY(p) + nU
                lJcol(lcount) = -1d+0
                lcount = lcount + 1
             end if

          else ! Linear Control
             lindJ(lcount) = (k-nY) + lrwind + 2*nlY(p)
             lJcol(lcount) = 0.5d+0
             lcount = lcount + 1
          end if

          ! Objective
          if (k == objL(p) .or. k+nY == objL(p)) then
             lindJ(lcount) = m
             lJcol(lcount) = 1.0d+0
             lcount = lcount + 1
          end if
       end do

       colind = colind + nVar
       nrwind = nrwind + 2*nnY(p) + nnC(p)
       lrwind = lrwind + nU + 2*nlY(p) + nlC(p)
       deallocate ( F, usrrow, usrval )
    end do

    ! Parameters
    nrwind = 0
    lrwind = nnCon
    do k = 1, nP
       nlocJ(k+colind) = ncount
       llocJ(k+colind) = lcount

       do p = 1, nPhs

          nJ  = neJ(p)
          ndS = ndPtr(p)
          ndE = ndPtr(p+1)
          nNodes = ndE - ndS + 1

          ! Get user-defined state equations
          allocate ( F(nY,ndS:ndE) )
          allocate ( usrrow(nJ), usrval(nJ,ndS:ndE) )
          call odecon( Status, usrode, p, nPhs, nY, nU, nP, nNodes, dFlag, &
                       F, usrrow, usrval, usrcol, nJ, dvar, pvar, 0_ip, 1_ip, &
                       cu, lencu, iu, leniu, ru, lenru, &
                       cw, lencw, iw, leniw, rw, lenrw )

          nrwind = nrwind + nnC(p)
          lrwind = lrwind + nlC(p)

          do j = ndS, ndE-2, 2
             hstep = step( (j+1)/2 )

             ! ODE Constraint 1 and 2
             do q = usrcol(k+nVar), usrcol(k+nVar+1)-1
                i = usrrow(q)

                if ( ytype(i,p) == 0 ) then
                   nindJ(ncount) = yrows(i,p) + nrwind
                   ncount = ncount + 1

                   nindJ(ncount) = yrows(i,p) + nrwind + nnY(p)
                   ncount = ncount + 1
                else
                   lindJ(lcount) = yrows(i,p) + lrwind
                   lJcol(lcount) = hstep*usrval(q,j)
                   lcount = lcount + 1

                   lindJ(lcount) = yrows(i,p) + lrwind + nlY(p)
                   lJcol(lcount) = 0.0
                   lcount = lcount + 1
                end if
             end do
             nrwind = nrwind + 2*nnY(p) + nnC(p)
             lrwind = lrwind + 2*nlY(p) + nlC(p) + nU
          end do

          lrwind = lrwind + nY
          deallocate ( F, usrrow, usrval )

          ! Objective
          if ( k+nY+nU == objL(p) ) then
             lindJ(lcount) = m
             lJcol(lcount) = 1.0
             lcount = lcount + 1
          end if
       end do
    end do

    deallocate ( usrcol )

    nlocJ(n+1) = ncount
    llocJ(n+1) = lcount

  end subroutine odeSetupHS

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine algSetupHS( nPhs, nY, nU, nP, nC, ctype, nnY, nlY, &
                         nnC, nlC, neG, nnCon, ndPtr, n, lenJ, &
                         nindJ, nlocJ, lJcol, lindJ, llocJ, &
                         algcon, usralg, &
                         cu, lencu, iu, leniu, ru, lenru ,&
                         cw, lencw, iw, leniw, rw, lenrw )
    integer(ip),  intent(in)     :: nPhs, nY, nU, nP, nC, ctype(nC,nPhs), &
                                    nnY(nPhs), nlY(nPhs), nnC(nPhs), nlC(nPhs), &
                                    neG(nPhs), nnCon, ndPtr(nPhs+1), n, lenJ, &
                                    lencu, leniu, lenru, lencw, leniw, lenrw

    integer(ip),  intent(inout) :: iu(leniu), iw(leniw)
    character(8), intent(inout) :: cu(lencu), cw(lencw)
    real(rp),     intent(inout) :: ru(lenru), rw(lenrw)

    integer(ip),  intent(out)    :: nindJ(lenJ), nlocJ(n+1), &
                                    lindJ(lenJ), llocJ(n+1)
    real(rp),     intent(out)    :: lJcol(lenJ)
    procedure(ictralg)           :: algcon
    external                     :: usralg

    !===========================================================================
    ! Sets up the SNOPT structures for the algebraic constraint
    ! derivatives.
    ! Differentiates linear and nonlinear constraints.
    !
    !         y,u(t0) | . | y,u(t1) | .... | . | y,u(tf) | p
    ! C(t0)      x                                         x
    ! C(t1)                   x                            x
    !   ...
    ! C(tf)                                         x      x
    !
    !
    ! 29 Dec 2008: First version of algSetupHS for snctrl v4.
    ! 19 Jan 2009: Sparse structs implemented.
    !===========================================================================
    integer(ip) :: Status, nVar, colind, ncount, lcount, nrwind, lrwind, &
                   crows(nC,nPhs), nNodes, ndS, ndE, nG, &
                   p, i, j, k, q, in, il
    real(rp)    :: dvar(n-nP), pvar(nP)
    integer(ip), allocatable :: usrrow(:), usrcol(:)
    real(rp),    allocatable :: C(:,:), usrval(:,:)


    Status = 0
    dvar   = 0.0
    pvar   = 0.0
    nVar   = nY+nU

    ! Set up row indices for nonlinear/linear state eqns.
    do p = 1, nPhs
       in = 0
       il = 0
       do i = 1, nC
          if ( ctype(i,p) == 0 ) then
             in = in + 1
             crows(i,p) = in
          else
             il = il + 1
             crows(i,p) = il
          end if
       end do
    end do

    ! Initialize variables
    nindJ = 0
    nlocJ = 0

    lJcol = 0
    lindJ = 0
    llocJ = 0

    nVar   = nY+nU
    ncount = 1
    lcount = 1

    colind = 0
    nrwind = 0
    lrwind = nnCon

    ! Go...
    allocate ( usrcol(nY+nU+nP+1) )

    do p = 1, nPhs

       nG  = neG(p)
       ndS = ndPtr(p)
       ndE = ndPtr(p+1)
       nNodes = ndE - ndS + 1

       ! Get user-defined Jacobian elements
       allocate( C(nC,ndS:ndE) )
       allocate( usrrow(nG), usrval(nG,ndS:ndE) )
       call algcon( Status, usralg, p, nPhs, nC, nY, nU, nP, nNodes, &
                    C, usrrow, usrval, usrcol, nG, dvar, pvar, 0_ip, 1_ip, &
                    cu, lencu, iu, leniu, ru, lenru, &
                    cw, lencw, iw, leniw, rw, lenrw )

       do j = ndS, ndE, 2
          do k = 1, nY+nU

             nlocJ(k+colind) = ncount
             llocJ(k+colind) = lcount

             do q = usrcol(k), usrcol(k+1)-1
                i = usrrow(q)

                if ( ctype(i,p) == 0 ) then
                   nindJ(ncount) = crows(i,p) + nrwind
                   ncount = ncount + 1
                else
                   lindJ(lcount) = crows(i,p) + lrwind
                   lJcol(lcount) = usrval(q,j)
                   lcount = lcount + 1
                end if
             end do
          end do
          colind = colind + nVar

          if ( j /= ndE ) then
             do k = 1, nY+nU
                nlocJ(k+colind) = ncount
                llocJ(k+colind) = lcount
             end do
             colind = colind + nVar
          end if

          nrwind = nrwind + nnC(p) + 2*nnY(p)
          lrwind = lrwind + nlC(p) + 2*nlY(p) + nU
       end do

       lrwind = lrwind + nY
       deallocate ( C, usrrow, usrval )
    end do

    ! Parameters
    nrwind = 0
    lrwind = nnCon

    do k = 1, nP
       nlocJ(k+colind) = ncount
       llocJ(k+colind) = lcount

       do p = 1, nPhs

          nG  = neG(p)
          ndS = ndPtr(p)
          ndE = ndPtr(p+1)
          nNodes = ndE - ndS + 1

          ! Get user-defined Jacobian elements
          allocate( C(nC,ndS:ndE) )
          allocate( usrrow(nG), usrval(nG,ndS:ndE) )
          call algcon( Status, usralg, p, nPhs, nC, nY, nU, nP, nNodes, &
                       C, usrrow, usrval, usrcol, nG, dvar, pvar, 0_ip, 1_ip, &
                       cu, lencu, iu, leniu, ru, lenru, &
                       cw, lencw, iw, leniw, rw, lenrw )

          do j = ndS, ndE, 2
             do q = usrcol(k+nVar), usrcol(k+nVar+1)-1
                i = usrrow(q)

                if (ctype(i,p) == 0) then
                   nindJ(ncount) = crows(i,p) + nrwind
                   ncount = ncount + 1
                else
                   lindJ(lcount) = crows(i,p) + lrwind
                   lJcol(lcount) = usrval(q,j)
                   lcount = lcount + 1
                end if
             end do
             nrwind = nrwind + nnC(p) + 2*nnY(p)
             lrwind = lrwind + nlC(p) + 2*nlY(p) + nU
          end do

          lrwind = lrwind + nY
          deallocate ( C, usrrow, usrval )
       end do
    end do

    deallocate ( usrcol )

    nlocJ(n+1) = ncount
    llocJ(n+1) = lcount

  end subroutine algSetupHS

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module ct28derv
