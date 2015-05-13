!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File  pendulumS.f90
!
!  This is a parameter estimation problem.
!
!  The original pendulum problem:
!      min            (1/2) * ( (y1(tf) - xf)^2 + (y2(tf) - yf)^2 )
!      subject to:
!                     y1' = y3
!                     y2' = y4
!                     y3' = u1 * y1 / m
!                     y4' = u1 * y2 / m - p1
!                     y1(0) = 0.4
!                     y2(0) = -0.3
!                     y3(0) = 0
!                     y4(0) = 0
!
!                     y3^2 + y4^2 + u1 * l^2/m - y2*p1 = 0
!
!                     .01 <= p1 <= 100
!
!      where xf, yf, l and m are constants.
!      l and m the length and mass of the pendulum respectively.
!      For this problem we set l = 0.5, m = .3, xf = -.231625, and
!      yf = -.443109.
!
!      p1 is the parameter representing the gravitational constant.
!      y1 and y2 represent the coordiantes of the pendulum.
!      The constraint y1^2 + y2^2 = l^2 is satisfied implicitly.
!
!  The problem input to SNCTRL:
!      min            u2(2)
!      subject to:
!                     y1' = y3
!                     y2' = y4
!                     y3' = u1 * y1 / m
!                     y4' = u1 * y2 / m - p1
!                     y1(0) = 0.4
!                     y2(0) = -0.3
!                     y3(0) = 0
!                     y4(0) = 0
!
!                     y3^2 + y4^2 + u1 * l^2/m - y2*p1 = 0
!
!            0 = u2(tf) - (1/2) * ( (y1(tf) - xf)^2 + (y2(tf) - yf)^2 )
!
!                     .01 <= p1 <= 100
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program pendulumS
  use SNCTRL
  implicit none

  !=====================================================================

  ! SNOPT variables and workspaces
  character(20) :: lfile
  real(rp)      :: sInf
  integer(ip)   :: Start, INFO, iPrint, iSpecs, nOut, &
                   iSumm, mincw, miniw, minrw, nS, nInf, &
                   lencu, leniu, lenru, &
                   lencw, leniw, lenrw
  parameter ( lenrw = 700000, leniw = 80000, lencw = 500, &
              lenru = 1,      leniu = 1,     lencu = 1)
  integer(ip)  :: iu(leniu), iw(leniw)
  real(rp)     :: ru(lenru), rw(lenrw)
  character(8) :: cu(lencu), cw(lencw)

  external     :: odecon, algcon, varbds
  type(ctProb) :: prob

  !=====================================================================

  iSpecs =  4  ! equivalenced to .spc
  iPrint =  9  ! equivalenced to .out
  iSumm  =  6  ! summary file goes to standard output...
  nOut   =  6  ! ... as do messages from this program.

  ! Open the Specs and print files.
  lfile = 'pendulumS.spc'
  open( iSpecs, file=lfile, status='OLD',     err=800 )
  lfile = 'pendulumS.out'
  open( iPrint, file=lfile, status='UNKNOWN', err=800 )


  ! Set options to default values.
  call sncInit( iPrint, iSumm, prob, cw, lencw, iw, leniw, rw, lenrw )


  ! Read a Specs file.
  call sncSpec( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

  if (INFO /= 101  .AND.  INFO /= 107) go to 910


  ! Setup the problem.
  prob%probname = 'pendS'
  prob%nY   = 4     ! number of states
  prob%nU   = 2     ! number of controls
  prob%nP   = 1     ! number of parameters
  prob%nC   = 2     ! number of algebraic constraints
  prob%nPhs = 1     ! number of phases

  ! Derivatives
  allocate ( prob%neJ(1), prob%neG(1) )
  prob%neJ(1) = 7
  prob%neG(1) = 8

  ! Objective
  allocate ( prob%objL(1) )
  prob%objL = prob%nY+2

  ! Interval
  allocate ( prob%npInt(1), prob%phsPt(2) )
  prob%npInt(1) = 25

  prob%phsPt(1) = 0
  prob%phsPt(2) = 2.0d+0

  ! Linearity
  allocate ( prob%ytype(4,1), prob%ctype(2,1) )
  prob%ytype(1,1) = 1
  prob%ytype(2,1) = 1
  prob%ytype(3,1) = 0
  prob%ytype(4,1) = 0

  prob%ctype(1,1) = 0
  prob%ctype(2,1) = 0


  !---------------------------------------------------------------------
  ! Solve the problem.
  !---------------------------------------------------------------------
  Start = 0
  call snctrlS( Start, prob, odecon, algcon, varbds, &
                INFO, mincw, miniw, minrw, nS, nInf, sInf, &
                cu, lencu, iu, leniu, ru, lenru, &
                cw, lencw, iw, leniw, rw, lenrw )

  if (INFO == 82 .OR. INFO == 83 .OR. INFO == 84) then
     go to 910
  end if

  write(nOut, *) ' '
  write(nOut, *) 'snctrl finished.'
  write(nOut, *) 'INFO          =', INFO
  write(nOut, *) 'nInf          =', nInf
  write(nOut, *) 'sInf          =', sInf

  call sncEnd(prob)

  if (INFO >= 30) go to 900

  stop

  !---------------------------------------------------------------------
  ! Error exit.
  !---------------------------------------------------------------------
800 write(nOut, 4000) 'Error while opening file', lfile
  stop

900 write(nOut, *) ' '
  write(nOut, *) 'STOPPING because of error condition'
910 stop

1000 format(f10.8)
4000 format(/  a, 2x, a  )

end program pendulumS

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine varbds ( curPhs, nPhs, nY, nU, nP, nC, nNodes, lbds, ubds, x, &
                    plbds, pubds, p, clbds, cubds )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in)  :: curPhs, nPhs, nY, nU, nP, nC, nNodes
  real(rp),    intent(out) :: lbds(nY+nU,nNodes), ubds(nY+nU,nNodes), &
                              x(nY+nU,nNodes), plbds(nP), pubds(nP), p(nP), &
                              clbds(nC), cubds(nC)

  !=============================================================================
  integer(ip) :: j
  real(rp)    :: xf, yf, l, bplus, bminus
  parameter    ( bplus  = 1.0d+20, bminus = -bplus )

  ! Some constants.
  xf = -.231625d+0
  yf = -.443101d+0
  l  = .5d+0

  ! Default set everything to free bounds and initial point is 0.
  lbds = bminus
  ubds = bplus
  x    = 0.0
  p    = 0.0

  ! Set up initial point (x)
  do j = 1, nNodes
     x(1,j) = .4d+0
     x(2,j) = -.3d+0
     x(3,j) = 0.0
     x(4,j) = 0.0

     x(5,j) = -5.0d+0
     x(6,j) = 0.0
  end do
  x(6,nNodes) = 0.5d+0*((x(1,nNodes)-xf)**2 + (x(2,nNodes)-yf)**2)

  ! At t0
  lbds(1,1) = .4d+0
  ubds(1,1) = .4d+0

  lbds(2,1) = -.3d+0
  ubds(2,1) = -.3d+0

  lbds(3,1) = 0.0
  ubds(3,1) = 0.0

  lbds(4,1) = 0.0
  ubds(4,1) = 0.0

  ! Parameters
  plbds(1) = 1d-2
  pubds(1) = 1d+2
  p(1) = 20.0

  ! Algebraic constraint bds
  clbds(1) = 0.0
  cubds(1) = 0.0
  clbds(2) = 0.0
  cubds(2) = 0.0

end subroutine varbds

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine odecon ( snStat, curPhs, nPhs, nY, nU, nP, nNodes, F, &
                    Jrow, Jval, Jcol, lenJ, dvar, pvar, needF, needJ, &
                    cu, lencu, iu, leniu, ru, lenru )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in) :: snStat, curPhs, nPhs, nY, nU, nP, nNodes, lenJ, &
                             needF, needJ, lencu, leniu, lenru
  real(rp),    intent(in) :: dvar(nY+nU,nNodes), pvar(nP)

  integer(ip), intent(inout) :: iu(leniu)
  real(rp),    intent(inout) :: ru(lenru)
  character(8),intent(inout) :: cu(lencu)

  integer(ip), intent(out) :: Jrow(lenJ), Jcol(nY+nU+nP)
  real(rp),    intent(out) :: F(nY,nNodes), Jval(lenJ,nNodes)

  !=============================================================================
  integer(ip) :: jt, neJ
  real(rp)    :: m, l

  m  = .3d+0
  l  = .5d+0

  if ( needF > 0 ) then
     do jt = 1, nNodes
        F(1,jt) = dvar(3,jt)
        F(2,jt) = dvar(4,jt)
        F(3,jt) = dvar(5,jt)*dvar(1,jt) / m
        F(4,jt) = dvar(5,jt)*dvar(2,jt) / m - pvar(1)
     end do
  end if

  if ( needJ > 0 ) then
     do jt = 1, nNodes
        neJ = 1
        ! Column 1
        Jcol(1) = neJ

        Jrow(neJ) = 3
        Jval(neJ,jt) = dvar(5,jt) / m
        neJ = neJ + 1

        ! Column 2
        Jcol(2) = neJ

        Jrow(neJ) = 4
        Jval(neJ,jt) = dvar(5,jt) / m
        neJ = neJ + 1

        ! Column 3
        Jcol(3) = neJ

        Jrow(neJ) = 1
        Jval(neJ,jt) = 1.0
        neJ = neJ + 1

        ! Column 4
        Jcol(4) = neJ

        Jrow(neJ) = 2
        Jval(neJ,jt) = 1.0
        neJ = neJ + 1

        ! Column 5
        Jcol(5) = neJ

        Jrow(neJ) = 3
        Jval(neJ,jt) = dvar(1,jt) / m
        neJ = neJ + 1

        Jrow(neJ) = 4
        Jval(neJ,jt) = dvar(2,jt) / m
        neJ = neJ + 1

        ! Column 6
        Jcol(6) = neJ

        ! Column 7
        Jcol(7) = neJ

        Jrow(neJ) = 4
        Jval(neJ,jt) = -1.0
        neJ = neJ + 1

        ! Finish off ptrs
        Jcol(8) = neJ
        neJ = neJ - 1
     end do
  end if

end subroutine odecon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine algcon ( snStat, curPhs, nPhs, nC, nY, nU, nP, nNodes, C, &
                    Grow, Gval, Gcol, lenG, dvar, pvar, &
                    needC, needG, cu, lencu, iu, leniu, ru, lenru )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in) :: snStat, curPhs, nPhs, nC, nY, nU, nP, nNodes, &
                             lenG, needC, needG, lencu, leniu, lenru
  real(rp),    intent(in) :: dvar(nY+nU,nNodes), pvar(nP)

  integer(ip), intent(inout) :: iu(leniu)
  real(rp),    intent(inout) :: ru(lenru)
  character(8),intent(inout) :: cu(lencu)

  integer(ip), intent(out) :: Grow(lenG), Gcol(nY+nU+nP+1)
  real(rp),    intent(out) :: C(nC,nNodes), Gval(lenG,nNodes)

  !=============================================================================
  integer(ip) :: jt, neG
  real(rp)    :: m, l, xf, yf

  m  = .3d+0
  l  = .5d+0
  xf = -.231625d+0
  yf = -.443101d+0

  if ( needC > 0 ) then
     do jt = 1, nNodes
        C(1,jt) = dvar(3,jt)**2 + dvar(4,jt)**2 + &
                  dvar(5,jt)*(l**2) / m - dvar(2,jt)*pvar(1)
        C(2,jt) = dvar(6,jt) - 0.5*( (dvar(1,jt) - xf)**2 + (dvar(2,jt) - yf)**2 )
     end do
  end if

  if ( needG > 0 ) then
     do jt = 1, nNodes
        neG = 1

        ! Column 1
        Gcol(1) = neG

        Grow(neG) = 2
        Gval(neG,jt) = -(dvar(1,jt) - xf)
        neG = neG + 1

        ! Column 2
        Gcol(2) = neG

        Grow(neG) = 1
        Gval(neG,jt) = -pvar(1)
        neG = neG + 1

        Grow(neG) = 2
        Gval(neG,jt) = -(dvar(2,jt) - yf)
        neG = neG + 1

        ! Column 3
        Gcol(3) = neG

        Grow(neG) = 1
        Gval(neG,jt) = 2*dvar(3,jt)
        neG = neG + 1

        ! Column 4
        Gcol(4) = neG

        Grow(neG) = 1
        Gval(neG,jt) = 2*dvar(4,jt)
        neG = neG + 1

        ! Column 5
        Gcol(5) = neG

        Grow(neG) = 1
        Gval(neG,jt) = l**2/m
        neG = neG + 1

        ! Column 6
        Gcol(6) = neG

        Grow(neG) = 2
        Gval(neG,jt) = 1d+0
        neG = neG + 1

        ! Column 7
        Gcol(7) = neG

        Grow(neG) = 1
        Gval(neG,jt) = -dvar(2,jt)
        neG = neG + 1

        ! Finish off ptrs
        Gcol(8) = neG
        neG = neG - 1
     end do
  end if

end subroutine algcon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
