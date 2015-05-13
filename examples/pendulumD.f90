!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File  pendulumD.f90
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

program pendulumD
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
  parameter ( lenrw = 50000, leniw = 50000, lencw = 500, &
              lenru = 1,     leniu = 1,     lencu = 1)
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
  lfile = 'pendulumD.spc'
  open( iSpecs, file=lfile, status='OLD',     err=800 )
  lfile = 'pendulumD.out'
  open( iPrint, file=lfile, status='UNKNOWN', err=800 )


  ! Set options to default values.
  call sncInit( iPrint, iSumm, prob, cw, lencw, iw, leniw, rw, lenrw )


  ! Read a Specs file.
  call sncSpec( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

  if (INFO /= 101  .AND.  INFO /= 107) go to 910


  ! Setup the problem.
  prob%probname(1:8) = 'pendD'
  prob%nY   = 4     ! number of states
  prob%nU   = 2     ! number of controls
  prob%nP   = 1     ! number of parameters
  prob%nC   = 2     ! number of algebraic constraints
  prob%nPhs = 1     ! number of phases

  ! Objective
  allocate ( prob%objL(1) )
  prob%objL(1) = prob%nY+2

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
  call snctrlD ( Start, prob, odecon, algcon, varbds, &
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

end program pendulumD

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine varbds ( curPhs, nPhs, nY, nU, nP, nC, nNodes, lbds, ubds, &
                    x, plbds, pubds, p, clbds, cubds )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in)  :: curPhs, nPhs, nY, nU, nP, nC, nNodes
  real(rp),    intent(out) :: lbds(nY+nU,nNodes), ubds(nY+nU,nNodes), &
                              x(nY+nU,nNodes), plbds(nP), pubds(nP), p(nP), &
                              clbds(nC), cubds(nC)

  !=============================================================================
  integer(ip) :: j
  real(rp)    :: xf, yf, l, bplus, bminus
  parameter ( bplus  = 1.0d+20, bminus = -bplus )

  ! Some constants.
  xf = -.231625
  yf = -.443101
  l  = 0.5

  ! Default set everything to free bounds and initial point is 0.
  lbds = bminus
  ubds = bplus
  x    = 0.0
  p    = 0.0

  ! Set up initial point (x)
  do j = 1, nNodes
     x(1,j) = .4
     x(2,j) = -.3
     x(3,j) = 0.0
     x(4,j) = 0.0

     x(5,j) = -5.0
     x(6,j) = 0.0
  end do
  x(6,nNodes) = 0.5*((x(1,nNodes)-xf)**2 + (x(2,nNodes)-yf)**2)

  ! At t0
  lbds(1,1) = .40
  ubds(1,1) = .40

  lbds(2,1) = -.30
  ubds(2,1) = -.30

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

subroutine odecon ( snStat, curPhs, nPhs, nY, nU, nP, nNodes, F, J, &
                    dvar, pvar, needF, needJ, cu, lencu, iu, leniu, ru, lenru )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in) :: snStat, curPhs, nPhs, nY, nU, nP, nNodes, &
                             needF, needJ, lencu, leniu, lenru
  real(rp),    intent(in) :: dvar(nY+nU,nNodes), pvar(nP)

  integer(ip), intent(inout) :: iu(leniu)
  real(rp),    intent(inout) :: ru(lenru)
  character(8),intent(inout) :: cu(lencu)

  real(rp),    intent(out) :: F(nY,nNodes), J(nY,nY+nU+nP,nNodes)

  !=============================================================================
  integer(ip) :: jt
  real(rp)    :: m, l

  m  = .3
  l  = .5

  if (needF > 0) then
     do jt = 1, nNodes
        F(1,jt) = dvar(3,jt)
        F(2,jt) = dvar(4,jt)
        F(3,jt) = dvar(5,jt)*dvar(1,jt) / m
        F(4,jt) = dvar(5,jt)*dvar(2,jt) / m - pvar(1)
     end do
  end if

  if (needJ > 0) then
     do jt = 1, nNodes
        J(1,1,jt) = 0.0
        J(1,2,jt) = 0.0
        J(1,3,jt) = 1.0
        J(1,4,jt) = 0.0
        J(1,5,jt) = 0.0
        J(1,6,jt) = 0.0
        J(1,7,jt) = 0.0

        J(2,1,jt) = 0.0
        J(2,2,jt) = 0.0
        J(2,3,jt) = 0.0
        J(2,4,jt) = 1.0
        J(2,5,jt) = 0.0
        J(2,6,jt) = 0.0
        J(2,7,jt) = 0.0

        J(3,1,jt) = dvar(5,jt) / m
        J(3,2,jt) = 0.0
        J(3,3,jt) = 0.0
        J(3,4,jt) = 0.0
        J(3,5,jt) = dvar(1,jt) / m
        J(3,6,jt) = 0.0
        J(3,7,jt) = 0.0

        J(4,1,jt) = 0.0
        J(4,2,jt) = dvar(5,jt) / m
        J(4,3,jt) = 0.0
        J(4,4,jt) = 0.0
        J(4,5,jt) = dvar(2,jt) / m
        J(4,6,jt) = 0.0
        J(4,7,jt) = -1.0
     end do
  end if

end subroutine odecon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine algcon ( snStat, curPhs, nPhs, nC, nY, nU, nP, nNodes, C, G, &
                    dvar, pvar, needC, needG, cu, lencu, iu, leniu, ru, lenru )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in) :: snStat, curPhs, nPhs, nC, nY, nU, nP, nNodes, &
                             needC, needG, lencu, leniu, lenru
  real(rp),    intent(in) :: dvar(nY+nU,nNodes), pvar(nP)

  integer(ip), intent(inout) :: iu(leniu)
  real(rp),    intent(inout) :: ru(lenru)
  character(8),intent(inout) :: cu(lencu)

  real(rp),    intent(out) :: C(nC,nNodes), G(nC,nY+nU+nP,nNodes)

  !=============================================================================
  integer(ip) :: jt
  real(rp)    :: m, l, xf, yf

  m  = .3
  l  = .5
  xf = -.231625
  yf = -.443101

  if ( needC > 0 ) then
     do jt = 1, nNodes
        C(1,jt) = dvar(3,jt)**2 + dvar(4,jt)**2 + &
                  dvar(5,jt)*(l**2) / m - dvar(2,jt)*pvar(1)
!        C(2,jt) = 0d+0
        C(2,jt) = dvar(6,jt) - .5*( (dvar(1,jt) - xf)**2 + (dvar(2,jt) - yf)**2 )
     end do
!     jt = nNodes
!     C(2,jt) = dvar(6,jt) - .5*( (dvar(1,jt) - xf)**2 + (dvar(2,jt) - yf)**2 )

  end if

  if ( needG > 0 ) then
     do jt = 1, nNodes
        G(1,1,jt) = 0.0
        G(1,2,jt) = -pvar(1)
        G(1,3,jt) = 2*dvar(3,jt)
        G(1,4,jt) = 2*dvar(4,jt)
        G(1,5,jt) = l**2/m
        G(1,6,jt) = 0.0
        G(1,7,jt) = -dvar(2,jt)

     G(2,1,jt) = -(dvar(1,jt) - xf)
     G(2,2,jt) = -(dvar(2,jt) - yf)
     G(2,3,jt) = 0.0
     G(2,4,jt) = 0.0
     G(2,5,jt) = 0.0
     G(2,6,jt) = 1.0
     G(2,7,jt) = 0.0
!!$        G(2,1,jt) = 0.0
!!$        G(2,2,jt) = 0.0
!!$        G(2,3,jt) = 0.0
!!$        G(2,4,jt) = 0.0
!!$        G(2,5,jt) = 0.0
!!$        G(2,6,jt) = 0.0
!!$        G(2,7,jt) = 0.0
     end do

     jt = nNodes

  end if

end subroutine algcon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
