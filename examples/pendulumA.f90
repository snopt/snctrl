!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File  pendulumA.f90
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

program pendulumA
  use SNCTRL
  implicit none

  !=====================================================================

  ! SNOPT variables and workspaces
  character(20) :: lfile
  real(rp)      :: sInf
  integer(ip)   :: Start, INFO, iPrint, iSpecs, nOut, Errors, &
                   iSumm, mincw, miniw, minrw, nS, nInf, &
                   lencu, leniu, lenru, &
                   lencw, leniw, lenrw
  parameter ( lenrw = 50000, leniw = 50000, lencw = 500, &
              lenru = 1,     leniu = 1,     lencu = 1)
  integer(ip)  :: iu(leniu), iw(leniw)
  real(rp)     :: ru(lenru), rw(lenrw)
  character(8) :: cu(lencu), cw(lencw)

  integer(ip)  :: ind, j
  real(rp)     :: xf, yf

  external     :: odecon, algcon, varbds
  type(ctProb) :: prob

  !=====================================================================

  iSpecs =  4  ! equivalenced to .spc
  iPrint =  9  ! equivalenced to .out
  iSumm  =  6  ! summary file goes to standard output...
  nOut   =  6  ! ... as do messages from this program.

  ! Open the Specs and print files.
  lfile = 'pendulumA.spc'
  open( iSpecs, file=lfile, status='OLD',     err=800 )
  lfile = 'pendulumA.out'
  open( iPrint, file=lfile, status='UNKNOWN', err=800 )


  ! Set options to default values.
  call sncInit( iPrint, iSumm, prob, cw, lencw, iw, leniw, rw, lenrw )


  ! Read a Specs file.
  call sncSpec( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

  if (INFO /= 101  .AND.  INFO /= 107) go to 910


  ! Setup the problem.
  prob%probname(1:8) = 'pendulmA'
  prob%nY   = 4     ! number of states
  prob%nU   = 2     ! number of controls
  prob%nP   = 1     ! number of parameters
  prob%nC   = 2     ! number of algebraic constraints
  prob%nPhs = 1     ! number of phases

  ! Objective
  allocate ( prob%objL(1) )
  prob%objL = prob%nY+2

  ! Interval
  allocate ( prob%npInt(1), prob%phsPt(2) )
  prob%npInt(1) = 25

  prob%phsPt(1) = 0
  prob%phsPt(2) = 2.0d+0

  ! Linearity
  allocate ( prob%ytype(4,1) )
  prob%ytype(1,1) = 1
  prob%ytype(2,1) = 1
  prob%ytype(3,1) = 0
  prob%ytype(4,1) = 0

  allocate ( prob%ctype(2,1) )
  prob%ctype(1,1) = 0
  prob%ctype(2,1) = 0


  ! Initial point
  ! (TR problem -> number of nodes = number of intervals + 1)
  call sncSet( 'Discretization TR', iPrint, iSumm, Errors, &
               cw, lencw, iw, leniw, rw, lenrw )

  xf = -.231625d+0
  yf = -.443101d+0
  allocate ( prob%x ((4+2)*26 + 1))
  allocate ( prob%hs((4+2)*26 + 1))
  prob%x  = 0.0
  prob%hs = 0

  ind = 0
  do j = 0, 25
     ind = j*(4+2)
     prob%x(1+ind) = .4
     prob%x(2+ind) = -.3
     prob%x(3+ind) = 0.0
     prob%x(4+ind) = 0.0

     prob%x(5+ind) = -5.0
     prob%x(6+ind) = 0.0
  end do
  ind = 25*6
  prob%x(6+ind) = 0.5*((prob%x(1+ind)-xf)**2 + (prob%x(2+ind)-yf)**2)
  prob%x(7+ind) = 20.0

  !---------------------------------------------------------------------
  ! Solve the problem.
  !---------------------------------------------------------------------
  Start = 0
  call snctrlA ( Start, prob, odecon, algcon, varbds, &
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

end program pendulumA

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine varbds ( curPhs, nPhs, nY, nU, nP, nC, y0low, y0upp, yflow, yfupp, &
                    ylow, yupp, ulow, uupp, plow, pupp, p, clow, cupp )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in)  :: nPhs, curPhs, nY, nU, nP, nC
  real(rp),    intent(out) :: y0low(nY), y0upp(nY), yflow(nY), yfupp(nY), &
                             ylow(nY),  yupp(nY),  ulow(nU),  uupp(nU), &
                             plow(nP), pupp(nP), p(nP), clow(nC), cupp(nC)

  !=============================================================================
  real(rp) :: bplus, bminus
  parameter ( bplus  = 1.0d+20, bminus = -bplus )

  ! Default (free bounds)
  y0low = bminus
  y0upp = bplus

  yflow = bminus
  yfupp = bplus

  ylow = bminus
  yupp = bplus

  ! At t0
  y0low(1) = .4
  y0upp(1) = .4

  y0low(2) = -.3
  y0upp(2) = -.3

  y0low(3) = 0.0
  y0upp(3) = 0.0

  y0low(4) = 0.0
  y0upp(4) = 0.0

  ! u1 free
  ulow(1) = bminus
  uupp(1) = bplus

  ulow(2) = bminus
  uupp(2) = bplus

  ! Parameter
  plow(1) = 1d-2
  pupp(1) = 1d+2
  p(1) = 20.0

  ! Algebraic constraints
  clow(1) = 0.0
  cupp(1) = 0.0

  clow(2) = 0.0
  cupp(2) = 0.0

end subroutine varbds

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine odecon ( snStat, curPhs, nPhs, nY, nU, nP, F, J, y, u, p, &
                    needF, needJ, cu, lencu, iu, leniu, ru, lenru )
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in)    :: snStat, curPhs, nPhs, nY, nU, nP, &
                               needF, needJ, lencu, leniu, lenru
  real(rp),    intent(in)    :: y(nY), u(nU), p(nP)

  integer(ip), intent(inout) :: iu(leniu)
  real(rp),    intent(inout) :: ru(lenru)
  character(8),intent(inout) :: cu(lencu)

  real(rp),    intent(out)   :: F(nY), J(nY, nY+nU+nP)

  !=============================================================================
  real(rp)    :: m, l

  m  = .3
  l  = .5

  if (needF > 0) then
     F(1) = y(3)
     F(2) = y(4)
     F(3) = u(1)*y(1) / m
     F(4) = u(1)*y(2) / m - p(1)
  end if

  if (needJ > 0) then
     J(1,1) = 0.0
     J(1,2) = 0.0
     J(1,3) = 1.0
     J(1,4) = 0.0
     J(1,5) = 0.0
     J(1,6) = 0.0
     J(1,7) = 0.0

     J(2,1) = 0.0
     J(2,2) = 0.0
     J(2,3) = 0.0
     J(2,4) = 1.0
     J(2,5) = 0.0
     J(2,6) = 0.0
     J(2,7) = 0.0

     J(3,1) = u(1) / m
     J(3,2) = 0.0
     J(3,3) = 0.0
     J(3,4) = 0.0
     J(3,5) = y(1) / m
     J(3,6) = 0.0
     J(3,7) = 0.0

     J(4,1) = 0.0
     J(4,2) = u(1) / m
     J(4,3) = 0.0
     J(4,4) = 0.0
     J(4,5) = y(2) / m
     J(4,6) = 0.0
     J(4,7) = -1.0
  end if

end subroutine odecon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

subroutine algcon ( snStat, curPhs, nPhs, nC, nY, nU, nP, C, G, y, u, p, &
                    needC, needG, cu, lencu, iu, leniu, ru, lenru)
  use precision, only : ip, rp
  implicit none
  integer(ip), intent(in)    :: snStat, curPhs, nPhs, nC, nY, nU, nP, &
                               needC, needG, lencu, leniu, lenru
  real(rp),    intent(in)    :: y(nY), u(nU), p(nP)

  integer(ip), intent(inout) :: iu(leniu)
  real(rp),    intent(inout) :: ru(lenru)
  character(8),intent(inout) :: cu(lencu)

  real(rp),    intent(out)   :: C(nC), G(nC, nY+nU+nP)

  !=============================================================================
  real(rp) :: m, l, xf, yf

  m  = .3d+0
  l  = .5d+0
  xf = -.231625d+0
  yf = -.443101d+0

  if (needC > 0) then
     C(1) = y(3)**2 + y(4)**2 + u(1)*(l**2) / m - y(2)*p(1)
     C(2) = u(2) - .50*( (y(1) - xf)**2 + (y(2) - yf)**2 )
  end if

  if (needG > 0) then
     G(1,1) = 0.0
     G(1,2) = -p(1)
     G(1,3) = 2*y(3)
     G(1,4) = 2*y(4)
     G(1,5) = l**2/m
     G(1,6) = 0.0
     G(1,7) = -y(2)

     G(2,1) = -(y(1) - xf)
     G(2,2) = -(y(2) - yf)
     G(2,3) = 0.0
     G(2,4) = 0.0
     G(2,5) = 0.0
     G(2,6) = 1.0
     G(2,7) = 0.0
  end if

end subroutine algcon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
