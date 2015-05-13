!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File rocketA.f90
!
!  The original problem formulation:
!        min            t1
!        subject to:
!                       y1'   = y2
!                       y2'   = 2*u1 - 1
!
!                       y1(0) = 0
!                       y2(0) = 0
!
!                       y1(tf) = 100
!                       u1 = 1 for t in [t0,t1]
!                       u1 = 0 for t in [t1,tf]
!
!  The problem input to SNCTRL:
!        min            p1
!        subject to:
!                       y1'   = y2*p1
!                       y2'   = (2*u1 - 1)*p1
!
!                       y1(0) = 0
!                       y2(0) = 0
!
!                       y1(tf) = 100
!                       u1 = 1 for t in [t0,t1]
!                       u1 = 0 for t in [t1,tf]
!                       p1 >= 0
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program rocketA
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

  integer(ip)  :: j

  external     :: odecon, algcon, varbds
  type(ctProb) :: prob

  !=====================================================================

  iSpecs =  4  ! equivalenced to .spc
  iPrint =  9  ! equivalenced to .out
  iSumm  =  6  ! summary file goes to standard output...
  nOut   =  6  ! ... as do messages from this program.

  ! Open the Specs and print files.
  lfile = 'rocketA.spc'
  open( iSpecs, file=lfile, status='OLD',     err=800 )
  lfile = 'rocketA.out'
  open( iPrint, file=lfile, status='UNKNOWN', err=800 )


  ! Set options to default values.
  call sncInit( iPrint, iSumm, prob, cw, lencw, iw, leniw, rw, lenrw )


  ! Read a Specs file.
  call sncSpec( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

  if (INFO /= 101  .AND.  INFO /= 107) go to 910

  ! Set up problem
  prob%probname(1:8) = 'rocketA'
  prob%nY = 2       ! number of states
  prob%nU = 1       ! number of controls
  prob%nP = 1       ! number of parameters
  prob%nC = 0       ! number of algebraic constraints
  prob%nPhs = 2     ! number of phases

  ! Objective
  allocate ( prob%objL(2) )
  prob%objL(1) = prob%nY+prob%nU+1
  prob%objL(2) = 0

  ! Interval
  allocate ( prob%phsPt(3), prob%npInt(2) )
  prob%phsPt(1) = 0.0
  prob%phsPt(2) = 0.5
  prob%phsPt(3) = 1.0

  prob%npInt(1) = 5
  prob%npInt(2) = 5

  ! Linearity
  allocate ( prob%ytype(2,2) )
  do j = 1, 2
     prob%ytype(1,j) = 0
     prob%ytype(2,j) = 0
  end do



  !---------------------------------------------------------------------
  ! Solve the problem, call snctrl.
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

4000 format(/  a, 2x, a  )

end program rocketA

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
  integer(ip) :: j
  real(rp)    :: bplus, bminus
  parameter    ( bplus  = 1.0d+20, bminus = -bplus )

  y0low = bminus
  y0upp = bplus

  yflow = bminus
  yfupp = bplus

  ylow = bminus
  yupp = bplus

  ulow = bminus
  uupp = bplus

  if ( curPhs == 1 ) then
     ! At t0
     y0low(1) = 0.0
     y0upp(1) = 0.0
     y0low(2) = 0.0
     y0upp(2) = 0.0

     ulow(1) = 1.0
     uupp(1) = 1.0
  end if

  if ( curPhs == nPhs ) then
     yflow(1) = 100.0
     yfupp(1) = 100.0

     ulow(1) = 0.0
     uupp(1) = 0.0
  end if

  ! Parameter
  plow(1) = 0.0
  pupp(1) = bplus
  p(1) = 1.0

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

  if ( needF > 0 ) then
     F(1) = y(2)*p(1)
     F(2) = (2*u(1) - 1) * p(1)
  end if

  if ( needJ > 0 ) then
     J(1,1) = 0.0
     J(1,2) = p(1)
     J(1,3) = 0.0
     J(1,4) = y(2)

     J(2,1) = 0.0
     J(2,2) = 0.0
     J(2,3) = 2*p(1)
     J(2,4) = 2*u(1) - 1
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

end subroutine algcon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
