!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File catmixS.f90
!
!  This problem has an algebraic objective.
!  The original problem:
!      min            y1 + y2 - 1
!      subject to:
!                     y1' = u1 * (10*y2 - y1)
!                     y2' = u1 * (y1 - 10*y2) - (1-u1)*y2
!
!                     y1(0) = 1
!                     y2(0) = 0
!                     0 <= u1 <= 1
!
!  The problem input to SNCTRL with an added control to represent
!  the algebraic objective:
!      min            u2
!      subject to:
!                     y1' = u1 * (10*y2 - y1)
!                     y2' = u1 * (y1 - 10*y2) - (1-u1)*y2
!
!                     y1(0) = 1
!                     y2(0) = 0
!                     0 <= u1 <= 1
!
!                     u2 = y1 + y2 - 1
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

program catmixS
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
  lfile = 'catmixS.spc'
  open( iSpecs, file=lfile, status='OLD',     err=800 )
  lfile = 'catmixS.out'
  open( iPrint, file=lfile, status='UNKNOWN', err=800 )


  ! Set options to default values.
  call sncInit( iPrint, iSumm, prob, cw, lencw, iw, leniw, rw, lenrw )


  ! Read a Specs file.
  call sncSpec( iSpecs, INFO, cw, lencw, iw, leniw, rw, lenrw )

  if (INFO /= 101  .AND.  INFO /= 107) go to 910


  ! Setup the problem.
  prob%probname(1:8) = ' catmixS'
  prob%nY   = 2     ! number of states
  prob%nU   = 2     ! number of controls
  prob%nP   = 0     ! number of parameters
  prob%nC   = 1     ! number of algebraic constraints
  prob%nPhs = 1     ! number of phases

  ! Objective
  allocate ( prob%objL(1) )
  prob%objL(1) = prob%nY+2

  ! Derivatives
  allocate ( prob%neJ(1), prob%neG(1) )
  prob%neJ(1) = 6
  prob%neG(1) = 3

  ! Interval
  allocate ( prob%npInt(1), prob%phsPt(2) )
  prob%npInt(1) = 10
  prob%phsPt(1) = 0_rp
  prob%phsPt(2) = 1_rp

  ! Linearity
  allocate ( prob%ytype(2,1), prob%ctype(1,1) )
  prob%ytype(1,1) = 0
  prob%ytype(2,1) = 0

  prob%ctype(1,1) = 1


  !---------------------------------------------------------------------
  ! Solve the problem
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

4000 format(/  a, 2x, a  )

end program catmixS

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
  real(rp)    :: bplus, bminus
  parameter    ( bplus  = 1.0d+20, bminus = -bplus )

  ! Default set everything to free bounds.
  lbds = bminus
  ubds = bplus
  x    = 0.0

  ! At t0
  lbds(1,1) = 1.0
  lbds(2,1) = 0.0

  ubds(1,1) = 1.0
  ubds(2,1) = 0.0

  ! 0 <= u1 <= 1
  do j = 1, nNodes
     lbds(1+nY,j) = 0.0
     ubds(1+nY,j) = 1.0
  end do

  ! algebraic constraint
  clbds(1) = 1.0
  cubds(1) = 1.0

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

  if ( needF > 0 ) then
     do jt = 1, nNodes
        F(1,jt) = dvar(1+nY,jt)*(10*dvar(2,jt) - dvar(1,jt))
        F(2,jt) = dvar(1+nY,jt)*(dvar(1,jt) - 10*dvar(2,jt)) - &
             (1-dvar(1+nY,jt))*dvar(2,jt)
     end do
  end if

  if ( needJ > 0 ) then
     do jt = 1, nNodes
        neJ = 1

        ! Column 1
        Jcol(1) = neJ

        Jrow(neJ) = 1
        Jval(neJ,jt) = -dvar(1+nY,jt)
        neJ = neJ + 1

        Jrow(neJ) = 2
        Jval(neJ,jt) = dvar(1+nY,jt)
        neJ = neJ + 1

        ! Column 2
        Jcol(2) = neJ

        Jrow(neJ) = 1
        Jval(neJ,jt) = 10.0*dvar(1+nY,jt)
        neJ = neJ + 1

        Jrow(neJ) = 2
        Jval(neJ,jt) = -10.0*dvar(1+nY,jt) - (1-dvar(1+nY,jt))
        neJ = neJ + 1

        ! Column 3
        Jcol(3) = neJ

        Jrow(neJ) = 1
        Jval(neJ,jt) = 10.0*dvar(2,jt) - dvar(1,jt)
        neJ = neJ + 1

        Jrow(neJ) = 2
        Jval(neJ,jt) = dvar(1,jt) - 10.0*dvar(2,jt) + dvar(2,jt)
        neJ = neJ + 1

        ! Column 4
        Jcol(4) = neJ

        ! Finish off pts
        Jcol(5) = neJ
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

  if ( needC > 0 ) then
     do jt = 1, nNodes
        C(1,jt) = dvar(1,jt) + dvar(2,jt) - dvar(2+nY,jt)
     end do
  end if

  if ( needG > 0 ) then
     do jt = 1, nNodes
        neG = 1

        ! Column 1
        Gcol(1) = neG

        Grow(neG) = 1
        Gval(neG,jt) = 1.0
        neG = neG + 1

        ! Column 2
        Gcol(2) = neG

        Grow(neG) = 1
        Gval(neG,jt) = 1.0
        neG = neG + 1

        ! Column 3
        Gcol(3) = neG

        ! Column 4
        Gcol(4) = neG

        Grow(neG) = 1
        Gval(neG,jt) = -1.0
        neG = neG + 1

        ! Finish off ptrs
        Gcol(5) = neG
        neG = neG - 1
     end do
  end if

end subroutine algcon

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
