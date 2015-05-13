!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File:  ct30ker.f90
! Kernel for the control interface.
!
! 06 Feb 2010: Current version.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module ct30ker
  use precision
  use ct02data, only : ctProb, ctWork, ctGrid, c0grid, c0gtrsh, c0wtrsh

  use ct15usr,  only : odeconS, odeconA, odeconD, algconS, algconA, algconD
  use ct16refn, only : c1refn
  use ct20eval, only : sncEvalTR, sncEvalHS
  use ct28derv, only : odeSetupTR, algSetupTR, odeSetupHS, algSetupHS

  implicit none

  private
  public  :: c3ker
  private :: s0fgctrl, ctrlObj, ctrlCon, &
             c2solv, c3err, c1prnt, c1psol, c2setup

  real(rp),    private, pointer :: gstep(:)
  integer(ip), private, pointer :: neJ(:), neG(:), gndPtr(:), &
                                   ytypeG(:,:), ctypeG(:,:)

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c3ker( Start, prob, usrode, usralg, usrbds, &
                    INFO, mincw, miniw, minrw, nS, nInf, sInf, &
                    cu, lencu, iu, leniu, ru, lenru, &
                    cw, lencw, iw, leniw, rw, lenrw )
    integer(ip),  intent(in)    :: Start, lencu, leniu, lenru, &
                                   lencw, leniw, lenrw

    integer(ip),  intent(inout) :: nS, iu(leniu), iw(leniw)
    real(rp),     intent(inout) :: rw(lenrw), ru(lenru)
    character(8), intent(inout) :: cw(lencw), cu(lencu)

    integer(ip),  intent(out)   :: INFO, mincw, miniw, minrw, nInf
    real(rp),     intent(out)   :: sInf

    external             :: usrode, usralg, usrbds
    type(ctProb), target :: prob

    !===========================================================================
    ! On entry, global option pointers have been set.
    ! 1. Sets up global problem pointers.
    ! 2. Create the uniform grid based on user data and set up the workspace.
    ! 3. Call c2solv to solve problem.
    ! 4. Print solution.
    ! 5. Check for refinement.
    !
    ! 26 Dec 2008: First version.
    ! 09 Feb 2010: v5.
    !===========================================================================
    integer(ip)    :: Errors, nItn, refLoop, nAdded, i, nY, nU, nC, nP, nPhs
    real(rp)       :: r, infBnd
    character(80)  :: str
    type(ctGrid)   :: grid
    type(ctWork)   :: work

    integer(ip), pointer     :: ctype(:), ytype(:), objL(:)
    integer(ip), allocatable :: nnY(:), nlY(:), nnC(:), nlC(:)

    ! Error checks
    call c3err( iw(iIntf), prob, Errors, iw, leniw )
    if ( Errors /= 0 ) then
       INFO = 99
       return
    end if


    ! If refining, set to TR method.
    if ( iw(iRefn) /= 0 ) iw(iDisc) = 0

    nItn    = 0
    refLoop = 0
    if ( iw(iRefL) < 0 ) &   ! Refine until no nodes are added.
         refLoop = 1

    infBnd = rw(70)
    if ( infBnd <= 0.0 ) infBnd = 1.0d+20

    nY        = prob%nY
    nU        = prob%nU
    nP        = prob%nP
    nC        = prob%nC
    nPhs      = prob%nPhs

    ! Save problem dimensions to workspace.
    iw(inY)   = nY
    iw(inU)   = nU
    iw(inP)   = nP
    iw(inC)   = nC
    iw(inPhs) = nPhs


    objL   => prob%objL
    ytypeG => prob%ytype
    if ( nC > 0 ) then
       ctypeG => prob%ctype
    else
       allocate( ctypeG(1,nPhs) )
       ctypeG = 0
    end if

    allocate( nnY(nPhs), nlY(nPhs), nnC(nPhs), nlC(nPhs) )
    nnY = 0
    nlY = 0
    nnC = 0
    nlC = 0

    nlY = sum( ytypeG, dim=1 )  ! count # linear ode constraints
    nnY = nY - nlY  ! # nonlinear ode constraints

    if ( nC > 0 ) then
       nlC = sum( ctypeG, dim=1 ) ! # linear alg constraints
       nnC = nC - nlC  ! # nonlinear alg constraints
    end if

    ! Derivatives in each phase.
    if ( iw(iIntf) == 0 ) then
       neJ => prob%neJ

       if ( nC > 0 ) then
          neG => prob%neG
       else
          allocate( neG(nPhs) )
          neG = 0
       end if

    else
       allocate( neJ(nPhs), neG(nPhs) )
       neJ = nY*( nY + nU + nP )
       neG = nC*( nY + nU + nP )
    end if


    ! Set the (uniform) grid.
    call c0grid( iw(iDisc), grid, nPhs, prob%phsPt, prob%npInt )

    ! Set the workspace.
    call c2setup( iw(iDisc), nPhs, nY, nU, nP, nC, objL, ytypeG, ctypeG, &
                  nnY, nlY, nnC, nlC, neJ, neG, prob%x, prob%hs, &
                  grid, work, usrode, usralg, usrbds, infBnd, &
                  cu, lencu, iu, leniu, ru, lenru, &
                  cw, lencw, iw, leniw, rw, lenrw )

    ! Solve the problem.
    call c2solv( Start, prob%probName, prob%x, prob%hs, work, grid, usrode, usralg, &
                 INFO, mincw, miniw, minrw, nS, nInf, sInf, &
                 cu, lencu, iu, leniu, ru, lenru, &
                 cw, lencw, iw, leniw, rw, lenrw )
    nItn = nItn + 1

    ! Print solution
    if ( iw(iCPrt) /= 0 ) call c1prnt( prob, grid, prob%x, work, iw, leniw, rw, lenrw )

    r = INFO / 10.0
    if ( iw(iRefn) == 0 .or. floor(r) /= 0 ) go to 900


    !---------------------------------------------------------------------------
    ! Adaptive Refinement
    ! sncRefn will create a new grid (replacing the old one) and calculate a
    ! new initial point for the refined problem.

100 call c1refn( iw(iIntf), rw(ctrefTol), nPhs, nY, nU, nP, neJ, work%n, &
                 prob%x, prob%hs, grid, usrode, nAdded, &
                 cu, lencu, iu, leniu, ru, lenru, &
                 cw, lencw, iw, leniw, rw, lenrw )


    ! Print a header...
    str = ''
    call snPRNT(13, str, iw, leniw)
    do i = 1, 80
       str = str(1:i) // '*'
    end do
    call snPRNT(13, str, iw, leniw)
    write(str,*) '   Refinement Run: ', nAdded, ' nodes added.'
    call snPRNT( 3, str, iw, leniw)

    if( nAdded == 0 ) go to 900

    ! Set up the optimal control problem.
    call c0wtrsh( work )
    call c2setup( iw(iDisc), nPhs, nY, nU, nP, nC, objL, ytypeG, ctypeG, &
                  nnY, nlY, nnC, nlC, neJ, neG, prob%x, prob%hs, &
                  grid, work, usrode, usralg, usrbds, infBnd, &
                  cu, lencu, iu, leniu, ru, lenru, &
                  cw, lencw, iw, leniw, rw, lenrw )

    ! Solve the problem.
    call c2solv( 0_ip, prob%probName, prob%x, prob%hs, work, grid, &
                 usrode, usralg, &
                 INFO, mincw, miniw, minrw, nS, nInf, sInf, &
                 cu, lencu, iu, leniu, ru, lenru, &
                 cw, lencw, iw, leniw, rw, lenrw )

    nItn = nItn + 1

    ! Print solution
    if ( iw(iCPrt) /= 0 ) call c1prnt( prob, grid, prob%x, work, iw, leniw, rw, lenrw )

    ! Error exit?
    r = INFO / 10.0
    if ( floor(r) /= 0 ) &
         go to 900

    if ( nItn <= iw(iRefL) .or. refLoop == 1 ) &
         go to 100


    ! Clean up.
900 call c0wtrsh( work )
    call c0gtrsh( grid )

    deallocate( nnY, nlY, nnC, nlC )

    if ( nC == 0 ) deallocate( ctypeG )

    if ( iw(iIntf) /= 0 ) then
       deallocate( neJ, neG )
    else
       if ( nC == 0 ) deallocate( neG )
    end if

  end subroutine c3ker

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c2solv ( Start, prob, x, hs, work, grid, usrode, usralg, &
                      INFO, mincw, miniw, minrw, nS, nInf, sInf, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: Start, lencu, leniu, lenru, lencw, leniw, lenrw
    character(8),intent(in)    :: prob
    integer(ip), intent(inout) :: nS, hs(:), iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: x(:), ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cw(lencw), cu(lencu)
    integer(ip), intent(out)   :: INFO, mincw, miniw, minrw, nInf
    real(rp),    intent(out)   :: sInf
    external     :: usrode, usralg
    type(ctWork), target :: work
    type(ctGrid), target :: grid

    !===========================================================================
    ! Calls the SNOPT(B) kernel.
    !
    ! 24 Dec 2008: Renamed sncSolve.
    ! 09 Feb 2010: v5.
    !===========================================================================
    character    :: CStart*4, Names(1)*8
    integer(ip)  :: m, n, neJ, nnObj, nnJac, nnCon, iObj
    real(rp)     :: ObjAdd, Obj

    integer(ip), parameter :: nName = 1
    integer(ip), pointer   :: indJ(:), locJ(:)
    real(rp),    pointer   :: Jcol(:), bl(:), bu(:), rc(:), pi(:)

    external     :: snLog, snLog2, sqLog, snSTOP


    ! Set global grid stuff (passed to ctrlCon below).
    gstep   => grid%step
    gndPtr  => grid%ndPtr

    ObjAdd = 0.0

    CStart = 'Cold'
    if ( Start == 2 ) then
       CStart = 'Warm'
    end if

    m     = work%m
    n     = work%n
    neJ   = work%lenJ
    nnObj = 0
    nnJac = work%n
    nnCon = work%nnCon
    iObj  = work%m

    Jcol => work%Jcol
    indJ => work%indJ
    locJ => work%locJ
    bl   => work%bl
    bu   => work%bu
    pi   => work%pi
    rc   => work%rc

    call snKerCT( CStart, m, n, neJ, nName, nnCon, nnObj, nnJac, &
                  iObj, ObjAdd, prob, s0fgctrl, usrode, usralg, &
                  snLog, snLog2, sqLog, snSTOP, &
                  Jcol, indJ, locJ, bl, bu, Names, hs, x, pi, rc, &
                  INFO, mincw, miniw, minrw, nS, nInf, sInf, Obj, &
                  cu, lencu, iu, leniu, ru, lenru, &
                  cw, lencw, iw, leniw, rw, lenrw )
    nullify ( gstep, gndPtr )

  end subroutine c2solv

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine s0fgctrl ( iExit, modefg, NonlinearCon, NonlinearObj, &
                        n, negCon, nnCon0, nnCon, nnJac, nnL, nnObj0, nnObj, &
                        usrode, usralg, dummyH, x, yCon, &
                        neJ, nlocJ, locJ, indJ, &
                        neH, nlocH, locH, indH, Hcol, &
                        fCon, fObj, gCon, gObj, &
                        cu, lencu, iu, leniu, ru, lenru, &
                        cw, lencw, iw, leniw, rw, lenrw )
    external :: &
         usrode, usralg, dummyH
    logical :: &
         NonlinearCon, NonlinearObj
    integer(ip) :: &
         iExit, lencu, leniu, lenru, lencw, leniw, lenrw, modefg, n, &
         neJ, neH, negCon, nlocJ, nlocH, nnCon0, nnCon, nnJac, nnL, &
         nnObj0, nnObj, indJ(neJ), locJ(nlocJ), &
         indH(neH), locH(nlocH), iu(leniu), iw(leniw)
    real(rp) :: &
         fObj, fCon(nnCon0), gObj(nnObj0), gCon(negCon), Hcol(neH), &
         x(n), yCon(nnCon0), ru(lenru), rw(lenrw)
    character(8) :: &
         cu(lencu), cw(lencw)
    ! ==========================================================================
    ! s0fgctrl  is an instance of fgwrap that calls the routines
    !   ctrlCon   and   ctrlObj   to evaluate the problem functions
    ! and possibly their gradients.
    ! usrode and usralg are user-written routines which provide
    ! the ODEs and algebraic constraints and the Jacobians.
    !
    ! Arguments  ctrlCon  and  ctrlObj  are called using modefg to
    ! control the gradients as follows:
    !
    ! modefg        Task
    ! ------        ----
    ! 2     Assign fCon, fObj and all known elements of gCon and gObj.
    ! 1     Assign all known elements of gCon and gObj.
    ! (fObj and fCon are ignored).
    ! 0     Assign fObj, fCon.  (gCon and gObj are ignored).
    !
    ! If s0fgctrl is called with minmax = 0 (feasible point only) then
    ! nnObj = max(nnJac,nnObj)  and the user objective is not used.
    !
    ! 06 Sep 2008: s0fgB adapted for control interface.
    ! 29 Apr 2015: Updated.
    ! ==========================================================================
    character &
         str*80
    external :: &
         ddot
    logical :: &
         FPonly, scaled
    integer(ip) :: &
         gotFD, gotG, gotGl, l, lG, lgConu, lgObju, lvlScl, lvlTim, &
         lvlDer, lAscal, lx0, lxScal, minmax, modeC, modeF, ngrad, &
         nfCon1, nfCon2, nfObj1, nfObj2, nnGlin, nnObjU, Status, StatusUser
    real(rp) :: &
         ddot, Gdummy, proxWeight
    ! ------------------------------------------------------------------
    parameter         (lvlDer =  70) ! = 0,1,2 or 3, deriv. level
    parameter         (gotFD  = 183) ! > 0 => some differences needed
    parameter         (gotG   = 184) ! > 0 => some exact derivs
    parameter         (gotGl  = 185) ! > 0 => constant Jacob elements
    parameter         (nfCon1 = 189) ! calls to fCon: mode = 0
    parameter         (nfCon2 = 190) ! calls to fCon  mode > 0
    parameter         (nfObj1 = 194) ! calls to fObj: mode = 0
    parameter         (nfObj2 = 195) ! calls to fObj: mode > 0

    real(rp) ::   half,            one
    parameter    (half   = 0.5d+0, one   = 1.0d+0)
    ! ------------------------------------------------------------------
    nnObjU    = iw( 22) ! # of objective variables
    lvlScl    = iw( 75) ! scale option
    minmax    = iw( 87) ! 1, 0, -1  => MIN, FP, MAX
    lvlTim    = iw(182) ! Timing level

    lAscal    = iw(296) ! Ascale(nb) = row and column scales
    lx0       = iw(298) ! x0(nnL)    = proximal point starting point
    lxScal    = iw(302) ! xScal(n)   = copy of scaled  x
    lgConu    = iw(319) ! record of unknown derivatives and constants
    lgObju    = iw(323) ! record of unknown derivatives

    Gdummy    = rw( 69) ! definition of an 'unset' value
    proxWeight= rw( 91) ! Proximal-point weight

    iExit     = 0

    FPonly    = minmax == 0
    scaled    = lvlScl >  2

    modeC     = modefg
    modeF     = modefg

    ! Determine the user-function call-status.
    call s8callStatus( Status, iw, leniw )

    if (Status == 1) then
       !--------------------------------------------------------------
       ! First evaluation of the problem functions in snOptB
       ! On entry, lvlScl = 0.
       !--------------------------------------------------------------
       iw(gotFD) =  0
       iw(gotG)  =  0
       iw(gotGl) =  0
       call snPRNT( 13, ' ', iw, leniw )
       call dload ( negCon, Gdummy, gCon, 1 )
       call dload ( nnObj , Gdummy, gObj, 1 )
    end if

    !-----------------------------------------------------------------
    ! Unscale x (never required for Status = 1)
    !-----------------------------------------------------------------
    if ( scaled ) then
       call dcopy( nnL, x         , 1, rw(lxScal), 1 )
       call ddscl( nnL, rw(lAscal), 1, x         , 1 )

       ! If the Jacobian has some constant elements, they are wrecked
       ! by the scaling.  Restore them from gConu.

       if ( NonlinearCon ) then
          if (modefg > 0  .AND.  iw(gotGl) > 0) then
             call dcopy ( negCon, rw(lgConu), 1, gCon, 1 )
          end if
       end if
    end if

    !=================================================================
    ! Compute the constraint functions.
    !=================================================================
    if ( NonlinearCon ) then
       statusUser = Status
       if (lvlTim >= 2) call s1time( 4, 0, iw, leniw, rw, lenrw )
       call ctrlCon &
            ( modeC, nnCon, nnJac, negCon, &
              x, fCon, gCon, statusUser, &
              usrode, usralg, &
              cu, lencu, iu, leniu, ru, lenru, &
              cw, lencw, iw, leniw, rw, lenrw )
       if (lvlTim >= 2) call s1time(-4, 0, iw, leniw, rw, lenrw )
       iw(nfCon1) = iw(nfCon1) + 1
       if (modefg > 0) &
            iw(nfCon2) = iw(nfCon2) + 1
    end if

    !=================================================================
    ! Compute the objective function.
    !=================================================================
    if (NonlinearObj  .AND.  modeC >= 0) then
       if ( FPonly ) then
          call dcopy( nnObj, x, 1, gObj, 1 )
          call daxpy( nnObj, (-one), rw(lx0), 1, gObj, 1 )
          fObj = half*proxWeight*ddot ( nnObj, gObj, 1, gObj, 1 )
          call dscal(nnObj, proxWeight, gObj, 1 )
       else ! nnObj = nnObj
          statusUser = Status
          if (lvlTim >= 2) call s1time( 5, 0, iw, leniw, rw, lenrw )
          call ctrlObj &
               ( modeF, nnObjU, x, fObj, gObj, statusUser, usrode, &
                 cu, lencu, iu, leniu, ru, lenru, &
                 cw, lencw, iw, leniw, rw, lenrw )
          if (lvlTim >= 2) call s1time(-5, 0, iw, leniw, rw, lenrw )
          iw(nfObj1) = iw(nfObj1) + 1
          if (modefg > 0) &
               iw(nfObj2) = iw(nfObj2) + 1
       end if
    end if

    !-----------------------------------------------------------------
    ! Scale  x and the derivatives.
    !-----------------------------------------------------------------
    if ( scaled ) then
       call dcopy ( nnL, rw(lxScal), 1, x, 1 )

       if ( NonlinearCon ) then
          call dddiv ( nnCon, rw(lAscal+n), 1, fCon, 1 )
          if (modefg > 0  .AND.  iw(gotG) > 0) then
             call s8scaleJ &
                  ( nnCon, nnJac, negCon, n, rw(lAscal), &
                    neJ, nlocJ, locJ, indJ, gCon, rw, lenrw )
          end if
       end if

       if (NonlinearObj  .AND.  modeC >= 0) then
          if (modefg > 0  .AND.  iw(gotG) > 0) then
             call s8scaleG &
                  ( nnObj, rw(lAscal), gObj, rw, lenrw )
          end if
       end if
    end if

    if (modeC < 0  .OR.  modeF < 0) then
       !--------------------------------------------------------------
       ! The user may be saying the function is undefined (mode = -1)
       ! or may just want to stop                         (mode < -1).
       !--------------------------------------------------------------
       if (modeC == -1  .OR.  modeF == -1) then
          iExit = -1
       else
          if (modeC < 0) then
             iExit = 72
          else
             iExit = 73
          end if
       end if
    end if

    !=================================================================
    ! Do some housekeeping on the first snOptB entry.
    !=================================================================
    if (Status == 1  .AND.  iExit == 0) then
       if ( NonlinearCon ) then
          !-----------------------------------------------------------
          ! Count how many Jacobian elements are provided.
          !-----------------------------------------------------------
          nnGlin = 0
          ngrad  = 0
          do l = 1, negCon
             if (gCon(l) /= Gdummy) ngrad  = ngrad + 1
          end do

          write(str, 1100) ngrad, negCon
          call snPRNT( 3, str, iw, leniw )

          if (ngrad < negCon) then

             ! Some Jacobian elements are missing.

             if (iw(lvlDer) >= 2) then
                !-----------------------------------------------------
                ! All the Jacobian is known.  Any undefined elements
                ! are assumed constant, and are restored from gConu.
                !-----------------------------------------------------
                call snPRNT( 3, &
                     ' ==>  Some constraint derivatives are missing, ' &
                     //' assumed constant.', iw, leniw )
                call snPRNT( 3, ' ', iw, leniw )

                lG  = lgConu
                do l  = 1, negCon
                   if (gCon(l) == Gdummy) then
                      gCon(l) = rw(lG)
                      nnGlin  = nnGlin + 1
                   end if
                   lG = lG + 1
                end do
             else
                !-----------------------------------------------------
                ! Save a permanent copy of gCon in gConu so that we know
                ! which derivatives must be estimated.
                !-----------------------------------------------------
                call dcopy ( negCon, gCon, 1, rw(lgConu), 1 )
             end if
          end if ! ngrad < negCon
          if (ngrad + nnGlin < negCon) iw(gotFD) = 1
          if (ngrad          >      0) iw(gotG ) = 1
          if (nnGlin         >      0) iw(gotGl) = 1
       end if

       if ( NonlinearObj ) then
          !-----------------------------------------------------------
          ! Count how many working gradient elements are known.
          ! (These may be the gradients of the FP objective.)
          !-----------------------------------------------------------
          if ( FPonly ) then
             ngrad    = nnObj
             iw(gotG) = 1

             write(str, 2010) nnObj
             call snPRNT( 3, str, iw, leniw )
          else
             ngrad = 0
             do l = 1, nnObjU
                if (gObj(l) /= Gdummy) ngrad = ngrad + 1
             end do
          end if

          write(str, 2000) ngrad, nnObj
          call snPRNT( 3, str, iw, leniw )

          if ( ngrad > 0 ) then
             iw(gotG) = 1
          end if

          if (ngrad < nnObjU) then

             ! Some objective gradients are missing.
             iw(gotFD) = 1

             if (iw(lvlDer) == 1  .OR.  iw(lvlDer) == 3) then
                !-----------------------------------------------------
                ! The objective gradient was meant to be known.
                !-----------------------------------------------------
                iw(lvlDer) = iw(lvlDer) - 1
                write(str, 2100) iw(lvlDer)
                call snPRNT( 3, str, iw, leniw )
             end if
          end if

          !--------------------------------------------------------
          ! Copy gObj into gObju.
          !--------------------------------------------------------
          call dcopy ( nnObj, gObj, 1, rw(lgObju), 1 )
       end if
    end if

    return

1100 format(' The user has defined', i8, '   out of', i8, &
         '   constraint gradients.')
2000 format(' The user has defined', i8, '   out of', i8, &
         '   objective  gradients.')
2010 format(' SnOptB  will define ', i8, '   gradients for the ', &
         ' FP objective.')
2100 format(' XXX  Some objective  derivatives are missing ---', &
         ' derivative level reduced to', i3)

  end subroutine s0fgctrl

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine ctrlCon( mode, nnCon, nnJac, neJac, x, fCon, gCon, Status, &
                      usrode, usralg, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: mode, nnCon, nnJac, neJac, Status, &
                                  lencu, leniu, lenru, lencw, leniw, lenrw
    real(rp),    intent(in)    :: x(nnJac)
    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)
    real(rp),    intent(out)   :: fCon(nnCon), gCon(neJac)
    external     :: usrode, usralg

    !===========================================================================
    ! Returns the constraint/derivative values.
    !===========================================================================
    integer(ip) :: s

    s  = size(gstep)

    ! SNCTRLA
    if ( iw(iDisc) == 0 ) then
       call sncEvalTR( Status, mode, &
                       iw(inY), iw(inU), iw(inP), iw(inC), iw(inPhs), &
                       gndPtr, gstep, s, ytypeG, ctypeG, neJ, neG, &
                       x, nnJac, fCon, nnCon, gCon, neJac, &
                       usrode, usralg, &
                       cu, lencu, iu, leniu, ru, lenru, &
                       cw, lencw, iw, leniw, rw, lenrw )
    else
       call sncEvalHS( Status, mode, &
                       iw(inY), iw(inU), iw(inP), iw(inC), iw(inPhs), &
                       gndPtr, gstep, s, ytypeG, ctypeG, neJ, neG, &
                       x, nnJac, fCon, nnCon, gCon, neJac, &
                       usrode, usralg, &
                       cu, lencu, iu, leniu, ru, lenru, &
                       cw, lencw, iw, leniw, rw, lenrw )
    end if

  end subroutine ctrlCon

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine ctrlObj ( mode, nnObj, x, fObj, gObj, Status, usrode, &
                       cu, lencu, iu, leniu, ru, lenru, &
                       cw, lencw, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: mode, nnObj, Status, &
                                  lencu, leniu, lenru, lencw, leniw, lenrw
    real(rp),    intent(in)    :: x(nnObj)
    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)
    real(rp),    intent(out)   :: fObj, gObj(nnObj)
    external     :: usrode

    !===========================================================================
    ! This is funObj for SNOPTB.
    ! Nothing to do!  Linear objective!
    !===========================================================================

  end subroutine ctrlObj

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c3err( iIntf, prob, Errors, iw, leniw )
    integer(ip), intent(in)    :: iIntf, leniw
    integer(ip), intent(inout) :: iw(leniw)
    integer(ip), intent(out)   :: Errors
    type(ctProb) :: prob

    !===========================================================================
    ! Very simple error checks on the ctProb structure.
    !
    ! 21 Aug 2008: First version of errCheck.
    ! 24 Dec 2008: Updated for v4 of snctrl.
    ! 09 Feb 2010: v5.
    !===========================================================================
    character :: str*80

    Errors = 0

    if ( prob%nY < 0 ) then
       Errors = Errors + 1
       write(str,1000) 'nY', prob%nY
       call snPRNT( 13, str, iw, leniw )
    end if

    if ( prob%nU < 0 ) then
       Errors = Errors + 1
       write(str,1000) 'nU', prob%nU
       call snPRNT( 13, str, iw, leniw )
    end if

    if ( prob%nC < 0 ) then
       Errors = Errors + 1
       write(str,1000) 'nC', prob%nC
       call snPRNT( 13, str, iw, leniw )
    end if

    if ( prob%nP < 0 ) then
       Errors = Errors + 1
       write(str,1000) 'nP', prob%nP
       call snPRNT( 13, str, iw, leniw )
    end if

    if ( prob%nPhs < 0 ) then
       Errors = Errors + 1
       write(str,1000) 'nPhs', prob%nPhs
       call snPRNT( 13, str, iw, leniw )
    end if

    if ( .not. allocated(prob%objL) ) then
       Errors = Errors + 1
       write(str,1100) 'objL'
       call snPRNT( 13, str, iw, leniw )
    end if

    if ( .not. allocated(prob%ytype) ) then
       Errors = Errors + 1
       write(str,1100) 'ytype'
       call snPRNT( 13, str, iw, leniw )
    end if

    if ( prob%nC > 0 ) then
       if ( .not. allocated(prob%ctype) ) then
          Errors = Errors + 1
          write(str,1100) 'ctype'
          call snPRNT( 13, str, iw, leniw )
       end if
    end if

    if ( .not. allocated(prob%npInt) ) then
       Errors = Errors + 1
       write(str,1100) 'npInt'
       call snPRNT( 13, str, iw, leniw )
    end if

    if ( .not. allocated(prob%phsPt) ) then
       Errors = Errors + 1
       write(str,1100) 'phsPt'
       call snPRNT( 13, str, iw, leniw )
    end if

    if ( iIntf == 0 ) then
       if ( .not. allocated(prob%neJ) ) then
          Errors = Errors + 1
          write(str,1100) 'neJ'
          call snPRNT( 13, str, iw, leniw )
       end if

       if ( prob%nC > 0 ) then
          if ( .not. allocated(prob%neG) ) then
             Errors = Errors + 1
             write(str,1100) 'neG'
             call snPRNT( 13, str, iw, leniw )
          end if
       end if
    end if

    return

1000 format('XXX snctrl: argument out of range: ', a6, ' = ', i6 )
1100 format('XXX snctrl: argument not allocated: ', a6 )

  end subroutine c3err

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c2setup ( iDisc, nPhs, nY, nU, nP, nC, objL, &
                       ytype, ctype, nnY, nlY, nnC, nlC, neJ, neG, x, hs, &
                       grid, work, usrode, usralg, vbounds, infBnd, &
                       cu, lencu, iu, leniu, ru, lenru, &
                       cw, lencw, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: iDisc, nPhs, nY, nU, nP, nC, &
                                  objL(nPhs), ytype(nY,nPhs), ctype(nC,nPhs), &
                                  nnY(nPhs), nlY(nPhs), nnC(nPhs), nlC(nPhs), &
                                  neJ(nPhs), neG(nPhs), &
                                  lencu, leniu, lenru, lencw, leniw, lenrw
    real(rp),    intent(in)    :: infBnd

    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)

    integer(ip), intent(inout), pointer :: hs(:)
    real(rp),    intent(inout), pointer :: x(:)

    external             :: usrode, usralg, vbounds
    type(ctGrid), target :: grid
    type(ctWork), target :: work

    !===========================================================================
    ! Grid has already been setup.
    ! Sets the workspace, bounds, derivative structures.
    !
    ! 09 Feb 2010: Current version.
    !===========================================================================
    integer(ip) :: nInt, nNodes, n, nDisCon, m, lenY, lenC, nnCon, nlCon, &
                   ldiff, ndiff, ndS, ndE, indS, indE, jt, p, s
    integer(ip), pointer :: ndPtr(:), intPtr(:), phs(:)
    real(rp),    pointer :: step(:), bl(:), bu(:), rc(:), pi(:), px(:)
    real(rp),    allocatable :: clbds(:,:), cubds(:,:)

    nInt   =  grid%nInt
    nNodes =  grid%nNodes
    ndPtr  => grid%ndPtr
    intPtr => grid%intPtr
    step   => grid%step

    allocate ( clbds(nC,nPhs), cubds(nC,nPhs) )

    !---------------------------------------------------------------------------
    ! Compute values for the workspace.
    !   nnCon is the number of nonlinear constraints.
    !   m    is the number of constraints.
    !   n    is the number of variables.
    !   lenY is the (max) size of the Jacobian for the state equations.
    !   lenC is the (max) size of the Jacobian for the algebraic
    !        constraints.

    if ( iDisc == 0 ) then
       ! Variables
       n = (nY+nU)*nNodes + nP

       ! Constraints
       nDisCon = nInt*nY + &  ! discretized state eqns
                 nNodes*nC    ! discretized algebraic constraints

       ! Constraints/Objective
       m = 1 + nDisCon + nY*(nPhs-1)

       ! Derivatives
       lenY = 2*(nY+nU) * nInt*nY + &   ! wrt states/controls
              nP * nInt*nY + &          ! wrt parameters
              2 * (nPhs-1)*nY + &       ! phase continuity
              nPhs                      ! objectives
       lenC = (nY+nU) * nNodes*nC + &   ! wrt states/controls
              nP * nNodes*nC            ! wrt parameters

    else
       ! Variables
       n = (nY+nU)*nNodes + nP

       ! Constraints
       nDisCon = nInt*nU + &    ! linear controls
                 2*nInt*nY + &  ! discretized state eqns
                 (nInt+1)*nC    ! discretized algebraic constraints

       ! Constraints/Objective
       m = 1 + nDisCon + nY*(nPhs-1)

       ! Derivatives
       lenY = 3*(nY+nU) * 2*nInt*nY + &   ! wrt states/controls
              nP * 2*nInt*nY + &          ! wrt parameters
              2 * (nPhs-1)*nY + &         ! phase continuity
              nPhs + &                    ! objectives
              3 * nInt*nU                 ! linear controls

       lenC = (nY+nU) * (nInt+1)*nC + &   ! wrt states/controls
              nP * (nInt+1)*nC            ! wrt parameters

    end if

    ! Nonlinear constraints
    !    Count the number of nonlinear constraints based on user's info on
    !    nonlinear/linear states.
    nnCon = 0
    nlCon = 0

    do p = 1, nPhs
       ldiff = nlY(p)
       ndiff = nnY(p)
       if ( iDisc == 1 ) then
          ldiff = 2*nlY(p) + nU
          ndiff = 2*nnY(p)
       end if

       jt = intPtr(p+1) - intPtr(p)
       nnCon = nnCon + jt*(nnC(p) + ndiff)
       nlCon = nlCon + jt*(nlC(p) + ldiff)
    end do

    nnCon = nnCon + nnC(nPhs)
    nlCon = nlCon + nlC(nPhs)
    nlCon = nlCon + (nPhs-1)*nY

    if (nnCon + nlCon + 1 /= m) then
       call snPRNT (13,'snctrl: INTERNAL setup error.', iw, leniw )
       return
    end if

    ! Set the workspace:
    s = 0
    if ( associated(x) ) then
       s   = size(x)
       px  => x
       phs => hs
       nullify ( x, hs )
    end if

    allocate ( x(n+m), hs(n+m) )
    allocate ( work%bl(n+m), work%bu(n+m), work%rc(n+m), work%pi(m) )
    bl => work%bl
    bu => work%bu
    rc => work%rc
    pi => work%pi

    x  = 0.0
    hs = 0
    bl = -infBnd
    bu =  infBnd
    pi =  0.0
    rc =  0.0

    ! Initialize bounds on discretized constraints.
    bl(1+n:m+n) = 0.0
    bu(1+n:m+n) = 0.0

    ! Set the objective (free) bounds.
    bl(m+n) = -infBnd
    bu(m+n) =  infBnd

    ! Initialize workspace.
    work%n     = n
    work%m     = m
    work%lenJ  = 0
    work%nnCon = nnCon
    work%lenY  = lenY
    work%lenC  = lenC


    !---------------------------------------------------------------------------
    ! Variable bounds:
    if ( iw(iIntF) == 0 .or. iw(iIntF) == 1 ) then   ! S/D versions
       indS = 1
       do p = 1, nPhs
          ndS    = ndPtr(p)
          ndE    = ndPtr(p+1)
          nNodes = ndE - ndS + 1

          indE = indS-1 + nNodes*(nY+nU)

          call vbounds ( p, nPhs, nY, nU, nP, nC, nNodes, &
                         bl(indS:indE), bu(indS:indE), x(indS:indE), &
                         bl(n-nP+1:n), bu(n-nP+1:n), &
                         x(n-nP+1:n), clbds(:,p), cubds(:,p) )
          indS = indE + 1
       end do

    else     ! A version
       call c2vbdsA ( nPhs, nY, nU, nP, nC, ndPtr, n, &
                       bl(1:n), bu(1:n), x(1:n), clbds, cubds, vbounds )

    end if

    ! Set up constraint bounds.
    if ( nC > 0 ) &
         call c2algbds( iDisc, nPhs, nC, nU, ctype, nnY, nlY, &
                        nnCon, intPtr, clbds, cubds, m, &
                        bl(1+n:m+n), bu(1+n:m+n) )


    ! Set up constraint derivative structures (work%lenJ, indJ, locJ, Jcol).
    call c2deriv ( iDisc, iw(iIntf), nPhs, nY, nU, nP, nC, objL, &
                    ytype, ctype, nnY, nlY, nnC, nlC, neJ, neG, &
                    nnCon, ndPtr, step, nInt, work, usrode, usralg, &
                    cu, lencu, iu, leniu, ru, lenru, &
                    cw, lencw, iw, leniw, rw, lenrw )

    ! Copy input x, hs.
    if ( s > 0 ) then
       x(1:s)  = px(1:s)
       hs(1:s) = phs(1:s)
       deallocate ( px, phs )
    end if

    deallocate( clbds, cubds )

  end subroutine c2setup

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c2vbdsA ( nPhs, nY, nU, nP, nC, ndPtr, n, bl, bu, x, &
                      clb, cub, vbounds )
    integer(ip), intent(in)  :: nPhs, nY, nU, nP, nC, n, ndPtr(nPhs+1)
    real(rp),    intent(out) :: bl(n), bu(n), x(n), clb(nC,nPhs), cub(nC,nPhs)
    external     :: vbounds

    !=========================================================================
    ! Set up bounds for states, controls, parameters for the A version.
    ! Return algebraic constraints.
    !
    ! 31 Jan 2009: First version of c2vbdsA.
    ! 06 Feb 2010: Current version.
    !=========================================================================
    integer(ip) :: j, p, indv, nVar, ndS, ndE, nNodes
    real(rp)    :: y0low(nY), y0upp(nY), yflow(nY), yfupp(nY), &
                   ylow(nY),  yupp(nY),  ulow(nU),  uupp(nU)

    nVar = nY+nU
    indv = 0

    do p = 1, nPhs
       ndS    = ndPtr(p)
       ndE    = ndPtr(p+1)
       nNodes = ndE - ndS + 1

       ! Get user-defined bounds
       call vbounds( p, nPhs, nY, nU, nP, nC, y0low, y0upp, &
                     yflow, yfupp, ylow, yupp, ulow, uupp, &
                     bl(n-nP+1:n), bu(n-nP+1:n), x(n-nP+1:n), &
                     clb(:,p), cub(:,p) )

       j = 1
       bl(1+indv:nY+indv) = y0low
       bu(1+indv:nY+indv) = y0upp
       bl(1+indv+nY:nU+indv+nY) = ulow
       bu(1+indv+nY:nU+indv+nY) = uupp
       indv = indv + nVar

       do j = 2, nNodes-1
          bl(1+indv:nY+indv) = ylow
          bu(1+indv:nY+indv) = yupp

          bl(1+indv+nY:nU+indv+nY) = ulow
          bu(1+indv+nY:nU+indv+nY) = uupp
          indv = indv + nVar
       end do

       j = nNodes
       bl(1+indv:nY+indv) = yflow
       bu(1+indv:nY+indv) = yfupp
       bl(1+indv+nY:nU+indv+nY) = ulow
       bu(1+indv+nY:nU+indv+nY) = uupp
       indv = indv + nVar

    end do

  end subroutine c2vbdsA

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c2deriv( iDisc, iIntf, nPhs, nY, nU, nP, nC, objL, &
                      ytype, ctype, nnY, nlY, nnC, nlC, neJ, neG, &
                      nnCon, ndPtr, step, s, work, usrode, usralg, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: iDisc, iIntf, nPhs, nY, nU, nP, nC, s, &
                                  objL(nPhs), ytype(nY,nPhs), ctype(nC,nPhs), &
                                  nnY(nPhs), nlY(nPhs), nnC(nPhs), nlC(nPhs), &
                                  neJ(nPhs), neG(nPhs), nnCon, ndPtr(nPhs+1), &
                                  lencu, leniu, lenru, lencw, leniw, lenrw
    real(rp),    intent(in)    :: step(s)

    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)

    external :: usrode, usralg
    type(ctWork), target :: work

    !=========================================================================
    ! Structures are set up in odeSetup and algSetup.  This routine calls them
    ! and combines the results.
    !
    ! 29 Dec 2008: Updated conSetupTR for v4 of snctrl.
    ! 16 Jan 2009: Implemented sparse data structures for user Jacobian.
    ! 19 Jan 2009: Combined TR/HS versions.
    ! 06 Feb 2010: Current version of c2deriv.
    !=========================================================================
    integer(ip) :: m, n, lenY, lenC, count, neY, neC, nJ, k, j
    integer(ip), allocatable :: ncindJ(:), lcindJ(:), nclocJ(:), lclocJ(:), &
                                nfindJ(:), lfindJ(:), nflocJ(:), lflocJ(:)
    real(rp),    allocatable :: lfJcol(:), lcJcol(:)
    integer(ip), pointer :: locJ(:), indJ(:)
    real(rp),    pointer :: Jcol(:)


    n    = work%n
    m    = work%m
    lenY = work%lenY

    !-------------------------------------------------------------------
    ! ODE constraints
    !-------------------------------------------------------------------
    allocate ( nfindJ(lenY), nflocJ(n+1) )
    allocate ( lfindJ(lenY), lflocJ(n+1), lfJcol(lenY) )

    nfindJ = 0
    nflocJ = 0

    lfindJ = 0
    lflocJ = 0
    lfJcol = 0

    if ( iIntf == 0 ) then   ! S version
       if ( iDisc == 0 ) then
          call odeSetupTR( nPhs, nY, nU, nP, objL, ytype, nnY, nlY, nnC, nlC, &
                           neJ, nnCon, ndPtr, step, s, n, m, &
                           lenY, nfindJ, nflocJ, lfJcol, lfindJ, lflocJ, &
                           odeconS, usrode, &
                           cu, lencu, iu, leniu, ru, lenru, &
                           cw, lencw, iw, leniw, rw, lenrw )

       else
          call odeSetupHS( nPhs, nY, nU, nP, objL, ytype, nnY, nlY, nnC, nlC, &
                           neJ, nnCon, ndPtr, step, s, n, m, &
                           lenY, nfindJ, nflocJ, lfJcol, lfindJ, lflocJ, &
                           odeconS, usrode, &
                           cu, lencu, iu, leniu, ru, lenru, &
                           cw, lencw, iw, leniw, rw, lenrw )

       end if
    else if ( iIntf == 1 ) then   ! D version
       if ( iDisc == 0 ) then
          call odeSetupTR( nPhs, nY, nU, nP, objL, ytype, nnY, nlY, nnC, nlC, &
                           neJ, nnCon, ndPtr, step, s, n, m, &
                           lenY, nfindJ, nflocJ, lfJcol, lfindJ, lflocJ, &
                           odeconD, usrode, &
                           cu, lencu, iu, leniu, ru, lenru, &
                           cw, lencw, iw, leniw, rw, lenrw )

       else
          call odeSetupHS( nPhs, nY, nU, nP, objL, ytype, nnY, nlY, nnC, nlC, &
                           neJ, nnCon, ndPtr, step, s, n, m, &
                           lenY, nfindJ, nflocJ, lfJcol, lfindJ, lflocJ, &
                           odeconD, usrode, &
                           cu, lencu, iu, leniu, ru, lenru, &
                           cw, lencw, iw, leniw, rw, lenrw )

       end if

    else if ( iIntf == 2 ) then   ! A version
       if ( iDisc == 0 ) then
          call odeSetupTR( nPhs, nY, nU, nP, objL, ytype, nnY, nlY, nnC, nlC, &
                           neJ, nnCon, ndPtr, step, s, n, m, &
                           lenY, nfindJ, nflocJ, lfJcol, lfindJ, lflocJ, &
                           odeconA, usrode, &
                           cu, lencu, iu, leniu, ru, lenru, &
                           cw, lencw, iw, leniw, rw, lenrw )

       else
          call odeSetupHS( nPhs, nY, nU, nP, objL, ytype, nnY, nlY, nnC, nlC, &
                           neJ, nnCon, ndPtr, step, s, n, m, &
                           lenY, nfindJ, nflocJ, lfJcol, lfindJ, lflocJ, &
                           odeconA, usrode, &
                           cu, lencu, iu, leniu, ru, lenru, &
                           cw, lencw, iw, leniw, rw, lenrw )

       end if
    end if

    if ( nC == 0 ) then
       neY = nflocJ(n+1) + lflocJ(n+1) - 2

       work%lenJ = neY
       allocate ( work%locJ(n+1), work%Jcol(neY), work%indJ(neY) )

       locJ => work%locJ
       Jcol => work%Jcol
       indJ => work%indJ

       count = 1
       do k = 1, n
          locJ(k) = count

          do j = nflocJ(k), nflocJ(k+1)-1
             indJ(count) = nfindJ(j)
             Jcol(count) = 0.0
             count       = count + 1
          end do

          do j = lflocJ(k), lflocJ(k+1)-1
             indJ(count) = lfindJ(j)
             Jcol(count) = lfJcol(j)
             count       = count + 1
          end do
       end do
       locJ(n+1) = count

       deallocate ( nfindJ, nflocJ, lfindJ, lflocJ, lfJcol )
       return
    end if


    !-------------------------------------------------------------------
    ! Algebraic constraints
    !    Same as above but with the algebraic constraints.
    !-------------------------------------------------------------------
    lenC = work%lenC
    allocate ( ncindJ(lenC), nclocJ(n+1) )
    allocate ( lcindJ(lenC), lclocJ(n+1), lcJcol(lenC) )

    ncindJ = 0
    nclocJ = 0

    lcindJ = 0
    lclocJ = 0
    lcJcol = 0

    if ( iIntf == 0 ) then   ! S version
       if ( iDisc == 0 ) then
          call algSetupTR( nPhs, nY, nU, nP, nC, ctype, nnY, nlY, nnC, nlC, &
                           neG, nnCon, ndPtr, n, lenC, &
                           ncindJ, nclocJ, lcJcol, lcindJ, lclocJ, &
                           algconS, usralg, &
                           cu, lencu, iu, leniu, ru, lenru, &
                           cw, lencw, iw, leniw, rw, lenrw )

       else
          call algSetupHS( nPhs, nY, nU, nP, nC, ctype, nnY, nlY, nnC, nlC, &
                           neG, nnCon, ndPtr, n, lenC, &
                           ncindJ, nclocJ, lcJcol, lcindJ, lclocJ, &
                           algconS, usralg, &
                           cu, lencu, iu, leniu, ru, lenru, &
                           cw, lencw, iw, leniw, rw, lenrw )

       end if
    else if ( iIntf == 1 ) then   ! D version
       if ( iDisc == 0 ) then
          call algSetupTR( nPhs, nY, nU, nP, nC, ctype, nnY, nlY, nnC, nlC, &
                           neG, nnCon, ndPtr, n, lenC, &
                           ncindJ, nclocJ, lcJcol, lcindJ, lclocJ, &
                           algconD, usralg, &
                           cu, lencu, iu, leniu, ru, lenru, &
                           cw, lencw, iw, leniw, rw, lenrw )

       else
          call algSetupHS( nPhs, nY, nU, nP, nC, ctype, nnY, nlY, nnC, nlC, &
                           neG, nnCon, ndPtr, n, lenC, &
                           ncindJ, nclocJ, lcJcol, lcindJ, lclocJ, &
                           algconD, usralg, &
                           cu, lencu, iu, leniu, ru, lenru, &
                           cw, lencw, iw, leniw, rw, lenrw )
       end if

    else if ( iIntf == 2 ) then   ! A version
       if ( iDisc == 0 ) then
          call algSetupTR( nPhs, nY, nU, nP, nC, ctype, nnY, nlY, nnC, nlC, &
                           neG, nnCon, ndPtr, n, lenC, &
                           ncindJ, nclocJ, lcJcol, lcindJ, lclocJ, &
                           algconA, usralg, &
                           cu, lencu, iu, leniu, ru, lenru, &
                           cw, lencw, iw, leniw, rw, lenrw )

       else
          call algSetupHS( nPhs, nY, nU, nP, nC, ctype, nnY, nlY, nnC, nlC, &
                           neG, nnCon, ndPtr, n, lenC, &
                           ncindJ, nclocJ, lcJcol, lcindJ, lclocJ, &
                           algconA, usralg, &
                           cu, lencu, iu, leniu, ru, lenru, &
                           cw, lencw, iw, leniw, rw, lenrw )

       end if
    end if


    !-------------------------------------------------------------------
    ! Combine the linear/nonlinear/ode/alg parts
    !-------------------------------------------------------------------
    neY = nflocJ(n+1) + lflocJ(n+1) - 2
    neC = nclocJ(n+1) + lclocJ(n+1) - 2
    nJ  = neY + neC

    work%lenJ = nJ
    allocate ( work%locJ(n+1), work%Jcol(nJ), work%indJ(nJ) )

    locJ => work%locJ
    Jcol => work%Jcol
    indJ => work%indJ

    count = 1
    do k = 1, n
       locJ(k) = count

       do j = nclocJ(k), nclocJ(k+1)-1
          indJ(count) = ncindJ(j)
          Jcol(count) = 0.0
          count       = count + 1
       end do

       do j = nflocJ(k), nflocJ(k+1)-1
          indJ(count) = nfindJ(j)
          Jcol(count) = 0.0
          count       = count + 1
       end do

       do j = lclocJ(k), lclocJ(k+1)-1
          indJ(count) = lcindJ(j)
          Jcol(count) = lcJcol(j)
          count       = count + 1
       end do

       do j = lflocJ(k), lflocJ(k+1)-1
          indJ(count) = lfindJ(j)
          Jcol(count) = lfJcol(j)
          count       = count + 1
       end do
    end do
    locJ(n+1) = count

    deallocate ( nfindJ, nflocJ, lfindJ, lflocJ, lfJcol )
    deallocate ( ncindJ, nclocJ, lcindJ, lclocJ, lcJcol )

  end subroutine c2deriv

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c2algbds( iDisc, nPhs, nC, nU, ctype, nnY, nlY, nnCon, intPtr, &
                       clbds, cubds, m, bl, bu )
    integer(ip), intent(in) :: iDisc, nPhs, nC, nU, ctype(nC,nPhs), &
                               nnY(nPhs), nlY(nPhs), nnCon, intPtr(nPhs), m
    real(rp),    intent(in) :: clbds(nC, nPhs), cubds(nC, nPhs)
    real(rp),   intent(out) :: bl(m), bu(m)

    !=========================================================================
    ! Set up algebraic constraint bounds.
    !
    ! 04 Jan 2009: First version of conSetBds.
    ! 06 Feb 2010: Current version of c2algbds.
    !=========================================================================
    integer(ip) :: indn, indl, ldiff, ndiff, p, i, j

    indn = 0
    indl = nnCon

    do p = 1, nPhs
       ldiff = nlY(p)
       ndiff = nnY(p)
       if ( iDisc == 1 ) then
          ldiff = 2*nlY(p) + nU
          ndiff = 2*nnY(p)
       end if

       ! alg constraint bds
       do j = intPtr(p), intPtr(p+1)-1
          do i = 1, nC
             if ( ctype(i,p) == 0 ) then ! nonlinear
                indn = indn + 1
                bl(indn) = clbds(i,p)
                bu(indn) = cubds(i,p)
             else ! linear
                indl = indl + 1
                bl(indl) = clbds(i,p)
                bu(indl) = cubds(i,p)
             end if
          end do
          indn = indn + ndiff
          indl = indl + ldiff
       end do
    end do

    ! at tf
    do i = 1, nC
       if ( ctype(i,nPhs) == 0 ) then ! nonlinear
          indn = indn + 1
          bl(indn) = clbds(i,nPhs)
          bu(indn) = cubds(i,nPhs)
       else ! linear
          indl = indl + 1
          bl(indl) = clbds(i,nPhs)
          bu(indl) = cubds(i,nPhs)
       end if
    end do

  end subroutine c2algbds

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c1prnt( prob, grid, x, work, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: leniw, lenrw
    real(rp),    intent(in)    :: x(:)

    integer(ip), intent(inout) :: iw(leniw)
    real(rp),    intent(inout) :: rw(lenrw)

    type(ctProb) :: prob
    type(ctGrid) :: grid
    type(ctWork) :: work

    !===========================================================================
    ! Control option to print control related output to Print file,
    ! using snPRNT.
    !
    ! 21 Aug 2008: Current version of sncPrnt.
    ! 19 Nov 2008: Added refinement limit option print.
    ! 24 Dec 2008: Updated for v4 of snctrl.
    ! 28 Jan 2009: Updated for new phase/sparse struct implementation.
    ! 09 Feb 2010: v5.
    !===========================================================================
    integer(ip) :: nY, nU, nP, nC, nPhs, nNodes
    character   :: str*120, out*3

    nY   = prob%nY
    nU   = prob%nU
    nP   = prob%nP
    nC   = prob%nC
    nPhs = prob%nPhs

    call snPRNT(13, ' -----------------------------', iw, leniw)
    call snPRNT( 3, '     SNCTRL Summary Output    ', iw, leniw)
    call snPRNT( 3, ' -----------------------------', iw, leniw)

    ! Discretization Opt
    if ( iw(iDisc) == 0 ) then
       out = ' TR'
    else
       out = ' HS'
    end if
    write(str, 1000) out
    call snPRNT(3, str, iw, leniw)

    ! Problem Data
    write(str, 1010) nY, nU
    call snPRNT(3, str, iw, leniw)
    write(str, 1020) nP, nC
    call snPRNT(3, str, iw, leniw)
    write(str, 1040) nPhs
    call snPRNT(3, str, iw, leniw)

    nNodes = grid%nNodes - (nPhs-1)
    write(str, 1030) grid%nInt, nNodes
    call snPRNT(3, str, iw, leniw)

    ! Refinement Opts
    if (iw(iRefn) == 0) then
       out = ' No'
    else
       out = 'Yes'
    end if
    call snPRNT(3, '', iw, leniw)
    write(str, 1050) out, rw(ctrefTol)
    call snPRNT(3, str, iw, leniw)
    write(str, 1060) iw(iRefL)
    call snPRNT(3, str, iw, leniw)
    call snPRNT(3, '', iw, leniw)

    if (iw(iIntf) == 0) then
       out = '  S'
    else if ( iw(iIntf) == 1 ) then
       out = '  D'
    else if ( iw(iIntf) == 2 ) then
       out = '  A'
    end if

    ! Internal stats
    write(str, 1100) out
    call snPRNT(3, str, iw, leniw)
    write(str, 1110) iw(ncOde), iw(ncAlg)
    call snPRNT(3, str, iw, leniw)

    call c1psol( nY, nU, nP, nPhs, grid, x, work%n, iw, leniw )

    return

1000 format ('  Discretization               ', a6)
1010 format ('  Number of states             ', i6, 2x, &
             '  Number of controls           ', i10 )
1020 format ('  Number of parameters         ', i6, 2x, &
             '  Number of add. constraints   ', i10 )
1030 format ('  Total number of intervals    ', i6, 2x, &
             '  Total number of nodes        ', i10)
1040 format ('  Number of phases             ', i6 )

1050 format ('  Refinement                   ', a6, 2x, &
             '  Refinement Tolerance       '  , e12.4)
1060 format ('  Refinement Limit             ', i6 )

1100 format ('  Control Interface            ', a6)
1110 format ('  Number of calls to usrode    ', i6, 2x, &
             '  Number of calls to usralg    ', i10 )

  end subroutine c1prnt

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c1psol( nY, nU, nP, nPhs, grid, x, n, iw, leniw )
    integer(ip),  intent(in)    :: nY, nU, nP, nPhs, n, leniw
    real(rp),     intent(in)    :: x(n)
    integer(ip),  intent(inout) :: iw(leniw)
    type(ctGrid), target        :: grid

    !===========================================================================
    ! Prints control interface output to the SNOPT output file.
    !
    ! Prints a column with the time nodes and then 5 columns of the
    ! state/control values.
    !
    ! 21 Aug 2008: Current version of c1psol.
    ! 24 Dec 2008: Updated for v4 of snctrl.
    ! 09 Feb 2010: v5.
    !===========================================================================
    integer(ip) :: nNodes, nVar, pind, tind, ind, vS, vE, which, &
                   i, j, p, c1, c2, ct
    integer(ip),  pointer :: intPtr(:), ndPtr(:)
    real(rp),     pointer :: step(:)
    real(rp), allocatable :: time(:)
    character(120)        :: str


    nNodes =  grid%nNodes
    intPtr => grid%intPtr
    ndPtr  => grid%ndPtr
    step   => grid%step

    allocate( time(nNodes+1) )
    nVar = nY+nU

    ! Set up time array
    tind = 1
    time(tind) = grid%t0

    if ( iw(iDisc) == 0 ) then
       do p = 1, nPhs
          do j = intPtr(p), intPtr(p+1)-1
             tind = tind + 1
             time(tind) = time(tind-1) + step(j)
          end do
          tind = tind + 1
          time(tind) = time(tind-1)
       end do
       time(nNodes+1) = grid%tf

    else
       do p = 1, nPhs
          do j = intPtr(p), intPtr(p+1)-1
             tind = tind + 1
             time(tind) = time(tind-1) + step(j) / 2.0

             tind = tind + 1
             time(tind) = time(tind-1) + step(j) / 2.0
          end do
          tind = tind + 1
          time(tind) = time(tind-1)
       end do
       time(nNodes+1) = grid%tf

    end if


    ! Start with states from the beginning...
    which = 0
    call snPRNT(11, ' State Variables', iw, leniw)
    call snPRNT( 1, ' ---------------', iw, leniw)

    vS  = 1
    vE  = min(5, nY)

100 tind = 0
    ind   = 0

    do p = 1, nPhs
       str = 's'

       do j = ndPtr(p), ndPtr(p+1)
          c1 = 2
          c2 = 16
          tind = tind + 1
          write(str(c1:c2), '(f12.5)') time(tind)

          ct = 0
          do i = vS, vE
             ct = ct + 1
             c1 = 16*ct
             c2 = 16*(ct+1)
             write(str(c1:c2), '(f12.5)') x(i + ind)
          end do
          ind = ind + nVar
          call snPRNT(1, str, iw, leniw)
          str = ' '
       end do

       str = ''
       call snPRNT(1, str, iw, leniw)
    end do


    ! Check where we're at...
    if ( which == 0 ) then
       if ( vE < nY ) then
          call snPRNT(11, ' State Variables', iw, leniw)
          call snPRNT( 1, ' ---------------', iw, leniw)
          vS = vE + 1
          vE = min(5, nY-vE) + vS - 1
          go to 100
       else
          str = ' '
          call snPRNT(11, ' Control Variables', iw, leniw)
          call snPRNT( 1, ' -----------------', iw, leniw)
          which = 1
          vS = 1 + nY
          vE = min(5, nU) + vS - 1
          go to 100
       end if
    end if

    if ( which == 1 ) then
       if ( vE < nU+nY ) then
          call snPRNT(11, ' Control Variables', iw, leniw)
          call snPRNT( 1, ' -----------------', iw, leniw)
          vS = vE + 1
          vE = min(5, nY+nU-vE) + vS - 1
          go to 100
       end if
    end if

   ! Print out parameters if any
   if ( nP > 0 ) then
      pind = n - nP
      call snPRNT(11, ' Parameters     ', iw, leniw)
      call snPRNT( 1, ' ---------------', iw, leniw)
      do i = 1, nP
          write(str, '(f12.5)') x(pind + i)
          call snPRNT(1, str, iw, leniw)
       end do
    end if

    deallocate ( time )

  end subroutine c1psol

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module ct30ker
