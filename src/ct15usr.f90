!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File:   ct15usr.f90
!   Wrappers for calling user-defined subroutines for objectives and constraints
!
! 30 Dec 2008: First version for v4 of snctrl.
! 28 Jan 2009: Sparse structures for user Jacobians implemented.
! 09 Feb 2010: v5.
! 30 Apr 2015: Updated.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module ct15usr
  use precision

  implicit none

  private
  public :: odeconS, algconS, odeconD, algconD, odeconA, algconA
  public :: iusrodeS, iusralgS, iusrodeA, iusralgA, iusrodeD, iusralgD

  interface

     !--------------------------------------------------------------------------

     subroutine iusrodeA( snStat, curPhs, nPhs, nY, nU, nP, F, J, y, u, p, &
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
     end subroutine iusrodeA

     subroutine iusralgA( snStat, curPhs, nPhs, nC, nY, nU, nP, C, G, y, u, p, &
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
     end subroutine iusralgA

     !--------------------------------------------------------------------------

     subroutine iusrodeD( snStat, curPhs, nPhs, nY, nU, nP, nNodes, F, J, &
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
     end subroutine iusrodeD

     subroutine iusralgD( snStat, curPhs, nPhs, nC, nY, nU, nP, nNodes, C, G, &
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

     end subroutine iusralgD

     !--------------------------------------------------------------------------

     subroutine iusrodeS( snStat, curPhs, nPhs, nY, nU, nP, nNodes, F, &
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
     end subroutine iusrodeS

     subroutine iusralgS( snStat, curPhs, nPhs, nC, nY, nU, nP, nNodes, C, &
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
     end subroutine iusralgS

     !--------------------------------------------------------------------------
  end interface

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine odeconS( Status, usrode, curPhs, nPhs, nY, nU, nP, nNodes, &
                      dFlag, F, Jrow, Jval, Jcol, lenJ, dvar, pvar, &
                      needF, needJ, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: Status, curPhs, nPhs, nY, nU, nP, nNodes, &
                                  lenJ, needF, needJ, &
                                  lencu, leniu, lenru, &
                                  lencw, leniw, lenrw
    real(rp),    intent(in)    :: dvar(nY+nU,nNodes), pvar(nP)

    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)

    integer(ip), intent(out)   :: Jrow(lenJ), Jcol(nY+nU+nP+1), dFlag(nY)
    real(rp),    intent(out)   :: F(nY,nNodes), Jval(lenJ,nNodes)
    procedure(iusrodeS)        :: usrode

    !===========================================================================
    ! Sparse version.
    ! This routine just calls the user-defined subroutine returning the state
    ! equations and Jacobian.
    !
    ! 28 Jan 2009: First version of usrodeS.
    ! 09 Feb 2010: v5.
    !===========================================================================
    integer(ip) :: i, k, q

    call usrode( Status, curPhs, nPhs, nY, nU, nP, nNodes, &
                 F, Jrow, Jval, Jcol, lenJ, dvar, pvar, &
                 needF, needJ, cu, lencu, iu, leniu, ru, lenru )
    iw(ncOde) = iw(ncOde) + 1

    ! Check diagonals
    if ( needJ > 0 ) then
       dFlag = 0
       do k = 1, nY
          do q = Jcol(k), Jcol(k+1)-1
             i = Jrow(q)
             if ( i == k ) &
                  dFlag(k) = 1
          end do
       end do
    end if

  end subroutine odeconS

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine algconS( Status, usralg, curPhs, nPhs, nC, nY, nU, nP, nNodes, &
                      C, Grow, Gval, Gcol, lenG, dvar, pvar, &
                      needC, needG, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: Status, curPhs, nPhs, nY, nU, nP, nC, nNodes, &
                                  lenG, needC, needG, &
                                  lencu, leniu, lenru, &
                                  lencw, leniw, lenrw
    real(rp),    intent(in)    :: dvar(nY+nU,nNodes), pvar(nP)

    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)

    integer(ip), intent(out)   :: Grow(lenG), Gcol(nY+nU+nP+1)
    real(rp),    intent(out)   :: C(nC,nNodes), Gval(lenG,nNodes)
    procedure(iusralgS)        :: usralg

    !===========================================================================
    ! Sparse version.
    ! This routine just calls the user-defined subroutine returning the
    ! algebraic constraints and Jacobian.
    !
    ! 28 Jan 2009: First version of usralgS.
    ! 09 Feb 2010: v5.
    !---------------------------------------------------------------------------
    call usralg( Status, curPhs, nPhs, nC, nY, nU, nP, nNodes, &
                 C, Grow, Gval, Gcol, lenG, dvar, pvar, &
                 needC, needG, cu, lencu, iu, leniu, ru, lenru )
    iw(ncAlg) = iw(ncAlg) + 1

  end subroutine algconS

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine odeconD( Status, usrode, curPhs, nPhs, nY, nU, nP, nNodes, &
                      dFlag, F, Jrow, Jval, Jcol, lenJ, dvar, pvar, &
                      needF, needJ, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: Status, curPhs, nPhs, nY, nU, nP, nNodes, &
                                  lenJ, needF, needJ, &
                                  lencu, leniu, lenru, lencw, leniw, lenrw
    real(rp),    intent(in)    :: dvar(nY+nU,nNodes), pvar(nP)

    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)

    integer(ip), intent(out)   :: Jrow(lenJ), Jcol(nY+nU+nP+1), dFlag(nY)
    real(rp),    intent(out)   :: F(nY,nNodes), Jval(lenJ,nNodes)
    procedure(iusrodeD)        :: usrode

    !===========================================================================
    ! Dense version.  Calls user-defined routine and sets up the sparse
    ! structures based on the dense structures provided by the user.
    !
    ! 28 Jan 2009: First version of usrodeD.
    ! 09 Feb 2010: v5.
    !===========================================================================
    integer(ip) :: i, k, jt, neJ
    real(rp)    :: J(nY,nY+nU+nP,nNodes)

    call usrode( Status, curPhs, nPhs, nY, nU, nP, nNodes, F, J, dvar, pvar, &
                 needF, needJ, cu, lencu, iu, leniu, ru, lenru )
    iw(ncOde) = iw(ncOde) + 1

    ! Dense to sparse
    do jt = 1, nNodes
       neJ = 1
       do k = 1, nY+nU+nP
          Jcol(k) = neJ

          do i = 1, nY
             Jrow(neJ)    = i
             Jval(neJ,jt) = J(i,k,jt)
             neJ = neJ + 1
          end do
       end do
    end do
    Jcol(nY+nU+nP+1) = neJ

    ! Diagonal flags
    dFlag = 1

  end subroutine odeconD

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine algconD( Status, usralg, curPhs, nPhs, nC, nY, nU, nP, nNodes, &
                      C, Grow, Gval, Gcol, lenG, dvar, pvar, &
                      needC, needG, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: Status, curPhs, nPhs, nY, nU, nP, nC, nNodes, &
                                  lenG, needC, needG, &
                                  lencu, leniu, lenru, lencw, leniw, lenrw
    real(rp),    intent(in)    :: dvar(nY+nU,nNodes), pvar(nP)

    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)

    integer(ip), intent(out)   :: Grow(lenG), Gcol(nY+nU+nP+1)
    real(rp),    intent(out)   :: C(nC,nNodes), Gval(lenG,nNodes)
    procedure(iusralgD)        :: usralg

    !===========================================================================
    ! Dense version.  Calls user-defined routine and sets up the sparse
    ! structures based on the dense structures provided by the user.
    !
    ! 28 Jan 2009: First version of usralgD.
    ! 09 Feb 2010: v5.
    !===========================================================================
    integer(ip) :: i, k, jt, neG
    real(rp)    :: G(nC,nY+nU+nP,nNodes)

    call usralg( Status, curPhs, nPhs, nC, nY, nU, nP, nNodes, C, G, &
                 dvar, pvar, needC, needG, cu, lencu, iu, leniu, ru, lenru )
    iw(ncAlg) = iw(ncAlg) + 1

    do jt = 1, nNodes
       neG = 1
       do k = 1, nY+nU+nP
          Gcol(k) = neG

          do i = 1, nC
             Grow(neG)    = i
             Gval(neG,jt) = G(i,k,jt)
             neG = neG + 1
          end do
       end do
    end do
    Gcol(nY+nU+nP+1) = neG

  end subroutine algconD

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine odeconA( Status, usrode, curPhs, nPhs, nY, nU, nP, nNodes, &
                      dFlag, F, Jrow, Jval, Jcol, lenJ, dvar, pvar, &
                      needF, needJ, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: Status, curPhs, nPhs, nY, nU, nP, nNodes, &
                                  lenJ, needF, needJ, &
                                  lencu, leniu, lenru, lencw, leniw, lenrw
    real(rp),    intent(in)    :: dvar(nY+nU,nNodes), pvar(nP)

    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)

    integer(ip), intent(out)   :: Jrow(lenJ), Jcol(nY+nU+nP+1), dFlag(nY)
    real(rp),    intent(out)   :: F(nY,nNodes), Jval(lenJ,nNodes)
    procedure(iusrodeA)        :: usrode

    !===========================================================================
    ! A version.  User-defined Jacobians are in dense format and are computed
    ! at a single node.
    !
    ! 31 Jan 2009: First version of usrodeA.
    ! 09 Feb 2010: v5.
    !===========================================================================
    integer(ip) :: i, k, jt, neJ, ind, nVar
    real(rp)    :: J(nY,nY+nU+nP)

    nVar = nY+nU
    ind  = 0

    do jt = 1, nNodes

       call usrode( Status, curPhs, nPhs, nY, nU, nP, F(:,jt), J, &
                    dvar(1:nY,jt), dvar(1+nY:1+nY+nU,jt), pvar, &
                    needF, needJ, cu, lencu, iu, leniu, ru, lenru )
       ind = ind + nVar
       iw(ncOde) = iw(ncOde) + 1

       neJ = 1
       do k = 1, nY+nU+nP
          Jcol(k) = neJ

          do i = 1, nY
             Jrow(neJ)    = i
             Jval(neJ,jt) = J(i,k)
             neJ = neJ + 1
          end do
       end do
    end do
    Jcol(nY+nU+nP+1) = neJ

    ! Diagonal flags
    dFlag = 1

  end subroutine odeconA

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine algconA( Status, usralg, curPhs, nPhs, nC, nY, nU, nP, nNodes, &
                      C, Grow, Gval, Gcol, lenG, dvar, pvar, &
                      needC, needG, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )
    integer(ip), intent(in)    :: Status, curPhs, nPhs, nY, nU, nP, nC, nNodes, &
                                  lenG, needC, needG, &
                                  lencu, leniu, lenru, lencw, leniw, lenrw
    real(rp),    intent(in)    :: dvar(nY+nU,nNodes), pvar(nP)

    integer(ip), intent(inout) :: iu(leniu), iw(leniw)
    real(rp),    intent(inout) :: ru(lenru), rw(lenrw)
    character(8),intent(inout) :: cu(lencu), cw(lencw)

    integer(ip), intent(out)   :: Grow(lenG), Gcol(nY+nU+nP+1)
    real(rp),    intent(out)   :: C(nC,nNodes), Gval(lenG,nNodes)
    procedure(iusralgA)        :: usralg

    !===========================================================================
    ! A version.  User-defined Jacobians are in dense format and are computed
    ! at a single node.
    !
    ! 31 Jan 2009: First version of usralgA.
    ! 09 Feb 2010: v5.
    !===========================================================================
    integer(ip) :: i, k, jt, neG, ind, nVar
    real(rp)    :: G(nC,nY+nU+nP)

    nVar = nY+nU
    ind  = 0

    do jt = 1, nNodes

       call usralg( Status, curPhs, nPhs, nC, nY, nU, nP, C(:,jt), G, &
                    dvar(1:nY,jt), dvar(1+nY:1+nY+nU,jt), pvar, &
                    needC, needG, cu, lencu, iu, leniu, ru, lenru )
       ind = ind + nVar
       iw(ncAlg) = iw(ncAlg) + 1

       neG = 1
       do k = 1, nY+nU+nP
          Gcol(k) = neG

          do i = 1, nC
             Grow(neG)    = i
             Gval(neG,jt) = G(i,k)
             neG = neG + 1
          end do
       end do
    end do
    Gcol(nY+nU+nP+1) = neG

  end subroutine algconA

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module ct15usr
