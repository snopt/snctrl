!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File:  ct02data.f90
! Type definition for structure.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module ct02data
  use precision,    only : ip, rp
  implicit none

  private
  public  :: ctProb
  public  :: ctWork, ctGrid, c0grid, c0gtrsh, c0wtrsh

  ! Control problem data
  type ctProb

     ! General problem data:
     character(8) :: probName

     integer(ip) :: nY,   &  ! number of states
                    nU,   &  ! number of controls
                    nP,   &  ! number of parameters
                    nC,   &  ! number of additional constraints
                    nPhs     ! number of phases

     integer(ip), allocatable :: objL(:),    &  ! location of objective
                                 ytype(:,:), &  ! nonlinear/linear state eqns
                                 ctype(:,:), &  ! nonlinear/linear algebraic constraints
                                 npInt(:)       ! number of intervals in each phase
     real(rp),    allocatable :: phsPt(:)       ! start/end points for each phase

     integer(ip), pointer     :: hs(:)
     real(rp),    pointer     :: x(:)

     ! For SPARSE structures:
     ! ( must be specified if using the sparse version )
     integer(ip), allocatable :: neJ(:), &  ! number of Jacobian elements for odecon
                                 neG(:)     ! number of Jacobian elements for algcon

   contains
     procedure, public, pass :: trash => ctTrash
  end type ctProb

  ! Control/SNOPT workspace
  type ctWork
     integer(ip) :: n, m, nnCon, lenC, lenY, lenJ

     integer(ip), allocatable :: indJ(:), locJ(:)
     real(rp),    allocatable :: Jcol(:), bl(:), bu(:), pi(:), rc(:)
  end type ctWork


  ! Grid structure
  type ctGrid
     integer(ip) :: nInt, nNodes
     real(rp)    :: t0, tf
     integer(ip), allocatable :: ndPtr(:), intPtr(:)
     real(rp),    allocatable :: step(:)
  end type ctGrid

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine ctTrash( prob )
    class(ctProb) :: prob

    prob%probName = ''
    prob%nY   = 0
    prob%nU   = 0
    prob%nP   = 0
    prob%nC   = 0
    prob%nPhs = 0

    if ( associated(prob%hs) ) deallocate(prob%hs)
    if ( associated(prob%x) )  deallocate(prob%x)
    prob%hs => null()
    prob%x  => null()

    if ( allocated(prob%objL) )  deallocate(prob%objL)
    if ( allocated(prob%ytype) ) deallocate(prob%ytype)
    if ( allocated(prob%ctype) ) deallocate(prob%ctype)
    if ( allocated(prob%npInt) ) deallocate(prob%npInt)
    if ( allocated(prob%phsPt) ) deallocate(prob%phsPt)
    if ( allocated(prob%neJ) )   deallocate(prob%neJ)
    if ( allocated(prob%neG) )   deallocate(prob%neG)
  end subroutine ctTrash

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c0wtrsh ( work )
    type(ctWork) :: work

    integer(ip) :: info

    if ( allocated(work%indJ) ) &
         deallocate ( work%indJ, stat=info )

    if ( allocated(work%locJ) ) &
         deallocate ( work%locJ, stat=info )

    if ( allocated(work%Jcol) ) &
         deallocate ( work%Jcol, stat=info )

    if ( allocated(work%bl) ) &
         deallocate ( work%bl, stat=info )

    if ( allocated(work%bu) ) &
         deallocate ( work%bu, stat=info )

    if ( allocated(work%pi) ) &
         deallocate ( work%pi, stat=info )

    if ( allocated(work%rc) ) &
         deallocate ( work%rc, stat=info )

  end subroutine c0wtrsh

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c0gtrsh ( grid )
    type(ctGrid) :: grid

    integer(ip) :: info

    if ( allocated(grid%ndPtr) ) &
         deallocate ( grid%ndPtr,  stat=info )

    if ( allocated(grid%intPtr) ) &
         deallocate ( grid%intPtr, stat=info )

    if ( allocated(grid%step) ) &
         deallocate ( grid%step,   stat=info )

  end subroutine c0gtrsh

   !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine c0grid ( iDisc, grid, nPhs, phsPt, npInt )
    integer(ip),  intent(in) :: iDisc, nPhs, npInt(nPhs)
    real(rp),     intent(in) :: phsPt(nPhs+1)
    type(ctGrid), target     :: grid

    !===========================================================================
    ! Given the number of phases (nPhs), phsPt, and npInt (which are user-
    ! defined), this subroutine creates a new uniform grid.
    !
    ! 08 Feb 2010: Current version.
    !===========================================================================
    integer(ip) :: nInt, nNodes, ind, i, j
    real(rp)    :: hstep
    integer(ip), pointer :: ndPtr(:), intPtr(:)
    real(rp),    pointer :: step(:)

    nInt = sum ( npInt )
    grid%nInt = nInt
    grid%t0   = phsPt(1)
    grid%tf   = phsPt(nPhs+1)

    allocate ( grid%step(nInt), grid%intPtr(nPhs+1), grid%ndPtr(nPhs+1) )
    step   => grid%step
    intPtr => grid%intPtr
    ndPtr  => grid%ndPtr

    ind = 0
    do i = 1, nPhs
       hstep = ( phsPt(i+1) - phsPt(i) ) / npInt(i)
       intPtr(i) = ind + 1
       do j = 1, npInt(i)
          ind = ind + 1
          step(ind) = hstep
       end do
    end do
    intPtr(nPhs+1) = ind + 1

    ! Set up node stuff.
    if ( iDisc == 0 ) then
       nNodes = nInt + 1 + (nPhs-1)
       ndPtr  = intPtr

    else
       nNodes = 2*nInt + 1 + (nPhs-1)
       ndPtr  = 2*intPtr - 1

    end if

    grid%nNodes = nNodes

  end subroutine c0grid

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module ct02data
