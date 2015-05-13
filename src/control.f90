!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File:  control.f90
! Contains general control subroutines for all version of snctrl.
! This is the module which the user includes.
!
! 26 Dec 2008: v4 of the control interface.
! 09 Feb 2010: v5.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module control
  use precision
  use ct02data, only : ctProb

  use ct15usr,  only: iusrodeS, iusralgS, &
                      iusrodeA, iusralgA, &
                      iusrodeD, iusralgD
  use ct30ker,  only : c3ker
  implicit none

  private
  public :: sncTitle, sncInit, sncSpec, snctrlS, snctrlD, snctrlA, sncEnd
  public :: sncSet, sncSetI, sncSetR, sncGet, sncGetC, sncGetI, sncGetR

contains

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sncTitle( title )
    character(30), intent(out) :: title

    !===========================================================================
    title  = 'S N C T R L  5      (Apr 2015)'
    !         123456789|123456789|123456789|

  end subroutine sncTitle

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sncInit( iPrint, iSumm, prob, cw, lencw, iw, leniw, rw, lenrw )
    integer(ip),  intent(in)    :: iPrint, iSumm, lencw, leniw, lenrw
    integer(ip),  intent(inout) :: iw(leniw)
    real(rp),     intent(inout) :: rw(lenrw)
    character(8), intent(inout) :: cw(lencw)
    type(ctProb), target        :: prob

    !===========================================================================
    ! Sets control options to default.  Calls snInit.
    !
    ! 16 Aug 2008: Current version of sncInit.
    ! 18 Nov 2008: Added snctrl title/version print out.
    ! 24 Dec 2008: Updated for v4 of snctrl.
    !===========================================================================
    character(30) :: ttl
    character(30), parameter :: str = '-------------------------------'

    call prob%trash
    call snInit( iPrint, iSumm, cw, lencw, iw, leniw, rw, lenrw )

    ! Default settings.
    iw(iIntf)    =  0
    iw(iDisc)    =  1
    iw(iRefn)    =  1
    iw(iRefL)    = -1
    iw(iCPrt)    =  0
    rw(ctRefTol) =  1.0d-4

    iw(ncOde)    =  0
    iw(ncAlg)    =  0


    ! Print out a title for SNCTRL !
    call sncTitle( ttl )
    call snPRNT( 13, ' '//str, iw, leniw )
    call snPRNT(  3, ' '//ttl, iw, leniw )
    call snPRNT(  3, ' '//str, iw, leniw )

  end subroutine sncInit

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sncSpec( iSpecs, iExit, cw, lencw, iw, leniw, rw, lenrw )
    integer(ip),  intent(in)    :: iSpecs, lencw, leniw, lenrw
    integer(ip),  intent(inout) :: iw(leniw)
    real(rp),     intent(inout) :: rw(lenrw)
    character(8), intent(inout) :: cw(lencw)
    integer(ip),  intent(out)   :: iExit

    !===========================================================================
    ! 09 May 2015: Updated to wrap SNOPT function.
    !===========================================================================

    call snSpec( iSpecs, iExit, cw, lencw, iw, leniw, rw, lenrw )

  end subroutine sncSpec

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sncSet ( buffer, iPrint, iSumm, Errors, &
                      cw, lencw, iw, leniw, rw, lenrw )
    character*(*), intent(in)    :: buffer
    integer(ip),   intent(in)    :: iPrint, iSumm, lencw, leniw, lenrw
    integer(ip),   intent(inout) :: Errors, iw(leniw)
    real(rp),      intent(inout) :: rw(lenrw)
    character(8),  intent(inout) :: cw(lencw)

    !===========================================================================
    ! sncSet  decodes the option contained in  buffer.
    !
    ! The buffer is output to file iPrint, minus trailing blanks.
    ! Error messages are output to files iPrint and iSumm.
    ! Buffer is echoed to iPrint but normally not to iSumm.
    ! It is echoed to iSumm before any error msg.
    !
    ! On entry,
    ! iPrint is the print   file.  no output occurs if iPrint .le 0.
    ! iSumm  is the Summary file.  no output occurs if iSumm  .le 0.
    ! Errors is the number of errors so far.
    !
    ! On exit,
    ! Errors is the number of errors so far.
    !
    ! 27 Jan 2008: Copy of snSet adapted for control interface.
    ! 23 Jul 2008: F90 version.
    ! 09 May 2015: Updated to wrap SNOPT function.
    !===========================================================================

    call snSet( buffer, iPrint, iSumm, Errors, cw, lencw, iw, leniw, rw, lenrw )

  end subroutine sncSet

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sncSeti( buffer, ivalue, iPrint, iSumm, Errors, &
                      cw, lencw, iw, leniw, rw, lenrw )
    character*(*), intent(in)    :: buffer
    integer(ip),   intent(in)    :: ivalue, iPrint, iSumm, lencw, leniw, lenrw
    integer(ip),   intent(inout) :: Errors, iw(leniw)
    real(rp),      intent(inout) :: rw(lenrw)
    character(8),  intent(inout) :: cw(lencw)

    !===================================================================
    ! sncSeti decodes the option contained in  buffer // ivalue.
    ! The parameters other than ivalue are as in snSet.
    !
    ! 27 Jan 2008: Copy of snSeti adapted for control interface.
    ! 23 Jul 2008: F90 version.
    ! 09 May 2015: Updated to wrap SNOPT function.
    !===================================================================

    call snSeti( buffer, ivalue, iPrint, iSumm, Errors, &
                 cw, lencw, iw, leniw, rw, lenrw )

  end subroutine sncSeti

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sncSetr( buffer, rvalue, iPrint, iSumm, Errors, &
                      cw, lencw, iw, leniw, rw, lenrw )
    character*(*), intent(in)    :: buffer
    integer(ip),   intent(in)    :: iPrint, iSumm, lencw, leniw, lenrw
    real(rp),      intent(in)    :: rvalue
    integer(ip),   intent(inout) :: Errors, iw(leniw)
    real(rp),      intent(inout) :: rw(lenrw)
    character(8),  intent(inout) :: cw(lencw)

    !===========================================================================
    ! sncSetr decodes the option contained in  buffer // rvalue.
    ! The parameters other than rvalue are as in snSet.
    !
    ! 27 Jan 2008: Copy of snSetr adapted for control interface.
    ! 23 Jul 2008: F90 version.
    ! 09 May 2015: Updated to wrap SNOPT function.
    !===========================================================================

    call snSetR( buffer, rvalue, iPrint, iSumm, Errors, &
                 cw, lencw, iw, leniw, rw, lenrw )

  end subroutine sncSetr

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  integer function sncGet( buffer, Errors, cw, lencw, iw, leniw, rw, lenrw )
    character*(*), intent(in)    :: buffer
    integer(ip),   intent(in)    :: lencw, leniw, lenrw
    integer(ip),   intent(inout) :: Errors, iw(leniw)
    real(rp),      intent(inout) :: rw(lenrw)
    character(8),  intent(inout) :: cw(lencw)

    !===========================================================================
    ! sncGet  decodes the option contained in  buffer
    ! and returns 1 if the option has previously been set, else 0.
    ! For example,
    ! i = snGet ( 'Maximize', Errors, cw, lencw, iw, leniw, rw, lenrw )
    !
    ! 27 Jan 2008: Copy of snGet adapted for control interface.
    ! 23 Jul 2008: F90 version.
    ! 09 May 2015: Updated to wrap SNOPT function.
    !===========================================================================
    integer :: snGet

    sncGet = snGet(buffer, Errors, cw, lencw, iw, leniw, rw, lenrw )

  end function sncGet

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sncGetc( buffer, cvalue, Errors, cw, lencw, iw, leniw, rw, lenrw )
    character*(*), intent(in)    :: buffer
    integer(ip),   intent(in)    :: lencw, leniw, lenrw
    integer(ip),   intent(inout) :: Errors, iw(leniw)
    real(rp),      intent(inout) :: rw(lenrw)
    character(8),  intent(inout) :: cw(lencw)
    character*(*), intent(out)   :: cvalue

    !===========================================================================
    ! sncGetc gets the value of the option contained in  buffer.
    ! The parameters other than cvalue are as in snSet.
    !
    ! 27 Jan 2008: Copy of snGetc adapted for control interface.
    ! 23 Jul 2008: F90 version.
    ! 09 May 2015: Updated to wrap SNOPT function.
    !===========================================================================
    call snGetC( buffer, cvalue, Errors, cw, lencw, iw, leniw, rw, lenrw )

  end subroutine sncGetc

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sncGeti( buffer, ivalue, Errors, cw, lencw, iw, leniw, rw, lenrw )
    character*(*), intent(in)    :: buffer
    integer(ip),   intent(in)    :: lencw, leniw, lenrw
    integer(ip),   intent(inout) :: Errors, iw(leniw)
    real(rp),      intent(inout) :: rw(lenrw)
    character(8),  intent(inout) :: cw(lencw)
    integer(ip),   intent(out)   :: ivalue

    !===========================================================================
    ! sncGeti gets the value of the option contained in  buffer.
    ! The parameters other than ivalue are as in snSet.
    !
    ! 27 Jan 2008: Copy of snGeti adapted for control interface.
    ! 23 Jul 2008: F90 version.
    ! 09 May 2015: Updated to wrap SNOPT function.
    !===========================================================================
    call snGetI( buffer, ivalue, Errors, cw, lencw, iw, leniw, rw, lenrw )

  end subroutine sncGeti

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sncGetr( buffer, rvalue, Errors, cw, lencw, iw, leniw, rw, lenrw )
    character*(*), intent(in)    :: buffer
    integer(ip),   intent(in)    :: lencw, leniw, lenrw
    integer(ip),   intent(inout) :: Errors, iw(leniw)
    real(rp),      intent(inout) :: rw(lenrw)
    character(8),  intent(inout) :: cw(lencw)
    real(rp),      intent(out)   :: rvalue

    !===========================================================================
    ! sncGetr gets the value of the option contained in  buffer.
    ! The parameters other than rvalue are as in snSet.
    !
    ! 27 Jan 2008: Copy of snGetr adapted for control interface.
    ! 23 Jul 2008: F90 version.
    ! 09 May 2015: Updated to wrap SNOPT function.
    !===========================================================================

    call snGetR( buffer, rvalue, Errors, cw, lencw, iw, leniw, rw, lenrw )

  end subroutine sncGetr

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine snctrlA( Start, prob, usrode, usralg, usrbds, &
                      INFO, mincw, miniw, minrw, nS, nInf, sInf, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )

    integer(ip),  intent(in)    :: Start, lencu, leniu, lenru, &
                                   lencw, leniw, lenrw

    integer(ip),  intent(inout) :: nS, nInf, iu(leniu), iw(leniw)
    real(rp),     intent(inout) :: sInf, rw(lenrw), ru(lenru)
    character(8), intent(inout) :: cw(lencw), cu(lencu)

    integer(ip),  intent(out)   :: INFO, mincw, miniw, minrw

    procedure(iusrodeA) :: usrode
    procedure(iusralgA) :: usralg
    !procedure(iusrvarA) :: usrbds
    external            :: usrbds
    type(ctProb)        :: prob

    !===========================================================================
    ! snctrlA, A version of the optimal control interface.
    !   ( dense version, 1 node )
    !
    ! 31 Jan 2009: Current version snctrlA.
    !===========================================================================
    iw(iIntf) = 2
    call c3ker( Start, prob, usrode, usralg, usrbds, &
                INFO, mincw, miniw, minrw, nS, nInf, sInf, &
                cu, lencu, iu, leniu, ru, lenru, &
                cw, lencw, iw, leniw, rw, lenrw )

  end subroutine snctrlA

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine snctrlD( Start, prob, usrode, usralg, usrbds, &
                      INFO, mincw, miniw, minrw, nS, nInf, sInf, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )

    integer(ip),  intent(in)    :: Start, lencu, leniu, lenru, &
                                   lencw, leniw, lenrw

    integer(ip),  intent(inout) :: nS, nInf, iu(leniu), iw(leniw)
    real(rp),     intent(inout) :: sInf, rw(lenrw), ru(lenru)
    character(8), intent(inout) :: cw(lencw), cu(lencu)

    integer(ip),  intent(out)   :: INFO, mincw, miniw, minrw

    procedure(iusrodeD) :: usrode
    procedure(iusralgD) :: usralg
    !procedure(iusrvarD) :: usrbds
    external            :: usrbds
    type(ctProb)        :: prob

    !===========================================================================
    ! snctrlD, D version of the optimal control interface.
    !   ( dense version, all nodes )
    !
    ! 29 Jan 2008: Current version snctrlD.
    !===========================================================================
    iw(iIntf) = 1
    call c3ker( Start, prob, usrode, usralg, usrbds, &
                INFO, mincw, miniw, minrw, nS, nInf, sInf, &
                cu, lencu, iu, leniu, ru, lenru, &
                cw, lencw, iw, leniw, rw, lenrw )

  end subroutine snctrlD

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine snctrlS( Start, prob, usrode, usralg, usrbds, &
                      INFO, mincw, miniw, minrw, nS, nInf, sInf, &
                      cu, lencu, iu, leniu, ru, lenru, &
                      cw, lencw, iw, leniw, rw, lenrw )

    integer(ip),  intent(in)    :: Start, lencu, leniu, lenru, &
                                   lencw, leniw, lenrw

    integer(ip),  intent(inout) :: nS, nInf, iu(leniu), iw(leniw)
    real(rp),     intent(inout) :: sInf, rw(lenrw), ru(lenru)
    character(8), intent(inout) :: cw(lencw), cu(lencu)

    integer(ip),  intent(out)   :: INFO, mincw, miniw, minrw

    procedure(iusrodeS) :: usrode
    procedure(iusralgS) :: usralg
    !procedure(iusrvarS) :: usrbds
    external            :: usrbds
    type(ctProb)        :: prob

    !===========================================================================
    ! snctrlS, S version of the optimal control interface.
    ! ( sparse structures, user defines elements at every node)
    !
    ! 02 Jan 2008: Current version snctrlS.
    ! 29 Jan 2009: Moved some stuff from wrp here.
    !===========================================================================
    iw(iIntf) = 0
    call c3ker( Start, prob, usrode, usralg, usrbds, &
                INFO, mincw, miniw, minrw, nS, nInf, sInf, &
                cu, lencu, iu, leniu, ru, lenru, &
                cw, lencw, iw, leniw, rw, lenrw )

  end subroutine snctrlS

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  subroutine sncEnd( prob )
    type(ctProb) :: prob

    call prob%trash

  end subroutine sncEnd

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

end module control
