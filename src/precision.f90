!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File:  precision.f90
! Definitions for integer and real kinds.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module precision
  use  iso_fortran_env
  implicit none
  public

  integer,   parameter :: ip = int32, rp = real64

  integer(ip), parameter :: iIntf    = 435, &
                            iDisc    = 436, &
                            iCPrt    = 437, &
                            iRefn    = 438, &
                            iRefL    = 439, &
                            ncOde    = 440, &
                            ncAlg    = 441, &
                            ctRefTol = 441
  integer(ip), parameter :: inY      = 442, &
                            inU      = 443, &
                            inP      = 444, &
                            inC      = 445, &
                            inPhs    = 446
end module precision
