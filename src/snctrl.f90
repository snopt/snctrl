!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File:  snctrl.f90
! This is the module that the user includes.
!
! 30 Apr 2015: Added for user to include.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module SNCTRL
  use precision, only : ip, rp
  use ct02data,  only : ctProb
  use control,   only : sncInit, sncSpec, snctrlS, snctrlD, snctrlA, &
                        sncGet, sncGetI, sncGetR, sncSet, sncSetI, sncSetR, &
                        sncEnd
end module SNCTRL
