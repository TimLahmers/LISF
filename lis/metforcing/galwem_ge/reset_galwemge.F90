!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.3
!
! Copyright (c) 2020 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
#include "LIS_misc.h"
!BOP
! !MODULE: reset_galwemge
!  \label{reset_galwemge}
!
! !REVISION HISTORY:
! 09 May 2022; Yeosang Yoon, Initial Code
!
! !INTERFACE:
subroutine reset_galwemge()
! !USES:
  use LIS_coreMod,       only : LIS_rc
  use galwemge_forcingMod, only : galwemge_struc
!
! !DESCRIPTION:
!  Routine to reset GALWEM-GE forcing related memory allocations.
!
!EOP
  implicit none

  integer   :: n
  integer   :: findex

  do n=1,LIS_rc%nnest
     galwemge_struc(n)%fcsttime1 = 3000.0
     galwemge_struc(n)%fcsttime2 = 0.0
  enddo

end subroutine reset_galwemge
