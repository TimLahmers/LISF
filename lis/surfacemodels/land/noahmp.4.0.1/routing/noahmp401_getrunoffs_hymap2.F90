!-----------------------BEGIN NOTICE -- DO NOT EDIT-----------------------
! NASA Goddard Space Flight Center
! Land Information System Framework (LISF)
! Version 7.4
!
! Copyright (c) 2022 United States Government as represented by the
! Administrator of the National Aeronautics and Space Administration.
! All Rights Reserved.
!-------------------------END NOTICE -- DO NOT EDIT-----------------------
!BOP
! !ROUTINE: noahmp401_getrunoffs_hymap2
!  \label{noahmp401_getrunoffs_hymap2}
!
! !REVISION HISTORY:
!  6 May 2011: Sujay Kumar; Initial Specification
!
! !INTERFACE:
subroutine noahmp401_getrunoffs_hymap2(n)
! !USES:
  use ESMF
  use LIS_coreMod, only : LIS_rc, LIS_masterproc
  use LIS_routingMod, only : LIS_runoff_state
  use LIS_logMod
  use LIS_historyMod
  use noahmp401_lsmMod, only : noahmp401_struc

  implicit none
! !ARGUMENTS: 
  integer,  intent(in)   :: n 
!
! !DESCRIPTION:
!  
!
! 
!EOP
  type(ESMF_Field)       :: sfrunoff_field
  type(ESMF_Field)       :: qrf_field
  type(ESMF_Field)       :: baseflow_field
  real, pointer          :: sfrunoff(:)
  real, pointer          :: qrf(:)
  real, pointer          :: baseflow(:)
  integer                :: t
  integer                :: c,r
  integer                :: status
  real, allocatable          :: runoff1(:)
  real, allocatable          :: runoff2(:)
  real, allocatable          :: runoff3(:)
  real, allocatable          :: runoff1_t(:)
  real, allocatable          :: runoff2_t(:)
  real, allocatable          :: runoff3_t(:)  

  real, allocatable          :: runoff_lsm(:)

  !ag (25Apr2017)
  type(ESMF_Field)       :: evapotranspiration_Field
  real,pointer           :: evapotranspiration(:)
  real, allocatable      :: evapotranspiration1(:)
  real, allocatable      :: evapotranspiration1_t(:)
  integer                :: evapflag

  character(30)      :: yfile
  character(4)       :: yr
  character(2)       :: mo,da,hr,mn
  character(1)       :: temp

  allocate(runoff1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(runoff2(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(runoff3(LIS_rc%npatch(n,LIS_rc%lsm_index)))
  allocate(runoff1_t(LIS_rc%ntiles(n)))
  allocate(runoff2_t(LIS_rc%ntiles(n)))
  allocate(runoff3_t(LIS_rc%ntiles(n)))

  allocate(runoff_lsm(LIS_rc%ntiles(n)))

  runoff1_t = -9999.0
  runoff2_t = -9999.0
  runoff3_t = -9999.0

  call ESMF_AttributeGet(LIS_runoff_state(n),"Routing model evaporation option",&
       evapflag, rc=status)
!if option is not defined, then assume that no evap calculations will be done
  if(status.ne.0)then 
     evapflag = 0
  endif
  
  call ESMF_StateGet(LIS_runoff_state(n),"Surface Runoff",&
       sfrunoff_field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for Surface Runoff')

  call ESMF_StateGet(LIS_runoff_state(n),"Groundwater River Water Flux",&
       qrf_field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for Groundwater River Water Flux')

  call ESMF_StateGet(LIS_runoff_state(n),"Subsurface Runoff",&
       baseflow_field,rc=status)
  call LIS_verify(status,'ESMF_StateGet failed for Subsurface Runoff')
  
  call ESMF_FieldGet(sfrunoff_field,localDE=0,farrayPtr=sfrunoff,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for Surface Runoff')

  call ESMF_FieldGet(qrf_field,localDE=0,farrayPtr=qrf,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for Groundwater River Water Flux')
  
  call ESMF_FieldGet(baseflow_field,localDE=0,farrayPtr=baseflow,rc=status)
  call LIS_verify(status,'ESMF_FieldGet failed for Subsurface Runoff')

  !remove runoff 3, use if-statement to use qsb or qrf...
  do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)  !units?
     runoff1(t) = NOAHMP401_struc(n)%noahmp401(t)%runsf
     runoff2(t) = NOAHMP401_struc(n)%noahmp401(t)%runsb
     if(NOAHMP401_struc(n)%run_opt .eq. 5) then
         runoff2(t) = 0.0
         runoff3(t) = noahmp401_struc(n)%noahmp401(t)%qrf 
     else
         runoff3(t) = 0.0
     endif
  enddo

! TML: add QRF to NoahMP401 structure and export here (like runsf or runsb)

  runoff1_t = LIS_rc%udef
  runoff2_t = LIS_rc%udef
  runoff3_t = LIS_rc%udef

  call LIS_patch2tile(n,1,runoff1_t, runoff1)
  call LIS_patch2tile(n,1,runoff2_t, runoff2)
  call LIS_patch2tile(n,1,runoff3_t, runoff3)

  sfrunoff = runoff1_t
  baseflow = runoff2_t
  qrf = runoff3_t

! Print HyMAP 2-way-cpl variables to ensure correct dimensions...
#if (defined SPMD)
        if ( LIS_masterproc ) then
        write(yr,'(i4)') LIS_rc%yr
        if (LIS_rc%mo.ge.10) then
            write(mo,'(i2)') LIS_rc%mo
        else
            write(temp,'(i1)') LIS_rc%mo
            mo = "0"//temp
        endif
        if (LIS_rc%da.ge.10) then
            write(da,'(i2)') LIS_rc%da
        else
            write(temp,'(i1)') LIS_rc%da
            da = "0"//temp
        endif
        if (LIS_rc%hr.ge.10) then
            write(hr,'(i2)') LIS_rc%hr
        else
            write(temp,'(i1)') LIS_rc%hr
            hr = "0"//temp
        endif
        if (LIS_rc%mn.ge.10) then
            write(mn,'(i2)') LIS_rc%mn
        else
            write(temp,'(i1)') LIS_rc%mn
            mn = "0"//temp
        endif
        yfile = yr//mo//da//hr//mn//"runoff_lsm_dmp.txt"
        !open(1,file=yfile, form='unformatted')
        !write(1) runoff_lsm
        !close(1)
        endif
#else
        write(yr,'(i4)') LIS_rc%yr
        if (LIS_rc%mo.ge.10) then
            write(mo,'(i2)') LIS_rc%mo
        else
            write(temp,'(i1)') LIS_rc%mo
            mo = "0"//temp
        endif
        if (LIS_rc%da.ge.10) then
            write(da,'(i2)') LIS_rc%da
        else
            write(temp,'(i1)') LIS_rc%da
            da = "0"//temp
        endif
        if (LIS_rc%hr.ge.10) then
            write(hr,'(i2)') LIS_rc%hr
        else
            write(temp,'(i1)') LIS_rc%hr
            hr = "0"//temp
        endif
        if (LIS_rc%mn.ge.10) then
            write(mn,'(i2)') LIS_rc%mn
        else
            write(temp,'(i1)') LIS_rc%mn
            mn = "0"//temp
        endif
        yfile = yr//mo//da//hr//mn//"runoff_lsm_ser.txt"
        !open(1,file=yfile, form='unformatted')
        !write(1) runoff_lsm
        !close(1)
#endif


  deallocate(runoff1)
  deallocate(runoff2)
  deallocate(runoff3)
  deallocate(runoff1_t)
  deallocate(runoff2_t)
  deallocate(runoff3_t)

  deallocate(runoff_lsm)

  !ag (05Jun2017)
  !Including meteorological forcings + evapotranspiration for computing evaporation from open waters in HyMAP2)
  if(evapflag.ne.0)then
    allocate(evapotranspiration1(LIS_rc%npatch(n,LIS_rc%lsm_index)))
    allocate(evapotranspiration1_t(LIS_rc%ntiles(n)))
    
    call ESMF_StateGet(LIS_runoff_state(n),"Total Evapotranspiration",&
         evapotranspiration_Field, rc=status)
    call LIS_verify(status, "noahmp401_getrunoffs_hymap2: ESMF_StateGet failed for Total Evapotranspiration")
        
    call ESMF_FieldGet(evapotranspiration_Field,localDE=0,&
         farrayPtr=evapotranspiration,rc=status)
    call LIS_verify(status, "noahmp401_getrunoffs_hymap2: ESMF_FieldGet failed for Total Evapotranspiration")
    
    do t=1, LIS_rc%npatch(n,LIS_rc%lsm_index)  
       evapotranspiration1(t)  =  &
            NOAHMP401_struc(n)%noahmp401(t)%ecan + &
            NOAHMP401_struc(n)%noahmp401(t)%etran + &
            NOAHMP401_struc(n)%noahmp401(t)%edir 
    enddo

    call LIS_patch2tile(n,1,evapotranspiration1_t, evapotranspiration1)
  
    evapotranspiration = evapotranspiration1_t

    deallocate(evapotranspiration1)
    deallocate(evapotranspiration1_t)
  endif  

end subroutine noahmp401_getrunoffs_hymap2
