subroutine mmf_start(n)
    use LIS_coreMod
    use NoahMP401_lsmMod
    use module_sf_noahmpdrv_401
    use LIS_historyMod, only: LIS_gather_masterproc_2d_local_to_global, &
                              LIS_scatter_global_to_local_grid
    use LIS_mpiMod
    use LIS_logMod, only     : LIS_logunit
    implicit none 
    integer :: i,j
    integer :: n, row, col, t, ridx, cidx, ierr
    real :: wtddt 
   ! SW, MMF 
    integer, allocatable,dimension(:,:) :: isltyp, ivgtyp
    !real, allocatable :: fdepth(:,:)
    real, allocatable,dimension(:,:) ::  fdepth, topo , area, rechclim, rivercond, &
                                        wtd, riverbed, eqwtd, pexp, smcwtdxy, &
                                        deeprechxy, rechxy, qslatxy, qrfsxy, qspringsxy  
    real, allocatable,dimension(:,:,:) :: smois, sh2o, smoiseq, tslb
#if (defined SPMD)
    integer, allocatable, dimension(:,:) :: gisltyp, givgtyp
    real, allocatable, dimension(:,:) :: gfdepth, gtopo , garea, grechclim, grivercond, &
                                         gwtd, griverbed, geqwtd, gpexp, gsmcwtdxy, &
                                         gdeeprechxy, grechxy, gqslatxy, gqrfsxy, gqspringsxy
    real, allocatable,dimension(:,:,:) :: gsmois, gsh2o, gsmoiseq, gtslb
#endif
    wtddt = int(LIS_rc%ts/60) ! in minutes? 

    allocate(isltyp(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(ivgtyp(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(fdepth(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(topo(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(area(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(rechclim(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(rivercond(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(wtd(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(riverbed(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(eqwtd(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(pexp(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(smcwtdxy(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(deeprechxy(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(rechxy(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(qslatxy(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(qrfsxy(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(qspringsxy(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(smois(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, 1:NOAHMP401_struc(n)%nsoil, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(tslb(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, 1:NOAHMP401_struc(n)%nsoil, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(sh2o(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, 1:NOAHMP401_struc(n)%nsoil, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    allocate(smoiseq(NOAHMP401_struc(n)%col_min:NOAHMP401_struc(n)%col_max, 1:NOAHMP401_struc(n)%nsoil, NOAHMP401_struc(n)%row_min:NOAHMP401_struc(n)%row_max))
    
    do row=NOAHMP401_struc(n)%row_min, NOAHMP401_struc(n)%row_max
        do col=NOAHMP401_struc(n)%col_min, NOAHMP401_struc(n)%col_max
            t = NOAHMP401_struc(n)%rct_idx(col,row)  ! rct_idx is col x row TML
            if(t .ne. LIS_rc%udef) then
                isltyp(col,row) = NOAHMP401_struc(n)%noahmp401(t)%soiltype
                ivgtyp(col,row) = NOAHMP401_struc(n)%noahmp401(t)%vegetype
            else
                isltyp(col,row) = NOAHMP401_struc(n)%soil2d(col,row) ! soil2d is col x row TML
                ivgtyp(col,row) = NOAHMP401_struc(n)%vege2d(col,row) ! vege2d is col x row TML
            endif
            smois(col,:,row) = NOAHMP401_struc(n)%init_smc(:)
            sh2o(col,:,row) = NOAHMP401_struc(n)%init_smc(:)
            tslb(col,:,row) = NOAHMP401_struc(n)%init_tslb(:)
            smoiseq(col,:,row) = 0.0
        enddo
    enddo
    !!! 2-D, MMF, SW
    do row=NOAHMP401_struc(n)%row_min, NOAHMP401_struc(n)%row_max
        do col=NOAHMP401_struc(n)%col_min, NOAHMP401_struc(n)%col_max
            ridx = row - NOAHMP401_struc(n)%row_min + 1
            cidx = col - NOAHMP401_struc(n)%col_min + 1
            fdepth(col,row)    = NOAHMP401_struc(n)%fdepth(cidx, ridx)
            topo(col,row)      = NOAHMP401_struc(n)%topo(cidx, ridx)
            area(col,row)      = NOAHMP401_struc(n)%area(cidx, ridx)
            riverbed(col,row)  = NOAHMP401_struc(n)%riverbed(cidx, ridx)
            eqwtd(col,row)     = NOAHMP401_struc(n)%eqwtd(cidx, ridx)
            rivercond(col,row) = NOAHMP401_struc(n)%rivercond(cidx, ridx)
            rechclim(col,row)  = NOAHMP401_struc(n)%rechclim(cidx, ridx)
            pexp(col,row)      = 1.0
        enddo
    enddo

    !print*, 'mmf_start INIT VARIABLES: '
    !print*, 'SMC-1 = ',smois(18,1,12)
    !print*, 'SMC-2 = ',smois(18,2,12)
    !print*, 'SMC-3 = ',smois(18,3,12)
    !print*, 'SMC-4 = ',smois(18,4,12)
    !print*, 'SMCWTD = ',smcwtdxy(18,12)
    !print*, 'EQWTD = ',eqwtd(18,12)
    !print*, 'WTD = ',wtd(18,12)

#if (defined SPMD)
    if ( LIS_masterproc ) then
        ! Allocate global in and inout variables
        write(LIS_logunit,*) "[INFO] Allocating Arrays for MMF Init."
        allocate(gisltyp(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(givgtyp(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(gfdepth(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(gtopo(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(griverbed(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(geqwtd(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(gpexp(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(garea(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(gwtd(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(gsmcwtdxy(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(gdeeprechxy(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(grechxy(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(gqslatxy(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(gqrfsxy(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(gqspringsxy(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(grechclim(LIS_rc%gnc(n), LIS_rc%gnr(n)))
        allocate(gsmois(LIS_rc%gnc(n), NOAHMP401_struc(n)%nsoil, LIS_rc%gnr(n)))
        allocate(gtslb(LIS_rc%gnc(n), NOAHMP401_struc(n)%nsoil, LIS_rc%gnr(n)))
        allocate(gsh2o(LIS_rc%gnc(n), NOAHMP401_struc(n)%nsoil, LIS_rc%gnr(n)))
        allocate(gsmoiseq(LIS_rc%gnc(n), NOAHMP401_struc(n)%nsoil, LIS_rc%gnr(n)))
    else
        ! Allocate dummy "global" in and inout variables
        write(LIS_logunit,*) "[INFO] Allocating Dummy Arrays for MMF Init."
        allocate(gisltyp(1,1))
        allocate(givgtyp(1,1))
        allocate(gfdepth(1,1))
        allocate(gtopo(1,1))
        allocate(griverbed(1,1))
        allocate(geqwtd(1,1))
        allocate(gpexp(1,1))
        allocate(garea(1,1))
        allocate(gwtd(1,1))
        allocate(gsmcwtdxy(1,1))
        allocate(gdeeprechxy(1,1))
        allocate(grechxy(1,1))
        allocate(gqslatxy(1,1))
        allocate(gqrfsxy(1,1))
        allocate(gqspringsxy(1,1))
        allocate(grechclim(1,1))
        allocate(gsmois(1,1,1))
        allocate(gtslb(1,1,1))
        allocate(gsh2o(1,1,1))
        allocate(gsmoiseq(1,1,1))
    endif

    !TML ADD GATHER CALL HERE...
    ! Gather in and inout variables
    write(LIS_logunit,*) "[INFO] Gathering Arrays for MMF Init."
    write(LIS_logunit,*) "[INFO] Gathering Array: isltyp"
    call LIS_gather_masterproc_2d_local_to_global(n, isltyp, gisltyp)
    write(LIS_logunit,*) "[INFO] Gathering Array: ivgtyp"
    call LIS_gather_masterproc_2d_local_to_global(n, ivgtyp, givgtyp)
    write(LIS_logunit,*) "[INFO] Gathering Array: fdepth"
    call LIS_gather_masterproc_2d_local_to_global(n, fdepth, gfdepth)
    write(LIS_logunit,*) "[INFO] Gathering Array: topo"
    call LIS_gather_masterproc_2d_local_to_global(n, topo, gtopo)
    call LIS_gather_masterproc_2d_local_to_global(n, riverbed, griverbed)
    call LIS_gather_masterproc_2d_local_to_global(n, eqwtd, geqwtd)
    call LIS_gather_masterproc_2d_local_to_global(n, pexp, gpexp)
    call LIS_gather_masterproc_2d_local_to_global(n, area, garea)
    call LIS_gather_masterproc_2d_local_to_global(n, wtd, gwtd)
    call LIS_gather_masterproc_2d_local_to_global(n, smcwtdxy, gsmcwtdxy)
    call LIS_gather_masterproc_2d_local_to_global(n, deeprechxy, gdeeprechxy)
    call LIS_gather_masterproc_2d_local_to_global(n, rechxy, grechxy)
    call LIS_gather_masterproc_2d_local_to_global(n, qslatxy, gqslatxy)
    call LIS_gather_masterproc_2d_local_to_global(n, qrfsxy, gqrfsxy)
    call LIS_gather_masterproc_2d_local_to_global(n, qspringsxy, gqspringsxy)
    call LIS_gather_masterproc_2d_local_to_global(n, rechclim, grechclim)
    call LIS_gather_masterproc_2d_local_to_global(n, smois(:,1,:), gsmois(:,1,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smois(:,2,:), gsmois(:,2,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smois(:,3,:), gsmois(:,3,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smois(:,4,:), gsmois(:,4,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smois(:,5,:), gsmois(:,5,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smois(:,6,:), gsmois(:,6,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smois(:,7,:), gsmois(:,7,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smois(:,8,:), gsmois(:,8,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smois(:,9,:), gsmois(:,9,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smois(:,10,:), gsmois(:,10,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smois(:,11,:), gsmois(:,11,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smois(:,12,:), gsmois(:,12,:))
    call LIS_gather_masterproc_2d_local_to_global(n, sh2o(:,1,:), gsh2o(:,1,:))
    call LIS_gather_masterproc_2d_local_to_global(n, sh2o(:,2,:), gsh2o(:,2,:))
    call LIS_gather_masterproc_2d_local_to_global(n, sh2o(:,3,:), gsh2o(:,3,:))
    call LIS_gather_masterproc_2d_local_to_global(n, sh2o(:,4,:), gsh2o(:,4,:))
    call LIS_gather_masterproc_2d_local_to_global(n, sh2o(:,5,:), gsh2o(:,5,:))
    call LIS_gather_masterproc_2d_local_to_global(n, sh2o(:,6,:), gsh2o(:,6,:))
    call LIS_gather_masterproc_2d_local_to_global(n, sh2o(:,7,:), gsh2o(:,7,:))
    call LIS_gather_masterproc_2d_local_to_global(n, sh2o(:,8,:), gsh2o(:,8,:))
    call LIS_gather_masterproc_2d_local_to_global(n, sh2o(:,9,:), gsh2o(:,9,:))
    call LIS_gather_masterproc_2d_local_to_global(n, sh2o(:,10,:), gsh2o(:,10,:))
    call LIS_gather_masterproc_2d_local_to_global(n, sh2o(:,11,:), gsh2o(:,11,:))
    call LIS_gather_masterproc_2d_local_to_global(n, sh2o(:,12,:), gsh2o(:,12,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smoiseq(:,1,:), gsmoiseq(:,1,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smoiseq(:,2,:), gsmoiseq(:,2,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smoiseq(:,3,:), gsmoiseq(:,3,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smoiseq(:,4,:), gsmoiseq(:,4,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smoiseq(:,5,:), gsmoiseq(:,5,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smoiseq(:,6,:), gsmoiseq(:,6,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smoiseq(:,7,:), gsmoiseq(:,7,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smoiseq(:,8,:), gsmoiseq(:,8,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smoiseq(:,9,:), gsmoiseq(:,9,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smoiseq(:,10,:), gsmoiseq(:,10,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smoiseq(:,11,:), gsmoiseq(:,11,:))
    call LIS_gather_masterproc_2d_local_to_global(n, smoiseq(:,12,:), gsmoiseq(:,12,:))
    call LIS_gather_masterproc_2d_local_to_global(n, tslb(:,1,:), gtslb(:,1,:))
    call LIS_gather_masterproc_2d_local_to_global(n, tslb(:,2,:), gtslb(:,2,:))
    call LIS_gather_masterproc_2d_local_to_global(n, tslb(:,3,:), gtslb(:,3,:))
    call LIS_gather_masterproc_2d_local_to_global(n, tslb(:,4,:), gtslb(:,4,:))
    call LIS_gather_masterproc_2d_local_to_global(n, tslb(:,5,:), gtslb(:,5,:))
    call LIS_gather_masterproc_2d_local_to_global(n, tslb(:,6,:), gtslb(:,6,:))
    call LIS_gather_masterproc_2d_local_to_global(n, tslb(:,7,:), gtslb(:,7,:))
    call LIS_gather_masterproc_2d_local_to_global(n, tslb(:,8,:), gtslb(:,8,:))
    call LIS_gather_masterproc_2d_local_to_global(n, tslb(:,9,:), gtslb(:,9,:))
    call LIS_gather_masterproc_2d_local_to_global(n, tslb(:,10,:), gtslb(:,10,:))
    call LIS_gather_masterproc_2d_local_to_global(n, tslb(:,11,:), gtslb(:,11,:))
    call LIS_gather_masterproc_2d_local_to_global(n, tslb(:,12,:), gtslb(:,12,:))

    if ( LIS_masterproc ) then
        ! Allocate out variables
        write(LIS_logunit,*) "[INFO] Allocating Output Arrays for MMF Init."
        allocate(grivercond(LIS_rc%gnc(n), LIS_rc%gnr(n)))
    else
        ! Allocate dummy "global" output variables
        write(LIS_logunit,*) "[INFO] Allocating Dummy Output Arrays for MMF Init."
        allocate(grivercond(1,1))
    endif

    if ( LIS_masterproc ) then
    write(LIS_logunit,*) "[INFO] Initialization Groundwater for MMF."

    !write(*,*) "mmf_start GROUNDWATER VARIABLES: "
    !do i = 50,60
    !    write(*,*) "EQWTD(",i,",35) = ",geqwtd(i,35)
    !enddo
    !do j = 30,40
    !    write(*,*) "EQWTD(55,",j,") = ",geqwtd(55,j)
    !enddo
    !do i = 50,60
    !    write(*,*) "ISLTYP(",i,",35) = ",gisltyp(i,35)
    !enddo
    !do j = 30,40
    !    write(*,*) "ISLTYP(55,",j,") = ",gisltyp(55,j)
    !enddo

    call  groundwater_init (noahmp401_struc(n)%nsoil,  & !nsoil ,
                            noahmp401_struc(n)%sldpth, & !dzs,
                            gisltyp, givgtyp, wtddt ,    &
                            gfdepth, gtopo, griverbed, geqwtd, grivercond, gpexp , garea ,gwtd , &
                            gsmois,gsh2o, gsmoiseq, gsmcwtdxy, gdeeprechxy, grechxy ,  &
                            gqslatxy, gqrfsxy, gqspringsxy,                  &
                            grechclim  ,  gtslb ,  &
                            1,             & !ids,
                            LIS_rc%gnc(n), & !ide, +1 for test
                            1,             & !jds,
                            LIS_rc%gnr(n), & !jde,
                            1, 1,          & !kds,kde,
                            1,             & !ims,
                            LIS_rc%gnc(n), & !ime,
                            1,             & !jms,
                            LIS_rc%gnr(n), & !jme,
                            1, 1,          & !kms,kme,
                            1,             & !ips,
                            LIS_rc%gnc(n), & !ipe,
                            1,             & !jps,
                            LIS_rc%gnr(n), & !jpe,
                            1,1,           & !kps,kpe,
                            1,             & !its,
                            LIS_rc%gnc(n), & !ite,
                            1,             & !jts,
                            LIS_rc%gnr(n), & !jte,
                            1,1)             !kts,kte
    endif
    write(LIS_logunit,*) "[INFO] Waiting for MMF Processors."
    call MPI_Barrier(LIS_MPI_COMM, ierr)

    ! Scatter out and inout variables
    write(LIS_logunit,*) "[INFO] Scattering MMF Init. Arrays."
    call LIS_scatter_global_to_local_grid(n, gpexp, pexp)
    call LIS_scatter_global_to_local_grid(n, geqwtd, eqwtd)
    call LIS_scatter_global_to_local_grid(n, griverbed, riverbed)
    call LIS_scatter_global_to_local_grid(n, gwtd, wtd)
    call LIS_scatter_global_to_local_grid(n, gsmcwtdxy, smcwtdxy)
    call LIS_scatter_global_to_local_grid(n, gdeeprechxy, deeprechxy)
    call LIS_scatter_global_to_local_grid(n, grechxy, rechxy)
    call LIS_scatter_global_to_local_grid(n, gqslatxy, qslatxy)
    call LIS_scatter_global_to_local_grid(n, gqrfsxy, qrfsxy)
    call LIS_scatter_global_to_local_grid(n, gqspringsxy, qspringsxy)
    call LIS_scatter_global_to_local_grid(n, gsmois(:,1,:), smois(:,1,:))
    call LIS_scatter_global_to_local_grid(n, gsmois(:,2,:), smois(:,2,:))
    call LIS_scatter_global_to_local_grid(n, gsmois(:,3,:), smois(:,3,:))
    call LIS_scatter_global_to_local_grid(n, gsmois(:,4,:), smois(:,4,:))
    call LIS_scatter_global_to_local_grid(n, gsmois(:,5,:), smois(:,5,:))
    call LIS_scatter_global_to_local_grid(n, gsmois(:,6,:), smois(:,6,:))
    call LIS_scatter_global_to_local_grid(n, gsmois(:,7,:), smois(:,7,:))
    call LIS_scatter_global_to_local_grid(n, gsmois(:,8,:), smois(:,8,:))
    call LIS_scatter_global_to_local_grid(n, gsmois(:,9,:), smois(:,9,:))
    call LIS_scatter_global_to_local_grid(n, gsmois(:,10,:), smois(:,10,:))
    call LIS_scatter_global_to_local_grid(n, gsmois(:,11,:), smois(:,11,:))
    call LIS_scatter_global_to_local_grid(n, gsmois(:,12,:), smois(:,12,:))
    call LIS_scatter_global_to_local_grid(n, gsh2o(:,1,:), sh2o(:,1,:))
    call LIS_scatter_global_to_local_grid(n, gsh2o(:,2,:), sh2o(:,2,:))
    call LIS_scatter_global_to_local_grid(n, gsh2o(:,3,:), sh2o(:,3,:))
    call LIS_scatter_global_to_local_grid(n, gsh2o(:,4,:), sh2o(:,4,:))
    call LIS_scatter_global_to_local_grid(n, gsh2o(:,5,:), sh2o(:,5,:))
    call LIS_scatter_global_to_local_grid(n, gsh2o(:,6,:), sh2o(:,6,:))
    call LIS_scatter_global_to_local_grid(n, gsh2o(:,7,:), sh2o(:,7,:))
    call LIS_scatter_global_to_local_grid(n, gsh2o(:,8,:), sh2o(:,8,:))
    call LIS_scatter_global_to_local_grid(n, gsh2o(:,9,:), sh2o(:,9,:))
    call LIS_scatter_global_to_local_grid(n, gsh2o(:,10,:), sh2o(:,10,:))
    call LIS_scatter_global_to_local_grid(n, gsh2o(:,11,:), sh2o(:,11,:))
    call LIS_scatter_global_to_local_grid(n, gsh2o(:,12,:), sh2o(:,12,:))
    call LIS_scatter_global_to_local_grid(n, gsmoiseq(:,1,:), smoiseq(:,1,:))
    call LIS_scatter_global_to_local_grid(n, gsmoiseq(:,2,:), smoiseq(:,2,:))
    call LIS_scatter_global_to_local_grid(n, gsmoiseq(:,3,:), smoiseq(:,3,:))
    call LIS_scatter_global_to_local_grid(n, gsmoiseq(:,4,:), smoiseq(:,4,:))
    call LIS_scatter_global_to_local_grid(n, gsmoiseq(:,5,:), smoiseq(:,5,:))
    call LIS_scatter_global_to_local_grid(n, gsmoiseq(:,6,:), smoiseq(:,6,:))
    call LIS_scatter_global_to_local_grid(n, gsmoiseq(:,7,:), smoiseq(:,7,:))
    call LIS_scatter_global_to_local_grid(n, gsmoiseq(:,8,:), smoiseq(:,8,:))
    call LIS_scatter_global_to_local_grid(n, gsmoiseq(:,9,:), smoiseq(:,9,:))
    call LIS_scatter_global_to_local_grid(n, gsmoiseq(:,10,:), smoiseq(:,10,:))
    call LIS_scatter_global_to_local_grid(n, gsmoiseq(:,11,:), smoiseq(:,11,:))
    call LIS_scatter_global_to_local_grid(n, gsmoiseq(:,12,:), smoiseq(:,12,:))
    call LIS_scatter_global_to_local_grid(n, gtslb(:,1,:), tslb(:,1,:))
    call LIS_scatter_global_to_local_grid(n, gtslb(:,2,:), tslb(:,2,:))
    call LIS_scatter_global_to_local_grid(n, gtslb(:,3,:), tslb(:,3,:))
    call LIS_scatter_global_to_local_grid(n, gtslb(:,4,:), tslb(:,4,:))
    call LIS_scatter_global_to_local_grid(n, gtslb(:,5,:), tslb(:,5,:))
    call LIS_scatter_global_to_local_grid(n, gtslb(:,6,:), tslb(:,6,:))
    call LIS_scatter_global_to_local_grid(n, gtslb(:,7,:), tslb(:,7,:))
    call LIS_scatter_global_to_local_grid(n, gtslb(:,8,:), tslb(:,8,:))
    call LIS_scatter_global_to_local_grid(n, gtslb(:,9,:), tslb(:,9,:))
    call LIS_scatter_global_to_local_grid(n, gtslb(:,10,:), tslb(:,10,:))
    call LIS_scatter_global_to_local_grid(n, gtslb(:,11,:), tslb(:,11,:))
    call LIS_scatter_global_to_local_grid(n, gtslb(:,12,:), tslb(:,12,:))
    call LIS_scatter_global_to_local_grid(n, grivercond, rivercond)

    ! Deallocate temporary global variables
    deallocate(gisltyp)
    deallocate(givgtyp)
    deallocate(gfdepth)
    deallocate(gtopo)
    deallocate(griverbed)
    deallocate(geqwtd)
    deallocate(grivercond)
    deallocate(gpexp)
    deallocate(garea)
    deallocate(gwtd)
    deallocate(gsmcwtdxy)
    deallocate(gdeeprechxy)
    deallocate(grechxy)
    deallocate(gqslatxy)
    deallocate(gqrfsxy)
    deallocate(gqspringsxy)
    deallocate(grechclim)
    deallocate(gsmois)
    deallocate(gsh2o)
    deallocate(gsmoiseq)
    deallocate(gtslb)

#else
    
    !write(*,*) "mmf_start GROUNDWATER VARIABLES: "
    !do i = 50,60
    !    write(*,*) "EQWTD(",i,",35) = ",eqwtd(i,35)
    !enddo
    !do j = 30,40
    !    write(*,*) "EQWTD(55,",j,") = ",eqwtd(55,j)
    !enddo
    !do i = 50,60
    !    write(*,*) "ISLTYP(",i,",35) = ",isltyp(i,35)
    !enddo
    !do j = 30,40
    !    write(*,*) "ISLTYP(55,",j,") = ",isltyp(55,j)
    !enddo

    call  groundwater_init (noahmp401_struc(n)%nsoil,  & !nsoil ,
                            noahmp401_struc(n)%sldpth, & !dzs, 
                            isltyp, ivgtyp, wtddt ,    &
                            fdepth, topo, riverbed, eqwtd, rivercond, pexp , area ,wtd ,  &
                            smois,sh2o, smoiseq, smcwtdxy, deeprechxy, rechxy ,  &
                            qslatxy, qrfsxy, qspringsxy,                  &
                            rechclim  , tslb  ,    &
                            NOAHMP401_struc(n)%col_min, & !ids,
                            NOAHMP401_struc(n)%col_max, & !ide, +1 for test
                            NOAHMP401_struc(n)%row_min, & !jds,
                            NOAHMP401_struc(n)%row_max, & !jde, 
                            1, 1,                       & !kds,kde,
                            NOAHMP401_struc(n)%col_min, & !ims,
                            NOAHMP401_struc(n)%col_max, & !ime, 
                            NOAHMP401_struc(n)%row_min, & !jms,
                            NOAHMP401_struc(n)%row_max, & !jme, 
                            1, 1,                       & !kms,kme,
                            NOAHMP401_struc(n)%col_min, & !ips,
                            NOAHMP401_struc(n)%col_max, & !ipe, 
                            NOAHMP401_struc(n)%row_min, & !jps,
                            NOAHMP401_struc(n)%row_max, & !jpe, 
                            1,1,                        & !kps,kpe,
                            NOAHMP401_struc(n)%col_min, & !its,
                            NOAHMP401_struc(n)%col_max, & !ite, 
                            NOAHMP401_struc(n)%row_min, & !jts,
                            NOAHMP401_struc(n)%row_max, & !jte, 
                            1,1)                          !kts,kte
#endif
    !print*, 'mmf_start GROUNDWATER VARIABLES: '
    !print*, 'SMC-1 = ',smois(18,1,12)
    !print*, 'SMC-2 = ',smois(18,2,12)
    !print*, 'SMC-3 = ',smois(18,3,12)
    !print*, 'SMC-4 = ',smois(18,4,12)
    !print*, 'SMCWTD = ',smcwtdxy(18,12)
    !print*, 'EQWTD = ',eqwtd(18,12)
    !print*, 'WTD = ',wtd(18,12)

    do row=NOAHMP401_struc(n)%row_min, NOAHMP401_struc(n)%row_max
        do col=NOAHMP401_struc(n)%col_min, NOAHMP401_struc(n)%col_max
            t = NOAHMP401_struc(n)%rct_idx(col,row)
            !print*, col
            !print*, row
            !print*, wtd(col,row)
            !print*, t
            ! Added if statement to deal with no-data values. TML
            if(t .ne. LIS_rc%udef) then
                NOAHMP401_struc(n)%noahmp401(t)%wtd       = wtd(col,row) 
                NOAHMP401_struc(n)%noahmp401(t)%zwt       = wtd(col,row)  !!!! zwt should be the same as wtd 
                NOAHMP401_struc(n)%noahmp401(t)%rivercond = rivercond(col,row) 
                NOAHMP401_struc(n)%noahmp401(t)%smc(:)    = smois(col,:,row) ! smois
                NOAHMP401_struc(n)%noahmp401(t)%sh2o(:)   = sh2o(col,:,row)
                NOAHMP401_struc(n)%noahmp401(t)%smoiseq(:)= smoiseq(col,:,row) 
                NOAHMP401_struc(n)%noahmp401(t)%tslb(:)= tslb(col,:,row)
                if(isltyp(col,row) .eq. 14) then
                    NOAHMP401_struc(n)%noahmp401(t)%smcwtd = 1.0
                else
                    NOAHMP401_struc(n)%noahmp401(t)%smcwtd = smcwtdxy(col,row)
                endif
                NOAHMP401_struc(n)%noahmp401(t)%deeprech  = deeprechxy(col,row)
                NOAHMP401_struc(n)%noahmp401(t)%rech      = rechxy(col,row)
                NOAHMP401_struc(n)%noahmp401(t)%qslat     = qslatxy(col,row) 
                NOAHMP401_struc(n)%noahmp401(t)%qrfs      = qrfsxy(col,row) 
                NOAHMP401_struc(n)%noahmp401(t)%qsprings  = qspringsxy(col,row)
            endif 
            NOAHMP401_struc(n)%rivercond(col, row)    = rivercond(col,row) !!! make a copy to the 2D paramter data structure 
            NOAHMP401_struc(n)%riverbed(col, row)     = riverbed(col,row)  !!! make a copy to the 2D paramter data structure 
            NOAHMP401_struc(n)%eqwtd(col, row)        = eqwtd(col,row)     !!! make a copy 
        enddo
    enddo 
    deallocate(isltyp)
    deallocate(ivgtyp)
    deallocate(fdepth)
    deallocate(topo)
    deallocate(area)
    deallocate(rechclim)
    deallocate(rivercond)
    deallocate(wtd)
    deallocate(riverbed)
    deallocate(eqwtd)
    deallocate(pexp)
    deallocate(smcwtdxy)
    deallocate(deeprechxy)
    deallocate(rechxy)
    deallocate(qslatxy)
    deallocate(qrfsxy)
    deallocate(qspringsxy)
    deallocate(smois)
    deallocate(sh2o)
    deallocate(smoiseq)
    deallocate(tslb)
end subroutine mmf_start
