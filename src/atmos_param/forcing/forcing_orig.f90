module forcing_mod

!-----------------------------------------------------------------------
! Authors: Jeremiah Sjoberg
!          Luke Davis (lukelbd@gmail.com)

! This forcing module optionally allows for Polvani-Kushner setup instead
! of the conventional Held-Suarez default, and adds sponge layer to upper stratosphere.
!
! Also adds new 'pkmod' forcing mode, which implements the following changes:
! 1) Arbitrary "tropopause temperature" (i.e. temperature minimum around
!    15km). Accomplished by extending/contracting the region of constant
!    temperature, and contracting/extending the region of +2.8K/km lapse
!    rate. Note can only contract/extend so far before eliminating these
!    piecewise 'regions', so make sure not to offset temperature by too much.
! 2) Tropopause-conforming height coordinates, or arbitrary switch to the
!    US standard atmosphere/polar vortex stratosphere. Accomplished by
!    offsetting height 'z' in the US standard temperature formula by the
!    difference between tropopause height X and reference height Y.
! 3) Polar vortex adaptation.  Add another segment, extending for 12km in
!    the vertical, with negative lapse rate (probably 2K/km is reasonable),
!    then have 4km region of constant temperature and then proceed with
!    the ordinary standard atmosphere specification.
!-----------------------------------------------------------------------

use     constants_mod, only: seconds_per_day, kappa, cp_air, grav, pi
use   mpp_domains_mod, only: mpp_get_global_domain
use    transforms_mod, only: grid_domain
use           fms_mod, only: error_mesg, FATAL, file_exist, field_size, &
                             read_data, open_namelist_file, check_nml_error, &
                             mpp_pe, mpp_root_pe, close_file, &
                             write_version_number, stdlog, &
                             uppercase

use  time_manager_mod, only: time_type, get_time

use  diag_manager_mod, only: register_diag_field, send_data

use  field_manager_mod, only: MODEL_ATMOS, parse
use tracer_manager_mod, only: query_method, get_number_tracers

implicit none
private

!-----------------------------------------------------------------------
!     ----- Interfaces -----

public :: forcing, forcing_init

!-----------------------------------------------------------------------
!     ----- Namelist -----
real, parameter :: H = 7.0 ! for calculation below
character(len=64) :: teq_mode = 'hs', damp_mode = 'hs', strat_damp = 'constant'
logical :: no_forcing = .false.
logical :: conserve_energy = .true.
logical :: strat_sponge = .true., strat_vtx = .true.
logical :: surf_schneider = .false. ! HS94 or Schneider surface temp?
logical :: locked_heating = .false., ndamp_decomp = .false., rdamp_decomp = .false., sponge_decomp = .false.
integer :: exp_b = 4 ! exponent in cosine for boundary layer; Schneider uses 8, Held-Suarez 4
integer :: exp_h = 0 ! exponent describing deflection from default meridional gradient sin^2(lat) = 1 - cos^(lat)
real :: sigma_b  = 0.7
real :: t_zero = 315.0, t_mean = 300.0, t_strat = 200.0, t_min = 100.0
real :: delh = 60.0, delv = 10.0, eps = 0.0
real :: vtx_edge = 50.0, vtx_width = 10.0, vtx_gamma = 2.0
real :: p_ref = 1000.0, p_sponge = 0.5, p_logeval = 200.0
real :: z_ozone = 20.0, z_kdepth = 50.0 ! pkswitch at 16km evaluates to roughly 100mb with scale height 7km
real :: lat_ref = 0
real :: q0_tropical = 0.0, x0_tropical = 0.0, y0_tropical = 0.3, sx_tropical = 0.4, sy_tropical = 0.11 ! upper troposphere is default
real :: q0_vortex = 0.0, x0_vortex = 1.57, y0_vortex = 0.1, sx_vortex = 0.4, sy_vortex = 0.14          ! lower stratosphere is default
real :: q0_arctic = 0.0, x0_arctic = 1.57, y0_arctic = 1.0
real :: q0_global = 0.0, q0_surface = 0.0, q0_realistic = 0.0
real :: q0_lsp = 0.0, m_lsp = 1.0, p0_lsp = 800.0, pt_lsp = 200.0, lat0_lsp = 45.0, slat_lsp = 10.0
real :: trflux = 1.e-5 ! surface flux for optional tracer
real :: trsink = -4.   ! damping time for tracer
!     ----- Allocatable params -----
! This mimics how surf_geopot in spectral_dynamics is allocated during
! initialization.
real, allocatable, dimension(:,:,:) :: tdt_locked
!     ----- Account for optional mean, anomaly component damping -----
! Namelist can specify just *one* value, and it will be
! inserted in first position, default Held-Suarez used for second position
real, parameter :: fill_value = 999. ! use this as dummy
real, dimension(2) :: kbl     = fill_value
real, dimension(2) :: ktrop   = fill_value
real, dimension(2) :: kstrat  = fill_value
real, dimension(2) :: kmeso   = fill_value
real, dimension(2) :: kfric   = fill_value
real, dimension(2) :: ksponge = fill_value

!-----------------------------------------------------------------------
namelist /forcing_nml/  no_forcing, locked_heating, conserve_energy, &
  teq_mode, damp_mode, strat_sponge, strat_vtx, strat_damp, surf_schneider, &
  t_zero, t_mean, t_strat, exp_h, delh, delv, eps, &
  lat_ref, vtx_edge, vtx_width, vtx_gamma, &
  exp_b, sigma_b, p_sponge, p_logeval, z_ozone, z_kdepth, &
  ktrop, kbl, kstrat, kmeso, kfric, ksponge, &
  trflux, trsink, &
  q0_lsp, m_lsp, p0_lsp, pt_lsp, lat0_lsp, slat_lsp, &
  q0_tropical, q0_vortex, q0_arctic, q0_global, q0_surface, q0_realistic, &
  x0_tropical, y0_tropical, sx_tropical, sy_tropical, &
  x0_vortex, y0_vortex, sx_vortex, sy_vortex, &
  x0_arctic, y0_arctic

!-----------------------------------------------------------------------
character(len=14) :: mod_name = 'forcing'
character(len=128) :: version='$Id: new_forcing.f90 2019/04/01 00:00:00 fms Exp $'
character(len=128) :: tagname='$Name: latest $'
real :: trdamp
real, dimension(2) :: tktrop, tkbl, tkstrat, tkmeso, vkfric, vksponge
integer :: id_forcing, id_teq, &
  id_ndamp, id_rdamp, id_sdamp, id_tdt, id_heating, id_udt, id_vdt, &
  id_ndamp_mean, id_ndamp_anom, id_rdamp_mean, id_rdamp_anom, &
  id_tdt_mean, id_tdt_anom, id_udt_mean, id_udt_anom, id_vdt_mean, id_vdt_anom, &
  id_tdt_diss, id_heat_diss, &
  id_udt_sponge, id_vdt_sponge
real    :: missing_value = -1.e10
logical :: module_is_initialized = .false.

contains

subroutine forcing ( is, ie, js, je, dt, Time, lon, lat, p_half, p_full, z_full, &
                         u, v, t, r, um, vm, tm, rm, udt, vdt, tdt, &
                         rdt, mask, kbot )

  !-----------------------------------------------------------------------
  type(time_type), intent(in)                      :: Time
  integer, intent(in)                              :: is, ie, js, je
  real, intent(in)                                 :: dt
  real, intent(in),    dimension(:,:)              :: lon, lat
  real, intent(in),    dimension(:,:,:)            :: p_half, p_full, z_full
  real, intent(in),    dimension(:,:,:)            :: u, v, t, um, vm, tm
  real, intent(in),    dimension(:,:,:,:)          :: r, rm
  real, intent(inout), dimension(:,:,:)            :: udt, vdt, tdt
  real, intent(inout), dimension(:,:,:,:)          :: rdt
  real, intent(in),    dimension(:,:,:), optional  :: mask
  integer, intent(in),    dimension(:,:), optional :: kbot
  !-----------------------------------------------------------------------
  real, dimension(size(t,1),size(t,2))             :: ps, diss_heat
  real, dimension(size(t,1),size(t,2),size(t,3))   :: sigma, pmass ! sigma and mass
  real, dimension(size(t,1),size(t,2),size(t,3))   :: udt_sponge, vdt_sponge, sdamp
  real, dimension(size(t,1),size(t,2),size(t,3))   :: teq, tdt_force, tdt_diss, udt_diss, vdt_diss
  real, dimension(size(r,1),size(r,2),size(r,3))   :: rst, rtnd
  real, dimension(size(t,1),size(t,2),size(t,3),2) :: tdt_damp, tdamp, udt_damp, vdt_damp, rdamp
  integer :: i, j, k, kb, c, n, num_tracers
  logical :: used
  real    :: flux, sink, value
  character(len=128) :: scheme, params

  if (no_forcing) return
  if (.not.module_is_initialized) call error_mesg ('forcing','forcing_init has not been called', FATAL)

  !     Surface pressure and sigma coordinates
  if (present(kbot)) then
    do j=1,size(p_half,2)
      do i=1,size(p_half,1)
        kb = kbot(i,j)
        ps(i,j) = p_half(i,j,kb+1)
      enddo
    enddo
  else
    ps = p_half(:,:,size(p_half,3))
  endif
  do k=1,size(p_full,3)
    sigma(:,:,k) = p_full(:,:,k) / ps
  enddo

  !     Rayleigh damping of wind components near the surface
  udt_damp = 0.0
  vdt_damp = 0.0
  call friction_damping ( sigma, u, v, udt_damp, vdt_damp, rdamp, mask=mask )
  if (id_udt   > 0) used = send_data(id_udt,   sum(udt_damp,4),   Time, is, js)
  if (id_vdt   > 0) used = send_data(id_vdt,   sum(vdt_damp,4),   Time, is, js)
  if (id_rdamp > 0) used = send_data(id_rdamp, rdamp(:,:,:,1), Time, is, js)
  if (rdamp_decomp) then
    if (id_udt_mean   > 0) used = send_data(id_udt_mean,   udt_damp(:,:,:,1), Time, is, js)
    if (id_udt_anom   > 0) used = send_data(id_udt_anom,   udt_damp(:,:,:,2), Time, is, js)
    if (id_vdt_mean   > 0) used = send_data(id_vdt_mean,   vdt_damp(:,:,:,1), Time, is, js)
    if (id_vdt_anom   > 0) used = send_data(id_vdt_anom,   vdt_damp(:,:,:,2), Time, is, js)
    if (id_rdamp_mean > 0) used = send_data(id_rdamp_mean, rdamp(:,:,:,1), Time, is, js)
    if (id_rdamp_anom > 0) used = send_data(id_rdamp_anom, rdamp(:,:,:,2), Time, is, js)
  endif


  !     Sponge layer damping of wind components at the top
  udt_sponge = 0.0
  vdt_sponge = 0.0
  if (strat_sponge) then
    call sponge_damping ( p_full, u, v, udt_sponge, vdt_sponge, sdamp, mask=mask )
    if (id_udt_sponge > 0) used = send_data(id_udt_sponge, udt_sponge, Time, is, js)
    if (id_vdt_sponge > 0) used = send_data(id_vdt_sponge, vdt_sponge, Time, is, js)
    if (id_sdamp > 0) used = send_data(id_sdamp, sdamp, Time, is, js)
  endif

  !     Thermal damping for held & suarez (1994) benchmark calculation
  !     Alternatively load heating from file
  if (locked_heating) then
    tdt_damp(:,:,:,1) = tdt_locked
    tdt_damp(:,:,:,2) = 0.0
    if (id_heating > 0) used = send_data(id_heating, tdt_locked, Time, is, js)
  else
    call thermal_damping ( lat, sigma, p_full, z_full, t, tdt_damp, tdamp, teq, mask=mask )
  endif
  if (id_teq   > 0) used = send_data(id_teq,   teq,             Time, is, js)
  if (id_tdt   > 0) used = send_data(id_tdt,   sum(tdt_damp,4), Time, is, js)
  if (id_ndamp > 0) used = send_data(id_ndamp, tdamp(:,:,:,1),  Time, is, js)
  if (ndamp_decomp) then
    if (id_tdt_mean   > 0) used = send_data(id_tdt_mean,   tdt_damp(:,:,:,1), Time, is, js)
    if (id_tdt_anom   > 0) used = send_data(id_tdt_anom,   tdt_damp(:,:,:,2), Time, is, js)
    if (id_ndamp_mean > 0) used = send_data(id_ndamp_mean, tdamp(:,:,:,1),    Time, is, js)
    if (id_ndamp_anom > 0) used = send_data(id_ndamp_anom, tdamp(:,:,:,2),    Time, is, js)
  endif

  ! Energy conservation
  tdt_diss = 0.0
  if (conserve_energy) then
    ! Apply tempreature tendency due to kinetic energy dissipation
    udt_diss = sum(udt_damp,4)
    vdt_diss = sum(vdt_damp,4)
    tdt_diss = -((um + 0.5 * udt_diss * dt) * udt_diss + (vm + 0.5 * vdt_diss * dt) * vdt_diss) / cp_air
    ! Diagnostic output
    if (id_tdt_diss  > 0) used = send_data(id_tdt_diss, tdt_diss, Time, is, js)
    if (id_heat_diss > 0) then
      pmass = p_half(:,:,2:) - p_half(:,:,:size(p_half,3)-1)
      diss_heat = cp_air / grav * sum(tdt_diss * pmass, 3)
      used = send_data(id_heat_diss, diss_heat, Time, is, js)
    endif
  endif

  !     Butler et al. (2010) radiative forcing and other forcing
  !     Lindgren et al. (2018) wave heating
  !     Dissipative heating to conserve energy
  ! Heating terms
  tdt_force = 0.0
  if (q0_global .ne. 0) then ! global constant heating
    tdt_force = tdt_force + q0_global / seconds_per_day
  endif
  if (q0_realistic .ne. 0) then ! heating that mimics equilibrium surface temp change ! where (teq > t_strat)
    tdt_force = tdt_force + q0_realistic * ((p_full / 100000.0) ** kappa) / seconds_per_day
  endif
  if (q0_surface .ne. 0) then ! maximum strength at surface then decays to zero at boundary layer
    where (sigma > sigma_b)
      tdt_force = tdt_force + q0_surface * ((sigma - sigma_b) / (1 - sigma_b)) / seconds_per_day
    endwhere
  endif
  if ((q0_arctic .ne. 0) .or. (q0_tropical .ne. 0) .or. (q0_vortex .ne. 0)) then
    call butler_heating ( lat, sigma, tdt_force )
  endif
  if (q0_lsp .ne. 0) then
    call lsp_heating( lon, lat, sigma, tdt_force )
  endif
  if (present(mask)) then
    tdt_force = mask * tdt_force
  endif
  if (id_forcing > 0) used = send_data(id_forcing, tdt_force, Time, is, js)

  !     Apply forced wind and temperature tendencies
  udt = udt + udt_damp(:,:,:,1) + udt_damp(:,:,:,2) + udt_sponge
  vdt = vdt + vdt_damp(:,:,:,1) + vdt_damp(:,:,:,2) + vdt_sponge
  tdt = tdt + tdt_damp(:,:,:,1) + tdt_damp(:,:,:,2) + tdt_force + tdt_diss ! mean and anomaly components

  !     Tracers
  call get_number_tracers(MODEL_ATMOS, num_tracers=num_tracers)
  if (num_tracers == size(rdt,4)) then
    do n = 1, size(rdt,4)
      flux = trflux
      sink = trsink
      if (query_method('tracer_sms', MODEL_ATMOS, n, scheme, params)) then
        if (uppercase(trim(scheme)) == 'NONE') cycle
        if (uppercase(trim(scheme)) == 'OFF') then
          flux = 0.; sink = 0.
        else
          if (parse(params, 'flux', value) == 1) flux = value
          if (parse(params, 'sink', value) == 1) sink = value
        endif
      endif
      rst = rm(:,:,:,n) + dt * rdt(:,:,:,n)
      call tracer_source_sink ( flux, sink, p_half, rst, rtnd, kbot )
      rdt(:,:,:,n) = rdt(:,:,:,n) + rtnd
    enddo
  else if (num_tracers == 0 .and. size(rdt,4) == 1) then ! allow this as a getaround for a problem with the solo fv model
    flux = trflux
    sink = trsink
    rst = rm(:,:,:,1) + dt * rdt(:,:,:,1)
    call tracer_source_sink ( flux, sink, p_half, rst, rtnd, kbot )
    rdt(:,:,:,1) = rdt(:,:,:,1) + rtnd
  else
    call error_mesg('forcing','size(rdt,4) not equal to num_tracers', FATAL)
  endif

end subroutine forcing

subroutine forcing_init ( axes, is, ie, js, je, num_levels, Time )

  !-----------------------------------------------------------------------
  !
  !           Routine for initializing the model with an
  !           initial condition at rest (u & v = 0)
  !           Also loads heating data for locked heating runs
  !
  !-----------------------------------------------------------------------
  ! Notes
  ! * Starting point was topography.data.nc, which is loaded globally
  !   during call to atmosphere_init in atmos_model.f90, which calls
  !   spectral_dynamics_init, which calls read_restart_or_do_coldstart, which
  !   loads surface geopotential only if this is a cold start (because it is
  !   saved in restart files), which calls spectral_init_cond, which calls
  !   get_topography, which loads data on global grid.
  ! * Also input geopotential height instead of using scale height to
  !   approximate it. Easy in spectral, but hard for bgrid, because
  !   geopoential height is not normally passed throuth atmosphere.f90. Need
  !   to compute from compute_height.f90 in bgrid_vert.f90, by providing it
  !   with 'fisl' eta=1 geopotential from the Dynam structure and with
  !   temperature and stuff from the Var structure. Similar routine compute_geop_height
  !   is called in bgrid_horiz_adjust.f90 by press_grad, called in bgrid_driver.f90.
  ! * Seems that all init methods work with *global* grid, then subsequent
  !   calls to amosphere subroutine are given partial data.

  integer, intent(in) :: axes(4)
  integer, intent(in) :: is, ie, js, je, num_levels
  type(time_type), intent(in) :: Time
  integer :: global_num_lon, global_num_lat, siz(4)
  character(len=20) :: ctmp1='     by      by     ', ctmp2='     by      by     '

  !-----------------------------------------------------------------------
  integer unit, io, ierr, c

  !     ----- Read namelist -----

  if (file_exist('input.nml')) then
    unit = open_namelist_file ( )
    ierr=1; do while (ierr /= 0)
    read  (unit, nml=forcing_nml, iostat=io, end=10)
    ierr = check_nml_error (io, 'forcing_nml')
    enddo
  10     call close_file (unit)
  endif

  !     ----- Write version info and namelist to log file -----

  call write_version_number (version,tagname)
  if (mpp_pe() == mpp_root_pe()) write (stdlog(),nml=forcing_nml)

  if (no_forcing) return

  !     ----- Error checks -----

  if (vtx_edge <  0) then
    call error_mesg ('forcing_init','polar vortex edge must be in northern hemisphere', FATAL)
  end if
  if ((trim(damp_mode) .ne. 'hs') .and. (trim(damp_mode) .ne. 'pk') .and. (trim(damp_mode) .ne. 'pkmod')) then
    call error_mesg ('forcing_init', 'damp_mode="' // damp_mode // '" not recognized', FATAL)
  end if
  if ((trim(teq_mode) .ne. 'hs') .and. (trim(teq_mode) .ne. 'pk') .and. (trim(teq_mode) .ne. 'pkmod')) then
    call error_mesg ('forcing_init', 'teq_mode="' // teq_mode // '" not recognized', FATAL)
  end if
  if ((trim(strat_damp) .ne. 'constant') .and. (trim(strat_damp) .ne. 'linear')) then
    call error_mesg ('forcing_init', 'strat_damp="' // strat_damp // '" not recognized', FATAL)
  endif

  !     ----- Allocate and load damping -----
  if (locked_heating) then
    allocate(tdt_locked(is:ie, js:je, num_levels))
    if (file_exist('INPUT/heating.data.nc')) then
      call mpp_get_global_domain(grid_domain, xsize=global_num_lon, ysize=global_num_lat)
      call field_size('INPUT/heating.data.nc', 'tdt', siz)
      if (( siz(1) == global_num_lon .or. siz(2) == global_num_lat ) .and. siz(3) == num_levels) then
        call read_data('INPUT/heating.data.nc', 'tdt', tdt_locked, grid_domain)
      else
        write(ctmp1(1: 4),'(i4)') siz(1)
        write(ctmp1(9:12),'(i4)') siz(2)
        write(ctmp1(17:20),'(i4)') siz(3)
        write(ctmp2(1: 4),'(i4)') global_num_lon
        write(ctmp2(9:12),'(i4)') global_num_lat
        write(ctmp2(17:20),'(i4)') num_levels
        call error_mesg ('forcing_init','Topography file contains data on a '// &
                ctmp1//' grid, but atmos model grid is '//ctmp2, FATAL)
      endif
    else
      call error_mesg('forcing_init', 'locked_heating=.true. but INPUT/heating.data.nc does not exist', FATAL)
    endif
  endif

  !     ----- For tracers -----

  trdamp = 0.0
  if (trsink < 0) trsink = -seconds_per_day * trsink
  if (trsink > 0) trdamp = 1.0 / trsink

  !     ----- Apply defaults -----

  if (kbl(1) == fill_value)     kbl(1)     = -4.
  if (ktrop(1) == fill_value)   ktrop(1)   = -40.
  if (kstrat(1) == fill_value)  kstrat(1)  = -20.
  if (kmeso(1) == fill_value)   kmeso(1)   = -1.
  if (kfric(1) == fill_value)   kfric(1)   = -1.
  if (ksponge(1) == fill_value) ksponge(1) = -2.

  !     ----- If one number specified in namelist, use it for *both* components -----

  if (kbl(2) == fill_value) then
    kbl(2) = kbl(1)
  else
    ndamp_decomp = .true.
  end if
  if (ktrop(2) == fill_value) then
    ktrop(2) = ktrop(1)
  else
    ndamp_decomp = .true.
  end if
  if (kstrat(2) == fill_value) then
    kstrat(2) = kstrat(1)
  else
    ndamp_decomp = .true.
  end if
  if (kmeso(2) == fill_value) then
    kmeso(2) = kmeso(1)
  else
    ndamp_decomp = .true.
  end if
  if (kfric(2) == fill_value) then
    kfric(2) = kfric(1)
  else
    rdamp_decomp = .true.
  end if
  if (ksponge(2) == fill_value) then
    ksponge(2) = ksponge(1)
  else
    sponge_decomp = .true.
  end if

  !     ----- Compute coefficients -----

  tkbl     = 0.0
  tktrop   = 0.0
  vkfric   = 0.0
  vksponge = 0.0
  tkstrat  = 0.0
  tkmeso   = 0.0
  do c=1,2
    ! Convert timescales from days to seconds
    if (kbl(c)     < 0) kbl(c)     = -seconds_per_day * kbl(c)
    if (ktrop(c)   < 0) ktrop(c)   = -seconds_per_day * ktrop(c)
    if (kfric(c)   < 0) kfric(c)   = -seconds_per_day * kfric(c)
    if (ksponge(c) < 0) ksponge(c) = -seconds_per_day * ksponge(c)
    if (kstrat(c)  < 0) kstrat(c)  = -seconds_per_day * kstrat(c)
    if (kmeso(c)   < 0) kmeso(c)   = -seconds_per_day * kmeso(c)
    ! Get coefficients from timescales
    if (kbl(c)     > 0) tkbl(c)     = 1.0 / kbl(c)
    if (ktrop(c)   > 0) tktrop(c)   = 1.0 / ktrop(c)
    if (kfric(c)   > 0) vkfric(c)   = 1.0 / kfric(c)
    if (ksponge(c) > 0) vksponge(c) = 1.0 / ksponge(c)
    if (kstrat(c)  > 0) tkstrat(c)  = 1.0 / kstrat(c)
    if (kmeso(c)   > 0) tkmeso(c)   = 1.0 / kmeso(c)
  enddo

  !     ----- Register diagnostic fields -----

  ! Radiative forcing
  id_forcing = register_diag_field ( mod_name, 'forcing', axes(1:3), Time, &
    'radiative forcing', 'deg_K/sec', &
    missing_value=missing_value )

  ! Equilibrium temp
  id_teq = register_diag_field ( mod_name, 'teq', axes(1:3), Time, &
    'equilibrium temperature', 'deg_K', &
    missing_value=missing_value, range=(/0.,500./) )

  ! Rate of temp change due to damping and damping rate
  id_tdt = register_diag_field ( mod_name, 'tdt', axes(1:3), Time, &
    'heating', 'deg_K/sec', &
    missing_value=missing_value )

  if (locked_heating) then
    id_heating = register_diag_field (mod_name, 'heating', axes(1:3), Time, &
      'locked heating', 'deg_K/sec', &
      missing_value=missing_value )
  endif

  id_ndamp = register_diag_field ( mod_name, 'ndamp', axes(1:3), Time, &
    'thermal damping coefficient', 'sec-1', &
    missing_value=missing_value )

  if (ndamp_decomp) then
    id_tdt_mean = register_diag_field ( mod_name, 'tdt_mean', axes(1:3), Time, &
      'heating, mean component', 'deg_K/sec', &
      missing_value=missing_value )

    id_tdt_anom = register_diag_field ( mod_name, 'tdt_anom', axes(1:3), Time, &
      'heating, anomaly component', 'deg_K/sec', &
      missing_value=missing_value )

    id_ndamp_mean = register_diag_field ( mod_name, 'ndamp_mean', axes(1:3), Time, &
      'thermal damping coefficient, mean component', 'sec-1', &
      missing_value=missing_value )

    id_ndamp_anom = register_diag_field ( mod_name, 'ndamp_anom', axes(1:3), Time, &
      'thermal damping coefficient, anomaly component', 'sec-1', &
      missing_value=missing_value )
  endif

  ! Heating by dissipation
  if (conserve_energy) then
    id_tdt_diss = register_diag_field ( mod_name, 'tdt_diss_rdamp', axes(1:3), &
      Time, 'dissipative heating from frictional damping', 'deg_K/sec',&
      missing_value=missing_value )

    id_heat_diss = register_diag_field ( mod_name, 'diss_rdamp', axes(1:2), &
      Time, 'integrated dissipative heating for frictional damping', 'W/m2')
  endif

  ! Rate of wind change due to friction
  id_udt = register_diag_field ( mod_name, 'udt', axes(1:3), Time, &
    'frictional zonal wind damping', 'm/s2',       &
    missing_value=missing_value     )

  id_vdt = register_diag_field ( mod_name, 'vdt', axes(1:3), Time, &
    'frictional meridional wind damping', 'm/s2',  &
    missing_value=missing_value )

  id_rdamp = register_diag_field ( mod_name, 'rdamp', axes(1:3), Time, &
    'frictional damping coefficient', 'sec-1', &
    missing_value=missing_value )

  if (rdamp_decomp) then
    id_udt_mean = register_diag_field ( mod_name, 'udt_mean', axes(1:3), Time, &
      'frictional zonal wind damping, mean component', 'm/s2',       &
      missing_value=missing_value )

    id_udt_anom = register_diag_field ( mod_name, 'udt_anom', axes(1:3), Time, &
      'frictional zonal wind damping, anomaly component', 'm/s2', &
      missing_value=missing_value )

    id_vdt_mean = register_diag_field ( mod_name, 'vdt_mean', axes(1:3), Time, &
      'frictional meridional wind damping, mean component', 'm/s2', &
      missing_value=missing_value )

    id_vdt_anom = register_diag_field ( mod_name, 'vdt_anom', axes(1:3), Time, &
      'frictional meridional wind damping, anomaly component', 'm/s2', &
      missing_value=missing_value )

    id_rdamp_mean = register_diag_field ( mod_name, 'rdamp_mean', axes(1:3), Time, &
      'frictional damping coefficient, mean component', 'sec-1', &
      missing_value=missing_value )

    id_rdamp_anom = register_diag_field ( mod_name, 'rdamp_anom', axes(1:3), Time, &
      'frictional damping coefficient, anomaly component', 'sec-1', &
      missing_value=missing_value )
  endif

  ! Rate of wind change by sponge
  if (strat_sponge) then
    id_udt_sponge = register_diag_field ( mod_name, 'udt_sponge', axes(1:3), Time, &
      'sponge zonal wind damping', 'm/s2',       &
      missing_value=missing_value     )

    id_vdt_sponge = register_diag_field ( mod_name, 'vdt_sponge', axes(1:3), Time, &
      'sponge meridional wind damping', 'm/s2',       &
      missing_value=missing_value     )

    id_sdamp = register_diag_field ( mod_name, 'sdamp', axes(1:3), Time, &
      'sponge damping coefficient', 'sec-1', &
      missing_value=missing_value )
  endif

  module_is_initialized  = .true.

end subroutine forcing_init

subroutine forcing_end 

  module_is_initialized = .false.

end subroutine forcing_end

subroutine thermal_damping ( lat, sigma, p_full, z_full, t, tdt, tdamp, teq, mask )

!-----------------------------------------------------------------------
!
!   Routine to compute thermal forcing for several benchmarks
!
!-----------------------------------------------------------------------
! Forcing parameters teq and tdamp have longitude dimension, but why?
! Easier to multiply with temperature, and allows us to impose longitudinally
! varying forcing profile with ease.
real, intent(in),  dimension(:,:)   :: lat
real, intent(in),  dimension(:,:,:) :: sigma, p_full, z_full, t
real, intent(in),  dimension(:,:,:), optional :: mask
real, intent(out), dimension(:,:,:)   :: teq
real, intent(out), dimension(:,:,:,:) :: tdt, tdamp
!-----------------------------------------------------------------------
real, dimension(size(t,1),size(t,2)) :: &
  t_surf, t_hs, t_pk, t_summer, t_winter, &
  lat_hfact, lat_vfact, w_vtx, &
  z_vortex, z_offset, z_norm, p_trop, p_norm ! intermediate containers
real, dimension(size(t,1),size(t,2),2) :: &
  t_decomp, tkhi_decomp, tksurf_decomp ! mean-anomaly decomp containers
real, dimension(size(t,2)) :: &
  t_bar ! just holds the mean
real :: p0, pb, pfact
real :: z_vortex_ref, p_trop_ref, t_surf_ref
real :: lat_ref_r, lat_vfact_ref, lat_hfact_ref
real :: vtx_edge_r, vtx_width_r
integer :: i, j, k, c, seconds, days

!-----------------------------------------------------------------------
!     Constants
! Note that Lindgren et al 2018 (lsp) wrote sigma_phi = 0.175 * 360/2pi but that is only if latitudes are in degrees, not radians; here they are in radians!
p0            = p_ref * 100.0
pfact         = p_logeval * 100.0 / p0
lat_ref_r     = lat_ref * pi / 180.0 ! latitude for retrieving reference tropopause height
vtx_width_r   = vtx_width * pi / 180.0
vtx_edge_r    = vtx_edge * pi / 180.0

!-----------------------------------------------------------------------
!     Surface temp and damping
lat_vfact = cos(lat) ** 2 ! for use later on
lat_vfact_ref = cos(lat_ref_r) ** 2
if (exp_h .ge. 0) then
  ! Poleward shift or default configuration
  lat_hfact = sin(lat) ** (2 + exp_h)
  lat_hfact_ref = sin(lat_ref_r) ** (2 + exp_h)
else
  ! Equatorward shift
  lat_hfact = 1 - cos(lat) ** (2 - exp_h)
  lat_hfact_ref = 1 - cos(lat_ref_r) ** (2 - exp_h)
endif
! Surface damping
do c=1,2
  tksurf_decomp(:,:,c) = (tkbl(c) - tktrop(c)) * cos(lat) ** exp_b
enddo

!-----------------------------------------------------------------------
!     Use HS94 surface, or special Schneider version?
! NOTE: Schneider just uses int(cosx * (t0 + delh * (1/3 - sin(lat)^2))), ensures
! average due to the sine term cancels with the offset term. Easy to do for
! exp_h > 1, hard for cosine version. We don't try to make a general form.
if (surf_schneider) then
  t_surf = t_mean + delh * (1.0 / 3.0 - lat_hfact) - eps * sin(lat)
  t_surf_ref = t_mean + delh * (1.0 / 3.0 - lat_hfact_ref) - eps * sin(lat_ref_r)
else
  t_surf = t_zero - delh * lat_hfact - eps * sin(lat)
  t_surf_ref = t_zero - delh * lat_hfact_ref - eps * sin(lat_ref_r)
endif

!-----------------------------------------------------------------------
!     Vortex weighting
w_vtx = 0.0 ! standard atmosphere everywhere
if (strat_vtx) then
  w_vtx = 0.5 * (1.0 + tanh((lat - abs(vtx_edge_r)) / vtx_width_r)) ! vortex in northern hemisphere
endif

!-----------------------------------------------------------------------
!     Tropopause height estimate
! Dimensions are longitude by <processor latitudes> by height
! Run model with one core for print statements

if (trim(teq_mode) .eq. 'pk') then
  ! In PK atmosphere, polar vortex and stratospheric warming both begin at
  ! an arbitrary input height. NOTE: This is slightly different! Stratospheric
  ! warming starts at same level as polar vortex cooling!
  z_vortex = z_ozone
  z_offset = (z_ozone - 20.0) ! z_ozone is height at which we want warming begins, 20.0 is default case

else if (trim(teq_mode) .eq. 'pkmod') then
  ! In PKMOD atmosphere, polar vortex begins right at the tropopause, and
  ! stratospheric warming begins at 20km at reference latitude
  p_trop = p0 * ( t_strat / (t_surf - delv * log(pfact) * lat_vfact) ) ** (1.0 / kappa)
  p_trop_ref = p0 * ( t_strat / (t_surf_ref - delv * log(pfact) * lat_vfact_ref) ) ** (1.0 / kappa)
  z_vortex = -H * log(p_trop / p0)
  z_offset = (z_vortex + H * log(p_trop_ref / p0)) + (z_ozone - 20.0) ! includes offset to accomodate configurable standard atmosphere height
  ! print *, p_trop(0,size(t,2) / 8), p_trop(0,size(t,2) / 4), p_trop(0,size(t,2) / 2)
  ! print *, z_vortex(0,size(t,2) / 8), z_vortex(0,size(t,2) / 4), z_vortex(0,size(t,2) / 2), z_vortex_ref
  ! print *, z_offset(0,size(t,2) / 8), z_offset(0,size(t,2) / 4), z_offset(0,size(t,2) / 2)
endif

! Loop through pressure
do k=1,size(t,3)
  !-----------------------------------------------------------------------
  !
  !     Equilibrium temp calculations
  !
  !-----------------------------------------------------------------------
  ! WARNING: In assignment e.g. a = max(a, 0), end up with float NaNs
  ! when fortran tried to assign the scalar. Some weird fortran thing.
  ! We need two variables for this kind of statement!
  z_norm = z_full(:,:,k) / 1000. ! from m to km
  p_norm = p_full(:,:,k) / p0
  t_hs   = (t_surf - delv * log(p_norm) * lat_vfact) * (p_norm) ** kappa
  if (trim(teq_mode) == 'hs') then
    teq(:,:,k) = max( t_hs, t_strat )
  else
    call tstd_summer( t_strat, z_offset, z_norm, t_summer )
    call tstd_winter( t_strat, vtx_gamma, z_vortex, z_norm, t_winter )
    t_pk = (1.0 - w_vtx) * t_summer + w_vtx * t_winter
    where (z_norm >= z_vortex)
      teq(:,:,k) = max( t_pk, t_min )
    elsewhere
      teq(:,:,k) = max( t_hs, t_strat )
    endwhere
  endif

  do c=1,2
    !-----------------------------------------------------------------------
    !
    !     Damping rate calculations
    !     Allows different damping rates applied to mean and anomaly components
    !
    !-----------------------------------------------------------------------
    ! Alternatively use smooth version below: see Holton-Mass model, 1976 (Journal of Atmospheric Sciences)
    ! tkhi = tkstrat + (1.0 + tanh((z_norm-35.0) / 7.0)) * (tkmeso-tkstrat) / 2.0
    if (trim(damp_mode) == 'hs') then
      where (sigma(:,:,k) > sigma_b)
        tdamp(:,:,k,c) = tktrop(c) + tksurf_decomp(:,:,c) * (sigma(:,:,k) - sigma_b) / (1.0 - sigma_b)
      elsewhere
        tdamp(:,:,k,c) = tktrop(c)
      endwhere
    else
      if (trim(strat_damp) == 'constant') then ! *constant* above tropopause region
        tkhi_decomp(:,:,c) = tkstrat(c)
      else
        tkhi_decomp(:,:,c) = min(tkmeso(c), &
          tkstrat(c) + (tkmeso(c) - tkstrat(c)) * (z_norm - z_vortex - z_kdepth) / (50 - z_vortex - z_kdepth))
      endif
      where (sigma(:,:,k) > sigma_b)
        tdamp(:,:,k,c) = tktrop(c) + tksurf_decomp(:,:,c) * (sigma(:,:,k) - sigma_b) / (1.0 - sigma_b)
      elsewhere (z_norm < z_vortex) ! below tropopause, not some fixed height
        tdamp(:,:,k,c) = tktrop(c)
      elsewhere (z_norm >= z_vortex .and. z_norm < z_vortex + z_kdepth)
        tdamp(:,:,k,c) = tktrop(c) - (tktrop(c) - tkstrat(c)) * (z_norm - z_vortex) / z_kdepth
      elsewhere
        tdamp(:,:,k,c) = tkhi_decomp(:,:,c)
      endwhere

    endif

    !-----------------------------------------------------------------------
    !
    !     Increment tendency
    !
    !-----------------------------------------------------------------------
    if (.not. ndamp_decomp) then
      ! Damp the *full* temperature and exit
      tdt(:,:,k,1) = -tdamp(:,:,k,1) * (t(:,:,k) - teq(:,:,k))
      tdt(:,:,k,2) = 0.0
      tdamp(:,:,k,2) = tdamp(:,:,k,1)
      exit ! i.e. break from top-level do loop
    else if (c==1) then
      ! Damp the mean component
      t_bar = sum(t(:,:,k), 1) / size(t, 1) ! mean
      do i=1,size(t,1)
        t_decomp(i,:,1) = t_bar
      enddo
      tdt(:,:,k,1) = -tdamp(:,:,k,1) * (t_decomp(:,:,1) - teq(:,:,k))
    else if (c==2) then
      ! Damp the anomalous component
      t_decomp(:,:,2) = t(:,:,k) - t_decomp(:,:,1) ! anomaly relative to mean
      tdt(:,:,k,2) = -tdamp(:,:,k,2) * t_decomp(:,:,2)
    endif
  enddo
enddo

if (present(mask)) then
  teq = teq * mask
  do c=1,2
    tdt(:,:,:,c) = tdt(:,:,:,c) * mask
  enddo
endif

end subroutine thermal_damping

subroutine tstd_summer ( t_tp, z_off, z_km, teq )

  !----------------------------------------------------------------------------!
  ! US standard atmosphere directly from lapse rate formula
  ! Add minor change, with 2.8K/km warming
  !----------------------------------------------------------------------------!

  real, intent(in)                  :: t_tp
  real, intent(in), dimension(:,:)  :: z_off, z_km
  real, intent(out), dimension(:,:) :: teq
  real, dimension(size(z_km,1),size(z_km,2)) :: z_coord
  real :: t_1, t_sp, t_2
  real :: z_extra

  z_coord = z_km - z_off ! want to add offset to RHS of below comparisons, same as subtracting from LHS!
  z_extra = (216.65 - t_tp) / 2.8 ! if t_tp is colder than US std, offset height in equation down

  ! See: https://en.wikipedia.org/wiki/Barometric_formula#Source_code
  ! We truncate the lowest and highest legs -- no tropospheric lapse rate, and no
  ! secondary mesospheric lapse rate (because never need to model that high)
  t_1 = t_tp + 1.0 * (32.0 - 20.0)           ! beginning of +2.8K lapse rate
  t_sp = t_1 + 2.8 * (47.0 - 32.0 + z_extra) ! stratopause temperature
  t_2 = t_sp - 2.8 * (71.0 - 51.0)           ! end of -2.8K lapse rate
  where (z_coord <= 20) ! maybe has to be array here? try it out
    teq = t_tp
  elsewhere (z_coord <= 32) ! lapse rate 1K for 12km
    teq = t_tp + 1.0 * (z_coord - 20)
  elsewhere (z_coord <= 47 + z_extra)
    teq = t_1 + 2.8 * (z_coord - 32)
  elsewhere (z_coord <= 51 + z_extra)
    teq = t_sp
  elsewhere (z_coord <= 71 + z_extra)
    teq = t_sp - 2.8 * (z_coord - 51 - z_extra)
  elsewhere
    teq = t_2 - 2.0 * (z_coord - 71 - z_extra)
  endwhere

end subroutine tstd_summer

subroutine tstd_winter ( t_tp, vtx_gamma, z_vortex, z_km, teq )

  !----------------------------------------------------------------------------!
  ! Polar vortex adaptation
  ! Will preserve the original stratopause temperature by extending a 2.8K
  ! lapse rate region far enough to reach this temperature
  !----------------------------------------------------------------------------!

  real, intent(in)                  :: t_tp, vtx_gamma
  real, intent(in), dimension(:,:)  :: z_km, z_vortex
  real, intent(out), dimension(:,:) :: teq

  ! New version with simple vortex lapse rate
  where (z_km <= z_vortex)
    teq = t_tp
  elsewhere
    teq = t_tp - vtx_gamma * (z_km - z_vortex)
  endwhere

  ! Previous version with calculation done in pressure coords
  ! p_pkswitch = p0 * exp(-z_ozone / H)
  ! t_vtx = t_strat_usstd * (p_full(:,:,k) / p_pkswitch) ** pexp

  ! Previous version with US standard atmosphere-like transition zones
  ! z_offset     = (t_min - 216.65) / 2.8
  ! if (z_offset > 4) then
  !   call error_mesg ('forcing', 'Minimum stratosphere temperature too warm', FATAL)
  ! end if
  ! t_1 = t_min - vtx_gamma * (32 - z_tp_ref)
  ! z_vtx_recover = z_offset + (t_sp - t_1) / 2.8
  ! where (z_coord <= z_vortex)
  !   teq = t_min
  ! elsewhere (z_coord <= 32 + z_offset) ! lapse rate 1K for 12km
  !   teq = t_min - vtx_gamma * (z_coord - (z_tp_ref + z_offset))
  ! elsewhere (z_coord <= 36 + z_offset)
  !   teq = t_1
  ! elsewhere (z_coord <= 36 + z_vtx_recover)
  !   teq = t_1 + 2.8 * (z_coord - (36 + z_offset))
  ! elsewhere (z_coord <= 36 + z_vtx_recover + 4)
  !   teq = t_sp
  ! elsewhere
  !   teq = t_sp - 2.8 * (z_coord - (36 + z_vtx_recover + 4))
  ! endwhere

end subroutine tstd_winter

subroutine butler_heating ( x, y, tdt )

  !-----------------------------------------------------------------------
  !
  !           Impose Butler et al. and other forcing perturbations
  !           Input is 'x' (latitude in radians) and 'y' (sigma level)
  !
  !-----------------------------------------------------------------------
  real, intent(in),    dimension(:,:)   :: x
  real, intent(in),    dimension(:,:,:) :: y
  real, intent(inout), dimension(:,:,:) :: tdt
  !-----------------------------------------------------------------------
  integer :: k

  do k = 1, size(tdt,3)
    !     ----- Tropical forcing, mimics lapse rate feedback -----
    if (q0_tropical .ne. 0) then
      tdt(:,:,k) = tdt(:,:,k) + q0_tropical * exp( &
         -((x(:,:)   - x0_tropical) ** 2 / (2 * sx_tropical ** 2) + &
           (y(:,:,k) - y0_tropical) ** 2 / (2 * sy_tropical ** 2)) &
           ) / seconds_per_day
    endif

    !     ----- Polar vortex forcing, mimics ozone depletion -----
    if (q0_vortex .ne. 0) then
      tdt(:,:,k) = tdt(:,:,k) + q0_vortex * exp( &
        -2 * ((x(:,:)   - x0_vortex) ** 2 / (2 * sx_vortex ** 2) + &
            (y(:,:,k) - y0_vortex) ** 2 / (2 * sy_vortex ** 2)) &
            ) / seconds_per_day
    endif

    !     ----- Arctic forcing, mimics Arctic amplification -----
    ! Verify this works!
    if (q0_arctic .ne. 0) then
      where (x .gt. 0)
        tdt(:,:,k) = tdt(:,:,k) + q0_arctic &
          * cos(x - x0_arctic) ** 15 * exp(6 * (y(:,:,k) - y0_arctic)) &
          / seconds_per_day
      endwhere
    endif
  enddo

end subroutine butler_heating

subroutine lsp_heating ( lon, lat, p_full, tdt )
  !-----------------------------------------------------------------------
  !
  !           Impose Lindgren,Sheshadri,Plumb 2018 localised heating
  !           perturbation to excite planetary waves
  !           DOI: 10.1029/2018JD028537
  !
  !-----------------------------------------------------------------------
  real, intent(in),    dimension(:,:)   :: lon, lat
  real, intent(in),    dimension(:,:,:) :: p_full
  real, intent(inout), dimension(:,:,:) :: tdt
  real :: pbot, ptop, phi0_lsp, sphi_lsp
  !-----------------------------------------------------------------------
  ! Note model longitudes and latitudes are in radians!
  integer :: k

  phi0_lsp = lat0_lsp * pi / 180.0 ! latitude at which localised heating is centred
  sphi_lsp = slat_lsp * pi / 180.0 ! 0.175 radians default, presumably related to the width of heating in latitude (standard deviation in Gaussian distribution)
  pbot = p0_lsp * 100.0
  ptop = pt_lsp * 100.0
  do k=1,size(tdt,3)
    where ((p_full(:,:,k) <= pbot) .and. (p_full(:,:,k) >= ptop))
      tdt(:,:,k) = tdt(:,:,k) + &
        q0_lsp * sin(m_lsp * lon(:,:)) &
        * exp(-0.5 * ((lat - phi0_lsp) / sphi_lsp) ** 2) &
        * sin(pi * log(p_full(:,:,k) / pbot) / log(ptop / pbot)) &
        / seconds_per_day
    endwhere 
  enddo

end subroutine

subroutine friction_damping ( sigma, u, v, udt, vdt, rdamp, mask )

  !-----------------------------------------------------------------------
  !
  !           Rayleigh damping of wind components near surface
  !
  !-----------------------------------------------------------------------
  real, intent(in),  dimension(:,:,:) :: sigma, u, v
  real, intent(in),  dimension(:,:,:), optional :: mask
  real, intent(out), dimension(:,:,:,:) :: udt, vdt, rdamp
  !-----------------------------------------------------------------------
  real, dimension(size(u,2)) :: u_mean, v_mean
  real, dimension(size(u,1),size(u,2),2) :: u_decomp, v_decomp
  integer :: i, j, k, c
  !-----------------------------------------------------------------------

  do k=1,size(u,3)
    do c=1,2
      !    ---- Determine damping coeffs ----
      where (sigma(:,:,k) <= 1.0 .and. sigma(:,:,k) > sigma_b)
        rdamp(:,:,k,c) = vkfric(c) * (sigma(:,:,k) - sigma_b) / (1.0 - sigma_b)
      elsewhere
        rdamp(:,:,k,c) = 0.0
      endwhere

      !    ---- Apply damping ----
      if (.not. rdamp_decomp) then
        ! Damp the total winds, then exit
        udt(:,:,k,1)    = -rdamp(:,:,k,1) * u(:,:,k)
        vdt(:,:,k,1)    = -rdamp(:,:,k,1) * v(:,:,k)
        udt(:,:,k,2)    = 0.0
        vdt(:,:,k,2)    = 0.0
        rdamp(:,:,k,2) = rdamp(:,:,k,1)
        rdamp(:,:,k,2) = rdamp(:,:,k,1)
        exit ! i.e. break from top-level do loop
      else if (c==1) then
        ! Damp the means
        u_mean = sum(u(:,:,k),1) / size(u,1) ! mean
        v_mean = sum(v(:,:,k),1) / size(v,1)
        do i=1,size(u,1)
          u_decomp(i,:,1) = u_mean
          v_decomp(i,:,1) = v_mean
        enddo
        udt(:,:,k,1)    = -rdamp(:,:,k,1) * u_decomp(:,:,1)
        vdt(:,:,k,1)    = -rdamp(:,:,k,1) * v_decomp(:,:,1)
      else if (c==2) then
        ! Damp the anomaly
        u_decomp(:,:,2) = u(:,:,k) - u_decomp(:,:,1) ! anomaly
        v_decomp(:,:,2) = v(:,:,k) - v_decomp(:,:,1)
        udt(:,:,k,2)    = -rdamp(:,:,k,2) * u_decomp(:,:,2)
        vdt(:,:,k,2)    = -rdamp(:,:,k,2) * v_decomp(:,:,2)
      endif
    enddo
  enddo

  if (present(mask)) then
    do c=1,2
      udt(:,:,:,c) = udt(:,:,:,c) * mask
      vdt(:,:,:,c) = vdt(:,:,:,c) * mask
    enddo
  endif

end subroutine friction_damping

subroutine sponge_damping ( p_full, u, v, uspg, vspg, sdamp, mask )

  !-----------------------------------------------------------------------
  !
  !     Damp upper stratospheric winds with sponge
  !
  !-----------------------------------------------------------------------
  real, intent(in),  dimension(:,:,:) :: p_full, u, v
  real, intent(out), dimension(:,:,:) :: uspg, vspg, sdamp
  real, intent(in),  dimension(:,:,:), optional :: mask
  !-----------------------------------------------------------------------
  real, dimension(size(u,1),size(u,2)) :: sp_fact, spcoeff, ksp
  real    :: p_sp
  integer :: i, j, k
  !-----------------------------------------------------------------------

  p_sp = p_sponge * 100.0
  ksp  = -vksponge(1) ! for time being, no option to separate mean/anomaly sponge damping

  do k = 1, size(u,3)
    where (p_full(:,:,k) .lt. p_sp)
      sp_fact = (p_sp - p_full(:,:,k)) / p_sp
      spcoeff = ksp * sp_fact * sp_fact
      uspg(:,:,k)  = spcoeff * u(:,:,k) 
      vspg(:,:,k)  = spcoeff * v(:,:,k)
      sdamp(:,:,k) = spcoeff
    elsewhere
      uspg(:,:,k)  = 0.0
      vspg(:,:,k)  = 0.0
    endwhere
  enddo

  if (present(mask)) then
    uspg = uspg * mask
    vspg = vspg * mask
  endif

end subroutine sponge_damping

subroutine tracer_source_sink ( flux, damp, p_half, r, rdt, kbot )

  !-----------------------------------------------------------------------
  !
  !     Tracer control
  !     Simple surface source and global sink
  !
  !-----------------------------------------------------------------------
  real, intent(in)  :: flux, damp, p_half(:,:,:), r(:,:,:)
  real, intent(out) :: rdt(:,:,:)
  integer, intent(in), optional :: kbot(:,:)
  !-----------------------------------------------------------------------
  real, dimension(size(r,1),size(r,2),size(r,3)) :: source, sink
  real, dimension(size(r,1),size(r,2))           :: pmass
  integer :: i, j, kb
  real    :: rdamp
  !-----------------------------------------------------------------------

  rdamp = damp
  if (rdamp < 0.) rdamp = -seconds_per_day * rdamp   ! convert days to seconds
  if (rdamp > 0.) rdamp = 1.0 / rdamp

  source=0.0
  if (present(kbot)) then
    do j=1,size(r,2)
      do i=1,size(r,1)
        kb = kbot(i,j)
        pmass (i,j)    = p_half(i,j,kb+1) - p_half(i,j,kb)
        source(i,j,kb) = flux / pmass(i,j)
      enddo
    enddo
  else
    kb = size(r,3)
    pmass     = p_half(:,:,kb+1) - p_half(:,:,kb)
    source(:,:,kb) = flux / pmass
  endif

  sink = rdamp * r
  rdt = source-sink

end subroutine tracer_source_sink

end module forcing_mod
