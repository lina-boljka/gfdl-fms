module atmosphere_mod

use                  fms_mod, only: set_domain, write_version_number, field_size, file_exist, &
                                    mpp_pe, mpp_root_pe, error_mesg, FATAL, read_data, write_data, nullify_domain

use            constants_mod, only: grav, pi

use           transforms_mod, only: trans_grid_to_spherical, trans_spherical_to_grid, &
                                    get_deg_lat, get_grid_boundaries, grid_domain,    &
                                    spectral_domain, get_grid_domain, get_lon_max, &
                                    get_lat_max, get_deg_lon

use         time_manager_mod, only: time_type, set_time, get_time, &
                                    operator( + ), operator( - ), operator( < )

use     press_and_geopot_mod, only: compute_pressures_and_heights

use    spectral_dynamics_mod, only: spectral_dynamics_init, spectral_dynamics, spectral_dynamics_end, &
                                    get_num_levels, get_axis_id, spectral_diagnostics, get_initial_fields, &
                                    complete_robert_filter

use          tracer_type_mod, only: tracer_type

use           forcing_mod, only: forcing_init, forcing

use        field_manager_mod, only: MODEL_ATMOS

use       tracer_manager_mod, only: get_number_tracers

implicit none
private
!=================================================================================================================================

character(len=128) :: version= &
'$Id: atmosphere.f90,v 13.0 2006/03/28 21:17:28 fms Exp $'
      
character(len=128) :: tagname= &
'$Name: latest $'
character(len=10), parameter :: mod_name='atmosphere'

!=================================================================================================================================

public :: atmosphere_init, atmosphere, atmosphere_end

!=================================================================================================================================

integer, parameter :: num_time_levels = 2
integer :: is, ie, js, je, num_levels, num_tracers, nhum
logical :: dry_model

real, allocatable, dimension(:,:,:) :: p_half, p_full
real, allocatable, dimension(:,:,:) :: z_half, z_full

type(tracer_type), allocatable, dimension(:) :: tracer_attributes
real, allocatable, dimension(:,:,:,:,:) :: grid_tracers
real, allocatable, dimension(:,:,:    ) :: psg, wg_full
real, allocatable, dimension(:,:,:,:  ) :: ug, vg, tg

real, allocatable, dimension(:,:    ) :: dt_psg
real, allocatable, dimension(:,:,:  ) :: dt_ug, dt_vg, dt_tg
real, allocatable, dimension(:,:,:,:) :: dt_tracers

real, allocatable, dimension(:)   :: deg_lat, rad_lat, deg_lon
real, allocatable, dimension(:,:) :: rad_lat_2d, rad_lon_2d

integer :: previous, current, future
logical :: module_is_initialized =.false.
character(len=4) :: ch_tmp1, ch_tmp2

integer         :: dt_integer
real            :: dt_real
type(time_type) :: Time_step

integer, dimension(4) :: axis_id

!=================================================================================================================================
contains
!=================================================================================================================================

subroutine atmosphere_init(Time_init, Time, Time_step_in)

type (time_type), intent(in)  :: Time_init, Time, Time_step_in

integer :: seconds, days, lon_max, lat_max, ntr, nt, i, j
integer, dimension(4) :: siz
real, dimension(2) :: time_pointers
character(len=64) :: file, tr_name
character(len=256) :: message

if(module_is_initialized) return

call write_version_number(version, tagname)

Time_step = Time_step_in
call get_time(Time_step, seconds, days)
dt_integer   = 86400*days + seconds
dt_real      = float(dt_integer)

call get_number_tracers(MODEL_ATMOS, num_prog=num_tracers)
allocate (tracer_attributes(num_tracers))
call spectral_dynamics_init(Time, Time_step, tracer_attributes, dry_model, nhum)
call get_grid_domain(is, ie, js, je)
call get_num_levels(num_levels)

allocate (p_half       (is:ie, js:je, num_levels+1))
allocate (z_half       (is:ie, js:je, num_levels+1))
allocate (p_full       (is:ie, js:je, num_levels))
allocate (z_full       (is:ie, js:je, num_levels))
allocate (wg_full      (is:ie, js:je, num_levels))
allocate (psg          (is:ie, js:je, num_time_levels))
allocate (ug           (is:ie, js:je, num_levels, num_time_levels))
allocate (vg           (is:ie, js:je, num_levels, num_time_levels))
allocate (tg           (is:ie, js:je, num_levels, num_time_levels))
allocate (grid_tracers (is:ie, js:je, num_levels, num_time_levels, num_tracers ))

allocate (dt_psg     (is:ie, js:je))
allocate (dt_ug      (is:ie, js:je, num_levels))
allocate (dt_vg      (is:ie, js:je, num_levels))
allocate (dt_tg      (is:ie, js:je, num_levels))
allocate (dt_tracers (is:ie, js:je, num_levels, num_tracers))

allocate (deg_lat    (       js:je))
allocate (rad_lat    (       js:je))
allocate (rad_lat_2d (is:ie, js:je))
allocate (deg_lon    (is:ie       ))
allocate (rad_lon_2d (is:ie, js:je))

p_half = 0.; z_half = 0.; p_full = 0.; z_full = 0.
wg_full = 0.; psg = 0.; ug = 0.; vg = 0.; tg = 0.; grid_tracers = 0.
dt_psg = 0.; dt_ug  = 0.; dt_vg  = 0.; dt_tg  = 0.; dt_tracers = 0.
!--------------------------------------------------------------------------------------------------------------------------------
file = 'INPUT/atmosphere.res.nc'
if(file_exist(trim(file))) then
  call get_lon_max(lon_max)
  call get_lat_max(lat_max)
  call field_size(trim(file), 'ug', siz)
  if(lon_max /= siz(1) .or. lat_max /= siz(2)) then
    write(message,*) 'Resolution of restart data does not match resolution specified on namelist. Restart data: lon_max=', &
                     siz(1),', lat_max=',siz(2),'  Namelist: lon_max=',lon_max,', lat_max=',lat_max
    call error_mesg('atmosphere_init', message, FATAL)
  endif
  call nullify_domain()
! call read_data(trim(file), 'previous', previous)           ! No interface of read_data exists to read integer scalars
! call read_data(trim(file), 'current',  current)            ! No interface of read_data exists to read integer scalars
  call read_data(trim(file), 'time_pointers', time_pointers) ! Getaround for no interface to read integer scalars
  previous = int(time_pointers(1))                           ! Getaround for no interface to read integer scalars
  current  = int(time_pointers(2))                           ! Getaround for no interface to read integer scalars
  do nt=1,num_time_levels
    call read_data(trim(file), 'ug',   ug(:,:,:,nt), grid_domain, timelevel=nt)
    call read_data(trim(file), 'vg',   vg(:,:,:,nt), grid_domain, timelevel=nt)
    call read_data(trim(file), 'tg',   tg(:,:,:,nt), grid_domain, timelevel=nt)
    call read_data(trim(file), 'psg', psg(:,:,  nt), grid_domain, timelevel=nt)
    do ntr = 1,num_tracers
      tr_name = trim(tracer_attributes(ntr)%name)
      call read_data(trim(file), trim(tr_name), grid_tracers(:,:,:,nt,ntr), grid_domain, timelevel=nt)
    enddo ! end loop over tracers
  enddo ! end loop over time levels
  call read_data(trim(file), 'wg_full', wg_full, grid_domain)
else
  previous = 1; current = 1
  call get_initial_fields(ug(:,:,:,1), vg(:,:,:,1), tg(:,:,:,1), psg(:,:,1), grid_tracers(:,:,:,1,:))
endif
!--------------------------------------------------------------------------------------------------------------------------------
if(dry_model) then
  call compute_pressures_and_heights(tg(:,:,:,current), psg(:,:,current), z_full, z_half, p_full, p_half)
else
  call compute_pressures_and_heights( &
       tg(:,:,:,current), psg(:,:,current), z_full, z_half, p_full, p_half, grid_tracers(:,:,:,current,nhum))
endif

call get_deg_lat(deg_lat)
do j=js,je
  rad_lat_2d(:,j) = deg_lat(j)*pi/180.
enddo

call get_deg_lon(deg_lon)
do i=is, ie
  rad_lon_2d(i,:) = deg_lon(i)*pi/180.
enddo

call forcing_init(get_axis_id(), Time)

module_is_initialized = .true.

return
end subroutine atmosphere_init

!=================================================================================================================================

subroutine atmosphere(Time)
type(time_type), intent(in) :: Time

real    :: delta_t
type(time_type) :: Time_next

if(.not.module_is_initialized) then
  call error_mesg('atmosphere','atmosphere module is not initialized',FATAL)
endif

dt_ug  = 0.0
dt_vg  = 0.0
dt_tg  = 0.0
dt_psg = 0.0
dt_tracers = 0.0

if(current == previous) then
  delta_t = dt_real
else
  delta_t = 2*dt_real
endif

Time_next = Time + Time_step

call forcing(1, ie-is+1, 1, je-js+1, delta_t, Time_next, rad_lat_2d(:,:), rad_lon_2d(:,:), &
                p_half(:,:,:         ),       p_full(:,:,:           ), &
                z_full(:,:,:         ),                                 &
                    ug(:,:,:,previous),           vg(:,:,:,previous  ), &
                    tg(:,:,:,previous), grid_tracers(:,:,:,previous,:), &
                    ug(:,:,:,previous),           vg(:,:,:,previous  ), &
                    tg(:,:,:,previous), grid_tracers(:,:,:,previous,:), &
                 dt_ug(:,:,:         ),        dt_vg(:,:,:           ), &
                 dt_tg(:,:,:         ),   dt_tracers(:,:,:,:))

if(previous == current) then
  future = num_time_levels + 1 - current
else
  future = previous
endif

call spectral_dynamics(Time, psg(:,:,future), ug(:,:,:,future), vg(:,:,:,future), &
                       tg(:,:,:,future), tracer_attributes, grid_tracers(:,:,:,:,:), future, &
                       dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers, wg_full, p_full, p_half, z_full)

call complete_robert_filter(tracer_attributes)

! I think this is a good spot to include time varying surface geopotential

if(dry_model) then
  call compute_pressures_and_heights(tg(:,:,:,future), psg(:,:,future), z_full, z_half, p_full, p_half)
else
  call compute_pressures_and_heights( &
     tg(:,:,:,future), psg(:,:,future), z_full, z_half, p_full, p_half, grid_tracers(:,:,:,future,nhum))
endif

call spectral_diagnostics(Time_next, psg(:,:,future), ug(:,:,:,future), vg(:,:,:,future), &
                          tg(:,:,:,future), wg_full, grid_tracers(:,:,:,:,:), future)

previous = current
current  = future

return
end subroutine atmosphere

!=================================================================================================================================

subroutine atmosphere_end
integer :: ntr, nt
character(len=64) :: file, tr_name

if(.not.module_is_initialized) return

file='RESTART/atmosphere.res'
call nullify_domain()
!call write_data(trim(file), 'previous', previous) ! No interface exists to write a scalar
!call write_data(trim(file), 'current',  current)  ! No interface exists to write a scalar
call write_data(trim(file), 'time_pointers', (/real(previous),real(current)/)) ! getaround for no interface to write a scalar
do nt=1,num_time_levels
  call write_data(trim(file), 'ug',   ug(:,:,:,nt), grid_domain)
  call write_data(trim(file), 'vg',   vg(:,:,:,nt), grid_domain)
  call write_data(trim(file), 'tg',   tg(:,:,:,nt), grid_domain)
  call write_data(trim(file), 'psg', psg(:,:,  nt), grid_domain)
  do ntr = 1,num_tracers
    tr_name = trim(tracer_attributes(ntr)%name)
    call write_data(trim(file), tr_name, grid_tracers(:,:,:,nt,ntr), grid_domain)
  enddo
enddo
call write_data(trim(file), 'wg_full', wg_full, grid_domain)

deallocate (p_half, z_half, p_full, z_full, wg_full, psg, ug, vg, tg, grid_tracers)
deallocate (dt_psg, dt_ug, dt_vg, dt_tg, dt_tracers)
deallocate (deg_lat, rad_lat, rad_lat_2d)

call set_domain(grid_domain)
call spectral_dynamics_end(tracer_attributes)
deallocate(tracer_attributes)

module_is_initialized = .false.

end subroutine atmosphere_end

!=================================================================================================================================

end module atmosphere_mod
