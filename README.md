GFDL dynamical core mods
========================

This repository contains changes to the "dry" dynamical core spectral model from the
GFDL Flexible Modeling System. The original model code can be found
[here](https://github.com/lukelbd/gfdl-fms-orig).

Additions
---------

Various features were added to the model, described in detail below and in the namelist
file [input.nml](exp/input.nml). Note that the `hs_forcing_mod` module and
`hs_forcing.f90` file are renamed to `forcing_mod` and
[forcing.f90](src/atmos_param/forcing/forcing.f90), and the `gaussian_topog_mod` module
and `gaussian_topog.f90` file are renamed to `idealized_topog_mod` and
[idealized_topog.f90](src/shared/topography/idealized_topog.f90). This renaming reflects
the availability of new model features. Also, the `forcing` subroutine now accepts
longitudes as well as latitudes.

* [vert_coordinate.f90](src/atmos_spectral/init/vert_coordinate.f90) implements the
  Polvani and Kushner method for constructing levels with high stratospheric resolution
  directly. It can be requested using `vert_coord_option='pk_sigma'`, analogous to
  `vert_coord_option='even_sigma'`, in the `&spectral_dynamics_nml` namelist (see
  [input.nml](exp/input.nml)).
* ECMWF hybrid coordinate data is located in [levels](exp/levels), and implemented as a
  `&vert_coordinate_nml` namelist in the [levels_ecmwf_60.nml](exp/levels_ecmwf_60.nml)
  file. This namelist can be copied into [input.nml](exp/input.nml) and applied by
  setting `vert_coord_option='input'` in the `&spectral_dynamics_nml` namelist (see
  below).
* [forcing.f90](src/atmos_param/forcing/forcing.f90) adds options to apply static
  diabatic heating perturbations from Butler at el. (2010) and Lindgren et al. (2018)
  (various `&forcing_nml` options), and the option to load locked diabatic heating from
  a file instead of using Newtonian damping (`locked_heating=.true.`).
* [forcing.f90](src/atmos_param/forcing/forcing.f90) adds new forcing options other than
  Held and Suarez (1994): the Polvani and Kushner (2002) `'pk'` scheme
  (`teq_mode='pk'`), and a *modified* Polvani and Kushner scheme (`teq_mode='pkmod'`),
  with stratospheric equilibrium temperatures conforming to the height at which the
  minimum tropospheric equilibrium temperature is reached. It also includes an option
  to change the damping timescale in the stratosphere (`damp_mode='pk'` or
  `damp_mode='pkmod'`), beginning where the minimum tropospheric equilibrium
  temperature is reached.
* [forcing.f90](src/atmos_param/forcing/forcing.f90) adds various options for modifying
  the damping scheme: specifying surface equilibrium temperature from its global mean
  surface value rather than its equatorial surface value (`surf_shneider=.true.` and
  e.g. `t_mean=300` instead of `surf_schneider=.false.` and e.g. `t_0=315`), changing
  the exponent to which the cosine latitude is raised to control how boundary layer
  damping intensity changes with latitude (e.g. `exp_b=4`), and changing the latitude of
  the most intense meridional equilibrium temperature gradients (e.g. `exp_h=-2` for
  equatorward motion, `exp_h=2` for poleward motion).
* [forcing.f90](src/atmos_param/forcing/forcing.f90) adds options for imposing stronger 
  tropopause inversion layer (TIL) - Boljka and Birner (in prep.).  
* [idealized_topog.f90](src/shared/topography/idealized_topog.f90) adds a new sinusoidal
  topography option, configurable with the `&idealized_topog_nml` namelist.

Documentation
-------------

The full, original model documentation is published on
[github-pages](https://lukelbd.github.io/gfdl-fms). Use the links at the top of each
page to navigate between sections! Below are some important notes regarding the
documentation.

* [bgrid.pdf](https://lukelbd.github.io/gfdl-fms/doc/bgrid.pdf),
  [spectral_core.pdf](https://lukelbd.github.io/gfdl-fms/doc/spectral_core.pdf),
  [shallow.pdf](https://lukelbd.github.io/gfdl-fms/doc/shallow.pdf),
  [barotropic.pdf](https://lukelbd.github.io/gfdl-fms/doc/barotropic.pdf), contain
  background information on the model physical equations.
* [diag_table](exp/diag_table) is a template "diagnostic table". It specifies the model
  output NetCDF files, time-averaging behavior, and variable names. A brief description
  of this file format in found in the [diag_table](exp/diag_table) comments. Full lists
  of diagnostic table parameters are found on the following pages.
    * [barotropic_diagnostics.html](https://lukelbd.github.io/gfdl-fms/src/atmos_spectral_barotropic/barotropic_diagnostics.html)
    * [bgrid_diagnostics.html](https://lukelbd.github.io/gfdl-fms/src/atmos_bgrid/tools/bgrid_diagnostics.html)
    * [forcing.html](https://lukelbd.github.io/gfdl-fms/src/atmos_param/forcing/forcing.html)
    * [shallow_diagnostics.html](https://lukelbd.github.io/gfdl-fms/src/atmos_spectral_shallow/shallow_diagnostics.html)
    * [spectral_dynamics.html](https://lukelbd.github.io/gfdl-fms/src/atmos_spectral/model/spectral_dynamics.html)
* [input.nml](exp/input.nml) is a template namelist file. It allows you to specify
  various model settings and forcing parameters at run-time, rather than compile-time.
  Brief descriptions of some important parameters are found in the
  [input.nml](exp/input.nml) comments. Full tables and descriptions for all namelist
  parameters are found on the following pages:
    * [atmos_model.html](https://lukelbd.github.io/gfdl-fms/src/atmos_solo/atmos_model.html)
    * [atmos_bgrid/atmosphere.html](https://lukelbd.github.io/gfdl-fms/src/atmos_bgrid/driver/solo/atmosphere.html)
    * [atmos_bgrid_shallow/atmosphere.html](https://lukelbd.github.io/gfdl-fms/src/atmos_bgrid/driver/shallow/atmosphere.html)
    * [atmos_spectral_barotropic/atmosphere.html](https://lukelbd.github.io/gfdl-fms/src/atmos_spectral_barotropic/atmosphere.html)
    * [atmos_spectral_shallow/atmosphere.html](https://lukelbd.github.io/gfdl-fms/src/atmos_spectral_shallow/atmosphere.html)
    * [barotropic_dynamics.html](https://lukelbd.github.io/gfdl-fms/src/atmos_spectral_barotropic/barotropic_dynamics.html)
    * [barotropic_physics.html](https://lukelbd.github.io/gfdl-fms/src/atmos_spectral_barotropic/barotropic_physics.html)
    * [bgrid_advection.html](https://lukelbd.github.io/gfdl-fms/src/atmos_bgrid/model/bgrid_advection.html)
    * [bgrid_cold_start.html](https://lukelbd.github.io/gfdl-fms/src/atmos_bgrid/tools/bgrid_cold_start.html)
    * [bgrid_core_driver.html](https://lukelbd.github.io/gfdl-fms/src/atmos_bgrid/model/bgrid_core_driver.html)
    * [bgrid_horiz_diff.html](https://lukelbd.github.io/gfdl-fms/src/atmos_bgrid/model/bgrid_horiz_diff.html)
    * [bgrid_integrals.html](https://lukelbd.github.io/gfdl-fms/src/atmos_bgrid/tools/bgrid_integrals.html)
    * [bgrid_sponge.html](https://lukelbd.github.io/gfdl-fms/src/atmos_bgrid/model/bgrid_sponge.html)
    * [diag_manager.html](https://lukelbd.github.io/gfdl-fms/src/shared/diag_manager/diag_manager.html)
    * [fms.html](https://lukelbd.github.io/gfdl-fms/src/shared/fms/fms.html)
    * [forcing.html](https://lukelbd.github.io/gfdl-fms/src/atmos_param/forcing/forcing.html)
    * [horiz_interp_spherical.html](https://lukelbd.github.io/gfdl-fms/src/shared/horiz_interp/horiz_interp_spherical.html)
    * [idealized_topog.html](https://lukelbd.github.io/gfdl-fms/src/shared/topography/idealized_topog.html)
    * [shallow_dynamics.html](https://lukelbd.github.io/gfdl-fms/src/atmos_spectral_shallow/shallow_dynamics.html)
    * [shallow_physics.html](https://lukelbd.github.io/gfdl-fms/src/atmos_param/shallow_physics/shallow_physics.html)
    * [shallow_physics.html](https://lukelbd.github.io/gfdl-fms/src/atmos_spectral_shallow/shallow_physics.html)
    * [spectral_dynamics.html](https://lukelbd.github.io/gfdl-fms/src/atmos_spectral/model/spectral_dynamics.html)
    * [topography.html](https://lukelbd.github.io/gfdl-fms/src/shared/topography/topography.html)

Input
-----

The model can be given two types of input files in the `INPUT` directory.

* Topography file `topography.data.nc`, containing a `zsurf` variable with longitude and
  latitude dimensions.
* Locked heating file `heating.data.nc`, containing a `t` variable or `tdt` variable
  with longitude, latitude, and level dimensions. This heating replaces the thermal
  damping term.

Truncation
----------

The following is a table of valid truncation numbers with the corresponding recommended
number of latitudes and cores for parallelization. Remember that to avoid aliasing of
the harmonics, `N >= (3M + 1)/2`, where `N` is the number of latitudes and `M` is the
truncation number.

| Truncation | Number of latitudes | Recommended cores (latitudes per file) |
| ---        | ---                 | ---                                    |
| 42         | 64                  | 16 (4)                                 |
| 63         | 96                  | 24 (4)                                 |
| 85         | 128                 | 32 (4)                                 |
| 106        | 160                 | 40 (4), 80 (2)                         |
| 159        | 240                 | 60 (4)                                 |
| 170        | 256                 | 64 (4)                                 |
| 266        | 400                 | 100 (4)                                |

