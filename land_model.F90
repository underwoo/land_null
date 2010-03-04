module land_model_mod ! This is the null version

use  mpp_domains_mod, only : CYCLIC_GLOBAL_DOMAIN, mpp_get_compute_domain
use  mpp_domains_mod, only : domain2d, mpp_define_layout, mpp_define_domains
use  mpp_domains_mod, only : mpp_get_ntile_count, mpp_get_tile_id, mpp_get_current_ntile

use time_manager_mod, only : time_type

use diag_manager_mod, only : diag_axis_init

use tracer_manager_mod, only : register_tracers

use field_manager_mod , only : MODEL_LAND

use          fms_mod, only : write_version_number, error_mesg, FATAL, mpp_npes

use         grid_mod, only : get_grid_ntiles, get_grid_size, define_cube_mosaic
use         grid_mod, only : get_grid_cell_vertices, get_grid_cell_centers

implicit none
private

! ==== public interfaces =====================================================

public land_model_init          ! initialize the land model
public land_model_end           ! finish land model calculations
public land_model_restart       ! saves the land model restart(s)
public update_land_model_fast   ! time-step integration
public update_land_model_slow   ! time-step integration
public atmos_land_boundary_type ! data from coupler to land
public land_data_type           ! data from land to coupler
public :: Lnd_stock_pe          ! return stocks of conservative quantities

! ==== end of public interfaces ==============================================

character(len=*), parameter :: &
     version = '$Id: land_model.F90,v 18.0 2010/03/02 23:37:40 fms Exp $', &
     tagname = '$Name: riga $'

type :: atmos_land_boundary_type
   real, dimension(:,:,:), pointer :: & ! (lon, lat, tile)
        t_flux    => NULL(), &   ! sensible heat flux, W/m2
        lw_flux   => NULL(), &   ! net longwave radiation flux, W/m2
        lwdn_flux => NULL(), &   ! downward longwave radiation flux, W/m2
        sw_flux   => NULL(), &   ! net shortwave radiation flux, W/m2
        swdn_flux => NULL(), &   ! downward shortwave radiation flux, W/m2
        lprec     => NULL(), &   ! liquid precipitation rate, kg/(m2 s)
        fprec     => NULL(), &   ! frozen precipitation rate, kg/(m2 s)
        tprec     => NULL(), &   ! temperture of precipitation, degK
        sw_flux_down_vis_dir   => NULL(), & ! visible direct 
        sw_flux_down_total_dir => NULL(), & ! total direct
        sw_flux_down_vis_dif   => NULL(), & ! visible diffuse
        sw_flux_down_total_dif => NULL(), & ! total diffuse
        dhdt      => NULL(), &   ! sensible w.r.t. surface temperature
        dhdq      => NULL(), &   ! sensible w.r.t. surface humidity
        drdt      => NULL(), &   ! longwave w.r.t. surface radiative temperature 
        cd_m      => NULL(), &   ! drag coefficient for momentum, dimensionless
        cd_t      => NULL(), &   ! drag coefficient for tracers, dimensionless
        ustar     => NULL(), &   ! turbulent wind scale, m/s
        bstar     => NULL(), &   ! turbulent buoyancy scale, m/s
        wind      => NULL(), &   ! abs wind speed at the bottom of the atmos, m/s
        z_bot     => NULL(), &   ! height of the bottom atmospheric layer above the surface, m
        drag_q    => NULL(), &   ! product of cd_q by wind
        p_surf    => NULL()      ! surface pressure, Pa

   real, dimension(:,:,:,:), pointer :: & ! (lon, lat, tile, tracer)
        tr_flux => NULL(),   &   ! tracer flux, including water vapor flux
        dfdtr   => NULL()        ! derivative of the flux w.r.t. tracer surface value, 
                                 ! including evap over surface specific humidity

   integer :: xtype             !REGRID, REDIST or DIRECT
end type atmos_land_boundary_type


type :: land_data_type
   logical :: pe ! data presence indicator for stock calculations
   real, pointer, dimension(:,:,:)   :: &  ! (lon, lat, tile)
        tile_size      => NULL(),  & ! fractional coverage of cell by tile, dimensionless
        t_surf         => NULL(),  & ! ground surface temperature, degK
        t_ca           => NULL(),  & ! canopy air temperature, degK
        albedo         => NULL(),  & ! broadband land albedo [unused?]
        albedo_vis_dir => NULL(),  & ! albedo for direct visible radiation
        albedo_nir_dir => NULL(),  & ! albedo for direct NIR radiation 
        albedo_vis_dif => NULL(),  & ! albedo for diffuse visible radiation 
        albedo_nir_dif => NULL(),  & ! albedo for diffuse NIR radiation
        rough_mom      => NULL(),  & ! surface roughness length for momentum, m
        rough_heat     => NULL(),  & ! roughness length for tracers and heat, m
        rough_scale    => NULL()     ! topographic scaler for momentum drag, m

   real, pointer, dimension(:,:,:,:)   :: &  ! (lon, lat, tile, tracer)
        tr    => NULL()              ! tracers, including canopy air specific humidity

   real, pointer, dimension(:,:) :: &  ! (lon, lat)
     discharge           => NULL(),  & ! liquid water flux from land to ocean
     discharge_heat      => NULL(),  & ! sensible heat of discharge (0 C datum)
     discharge_snow      => NULL(),  & ! solid water flux from land to ocean
     discharge_snow_heat => NULL()     ! sensible heat of discharge_snow (0 C datum)

   logical, pointer, dimension(:,:,:):: &
        mask => NULL()               ! true if land

   integer :: axes(2)        ! IDs of diagnostic axes
   type(domain2d) :: domain  ! our computation domain
   logical, pointer :: maskmap(:,:) 
end type land_data_type

contains

! ============================================================================
subroutine land_model_init (cplr2land, land2cplr, time_init, time, dt_fast, dt_slow)

  type(atmos_land_boundary_type), intent(inout) :: cplr2land ! boundary data
  type(land_data_type)          , intent(inout) :: land2cplr ! boundary data
  type(time_type), intent(in) :: time_init ! initial time of simulation (?)
  type(time_type), intent(in) :: time      ! current time
  type(time_type), intent(in) :: dt_fast   ! fast time step
  type(time_type), intent(in) :: dt_slow   ! slow time step

 ! ---- local vars
  integer :: nlon, nlat ! size of global grid in lon and lat directions
  integer :: ntiles     ! number of tiles in the mosaic grid
  integer :: is,ie,js,je,id_lon,id_lat,i
  type(domain2d) :: domain
  real, allocatable, dimension(:,:)  :: glon, glat
  integer, allocatable, dimension(:) :: tile_ids
  integer :: ntracers, ntprog, ndiag
  integer :: layout(2) = (/0,0/)

  call write_version_number (version, tagname)

  ! define the processor layout information according to the global grid size 
  call get_grid_ntiles('LND',ntiles)
  call get_grid_size('LND',1,nlon,nlat)
  if( layout(1)==0 .AND. layout(2)==0 ) then
    call mpp_define_layout( (/1,nlon,1,nlat/), mpp_npes()/ntiles, layout )
  endif
  if( layout(1)/=0 .AND. layout(2)==0 )layout(2) = mpp_npes()/(layout(1)*ntiles)
  if( layout(1)==0 .AND. layout(2)/=0 )layout(1) = mpp_npes()/(layout(2)*ntiles)
  
  ! defne land model domain
  if (ntiles==1) then
     call mpp_define_domains ((/1,nlon, 1, nlat/), layout, domain, xhalo=1, yhalo=1,&
          xflags = CYCLIC_GLOBAL_DOMAIN, name = 'LAND MODEL')
  else
     call define_cube_mosaic ('LND', domain, layout, halo=1)
  endif
  land2cplr%domain = domain

  call mpp_get_compute_domain(domain, is,ie,js,je)
  allocate(land2cplr%tile_size(is:ie,js:je,1))
  land2cplr%tile_size = 0.0

  allocate(tile_ids(mpp_get_current_ntile(domain)))
  tile_ids = mpp_get_tile_id(domain)
  allocate(glon(nlon,nlat), glat(nlon,nlat))
  call get_grid_cell_centers ('LND', tile_ids(1), glon,  glat)

  if(mpp_get_ntile_count(domain)==1) then
     ! grid has just one tile, so we assume that the grid is regular lat-lon
     ! define longitude axes and its edges
     id_lon  = diag_axis_init('lon', glon(:,1), 'degrees_E', 'X', &
          long_name='longitude', set_name='land', domain2=domain )

     ! define latitude axes and its edges
     id_lat = diag_axis_init ('lat', glat(1,:), 'degrees_N', 'Y', &
          long_name='latitude',  set_name='land', domain2=domain )
  else
     id_lon = diag_axis_init ( 'grid_xt', (/(real(i),i=1,nlon)/), 'degrees_E', 'X', &
          long_name='T-cell longitude', set_name='land', domain2=domain, aux='geolon_t' )
     id_lat = diag_axis_init ( 'grid_yt', (/(real(i),i=1,nlat)/), 'degrees_N', 'Y', &
          long_name='T-cell latitude',  set_name='land', domain2=domain, aux='geolat_t' )
  endif
  land2cplr%axes = (/id_lon,id_lat/)

  call register_tracers(MODEL_LAND, ntracers, ntprog, ndiag)

  allocate(land2cplr%mask          (is:ie,js:je,1))
  allocate(land2cplr%t_surf        (is:ie,js:je,1))
  allocate(land2cplr%t_ca          (is:ie,js:je,1))
  allocate(land2cplr%tr            (is:ie,js:je,1,ntprog))
  allocate(land2cplr%albedo        (is:ie,js:je,1))
  allocate(land2cplr%albedo_vis_dir(is:ie,js:je,1))
  allocate(land2cplr%albedo_nir_dir(is:ie,js:je,1))
  allocate(land2cplr%albedo_vis_dif(is:ie,js:je,1))
  allocate(land2cplr%albedo_nir_dif(is:ie,js:je,1))
  allocate(land2cplr%rough_mom     (is:ie,js:je,1))
  allocate(land2cplr%rough_heat    (is:ie,js:je,1))
  allocate(land2cplr%rough_scale   (is:ie,js:je,1))
  allocate(land2cplr%discharge     (is:ie,js:je))
  allocate(land2cplr%discharge_heat(is:ie,js:je))
  allocate(land2cplr%discharge_snow(is:ie,js:je))
  allocate(land2cplr%discharge_snow_heat(is:ie,js:je))

  land2cplr%mask              = .FALSE.
  land2cplr%t_surf            = 280.0
  land2cplr%t_ca              = 280.0
  land2cplr%tr                = 0.0
  land2cplr%albedo            = 0.0
  land2cplr%albedo_vis_dir    = 0.0
  land2cplr%albedo_nir_dir    = 0.0
  land2cplr%albedo_vis_dif    = 0.0
  land2cplr%albedo_nir_dif    = 0.0
  land2cplr%rough_mom         = 0.0
  land2cplr%rough_heat        = 0.0
  land2cplr%rough_scale       = 0.0
  land2cplr%discharge         = 0.0
  land2cplr%discharge_heat    = 0.0
  land2cplr%discharge_snow    = 0.0
  land2cplr%discharge_snow_heat = 0.0

  allocate(cplr2land%t_flux(is:ie,js:je,1) )
  allocate(cplr2land%lw_flux(is:ie,js:je,1) )
  allocate(cplr2land%sw_flux(is:ie,js:je,1) )
  allocate(cplr2land%lprec(is:ie,js:je,1) )
  allocate(cplr2land%fprec(is:ie,js:je,1) )
  allocate(cplr2land%tprec(is:ie,js:je,1) )
  allocate(cplr2land%dhdt(is:ie,js:je,1) )
  allocate(cplr2land%dhdq(is:ie,js:je,1) )
  allocate(cplr2land%drdt(is:ie,js:je,1) )
  allocate(cplr2land%p_surf(is:ie,js:je,1) )
  allocate(cplr2land%tr_flux(is:ie,js:je,1,ntprog) )
  allocate(cplr2land%dfdtr(is:ie,js:je,1,ntprog) )
  allocate(cplr2land%lwdn_flux(is:ie,js:je,1) )
  allocate(cplr2land%swdn_flux(is:ie,js:je,1) )
  allocate(cplr2land%sw_flux_down_vis_dir(is:ie,js:je,1) )
  allocate(cplr2land%sw_flux_down_total_dir(is:ie,js:je,1) )
  allocate(cplr2land%sw_flux_down_vis_dif(is:ie,js:je,1) )
  allocate(cplr2land%sw_flux_down_total_dif(is:ie,js:je,1) )
  allocate(cplr2land%cd_t(is:ie,js:je,1) )
  allocate(cplr2land%cd_m(is:ie,js:je,1) )
  allocate(cplr2land%bstar(is:ie,js:je,1) )
  allocate(cplr2land%ustar(is:ie,js:je,1) )
  allocate(cplr2land%wind(is:ie,js:je,1) )
  allocate(cplr2land%z_bot(is:ie,js:je,1) )
  allocate(cplr2land%drag_q(is:ie,js:je,1) )

  cplr2land%t_flux                 = 0.0
  cplr2land%lw_flux                = 0.0
  cplr2land%sw_flux                = 0.0
  cplr2land%lprec                  = 0.0
  cplr2land%fprec                  = 0.0
  cplr2land%tprec                  = 0.0
  cplr2land%dhdt                   = 0.0
  cplr2land%dhdq                   = 0.0
  cplr2land%drdt                   = 0.0
  cplr2land%p_surf                 = 1.0e5
  cplr2land%tr_flux                = 0.0
  cplr2land%dfdtr                  = 0.0
  cplr2land%lwdn_flux              = 0.0
  cplr2land%swdn_flux              = 0.0
  cplr2land%sw_flux_down_vis_dir   = 0.0
  cplr2land%sw_flux_down_total_dir = 0.0
  cplr2land%sw_flux_down_vis_dif   = 0.0
  cplr2land%sw_flux_down_total_dif = 0.0
  cplr2land%cd_t                   = 0.0
  cplr2land%cd_m                   = 0.0
  cplr2land%bstar                  = 0.0
  cplr2land%ustar                  = 0.0
  cplr2land%wind                   = 0.0
  cplr2land%z_bot                  = 0.0
  cplr2land%drag_q                 = 0.0

  deallocate(glon, glat, tile_ids)

end subroutine land_model_init

! ============================================================================
subroutine land_model_end (cplr2land, land2cplr)
  type(atmos_land_boundary_type), intent(inout) :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr
  
end subroutine land_model_end

! ============================================================================
subroutine land_model_restart(timestamp)
  character(*), intent(in), optional :: timestamp
  
end subroutine land_model_restart

! ============================================================================
subroutine update_land_model_fast ( cplr2land, land2cplr )
  type(atmos_land_boundary_type), intent(in)    :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  call error_mesg('update_land_model_fast','Should not be calling null version of update_land_model_fast',FATAL)

end subroutine update_land_model_fast


! ============================================================================
subroutine update_land_model_slow ( cplr2land, land2cplr )
  type(atmos_land_boundary_type), intent(inout) :: cplr2land
  type(land_data_type)          , intent(inout) :: land2cplr

  call error_mesg('update_land_model_slow','Should not be calling null version of update_land_model_slow',FATAL)

end subroutine update_land_model_slow

! ============================================================================
subroutine Lnd_stock_pe(bnd,index,value)

type(land_data_type), intent(in)  :: bnd 
integer             , intent(in)  :: index
real                , intent(out) :: value ! Domain water (Kg) or heat (Joules)

value = 0.0
 
end subroutine Lnd_stock_pe
! ============================================================================

end module land_model_mod
