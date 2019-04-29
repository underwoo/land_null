program land_model_init_test
  use time_manager_mod
  use land_model_mod

  type(atmos_land_boundary_type) :: cplr2land
  type(land_data_type)           :: land2cplr
  type(time_type) :: time_init ! initial time of simulation (?)
  type(time_type) :: time      ! current time
  type(time_type) :: dt_fast   ! fast time step
  type(time_type) :: dt_slow   ! slow time step

  call land_model_init(cplr2land, land2cplr, time_init, time, dt_fast, dt_slow)
end program land_model_init_test
