module soil_mod

use time_manager_mod, only : time_type

IMPLICIT NONE

public send_averaged_data   ! sends tile-averaged data to diagnostics manager

! ---- interface to tile-averaged diagnostic routines ------------------------
interface send_averaged_data
   module procedure send_averaged_data2d
   module procedure send_averaged_data3d
end interface

contains

! ============================================================================
function send_averaged_data2d ( id, field, area, time, mask )
! ============================================================================
! average the data over tiles and send then to diagnostics

  ! --- return value ---------------------------------------------------------
  logical                      :: send_averaged_data2d
  ! --- arguments ------------------------------------------------------------
  integer, intent(in)          :: id             ! id od the diagnostic field 
  real,    intent(in)          :: field(:,:,:)   ! field to average and send
  real,    intent(in)          :: area (:,:,:)   ! area of tiles (== averaging 
                                                 ! weights), arbitrary units
  type(time_type), intent(in)  :: time           ! current time
  logical, intent(in),optional :: mask (:,:,:)   ! land mask

  ! --- local vars -----------------------------------------------------------

  send_averaged_data2d = .false.

end function send_averaged_data2d

! ============================================================================
function send_averaged_data3d( id, field, area, time, mask )
! ============================================================================
! average the data over tiles and send then to diagnostics

  ! --- return value ---------------------------------------------------------
  logical                      :: send_averaged_data3d
  ! --- arguments ------------------------------------------------------------
  integer, intent(in)          :: id              ! id of the diagnostic field
  real,    intent(in)          :: field(:,:,:,:)  ! (lon, lat, tile, lev) field 
                                                  ! to average and send
  real,    intent(in)          :: area (:,:,:)    ! (lon, lat, tile) tile areas 
                                                  ! ( == averaging weights), 
                                                  ! arbitrary units
  type(time_type), intent(in)  :: time            ! current time
  logical, intent(in),optional :: mask (:,:,:)    ! (lon, lat, tile) land mask

  ! --- local vars -----------------------------------------------------------

  send_averaged_data3d = .false.

end function send_averaged_data3d

end module soil_mod
