! Copyright 2017- LabTerra

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.)

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

module soil_pools

    use types

    implicit none
    ! CONSTANTS and POOLS FOR CARBON DECAY

    !litter I (1) litter II (2) soilI (3) soil II (4)
    real(r_4), dimension(1), public :: sp_available_p = 0.204299955    ! g m-2 Yang et al., 2013
    real(r_4), dimension(1), public :: sp_available_n = 0.30775999     ! g m-2 Xu et al. 2013 ?
    real(r_4), dimension(1), public :: sp_so_p = 0.0
    real(r_4), dimension(1), public :: sp_in_p = 0.0
    real(r_4), dimension(4), public :: sp_csoil = (/500.0, 500.0, 500.0, 500.0/)
    real(r_4), dimension(8), public :: sp_snr  = 0.0

 end module soil_pools
