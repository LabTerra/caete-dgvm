! Copyright 2017- LabTerra

!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.

!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.

!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.

module types
   implicit none


   integer,parameter,public :: l_1 = 4  ! standart Logical type
   integer,parameter,public :: i_2 = 2  ! 16 bits integer
   integer,parameter,public :: i_4 = 4  ! 32 bits integer
   integer,parameter,public :: r_4 = 4  ! 32 bits float
   integer,parameter,public :: r_8 = 8  ! 64 bits float

end module types


module global_par
   use types
   implicit none

   real(r_4),parameter,public :: h = 1.0                         ! soil layer thickness (meters)
   real(r_4),parameter,public :: diffu = 1.036800e14             ! soil thermal diffusivity (m2/mes)
   real(r_4),parameter,public :: tau = (h**2)/(2.0*diffu)        ! e-folding times (months)
   real(r_4),parameter,public :: rcmax = 2500.0                  ! ResistÊncia estomática máxima s/m 
   real(r_4),parameter,public :: rcmin = 100                     ! ResistÊncia estomática mínima s/m
   real(r_4),parameter,public :: cmin = 1e-4                     ! Minimum to survive kg m-2
   real(r_4),parameter,public :: ca = 363.0                      ! ppmv - atm[CO2]
   real(r_4),parameter,public :: wmax = 500.0                    ! Maximum water soil capacity (Kg m-2)
   real(r_4),parameter,public :: tsnow = -1.0          
   real(r_4),parameter,public :: tice  = -2.5      
   real(r_8),parameter,public :: csru = 0.5D0                    ! Root attribute
   real(r_8),parameter,public :: alfm = 1.391D0                  ! Root attribute
   real(r_8),parameter,public :: gm = 3.26D0 * 86400D0           ! (*86400 transform s/mm to dia/mm)
   real(r_8),parameter,public :: sapwood = 0.05D0                ! Fraction of wood tissues that are sapwood
   real(r_8),parameter,public :: rfrac_leaf = 0.55D0             ! Fraction of NUtriens resorbed at tissue death
   real(r_8),parameter,public :: rfrac_froot = 0.50D0
   real(r_8),parameter,public :: rfrac_wood = 0.45D0
   real(r_4),parameter,public :: nodata = -9999.0
   real(r_4),parameter,public :: vulnerability_curve = 4         ! Vulnerability curve of xylem
   real(r_4),parameter,public :: g = 9.8                         ! Gravity (m/s2)
   real(r_4),parameter,public :: rho = 997.0                     ! Density of water (Kg/m3)
   integer(i_4),parameter,public :: npls = 500                   ! Number of Plant Life Strategies-PLSs simulated (Defined at compile time)
   integer(i_4),parameter,public :: nt1 = 42369                  ! Number of days between 01/1901-12/1970(inclusive)
   integer(i_4),parameter,public :: ntraits = 15                 ! Number of traits for each PLS
   logical(l_1),parameter,public :: debug = .false.              ! Logical variables - For model debug and development
   logical(l_1),parameter,public :: text_ts = .true.


end module global_par


module photo_par
   use types, only : r_8
   implicit none


   real(r_8),public, parameter ::       &
        a   = 0.8300D0       ,&          !Photosynthesis co-limitation coefficient
        a2  = 0.930D0        ,&          !Photosynthesis co-limitation coefficient
        p3  = 21200.0D0      ,&          !Atmospheric oxygen concentration (Pa)
        p4  = 0.080D0        ,&          !Quantum efficiency (mol electrons/Ein)
        p5  = 0.150D0        ,&          !Light scattering rate
        p6  = 2.0D0          ,&          !Parameter for jl
        p7  = 0.50D0         ,&          !Ratio of light limited photosynthesis to Rubisco carboxylation
        p8  = 5200.0D0       ,&          !Photo-respiration compensation point
        p9  = 0.570D0        ,&          !Photosynthesis co-limitation coefficient
        p10 = 0.100D0        ,&          !Q10 function
        p11 = 25.0D0         ,&          !Q10 function reference temperature (oC)
        p12 = 30.0D0         ,&          !Michaelis-Menten constant for CO2 (Pa)
        p13 = 2.100D0        ,&          !Michaelis-Menten constant for CO2
        p14 = 30000.0D0      ,&          !Michaelis-Menten constant for O2 (Pa)
        p15 = 1.20D0         ,&          !Michaelis-Menten constant for O2
        p19 = 0.90D0         ,&          !Maximum ratio of internal to external CO2
        p20 = 0.10D0         ,&          !Critical humidity deficit (kg/kg)
        p22 = 2.0D0          ,&          !Rubisco carboxylation rate
        p23 = 0.30D0         ,&          !Rubisco carboxylation rate
        p24 = 36.0D0         ,&          !Rubisco carboxylation rate (oC)
        p25 = 0.000008D0     ,&          !Maximum gross photosynthesis rate (molCO2/m2/s)
        p26 = 0.50D0         ,&          !light extinction coefficient for IPAR/sun (0.5/sen90)
        p27 = 1.50D0         ,&          !light extinction coefficient for IPAR/shade (0.5/sen20)
        p28 = 0.500D0        ,&          !Soil moisture at field capacity
        p29 = 0.205D0        ,&          !Soil moisture at wilting point
        p30 = 0.015D0        ,&          !Ratio of respiration to Rubisco carboxylation rates
        p31 = 3.850D0        ,&          !Whole plant to leaf respiration ratio
        p32 = 2.00D0         ,&
        p33 = 0.10D0         ,&
        p34 = 0.30D0         ,&
        p35 = 0.05D0         ,&          !soil pool carbon turnover rate
        p36 = 0.25D0         ,&
        alphap = 0.0913D0    ,&          ! 0.0913 parameter for v4m. Hard to explain. See Chen et al. 1994
        vpm25 =  85.0D0      ,&          ! µmol m-2 s-1 PEPcarboxylase CO2 saturated rate of carboxilation at 25°C
        h_vpm = 185075.0D0   ,&          ! Arrhenius eq. constant
        s_vpm = 591.0D0      ,&          ! Arrhenius eq. constant
        r_vpm = 8.314D0      ,&          ! Arrhenius eq. constant
        e_vpm = 60592.0D0    ,&          ! Arrhenius eq. constant
        kp25 = 82.0D0                    ! µmol mol-1 (ppm)  MM constant PEPcase at
end module photo_par
