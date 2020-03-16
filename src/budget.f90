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

! contacts :: David Montenegro Lapola <lapoladm ( at ) gmail.com> 
!             João Paulo Darela Filho <darelafilho ( at ) gmail.com>


module budget
  implicit none
  private

  public :: daily_budget

contains


!,mineral_n,labile_p

  subroutine daily_budget(dt, w1, g1, s1, ts, temp, prec, p0, ipar, rh&
       &, mineral_n, labile_p, sto_budg, cl1_pft, ca1_pft, cf1_pft, dleaf, dwood&
       &, droot, w2, g2, s2, smavg, ruavg, evavg, epavg&
       &, phavg, aravg, nppavg, laiavg, rcavg, f5avg&
       &, rmavg, rgavg, cleafavg_pft, cawoodavg_pft&
       &, cfrootavg_pft, ocpavg, wueavg&
       &, cueavg, c_defavg, vcmax, specific_la&
       &, nupt, pupt, litter_l, cwd, litter_fr, lnr)


    use types
    use global_par
    use photo, only: pft_area_frac, allocation
    use water, only: evpot2, penman, available_energy, runoff
    use productivity

!     ----------------------------INPUTS-------------------------------
    real(r_4),dimension(ntraits,npls),intent(in) :: dt
    real(r_4),dimension(npls),intent(in) :: w1   !Initial (previous month last day) soil moisture storage (mm)
    real(r_4),dimension(npls),intent(in) :: g1   !Initial soil ice storage (mm)
    real(r_4),dimension(npls),intent(in) :: s1   !Initial overland snow storage (mm)
    real(r_4),intent(in) :: ts                   ! Soil temperature (oC)
    real(r_4),intent(in) :: temp                 ! Surface air temperature (oC)
    real(r_4),intent(in) :: prec                 ! Precipitation (mm/day)
    real(r_4),intent(in) :: p0                   ! Surface pressure (mb)
    real(r_4),intent(in) :: ipar                 ! Incident photosynthetic active radiation mol Photons m-2 s-1
    real(r_4),intent(in) :: rh                   ! Relative humidity
    real(r_8),intent(in) :: mineral_n
    real(r_8),intent(in) :: labile_p
! State variables INPUTS & OUTPUTS 
!NEW inout only for test pourposes. these are inputs
! real(r_8),intent(in) :: mineral_n   ! Mineral pools (SOIL) kg(Nutrient) m-2
! real(r_8),intent(in) :: labile_p

    real(r_8),dimension(3,npls),intent(inout)  :: sto_budg ! Rapid Storage Pool (C,N,P)
    real(r_8),dimension(npls),intent(inout) :: cl1_pft  ! initial BIOMASS cleaf compartment
    real(r_8),dimension(npls),intent(inout) :: cf1_pft  !                 froot
    real(r_8),dimension(npls),intent(inout) :: ca1_pft  !                 cawood
    real(r_8),dimension(npls),intent(inout) :: dleaf  ! CHANGE IN cVEG (DAILY BASIS) TO GROWTH RESP
    real(r_8),dimension(npls),intent(inout) :: droot
    real(r_8),dimension(npls),intent(inout) :: dwood 


!     ----------------------------OUTPUTS------------------------------
    real(r_4),intent(out) :: epavg                          !Maximum evapotranspiration (mm/day)
    real(r_4),intent(out),dimension(npls) :: w2             !Final (last day) soil moisture storage (mm)
    real(r_4),intent(out),dimension(npls) :: g2             !Final soil ice storage (mm)
    real(r_4),intent(out),dimension(npls) :: s2             !Final overland snow storage (mm)
    real(r_8),intent(out),dimension(npls) :: smavg          !Snowmelt Daily average (mm/day)
    real(r_8),intent(out),dimension(npls) :: ruavg          !Runoff Daily average (mm/day)
    real(r_8),intent(out),dimension(npls) :: evavg          !Actual evapotranspiration Daily average (mm/day)
    real(r_8),intent(out),dimension(npls) :: phavg          !Daily photosynthesis (Kg m-2 day-1)
    real(r_8),intent(out),dimension(npls) :: aravg          !Daily autotrophic respiration (Kg m-2 day-1)
    real(r_8),intent(out),dimension(npls) :: nppavg         !Daily NPP (average between PFTs)(Kg m-2 day-1)
    real(r_8),intent(out),dimension(npls) :: laiavg         !Daily leaf area Index m2m-2
    real(r_8),intent(out),dimension(npls) :: rcavg          !Daily canopy resistence s/m
    real(r_8),intent(out),dimension(npls) :: f5avg          !Daily canopy resistence s/m
    real(r_8),intent(out),dimension(npls) :: rmavg,rgavg    !maintenance/growth respiration (Kg m-2 day-1)
    real(r_8),intent(out),dimension(npls) :: cleafavg_pft   !Carbon in plant tissues (kg m-2)
    real(r_8),intent(out),dimension(npls) :: cawoodavg_pft  ! 
    real(r_8),intent(out),dimension(npls) :: cfrootavg_pft  !
    real(r_8),intent(out),dimension(npls) :: ocpavg         ! [0-1]
    real(r_8),intent(out),dimension(npls) :: wueavg         ! 
    real(r_8),intent(out),dimension(npls) :: cueavg         ! [0-1]
    real(r_8),intent(out),dimension(npls) :: c_defavg       ! kg(C) m-2
    real(r_8),intent(out),dimension(npls) :: vcmax          ! µmol m-2 s-1
    real(r_8),intent(out),dimension(npls) :: specific_la    ! m2 g(C)-1
    real(r_8),intent(out),dimension(npls) :: nupt           ! g m-2
    real(r_8),intent(out),dimension(npls) :: pupt           ! g m-2
    real(r_8),intent(out),dimension(npls) :: litter_l       ! g m-2
    real(r_8),intent(out),dimension(npls) :: cwd            ! g m-2
    real(r_8),intent(out),dimension(npls) :: litter_fr      ! g m-2
! Lnr variables         [(lln2c),(rln2c),(cwdn2c),(llp2c),(rlp2c),(cwdp2c)]
    real(r_8),intent(out),dimension(6,npls) :: lnr         ! g(N) g(C)-1

!     -----------------------Internal Variables------------------------
    integer(i_4) :: p
    logical :: end_pls = .false.
    logical :: no_cell = .false.
    real(r_4),dimension(ntraits) :: dt1 ! Store pls attributes array (1D)
!     RELATED WITH GRIDCELL OCUPATION

    real(r_8),dimension(npls) :: ocp_mm
    real(r_8),dimension(npls) :: ocp_coeffs !,ocp_coeffs2
    logical(l_1),dimension(npls) :: ocp_wood !, ocp_wood2 - 

!     WBM COMMUNICATION (water balance)
    real(r_4) :: psnow                !Snowfall (mm/day)
    real(r_4) :: prain                !Rainfall (mm/day)
    real(r_4) :: emax

    real(r_4),dimension(npls) :: rimelt               !Runoff due to soil ice melting
    real(r_4),dimension(npls) :: smelt                !Snowmelt (mm/day)
    real(r_4),dimension(npls) :: w                    !Daily soil moisture storage (mm)
    real(r_4),dimension(npls) :: g                    !Daily soil ice storage (mm)
    real(r_4),dimension(npls) :: s                    !Daily overland snow storage (mm)
    real(r_4),dimension(npls) :: ds  
    real(r_4),dimension(npls) :: dw  
    real(r_4),dimension(npls) :: roff                 !Total runoff
    real(r_4),dimension(npls) :: evap                !Actual evapotranspiration (mm/day)


!c     Carbon Cycle
    real(r_4),dimension(npls) ::  ph             !Canopy gross photosynthesis (kgC/m2/yr)
    real(r_4),dimension(npls) ::  ar             !Autotrophic respiration (kgC/m2/yr)
    real(r_4),dimension(npls) ::  nppa           !Net primary productivity / auxiliar
    real(r_4),dimension(npls) ::  laia           !Leaf area index (m2 leaf/m2 area)
    real(r_4),dimension(npls) ::  rc2            !Canopy resistence (s/m)
    real(r_4),dimension(npls) ::  f1             !
    real(r_4),dimension(npls) ::  f5             !Photosynthesis (mol/m2/s)
    real(r_4),dimension(npls) ::  vpd            !Vapor Pressure deficit
    real(r_4),dimension(npls) ::  rm             ! maintenance & growth a.resp
    real(r_4),dimension(npls) ::  rg
    real(r_4),dimension(npls) ::  wue, cue, c_def
    real(r_8),dimension(npls) ::  cl1,cf1,ca1 ! carbon pre-allocation 
    real(r_8),dimension(npls) ::  cl2,cf2,ca2 ! carbon pos-allocation
    real(r_8),dimension(3,npls) :: day_storage   ! g m-2




!     Precipitation
!     =============     
    psnow = 0.0
    prain = 0.0
    if (temp.lt.tsnow) then
       psnow = prec !* 0.0328549 ! converts mm/month in mm/day
    else
       prain = prec !* 0.0328549
    endif

!     Initialization
!     --------------
    epavg   = 0.0
    w       = w1     ! hidrological pools state vars  
    g       = g1
    s       = s1
    smavg   = 0.0    !  plss vectors (outputs)
    ruavg   = 0.0
    evavg   = 0.0
    rcavg   = 0.0
    laiavg  = 0.0
    phavg   = 0.0
    aravg   = 0.0
    nppavg  = 0.0
    rmavg   = 0.0
    rgavg   = 0.0
    ocpavg  = 0.0
    wueavg  = 0.0
    cueavg  = 0.0 
    ocp_mm  = 0.0
    emax  = 0.0

    cl1 = cl1_pft ! daily initial carbonVEG pools 
    ca1 = ca1_pft
    cf1 = cf1_pft

    nppa  = 0.0
    ph    = 0.0
    ar    = 0.0
    laia  = 0.0
    f5    = 0.0
    f1    = 0.0
    vpd   = 0.0
    rc2   = 0.0
    rm    = 0.0
    rg    = 0.0 
    wue   = 0.0
    cue   = 0.0
    rc2   = 0.0
    f1    = 0.0
    f5    = 0.0

!     Grid cell area fraction (%) ocp_coeffs(pft(1), pft(2), ...,pft(p))
!     =================================================================     
    call pft_area_frac(cl1, cf1, ca1, ocp_coeffs, ocp_wood) ! def in funcs.f90


    if(debug) then
       write(1234,*) '-----Message from budget-----------------------'
       write(1234,*) '-----------------------------------------------'
       write(1234,*) 'cl1        cf1         ca1'
       write(1234,*) cl1, cf1, ca1 
       write(1234,*) '-------ocp_coeffs------------------------------'
       write(1234,*) ocp_coeffs
       write(1234,*) ocp_wood
    endif

!     Maximum evapotranspiration   (emax)
!     =================================
    emax = evpot2(p0,temp,rh,available_energy(temp))

!     Productivity & Growth (ph, ALLOCATION, aresp, vpd, rc2 & etc.) for each PLS
!     =================================
    do p = 1,npls

       if(debug) then
          write(1234,*) 'PFT n ---->', p
       endif

       dt1 = dt(:,p) ! Pick up the pls functional attributes list

       end_pls = .false. 

! subroutine prod(dt,light_limit,temp,p0,w,ipar,rh,emax,cl1_prod,&
!      & ca1_prod,cf1_prod,beta_leaf,beta_awood,beta_froot,sto1,ph,ar,&
!      & nppa,laia,f5,vpd,rm,rg,rc,wue,c_defcit,vm_out,sla,sto2)

       call prod(dt1,ocp_wood(p),temp,p0,w(p),ipar,rh,emax,cl1(p)&
            &,ca1(p),cf1(p),dleaf(p),dwood(p),droot(p)&
            &,sto_budg(:,p),ph(p),ar(p),nppa(p),laia(p)&
            &,f5(p),vpd(p),rm(p),rg(p),rc2(p),wue(p),c_def(p)&
            &,vcmax(p),specific_la(p),day_storage(:,p))

       sto_budg(:,p) = day_storage(:,p)

       if(debug) then
          write(1234,*) 'sto_budg ---->',sto_budg(:,p)
          write(1234,*) 'dt1 --------->',dt1
       endif
!!

       ! TODO 

       ! the root potential to extract nutrients from soil is 
       ! calculated here

       ! f(fine_root_mass, fine_root_residence_time)

       ! Specific root area (fine roots)
       ! Specific root length

       ! upt_capacity <- root specific uptake capacity g (P or N) m-2(root) day-1

       ! A function like F5


!     Carbon/Nitrogen/Phosphorus allocation/deallocation
!     =====================================================
       if(debug) then
          write(1234,*) 'INPUTS PRECEDING ALLOCATION'
          write(1234,*) 'dt1', dt1
          write(1234,*) 'nppa(p)', nppa(p)
          write(1234,*) 'mineral_n', mineral_n
          write(1234,*) 'labile_p', labile_p 
          write(1234,*) 'cl1(p)', cl1(p)
          write(1234,*) 'cf1(p)', cf1(p)
          write(1234,*) 'ca1(p)', ca1(p)
          write(1234,*) 'stop(p)', sto_budg(:,p)
       endif

       call allocation (dt1,nppa(p),mineral_n,labile_p,cl1(p),ca1(p)&
            &,cf1(p),sto_budg(:,p),day_storage(:,p),cl2(p),ca2(p)&
            &,cf2(p),litter_l(p),cwd(p)&
            &,litter_fr(p),nupt(p),pupt(p),lnr(:,p),end_pls)

       sto_budg(:,p) = day_storage(:,p)

! Se o PFT nao tem carbono goto 666-> TUDO ZERO
       if(end_pls) then
          no_cell = .true.
          dleaf(p) = 0.0
          dwood(p) = 0.0
          droot(p) = 0.0
          goto 666 ! gt hell
       endif

       if(ph(p) .eq. 0.0 .or. nppa(p) .eq. 0.0) then
          cue(p) = 0.0
       else
          cue(p) = nppa(p)/ph(p)
       endif

       dleaf(p) = cl2(p) - cl1(p)  !kg m-2
       dwood(p) = ca2(p) - ca1(p)
       droot(p) = cf2(p) - cf1(p)

       !!!!====================================================
!     Snow budget
!     ===========     
       smelt(p) = 2.63 + 2.55*temp + 0.0912*temp*prain !Snowmelt (mm/day)
       smelt(p) = amax1(smelt(p),0.)
       smelt(p) = amin1(smelt(p),s(p)+psnow)
       ds(p) = psnow - smelt(p)
       s(p) = s(p) + ds(p)

!     Water budget
!     ============
       if (ts.le.tice) then !Frozen soil
          g(p) = g(p) + w(p) !Soil moisture freezes
          w(p) = 0.0
          roff(p) = smelt(p) + prain !mm/day
          evap(p) = 0.0
          ph(p) = 0.0
          ar(p) = 0.0
          nppa(p) = 0.0
          laia(p) = 0.0
          rc2(p) = rcmin
          rm(p) = 0.0
          rg(p) = 0.0
! colocar as outras vars da prod aqui??

       else                !Non-frozen soil
          w(p) = w(p) + g(p)
          g(p) = 0.0
          rimelt(p) = 0.0
          if (w(p).gt.wmax) then
             rimelt(p) = w(p) - wmax !Runoff due to soil ice melting
             w(p) = wmax
          endif

          roff(p) = runoff(w(p)/wmax)       !Soil moisture runoff (roff, mm/day)

          evap(p) = penman(p0,temp,rh,available_energy(temp),rc2(p)) !Actual evapotranspiration (evap, mm/day)
          dw(p) = prain + smelt(p) - evap(p) - roff(p)
          w(p) = w(p) + dw(p)
          if (w(p).gt.wmax) then
             roff(p) = roff(p) + (w(p) - wmax)
             w(p) = wmax
          endif
          if (w(p).lt.0.) w(p) = 0.
          roff(p) = roff(p) + rimelt(p) !Total runoff
       endif

       if (p .eq. 1) epavg = emax !mm/day
       w2(p) = real(w(p),r_4) !* ocp_coeffs(p), r_4)
       g2(p) = real(g(p),r_4) !* ocp_coeffs(p), r_4)
       s2(p) = real(s(p),r_4) !* ocp_coeffs(p), r_4)
       smavg(p) = smelt(p) !* ocp_coeffs(p)
       ruavg(p) = roff(p)  !!* ocp_coeffs(p)   ! mm day-1
       evavg(p) = evap(p)  !!* ocp_coeffs(p)   ! mm day-1
       phavg(p) = ph(p)    !!* ocp_coeffs(p)   !kgC/m2/day
       aravg(p) = ar(p)    !!* ocp_coeffs(p)   !kgC/m2/year
       nppavg(p) = nppa(p) !!* ocp_coeffs(p)   !kgC/m2/day
       laiavg(p) = laia(p) !* ocp_coeffs(p)
       rcavg(p) = rc2(p)   !!* ocp_coeffs(p)   ! s m -1
! F5 goes out!
       f5avg(p) = f5(p) !* ocp_coeffs(p)
       rmavg(p) = rm(p) !* ocp_coeffs(p)
       rgavg(p) = rg(p) !* ocp_coeffs(p)

! Vcmax and SLA ! These are CWM: Dependent on (respectively) Foliar N and P; tau_leaf (resid. time)
!  vcmax(p) = vcmax(p) !* ocp_coeffs(p) 
!  specific_la(p) = specific_la(p) !* ocp_coeffs(p)

!  !Rapid Storage Pool
!  sto_budg(:,p) = sto_budg(:,p) !* ocp_coeffs(p)

!  ! Variables used to carbon in soil
!  nupt(p) = nupt(p) !* ocp_coeffs(p)
!  pupt(p) = pupt(p) !* ocp_coeffs(p)
!  litter_l(p) = litter_l(p) !* ocp_coeffs(p)
!  cwd(p) = cwd(p) !* ocp_coeffs(p)
!  litter_fr(p) = litter_fr(p) !* ocp_coeffs(p)
!  lnr(:,p) = lnr(:,p) !* ocp_coeffs(p)

       cleafavg_pft(p)  =  cl2(p) !- ((c_def(p) * 0.00273791) / 3.0)
       cawoodavg_pft(p) =  ca2(p) !- ((c_def(p) * 0.00273791) / 3.0)
       cfrootavg_pft(p) =  cf2(p) !- ((c_def(p) * 0.00273791) / 3.0)
       wueavg(p) = wue(p)  !!* ocp_coeffs(p)
       cueavg(p) = cue(p)  !!* ocp_coeffs(p)
       c_defavg(p) = c_def(p) !* ocp_coeffs(p)

!  gridcell Occupation
       ocpavg(p) = ocp_coeffs(p)


       if(c_def(p) .gt. 0.0) then
          if(dt1(7) .gt. 0.0) then
             cl1_pft(p) = cl2(p) - ((c_def(p) * 1e-3) * 0.333333333)
             ca1_pft(p) = ca2(p) - ((c_def(p) * 1e-3) * 0.333333333)
             cf1_pft(p) = cf2(p) - ((c_def(p) * 1e-3) * 0.333333333)
          else
             cl1_pft(p) = cl2(p) - ((c_def(p) * 1e-3) * 0.5)
             ca1_pft(p) = 0.0
             cf1_pft(p) = cf2(p) - ((c_def(p) * 1e-3) * 0.5)
          endif
       else
          if(dt1(7) .gt. 0.0) then
             cl1_pft(p) = cl2(p) 
             ca1_pft(p) = ca2(p) 
             cf1_pft(p) = cf2(p) 
          else
             cl1_pft(p) = cl2(p)
             ca1_pft(p) = 0.0
             cf1_pft(p) = cf2(p)
          endif
       endif

       no_cell = .false.

!  if(cf1_pft(p).le.0.0 .and. cl1_pft(p).le.0.0) then
!     no_cell = .true.
!  else if(ca1_pft(p).lt. 0.0 .and. dt1(7) .gt. 0.0) then
!     no_cell = .true.
!  else
!     continue
!  endif
! :,p
! 
666    continue
       if(no_cell) then
! All*** outputs are set to zero except epavg
          if(p .eq. 1) epavg = emax
          smavg(p) = 0.0
          ruavg(p) = 0.0
          evavg(p) = 0.0
          rcavg(p) = 0.0
          phavg(p) = 0.0
          aravg(p) = 0.0
          nppavg(p) = 0.0
          wueavg(p) = 0.0
          cueavg(p) = 0.0
          c_defavg(p) = 0.0
          laiavg(p) = 0.0
          f5avg(p) = 0.0
          rmavg(p) = 0.0
          rgavg(p) = 0.0
          vcmax(p) = 0.0
          specific_la(p) = 0.0
          sto_budg(:,p) =  0.0
          nupt(p) = 0.0
          pupt(p) = 0.0 
          litter_l(p) = 0.0
          cwd(p) = 0.0
          litter_fr(p) = 0.0
          lnr(:,p) = 0.0
          ocpavg(p) = 0.0
          cleafavg_pft(p)  = 0.0
          cawoodavg_pft(p) = 0.0
          cfrootavg_pft(p) = 0.0
          cl1_pft(p) = 0.0
          ca1_pft(p) = 0.0
          cf1_pft(p) = 0.0
          w2(p) = 0.0
          g2(p) = 0.0
          s2(p) = 0.0
       endif

    enddo ! end pls_loop (p)

! call pft_area_frac(cl1_pft, cf1_pft, ca1_pft, ocp_coeffs2, ocp_wood2)
! ocpavg = (ocp_coeffs + ocp_coeffs2) / 2.0
    if(debug) then      
       write(1234,*) '-----END Message from budget-------------------'
       write(1234,*) '------------------------------------------------'
    endif

! CLEAN NANs OF some outputs
    do p = 1,npls
       if(isnan(cleafavg_pft(p))) cleafavg_pft(p) = 0.0
       if(isnan(cawoodavg_pft(p))) cawoodavg_pft(p) = 0.0
       if(isnan(cfrootavg_pft(p))) cfrootavg_pft(p) = 0.0
    end do

    return
  end subroutine daily_budget

end module budget
