!**********************************************
!*       timestep.f95       
!**********************************************
!
!> FluidEarth entry point, also called from RUN_SCAPE
!>
!> Advances the model by one timestep.
!>
!> If one year has passed, increment the year and take year-change actions 
!> including writing output to files if required.
!>
!> Update the sea level rise and smaple the wave conditions.
!>
!> Transform the wave conditions and calculate the cross-shore distribution of sediment.
!>
!> Calculate sediment transport.
!>
!> Update the cliff/platforms in the boundary regions and add in the contribution of the talus to the beach.
!>
!> Adjust the location of the beach and adjust the position of the representation of talus surface.
!>
!> Slump the cliff.
!
!**********************************************
!
!   Copyright (C) 2016  Michael Walkden & James Hall
!
!   This program is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation (Version 3).
!
!   This program is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!       Mike.walkden@gmail.com
!       WSP | Parsons Brinckerhoff    
!       Keble House
!       Southernhay Gardens
!       Exeter, EX1 1NT                
!
!   The code was developed with the financial support of:
!        The Engineering and Physical Science Research Council
!        The Tyndall Centre for Climate Change Research
!        The Natural Environment Research Council
!
!   The authors would also like to gratefully acknowledge the contributions made to the code by:
!   Mark Dickson - development of process representations
!   James Thomas - code optimisation
!   John Barnes - code rationalisation, optimisation, development of process representations, 
!       development of the horizontal grid architecture, OpenMI/FluidEarth engine wrapping code 
!       and support and other contributions.
!
!**********************************************

module timestep_module
    implicit none

contains


subroutine update(success_code)

    use setup_module
    use global_data
    use local_data
    use beach_location_module
    use xshore_dist_t_data
    use dump_module
    use wave_transformation_module
    use xshore_dist_t_module
    use slump_module
    use platform_erosion_module
    use sediment_transport_module
    use SLR_module
    use utilities_module
    implicit none

    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

    real(kind=double) talus_width, max_tvolume, min_tvolume, apparent_tvolume
    real(kind=double) :: tidal_amp_value !< single value as read from Waves and Tides file
        !< to be scaled and interpolated along the model
    integer :: seawall_element ! index of cliffy element that defines the position of a seawall
    real(kind=double) :: new_beachos
    character(len=256) :: cwd
#ifdef __GFORTRAN__
#else
    character(len=256) :: CURDIR@
#endif
    logical :: mask_h(nYsections) ! mask sections with inshore wave height >= minHs

    success_code = 0
	    !  write(6, fmt='(A,I7)') ' nTide ', ntide ! To screen Debug code
        if (ntide > 1 .and. mod(ntide - 1, tidesPerYear) == 0) then ! One year has passed, increment the year
            year = year + 1
                
            fines_volume = 0
        
            total_beach_addition = 0.0D0

            ann_av_bvol = 0.0D0

            call GET_PERIOD_FACTOR(year, periodFactor, havePeriodFactor)
            call GET_HEIGHT_FACTOR(year, heightFactor, haveHeightFactor)
            call GET_ANGLE_CHANGE(year, angleChange, haveAngleChange)
    
            if (includeSeawall) then
                ! Seawalls and revetments on the Y sections
                if (seawallMode == defineSeawallsByYear) then
                    seawall_element = seawallBaseLevel / dz
                end if

                do section = 1, nYsections 
                    if (year == seawall_construction(section)) then ! Make a seawall active at this location
                        seawall_active(section) = 1
                        
                        if (seawallMode == defineSeawallsByYear) then
                            seawall_os(section) = cliffy(seawall_element, section)
                        end if
                    end if
                       
                    if (year == seawall_removal(section)) then
                        seawall_active(section) = 0    ! Make the seawall inactive at this section
                    end if
                end do ! section = 1, nYsections
            end if

            QAnn = 0.0D0
            QpAnn = 0.0D0
            groyneEffectAnnualSum = 0.0D0

            sea_rise = sea_rise + sea_level_next_year - msl_OD
            msl_OD = sea_level_next_year
            call SL_for_year(year + 1, sea_level_next_year, success_code)
            if (success_code < 0) return

            slr_per_tide = (sea_level_next_year - msl_OD) / tidesPerYear
        else   ! End of the year increment clause
            sea_rise = sea_rise + slr_per_tide
            msl_OD = msl_OD + slr_per_tide
        end if 

        call dump_real("run_scape11698", "bvolume", .false., bvolume)

!***********************************************
!
!    New Sea level rise algorithm WkBk p90(8)
!
!***********************************************

!*** Allocate temp arrays

          allocate(tempslrc(nce,nYsections), stat = alloc_error)    
          if(alloc_error /= 0)then
            write(6,*)'*** Error : could not allocate space'
            write(6,*)'*** Error : for cliffy, beachy & cliffheights'
            stop
          end if    

      tempslrc=0

      ! for sea level rise
      
    do while(sea_rise >= dz)
        sea_rise = sea_rise - dz
        baseline_OD = baseline_OD + dz
        tempslrc(1:nce-1,:) = cliffy(2:nce,:)
        tempslrc(nce,:) = cliffy(nce,:)
        cliffy = tempslrc
    end do
      
      ! for sea level fall
      
    do while(sea_rise <= -dz)
        sea_rise = sea_rise + dz
        baseline_OD = baseline_OD - dz
        tempslrc(2:nce,:) = cliffy(1:nce-1,:)
        tempslrc(1,:) = cliffy(2,:) + (cliffy(2,:) - cliffy(3,:))
        cliffy = tempslrc
    end do
      
    call dump_cliff_real("run_scape1738", "cliffy", cliffy)

!*** Deallocate temp arrays

          deallocate(tempslrc, stat=alloc_error)
          if(alloc_error /= 0)then
            write(6,*)'*** Unexpected deallocation error'
            write(6,*)'*** for tempslrc'
          end if   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Sample the wave conditions
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ierr=0
        period = 0.0D0

        do i = 1, wave_count
            read(wavesAndTidesFileUnit, *, iostat=ierr) tidal_amp_value, period, (hwp(wp), anglewp(wp), wp = 1, nWavePoints)
            
            if(ierr.ne.0)then    ! Prog has come to the end of the wave file
                rewind(wavesAndTidesFileUnit, iostat=ierr)        ! Rewind it
                read(wavesAndTidesFileUnit, *, iostat=ierr) tidal_amp_value, period, (hwp(wp), anglewp(wp), wp = 1, nWavePoints)
                success_code = SUCCESS_CODE_FROM_IOSTAT(wavesAndTidesFileUnit, ierr, 'Reading waves and tides')
                if (success_code < 0) return
            end if
        end do

        ! Adjust wave values for year factors
        if (havePeriodFactor) period = period * periodFactor
        if (haveHeightFactor) hwp = hwp * heightFactor
        if (haveAngleChange) then
            anglewp = anglewp + angleChange
            if (angleChange < 0.0D0) then
                where (anglewp < 0.0D0) anglewp = anglewp + 360.0D0
            else
                where (anglewp >= 360.0D0) anglewp = anglewp - 360.0D0
            end if
        end if

        wave_count = 1

        if (tidalAmpSection > 1) then
            call INTERPOLATE_VALUES(tidal_amp(1: tidalAmpSection), &
                tidal_amp_value * tidalAmpScaleRight, tidal_amp_value)
        end if

        if (tidalAmpSection < nYsections) then
            call INTERPOLATE_VALUES(tidal_amp(tidalAmpSection: nYsections), &
                tidal_amp_value, tidal_amp_value * tidalAmpScaleLeft)
        end if

        if (nYsections == 1) then
          tidal_amp(1) = tidal_amp_value
        end if

        anglewp = anglewp * radcon    ! Convert to radians
        
        ! Check to see whether to proceed into the main loop

        proceed_switch = 1

        if (maxval(hwp) < minHs)then
          proceed_switch = 0    
        end if
        if(period < 2)then
          proceed_switch = 0    
        end if
    
    pot_left_sediment_flux = 0.0d0
    pot_right_sediment_flux = 0.0d0
    left_sediment_flux = 0.0d0
    right_sediment_flux = 0.0d0

    !*******************************************
    !
    !    Check whether this is an event, if so then proceed with 
    !    the calculations of sediment transport and erosion
    !
    !*******************************************

    ann_av_bvol = ann_av_bvol + bvolume
    call dump_real("run_scape1851", "ann_av_bvol", .false., ann_av_bvol)

    if(proceed_switch == 1)then
            
        eventcounter = eventcounter + 1 ! An event is a tide during which 
        ! sediment transport and erosion occurs. 

        heightsurge = 0.0D0 ! This has never been used, or properly tested
        ! so far it has been assumed that the surge is accounted for in 
        ! the water level file

    call WAVE_TRANSFORMATION(nQpoints, nYsections, nce, &
            hwp, period, anglewp, depthWPMsl, tidal_amp, &  ! MD 300304
            heightsurge, angleosc, periodValues, depthValues, clut,    &
            cglut,breakrat, beachy,cliffy,        &
            msl_e,setup,setdown,h,anglewob,cgroup, theta, &
            setup_qpoint,setdown_qpoint,h_qpoint,anglewob_qpoint,cgroup_qpoint, &
            theta_qpoint, wpr)

    call dump_real("run_scape1877", "setup", .false., setup)
    call dump_real("run_scape1880", "setdown", .false., setdown)
    call dump_real("run_scape1881", "h", .false., h)
    call dump_real("run_scape1882", "anglewob", .false., anglewob)
    call dump_real("run_scape1883", "cgroup", .false., cgroup)
    call dump_real("run_scape1884", "setup", .true., setup_qpoint)
    call dump_real("run_scape1885", "setdown", .true., setdown_qpoint)
    call dump_real("run_scape1886", "h", .true., h_qpoint)
    call dump_real("run_scape1887", "anglewob", .true., anglewob_qpoint)
    call dump_real("run_scape1888", "cgroup", .true., cgroup_qpoint)
    call dump_real("run_scape1889", "theta", .true., theta_qpoint)
        
        ! mask out very small inshore waves
        mask_h = h > minHs
        
        do section = 1, nYsections
            dbreak(section) = h(section) / breakrat    ! Wave breaker depth for each section
            half_tide(section) = nint(tidal_amp(section) * 20.0d0 / dbreak(section))
            ! Half tide, this is a complicated one.
            ! Ultimately it is for the calculation of the tide-integrated cross-shore distributions 
            ! of drift and erosion. These use the static distribution of drift or erosion. These 
            ! are shape functions (sf) sf_drift & sf_erode. These have abscissa increments (dx) of
            ! 1/40th of the breaker depth. Half tide is the tidal amplitude measured 
            ! in these increments
            if (half_tide(section) < 1.0D0) then
                half_tide(section) = 1.0D0
            end if
        end do
        
        if(nYsections > 1)then    ! More than one section
                ! call the subroutine to calculate the TIDAL cross shore distribution of 
                ! longshore drift.

                if (tidalRangeVariation) then
                    do section = 1, nYsections
                        if (mask_h(section)) then

                          

                          
                            call xshore_dist_t(nint(half_tide(section)), sf_drift,    &
                              nce, dbreak(section), dsf_drift,        &
                              setup(section),setdown(section), heightsurge, np_sf_drift,        &
                              beachCrestLevelE - beachHeightE, beachCrestLevelE,        &
                              cliffheights, msl_m,                &
                              sect_drift(:, section), bot_drift(section), top_drift(section), 1, success_code)
                            if (success_code < 0) return
                        else ! (mask_h(section))
                            sect_drift(:, section) = 0.0D0
                            bot_drift(section) = beachCrestLevelE - beachHeightE
                            top_drift(section) = beachCrestLevelE - 1
                        end if ! (mask_h(section))
                    end do ! section = 1, nYsections
                else ! (tidalRangeVariation)
                    mean_dbreak = MASKED_AVERAGE_VALUE(dbreak, mask_h) 
                    mean_half_tide = MASKED_AVERAGE_VALUE(half_tide, mask_h)
                    mean_setup = MASKED_AVERAGE_VALUE(setup, mask_h) 
                    mean_setdown = MASKED_AVERAGE_VALUE(setdown, mask_h)
                     
                    if (COUNT(mask_h) > 0) then
                        call xshore_dist_t(nint(mean_half_tide), sf_drift,    &
                        nce, mean_dbreak, dsf_drift,        &
                        mean_setup, mean_setdown, heightsurge, np_sf_drift,        &
                        beachCrestLevelE - beachHeightE, beachCrestLevelE,        &
                        cliffheights, msl_m,                &
                        this_tide_drift, bottom, top, 1, success_code)
                        if (success_code < 0) return
                    else ! (COUNT(mask_h) > 0) 
                          
                             
                        bottom = beachCrestLevelE - beachHeightE
                        top = beachCrestLevelE - 1
                    end if ! (COUNT(mask_h) > 0) 
    
                    bot_drift = bottom    ! Sets the lower vertical limit of the drift distribution
                    top_drift = top ! Sets the upper vertical limit of the drift distribution
                    
                    do section = 1, nYsections
                        ! sect_drift is the cross-shore distribution of longshore drift
                        if (mask_h(section)) then
                            sect_drift(:,section) = this_tide_drift(:)
                        else
                            sect_drift(:, section) = 0.0D0
                        end if
                    end do ! section = 1, nYsections
                end if ! (tidalRangeVariation)

    !*******************************************
    !
    ! ----    Move sediment ----
    !
    !*******************************************

        call SEDIMENT_TRANSPORT(nce, nYsections, nQpoints, year, &
            offshore_st, bvolume, sediment_influx, obvolume, top_drift, bot_drift, h, &
            sect_drift, beachy, cliffy, QAnn, QPAnn, &
            period, tstep_secs, h_qpoint, cgroup_qpoint, theta_qpoint, cerck1)

        end if    !    if(nQpoints > 1)then

!*********************************************************************
!
!    Code developed for Norfolk to represent a storm beach profile
!    See  WkBk p155(10)
!
!*********************************************************************

        do section = firstActiveSection, lastActiveSection
            crossbeach_st = 0
            ! Check if storm waves
            if ((h(section) > minStormWaveHeight)) then ! (m) storm waves
                call CROSSBEACH_EXCHANGE(section, nYsections, nce, ubTop, ubBottom, lbTop, lbBottom, &
                    upBeachHeightE, lowbeachheightE, &
                    stormflag, crossbeach_st, &
                    upbvolume, upbruun, upbeach_os, upbeachy, &
                    lowbvolume, lowbruun, lowbeach_os, lowbeachy, &
                    cliffy, talusy, beachdepth, bvolume, h, offshore_st, &
                    period, tstep_secs)
            else
                stormflag(section) = 0
            endif
        enddo



        call PLATFORM_EROSION(nYsections, nce, &
            upbeachy, &
            lowbeachy, &
            beachelev, cliffy, beachy, talusy, this_erode, m, &
            temp5, bdepth, battenuation, &
            seawall_effect,  cliffheights, &
            beach_addition, fines_addition, h, &
            dbreak, half_tide, sand_fraction, resistance, tvolume, &
            seawall_os, seawall_active, setup, setdown, &
            sf_erode, np_sf_erode, period, tstep_secs, &
            block_size, dsf_erode, heightsurge, &
            msl_m, bottom, top, success_code)
        if (success_code < 0) return

!***************************
!    Update beach parameters
!***************************

      if (mod(eventcounter,100) == 0)then    ! Only do every 100 events
          do section = firstActiveSection, lastActiveSection
              ! Update the beach offset and berm width
              new_beachos = cliffy(beachCrestLevelE, section)    
              berm(section) = berm(section) + beachos(section) - new_beachos
              beachos(section) = new_beachos
          enddo
      end if

!****************************************************************
!
!       Update the cliff/platforms in the boundary regions
!
!****************************************************************

    ! Erosion is not calculated in the boundary regions to save time. The profile shapes
    ! in these regions are assumed to be the first eroding section adjacent to the boundary.
    ! same as others. The shoreline planshape angle is assumed to never change. Not a 
    ! good assumption! Mind you the boundary regions are not always used.

            if(mod(eventcounter,100) == 0)then ! Every 100 events
                do section = 1, firstActiveSection - 1    ! Update 'lower' boundary region
                    cliffy(:,section) = cliffy(:,firstActiveSection)-section_offsets(firstActiveSection)+section_offsets(section)
                end do
                do section = lastActiveSection + 1, nYsections ! Update 'higher' boundary region
                    cliffy(:,section) = cliffy(:,lastActiveSection)-section_offsets(lastActiveSection)+section_offsets(section)
                end do
            end if

!****************************************************************
!
!    Add in the contribution of the talus to the beach
!
!****************************************************************

    ! The material moving from the talus to the beach or from the platform to the beach
    ! is represented by 'beach_addition'
    ! Material lost as fines is represented by 'fines_addition'

        ! Assume some values for the boundary regions, where, of course, no talus or platform 
        ! erosion is happening
            do section = 1, firstActiveSection - 1
                beach_addition(section) = beach_addition(firstActiveSection)
                fines_addition(section) = fines_addition(firstActiveSection)
            end do
            do section = lastActiveSection + 1, nYsections
                beach_addition(section) = beach_addition(lastActiveSection)
                fines_addition(section) = fines_addition(lastActiveSection)
            end do

        ! Code added during creation of new talus representation (WkBk p136(11))
            do section = 1, nYsections
                bvolume(section) = bvolume(section) + beach_addition(section)
                total_beach_addition(section) = total_beach_addition(section) + beach_addition(section)
                fines_volume = fines_volume + fines_addition(section)
            end do
     call dump_real("timestep_607", "beach_addition", .false., beach_addition)
     call dump_real("timestep_607", "bvolume", .false., bvolume)
     call dump_real("timestep_607", "total_beach_addition", .false., total_beach_addition)

     if (year < startYear + beachZeroDelay) then
         bvolume = 0.0D0
     end if
            
!****************************************
!
!    Adjust the location of the beach
!
!****************************************
            !do section = 1, nYsections
                call BEACH_LOCATION(-1, nYsections, nce, beachHeightE, bTop, bBottom, &
                    bruun, beachos, beachy, &
                    cliffy, beachdepth, bvolume, &
                    cliffy, beach_regular, beachelev, berm)
            !end do

        if (vellingaActive) then
!****************************************
!
!    calculate a VELLINGA profile for this wave condition
!
!****************************************


            v_level_odn = 0.0D0

            ! Need to handle tidal_amp properly here
            vell_top_e = msl_e + nint(tidal_amp(1)/dz)
            vell_bott_e = msl_e - nint(hwp(1)*0.75/dz)

            ! define the horizontal positions of the Vellinga profile
            ! Nothing at the top
            ! Vellinga shape within a defined range
            do i = vell_top_e,vell_bott_e,-1
                vell_y = (vell_top_e - i)*dz
                vxs(i) = (((((7.6/hwp(1))*vell_y)/0.47)**2)-18)/(((7.6/hwp(1))**1.28)*0.962)
            end do

            ! Constant slope at the bottom
            do i = vell_bott_e,1,-1
                vxs(i) = vxs(i+1) + 12.5*dz
            end do



            ! This profile is defined by offshore conditions so is true for all
            ! sections. However the section location should vary.

!****************************************
!
!    Adjust the location of the VELLINGA profiles
!
!****************************************
            vellinga = 0
            v_ref = 0

            do i = 1, nvsect
                section = vellinga_sections(i)
                
                if (seawall_active(section) == 1) then
                    call BEACH_LOCATION(section, nYsections, nce, beachHeightE, bTop, bBottom, &
                        bruun, beachos, beachy, &
                        cliffy, beachdepth, bvolume, &
                        cliffy, beach_vellinga, beachelev, vell_berm, &
                        seawall_os, vellinga, vxs, vell_top_e)
            
                    ! find the level of the intersection of the seawall
                    v_ref(i) = vell_top_e
                    do while (vellinga(v_ref(i)) .lt. seawall_os(section))
                        v_ref(i)= v_ref(i)-1
                    end do
                    v_level_odn(i) = v_ref(i)*dz - msl_m + msl_OD
                end if ! (seawall_active(section) == 1)

                write(saveVellLevelsUnit,"(F7.3)")     v_level_odn(i)
            end do
        end if ! vellingaActive
        
    end if    ! End of erosion event check
    ! Code gets here every hour
!end do     ! End of tide hour loop



	! Output potential and actual sediment out flux at boundaries
    if (year > startYear + beachZeroDelay) then
         write(saveSedFluxLeftUnit,"(F14.2)") left_sediment_flux
         write(savePotSedFluxLeftUnit,"(F14.2)") pot_left_sediment_flux
         write(saveSedFluxRightUnit,"(F14.2)") right_sediment_flux
         write(savePotSedFluxRightUnit,"(F14.2)") pot_right_sediment_flux
     end if
    

!*****************************************************************
!
!    Adjust the position of the representation of talus surface
!
!*****************************************************************

if (mod(eventcounter, 1) == 0)then ! Do this every event
    do section = firstActiveSection, lastActiveSection
        apparent_tvolume = 0.0D0
        
        do i = 1, nce
            talus_width = talusy(i, section) - cliffy(i, section)
            
            if (talus_width > 0.0D0) apparent_tvolume = apparent_tvolume + talus_width
        end do
        
        apparent_tvolume = apparent_tvolume * dz * dx
        
        if (talusy(nce, section) > cliffy(nce,section))then
            ! The talus projects beyond the top of the model
!            call RAISE_EXCEPTION('The talus projects beyond the top of the model')
!            success_code = -1
!            return
            apparent_tvolume = apparent_tvolume +        &
                (0.5 * tan(talusSlope_radians) * (talusy(nce,section) - cliffy(nce,section)) *        &
                (talusy(nce,section) - cliffy(nce,section)))
        end if

        ! Only adjust if significantly different
        min_tvolume = tvolume(section) * 0.95D0
        max_tvolume = tvolume(section) * 1.05D0

        if(apparent_tvolume > max_tvolume) then
            talusy(:,section) = talusy(:,section) - 0.25D0
        elseif(apparent_tvolume < min_tvolume) then
            talusy(:,section) = talusy(:,section) + 0.25D0
        end if
        talus_offsets(section) = talusy(nce,section)
    end do
end if

    call dump_cliff_real("timestep749", "talus", talusy)
!**********************
!
!    Slump the cliff 
!
!**********************

    if(mod(eventcounter, slumpPeriod) == 0)then ! Do this every slumpPeriod events
        if(nQpoints == 1)then
            i=1
            call slump(i)
        else
            do i=firstActiveSection,lastActiveSection
                call slump(i)
            end do
        end if
    end if

    if (mod(ntide, tidesPerYear) == 0) then ! One year has passed, output the results for the year
#ifdef __GFORTRAN__
    call getcwd(cwd)
#else
    cwd = CURDIR@()
#endif    
        write(6, fmt='(A,A,I6)') trim(cwd), ' Year:', year ! To screen

        outputYearFlag = year == nextOutputYear

        if (outputYearFlag) then
            call SET_OUTPUT_YEAR(year)
        end if
        
        ! capture the next retreat
        retreat(:, year) = cliffy(nce, :) - initial_cliff_toe(:)

        if (outputYearFlag) then
            write(saveFinesVolUnit,"(F10.0)") fines_volume
        end if
        
        if (outputYearFlag .and. year >= firstYearTransportOutput) then
            ! Total annual beach addition
            write(saveBeachAddAnnUnit,"(F10.0)") total_beach_addition
        end if

        if (outputYearFlag) then
            if (year >= firstYearTransportOutput) then
                ann_av_bvol = ann_av_bvol / tidesPerYear
                write(saveAnnBvolUnit,"(F10.1)") ann_av_bvol 
            end if
        end if
        
        if (outputYearFlag) then
            ! Save some data on the position of the cliff toe
            ! i.e. the cliff toe position is saved once per year.
            write(saveRockContourUnit,"(F12.4)") cliffy(outputContourLevel_e, :)

            do section = 1, nYsections
                write(saveShoreContourUnit,"(F12.4)") &
                    max(cliffy(outputContourLevel_e, section), beachy(outputContourLevel_e, section))
            end do
            
            !Save info on where seawalls and revetments are
            write(saveSeawallActiveUnit,'(I1)') seawall_active

            !Save info on where groynes are
            if (year >= firstYearTransportOutput) then
                groyneEffectAnnualSum = groyneEffectAnnualSum / tidesPerYear
                write(saveGroyneEffectUnit,'(f6.5)') groyneEffectAnnualSum
            end if    

            if (year >= firstYearTransportOutput) then
                write(saveSedTransAnnUnit,"(F14.1)") QAnn
                write(savePotTransAnnUnit,"(F14.1)") QpAnn
            end if
            
            outputYearFlag = .false.
        end if

        if (firstYearProfileOutput /= -1000000 .and. &
            year >= firstYearProfileOutput .and. year <= lastYearProfileOutput .and. &
            mod(year - firstYearProfileOutput, profileOutputTimestep) == 0) then
            call output_profiles(saveRockProfilesUnit, &
                cliffy, 1, nce, firstSectionProfileOutput, lastSectionProfileOutput, success_code)
            call output_profiles(saveBeachProfilesUnit, &
                beachy, bBottom(1), bTop(1), firstSectionProfileOutput, lastSectionProfileOutput, success_code)
        end if

    end if  ! End of the year end results clause

    ntide = ntide + 1

end subroutine update


!> FluidEarth entry point
!>
!> return current time as a julian day
function get_current_time(success_code)
    use global_data, only: nTide
    use local_data, only: time_start
    implicit none
    
    integer, parameter :: double = kind(1D0)
    
    real(kind=double) :: get_current_time    
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning
    
    real(kind=double) :: daysElapsed
    real(kind=double) :: currentTime
    
    ! Tide period is 12.4164306 hours based on 706 tides per year
    daysElapsed = ((nTide - 1) * 12.4164306D0) / 24.0D0

    currentTime = time_start + daysElapsed

    get_current_time = ANINT (currentTime * 1.0D8) / 1.0D8
    success_code = 0
    
end function get_current_time

!> Output cliffy or beachy values for a range of section indices
subroutine output_profiles(unitNo, values, firstElement, lastElement, firstSection, lastSection, success_code)
    use global_data, only: nce, nYsections, msl_OD, msl_m
    use setup_module, only: dz, relative_to_OD
    use utilities_module
    implicit none
    
    integer, parameter :: double = kind(1D0)
    
    integer, intent(in) :: unitNo !< FORTRAN unit number for already open output file
    real(kind=double), intent(in), dimension(nce, nYsections) :: values !< cliffy or beachy values 
    integer, intent(in) :: firstElement !< first cliff element index
    integer, intent(in) :: lastElement !< last cliff element index
    integer, intent(in) :: firstSection !< first beach section index
    integer, intent(in) :: lastSection !< last beach section index
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

    integer :: section
    
    do section = firstSection, lastSection
        ! Elevation (to ODN) of the lowest element in the output
        write(unitNo,"(F12.4)") relative_to_OD(firstElement)
        write(unitNo,"(F12.4)") values(firstElement: lastElement, section)
    end do
    
    success_code = 0

end subroutine  output_profiles


end module timestep_module