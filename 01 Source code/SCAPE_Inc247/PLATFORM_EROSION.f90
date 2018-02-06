!**********************************************
!*
!*       PLATFORM_EROSION.f90       
!*
!**********************************************
!
!>    This SCAPE module describes the erosion of sections of platform, and lower cliff. 
!>
!>    Representation of 'storm' beach behaviour. 
!>    This was introduced to represent beach behaviour during a storm. It was 
!>    decided to do this because:
!>    (a) the regional beach is unresponsive to individual waves, so is more 
!>    stable than reality in that it is less likley to empty and reveal the 
!>    platform to wave attack, 
!>    (b) during storms material is moved offshore making the platform more 
!>    vulnerable, 
!>    (c) the storm beach profile is lower and flatter, again leaving the 
!>    platform more vulnerable. 
!>
!>    The algorithms here try to represent these things by dividing the beach 
!>    into two parts when a storm occurs. Each part (upper and lower) has its own 
!>    volume and its own discretised surface. The upper beach is steeper than the 
!>    usual beach, the lower beach is flatter. Material can move from the upper to
!>    the lower during a storm. The rates of movement are determined through a user defined look-up table. 
!>    Once the storm is over the beach returns to its 'normal' form, and its colume recovers using a behavioural rule.
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


module platform_erosion_module
    implicit none
    
contains
        
    !> Apply platform erosion to the whole beach.
    !>
    !> Calculate wave agression. Calculate the cross-shore distribution of erosion for this tide.
    !>
    !> Allow for storm waves, calling CROSSBEACH_EXCHANGE() and handling the upper and lower beach separately in a storm.
    !>
    !> Calculate the erosion, allowing for the presence of structures.
    subroutine PLATFORM_EROSION(nYsections, nce, &
        upbeachy, &
        lowbeachy, &
        beachelev, cliffy, beachy, talusy, this_erode, m, &
        temp5, bdepth, battenuation, &
        seawall_effect, cliffheights, &
        beach_addition, fines_addition, h, &
        dbreak, half_tide, sand_fraction, resistance, tvolume, &
        seawall_os, seawall_active, setup, setdown, &
        sf_erode, np_sf_erode, period, tstep_secs, &
        block_size, dsf_erode, heightsurge, &
        msl_m, bottom, top, success_code)
        
    use xshore_dist_t_module
    use setup_module
    use vertical_grid
    use utilities_module
    use global_data, only: sedimentContent
    use dump_module
    implicit none
    
    integer, intent(in) :: nYsections !< number of beach sections
    integer, intent(in) :: nce !< number of platform elements
    real(kind=double), intent(inout), dimension(nce,nYsections) :: upbeachy !< distance of surface of upper beach from baseline, active during storm conditions
    real(kind=double), intent(inout), dimension(nce, nYsections) :: lowbeachy !< distance of surface of lower beach from baseline, active during storm conditions
    integer, intent(in), dimension(nYsections) :: beachelev !< Elevation of the berm as cliff element number
    real(kind=double), intent(inout), dimension(nce, nYsections) :: cliffy !< Discretised cliff and platform y value
    real(kind=double), intent(inout), dimension(nce, nYsections) :: beachy !< Discretised platform y value (distance of beach surface from baseline)
    real(kind=double), intent(in), dimension(nce, nYsections) :: talusy !< Talus y value (horizontal distance of talus surface from baseline)
!    real(kind=double), intent(in), dimension(nce, nYsections) :: sect_erode
    real(kind=double), intent(inout), dimension(nce) :: this_erode !< Cross-shore distribution of erosion
    real(kind=double), intent(inout), dimension(nce) :: m !< slope of cliffy
    real(kind=double), intent(inout), dimension(nce) :: bdepth !< Beach depth
    real(kind=double), intent(inout), dimension(nce) :: temp5 !< Related to attenuation factor
    real(kind=double), intent(inout), dimension(nce) :: battenuation !< factor representing the protection of the platform/cliff by sufficiently thick beaches
    real(kind=double), intent(inout), dimension(nce) :: seawall_effect !< Erosion due to seawall
    real(kind=double), intent(inout), dimension(nce) :: cliffheights !< Heights of cliff elements
    real(kind=double), intent(inout), dimension(nYsections) :: beach_addition !< The material moving from the talus or from the platform to the beach
    real(kind=double), intent(inout), dimension(nYsections) :: fines_addition !< Material lost as fines
    real(kind=double), intent(in), dimension(nYsections) :: h !< Wave height
    real(kind=double), intent(in), dimension(nYsections) :: dbreak !< Wave breaker depth for each beach section
    real(kind=double), intent(in), dimension(nYsections) :: half_tide !< Tidal amplitude measured in shape function increments
    real(kind=double), intent(inout), dimension(nYsections) :: sand_fraction !< Sand fraction (sediment content) ratio for each beach section
    real(kind=double), intent(inout), dimension(nce, nYsections) :: resistance !< strength of the platform/cliff material
    real(kind=double), intent(inout), dimension(nYsections) :: tvolume !< talus volume for each beach section
    real(kind=double), intent(in), dimension(nYsections) :: seawall_os !< horizontal distance of the seawall for each beach section
    integer, intent(in), dimension(nYsections) :: seawall_active !< existance of a seawall for each beach section
!    integer, intent(in), dimension(nYsections) :: top_erode
!    integer, intent(in), dimension(nYsections) :: bot_erode
    real(kind=double), intent(in), dimension(nYsections) :: setup
    real(kind=double), intent(in), dimension(nYsections) :: setdown
    real(kind=double), intent(in), dimension(50,2) :: sf_erode !< Shape function - vertical distribution of erosion with breaking waves
        !< ALL SHAPE FUNCTIONS MUST HAVE THE SAME dX
    integer, intent(in) :: np_sf_erode !< Number of points in erosion shape function
    real(kind=double), intent(in) :: period !< Wave period
    real(kind=double), intent(in) :: tstep_secs !< Time step in seconds
    real(kind=double), intent(in) :: block_size !< Largest depth that can be eroded in one tide
    real(kind=double), intent(in) :: dsf_erode !< The distance between points in the shape function
    real(kind=double), intent(in) :: heightsurge !< Not used, left at zero
    real(kind=double), intent(in) :: msl_m !< Initial mean sea level from base of model (m)
    integer, intent(inout) :: bottom !< Bottom active cliff element
    integer, intent(inout) :: top !< Top active cliff element
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

!*** Local variables

    real(kind=double) :: wagression ! wave erosive potential (following Kamphuis), a function of wave height and period
    
    real(kind=double) :: sum_dtalus_minus ! material removed from the talus (through wave action)
    real(kind=double) :: sum_dtalus_plus
    real(kind=double) :: sum_dbeach_plus ! material added to the beach (through wave action on talus)
    real(kind=double) :: save_sum_dtalus_minus
    real(kind=double) :: sum_dtalus_minus_XSC ! components multiplied by sediment content
    real(kind=double) :: sum_dbeach_plus_XSC ! components multiplied by sediment content
    real(kind=double) :: sum_dtalus_minus_X1MSC ! components multiplied by (1 - sediment content)
    real(kind=double) :: sum_dbeach_plus_X1MSC ! components multiplied by (1 - sediment content)
    
    real(kind=double) :: talus_element_width ! horizontal width of the talus

    real(kind=double) :: toler  
    real(kind=double) :: dtalus !< change in talus volume

    integer :: j, index
    integer :: section
    integer :: yIndex
    real(kind=double) :: tstepRemaining ! number of seconds remaining in current timestep
    real(kind=double) :: recession ! recession amount (m)
    real(kind=double) :: recessionTime ! number of seconds taken by recession
    real(kind=double) :: totalRecession ! total recession amount for current timestep
    real(kind=double) :: saveRecession ! temporary storage
    real(kind=double) :: talusRecessionTime ! time taken to recede through the talus
    logical :: talusConsidered
    logical :: blocksizeLimitReached

    success_code = 0
    toler=1.0d-6

     call dump_cliff_real("platfor_erosion_235", "beachy", beachy)
     call dump_cliff_real("platfor_erosion_235", "upbeachy", upbeachy)
     call dump_cliff_real("platfor_erosion_235", "lowbeachy", lowbeachy)

    sections: do section = firstActiveSection, lastActiveSection
    
		 
  !  if (h(section) < minHs) then
   !     cycle ! skip section with small inshore waves
   ! end if
    
    this_erode = 0.0D0
    m = 0.0D0
    temp5 = 0.0D0
    bdepth = 0.0D0
    battenuation = 0.0D0


!******************************
!    Calculate wave agression
!******************************

! Wave agression is calculated using the equation of Kamphuis

    wagression = h(section)**(13.0d0/4.0d0) * period**(3.0d0/2.0d0)

    !*******************************************************
    !    Calculate the cross-shore distribution of erosion
    !    for this tide
    !*******************************************************
		
        
      call xshore_dist_t(nint(half_tide(section)),sf_erode,nce,    &
        dbreak(section),dsf_erode,setup(section),setdown(section),heightsurge, np_sf_erode,    &
        2,nce,cliffheights,msl_m,        &
        this_erode,bottom,top,2, success_code)
      if (success_code < 0) return

        this_erode = this_erode/dz
    
!******* Code for the new representation of talus WkBk p141(11) *******
    do j=bottom,top
        temp5(1) = max(cliffy(j-1,section),talusy(j-1,section)) - max(cliffy(j,section),talusy(j,section))
        if(abs(temp5(1)) > toler)then
            m(j) = abs(dz/temp5(1))
          else
            m(j) = 1.0d0/toler**2
          endif             
    end do
!******* End of new code for the new representation of talus (141(11)) *******


    temp5(1) = 0

    where (m > maxPlatformSlope_radians)
      m = maxPlatformSlope_radians   ! MW030304
    end where

!*****************************************************************
!    Calculate the beach depths over the platform, so its level 
!    of protection can be established.
!*****************************************************************

    if (h(section) < minStormWaveHeight)then ! not a storm, use the normal beach
        call BEACH_DEPTH(nYsections, nce, bottom, top, section, beachelev, cliffy, &
            h, temp5, bdepth, battenuation, beachy, lowbeachy, .false.)
    else
        call BEACH_DEPTH(nYsections, nce, bottom, top, section, beachelev, cliffy, &
            h, temp5, bdepth, battenuation, upbeachy, lowbeachy, .true.)
    endif
    

!*********************************************
!    Account for the presence of structures
!*********************************************

! Start of new 19 Sept 2002, WkBk p115(12)
    beach_addition(section) = 0
    fines_addition(section) = 0
    sum_dtalus_minus = 0
    sum_dtalus_plus = 0
    sum_dbeach_plus = 0
    sum_dtalus_minus_XSC = 0
    sum_dbeach_plus_XSC = 0
    sum_dtalus_minus_X1MSC = 0
    sum_dbeach_plus_X1MSC = 0
    talus_element_width = 0
    seawall_effect = 1.0d0        ! i.e. erosion allowed

! The effects of the revetments and seawalls are contained in vectors containing nce 
! numbers, one for each cliff element. If the seawall_effect at any elevation is 0
! then no erosion will occur there.

    if(seawall_active(section) == 1)then ! a seawall has been built at this section
        do j = bottom,top
            if(cliffy(j,section) <= seawall_os(section))then ! This elevation is protected
                cliffy(j,section) = seawall_os(section)
                seawall_effect(j) = 0        ! i.e. no erosion here
            endif
        enddo
    endif

!****************************
!    Calculate the erosion
!****************************

    do j = bottom, top    ! Work on elements within the range of erosion
        if (battenuation(j) > 0)then    ! The beach does not totally protect this element
            tstepRemaining = tstep_secs
            totalRecession = 0.0D0
            talusConsidered = .false.
            blocksizeLimitReached = .false.

            do while (tstepRemaining > 0.0D0)

                if (useVAgrid) then
                    if (cliffy(j, section) < dy) then
                        yIndex = 0
                    else
                        yIndex = int((cliffy(j, section) - toler) / dy)
                    end if

                    if (resistance(j, section) < 0.0D0) then ! get strength value from HA grid
                        do
                            call get_cell_values(relative_to_OD(j), yIndex, section, &
                                resistance(j, section), sedimentContent(j, section))
    
                            if (resistance(j, section) == 0.0D0) then ! cell empty
                                if (yIndex == 0) then
                                    cliffy(j, section) = 0.0D0 ! Set y to minimum allowed value
                                    resistance(j, section) = -1.0D0 ! Flag resistance unknown
                                    exit
                                else
                                    cliffy(j, section) = yIndex * dy ! recede to next cell
                                    yIndex = yIndex - 1
                                end if
                            else ! cell not empty
                                exit
                            end if
                        end do
                    end if
                end if ! useVAgrid
                
                if (cliffy(j, section) <= 0.0D0) then ! Do not allow y to go negative
                    exit
                end if

                recessionTime = 0.0D0
                talusRecessionTime = 0.0D0
                
                if (.not. talusConsidered .and. talusy(j, section) > cliffy(j, section)) then    
                    ! Talus is in front of the cliff at this element
                    recession = tstepRemaining * &
                        this_erode(j)*m(j)*battenuation(j)*seawall_effect(j)* wagression / &
                        (talusStrengthRatio * resistance(j, section))
                    talus_element_width = talusy(j,section) - cliffy(j,section)
                    
                    if (recession > talus_element_width) then
                        ! Recession penetrates cliff through talus
                        talusRecessionTime = tstepRemaining * talus_element_width / recession
                    
                        !******* MW030304 start *******
                        sum_dtalus_minus = sum_dtalus_minus + talus_element_width     ! Material added to the beach

                        if (useVAgrid) then
                            sum_dtalus_minus_XSC = &
                                sum_dtalus_minus_XSC + talus_element_width * sedimentContent(j, section)
                            sum_dtalus_minus_X1MSC = &
                                sum_dtalus_minus_X1MSC + talus_element_width * (1.0D0 - sedimentContent(j, section))
                        end if
                        
                        temp5(1) = cliffy(j-1,section) - cliffy(j,section)
                        
                        if(abs(temp5(1)) > toler)then
                          m(j) = abs(dz/temp5(1))
                        else
                          m(j) = 1.0d0/toler**2
                        endif             
    
                        temp5(1) = 0
                        if (m(j) > maxPlatformSlope_radians)then
                          m(j) = maxPlatformSlope_radians   ! MW030304
                        endif
                        !******* MW030304 end *******
                    
                        ! calculate recession into cliff
                        recession = (tstepRemaining - talusRecessionTime) * this_erode(j) * m(j) * &
                            battenuation(j) * seawall_effect(j) * wagression / resistance(j, section)
                    
                        if (recession > block_size) then
                            blocksizeLimitReached = .true.
                            saveRecession = recession
                            recession = block_size
                            recessionTime = (tstepRemaining - talusRecessionTime) * recession / saveRecession
                        else
                            recessionTime = tstepRemaining - talusRecessionTime
                        end if

                        if (useVAgrid .and. cliffy(j, section) - recession < yIndex * dy) then ! Recession to next cell
                            saveRecession = recession
                            recession = cliffy(j, section) - yIndex * dy
                            recessionTime = recessionTime * recession / saveRecession
                            resistance(j, section) = -1.0D0 ! Flag resistance unknown
                        end if
                        
                        sum_dbeach_plus = sum_dbeach_plus + recession    ! Material added to the beach

                        if (useVAgrid) then
                            sum_dbeach_plus_XSC = sum_dbeach_plus_XSC + recession * sedimentContent(j, section)
                            sum_dbeach_plus_X1MSC = sum_dbeach_plus_X1MSC + recession * (1.0D0 - sedimentContent(j, section))
                        end if
                        
                        cliffy(j,section) = cliffy(j,section) - recession    ! Update the profile
                    else
                        ! Recession does NOT penetrate cliff, only talus
                        talusRecessionTime = tstepRemaining
                        sum_dtalus_minus = sum_dtalus_minus + recession

                        if (useVAgrid) then
                            sum_dtalus_minus_XSC = sum_dtalus_minus_XSC + recession * sedimentContent(j, section)
                            sum_dtalus_minus_X1MSC = sum_dtalus_minus_X1MSC + recession * (1.0D0 - sedimentContent(j, section))
                        end if

                        recession = 0.0D0
                    endif
                else
                    ! Direct wave attack on the cliff, no talus
                    recession = tstepRemaining * &
                        this_erode(j) * m(j) * battenuation(j) * seawall_effect(j) * wagression / resistance(j, section)

                    if (totalRecession + recession > block_size) then
                        blocksizeLimitReached = .true.
                        saveRecession = recession
                        recession = block_size - totalRecession
                        recessionTime = tstepRemaining * recession / saveRecession
                    else
                        recessionTime = tstepRemaining
                    end if

                    if (useVAgrid .and. cliffy(j, section) - recession < yIndex * dy) then ! Recession to next cell
                        saveRecession = recession
                        recession = cliffy(j, section) - yIndex * dy
                        recessionTime = recessionTime * recession / saveRecession
                        resistance(j, section) = -1.0D0 ! Flag resistance unknown
                    end if
                    
                    cliffy(j,section) = cliffy(j,section) - recession    ! Update the profile
                    sum_dbeach_plus = sum_dbeach_plus + recession    ! Material heading for the beach

                    if (useVAgrid) then
                        sum_dbeach_plus_XSC = sum_dbeach_plus_XSC + recession * sedimentContent(j, section)
                        sum_dbeach_plus_X1MSC = sum_dbeach_plus_X1MSC + recession * (1.0D0 - sedimentContent(j, section))
                    end if
                endif
                
                talusConsidered = .true.

                if (blocksizeLimitReached) then
                    exit
                else
                    tstepRemaining = tstepRemaining - talusRecessionTime - recessionTime
                    totalRecession = totalRecession + recession
                end if
            end do ! while (tstepRemaining > 0.0D0)
        endif ! (battenuation(j) > 0)
    enddo ! j = bottom, top

    ! MD 170504 -- It is possible, under rare circumstances, for erosion to occur below a seawall
    ! with subsequent slumping adding sediment to beaches.  The code below prevents slumping
    ! by updating the profile immediately after erosion has occurred. 
  
    if(seawall_active(section) == 1)then ! a seawall has been built at this section
        do j = bottom,top
            if(cliffy(j,section) <= seawall_os(section))then ! This elevation is protected
                cliffy(j,section) = seawall_os(section)
                seawall_effect(j) = 0        ! i.e. no erosion here
            endif
        enddo
    endif

    dtalus = (sum_dtalus_plus - sum_dtalus_minus) * dz * dx
    
    if(tvolume(section) + dtalus < 0)then ! the talus empties, so
        save_sum_dtalus_minus = sum_dtalus_minus
        sum_dtalus_minus = tvolume(section)/(dz * dx) + sum_dtalus_plus ! reduce the loss

        if (useVAgrid .and. abs(save_sum_dtalus_minus) > toler) then
            sum_dtalus_minus_XSC = sum_dtalus_minus_XSC * sum_dtalus_minus / save_sum_dtalus_minus
            sum_dtalus_minus_X1MSC = sum_dtalus_minus_X1MSC * sum_dtalus_minus / save_sum_dtalus_minus
        end if
        
        tvolume(section) = 0
    else ! it still has some volume
        tvolume(section) = tvolume(section)  + dtalus
    end if
    
    ! Volume added to beach from the recession of the cliff and talus
    if (useVAgrid) then
        beach_addition(section) = (sum_dtalus_minus_XSC + sum_dbeach_plus_XSC) * dz * dx
        fines_addition(section) = (sum_dtalus_minus_X1MSC + sum_dbeach_plus_X1MSC) * dz * dx
    else
        beach_addition(section) = (sum_dtalus_minus + sum_dbeach_plus) * dz * dx * sand_fraction(section)
        fines_addition(section) = (sum_dtalus_minus + sum_dbeach_plus) * dz * dx * (1-sand_fraction(section))
    end if

! End of new 19 Sept 2002, p115(12)

        end do sections

        return
        end subroutine PLATFORM_EROSION                             


!> Calculate the beach depths over the platform, so its level of protection can be established.
subroutine BEACH_DEPTH(nYsections, nce, bottom, top, section, beachelev, cliffy, &
        h, temp5, bdepth, battenuation, tempBeachy, lowbeachy, storm)
        
    use setup_module
    implicit none
    
    integer, intent(in) :: nYsections !< number of beach sections
    integer, intent(in) :: nce !< number of platform elements
    integer, intent(in) :: bottom !< Bottom active cliff element
    integer, intent(in) :: top !< Top active cliff element
    integer, intent(in) :: section !< Beach Section index
    integer, intent(in), dimension(nYsections) :: beachelev !< Elevation of the berm as cliff element number at each Beach Section
    real(kind=double), intent(in), dimension(nce, nYsections) :: cliffy !< discretised cliff and platform y value
    real(kind=double), intent(in), dimension(nYsections) :: h !< Wave height
    real(kind=double), intent(out), dimension(nce) :: temp5
    real(kind=double), intent(out), dimension(nce) :: bdepth !< Beach depth
    real(kind=double), intent(out), dimension(nce) :: battenuation !< factor representing the protection of the platform/cliff by sufficiently thick beaches
    real(kind=double), intent(inout), dimension(nce, nYsections) :: tempBeachy !< beachy or upbeachy if storm
    real(kind=double), intent(in), dimension(nce, nYsections) :: lowbeachy !< needed if storm
    logical, intent(in) :: storm !< Storm conditions apply

    integer :: j
    integer :: k
    
    ! Corrected version Wk.Bk. 91(9) 127(9)
    bdepth = 0.D0    ! Beach depth
    
    do j = bottom, top        ! The range where erosion occurs
        if (EFFECTIVE_BEACHY(tempBeachy(j, section), lowbeachy(j, section), storm) > cliffy(j, section)) then    ! Beach is seaward of platform at this elevation
            k = j
           
            do while (EFFECTIVE_BEACHY(tempBeachy(k, section), lowbeachy(k, section), storm) > cliffy(j, section) &
                .and. k < (nce - 1))
                k = k + 1    ! Count up the number of elements of beach cover
                
                if (k >= beachelev(section) .and. bermSlope > 0.0D0)then    ! We are over the berm
                    tempBeachy(k, section) = tempBeachy(k - 1, section) - (dz * bermSlope)
                endif
                
                if((k - j - 1) * dz > beachDisturbanceRatio * h(section)) exit ! The beach is thick at this location and so is completely protective
            end do
            
            bdepth(j) = (k - j) * dz
        end if
    end do
    
    temp5 = bdepth / (beachDisturbanceRatio * h(section))
    where (temp5 > 1.0D0)
        temp5 = 1.0D0
    end where
    battenuation = 1.0D0 - temp5    ! This is the attenuation of the erosive capability of the 
    ! beach. 1 indicates no beach protection.

    if (tempBeachy(beachelev(section) + 1, section) > 0)then
        ! reset the top of the beach to a flat slope, return it to 'normal'
        ! following the modification of tempBeachy to represent the berm slope
        do j = beachelev(section) + 1, nce
            tempBeachy(j, section) = 0
        enddo
    endif
    
    return
end subroutine BEACH_DEPTH



!> Support function for subroutine BEACH_DEPTH
real(kind=double) function EFFECTIVE_BEACHY(tempBeachy, lowbeachy, storm)
    
    use setup_module, only: double
    implicit none
    
    real(kind=double), intent(in) :: tempBeachy !< beachy or upbeachy if storm
    real(kind=double), intent(in) :: lowbeachy !< required if storm
    logical, intent(in) :: storm !< Storm conditions apply
    
    if (storm) then
        EFFECTIVE_BEACHY = max(tempBeachy, lowbeachy)
    else
        EFFECTIVE_BEACHY = tempBeachy
    endif
    
end function EFFECTIVE_BEACHY


!> Calculate crossbeach exchange during storm waves.
subroutine CROSSBEACH_EXCHANGE(i, nYsections, nce, ubTop, ubBottom, lbTop, lbBottom, &
        upBeachHeightE, lowBeachHeightE, &
        stormflag, crossbeach_st, &
        upbvolume, upbruun, upbeach_os, upbeachy, &
        lowbvolume, lowbruun, lowbeach_os, lowbeachy, &
        cliffy, talusy, beachdepth, bvolume, h, offshore_st, &
        period, tstep_secs)
        
    use beach_location_module
    use setup_module
    implicit none
    
    integer, intent(in) :: i !< beach section index
    integer, intent(in) :: nYsections !< number of beach sections
    integer, intent(in) :: nce !< number of platform elements
    integer, intent(in), dimension(nYsections) :: ubTop !< indices of upper beach top elements
    integer, intent(in), dimension(nYsections) :: ubBottom !< indices of upper beach bottom elements
    integer, intent(in), dimension(nYsections) :: lbTop !< indices of lower beach top elements
    integer, intent(in), dimension(nYsections) :: lbBottom !< indices of lower beach bottom elements
    integer, intent(in) :: upBeachHeightE !< upper beach height in number of elements
    integer, intent(in) :: lowBeachHeightE !< lower beach height in number of elements
    integer, intent(inout), dimension(nYsections) :: stormflag !< 0 indicates first timestep of sequence of storm waves
    real(kind=double), intent(inout), dimension(nYsections) :: crossbeach_st !< Cross beach sediment transport
    real(kind=double), intent(inout), dimension(nYsections) :: upbvolume !< upper beach volume
    real(kind=double), intent(in), dimension(upBeachHeightE) :: upbruun !< upper beach 'Bruun' shape
    real(kind=double), intent(inout), dimension(nYsections) :: upbeach_os !< Upper beach offset at each beach section
    real(kind=double), intent(inout), dimension(nce,nYsections) :: upbeachy !< distance of surface of upper beach from baseline, active during storm conditions
    real(kind=double), intent(inout), dimension(nYsections) :: lowbvolume !< lower beach volume
    real(kind=double), intent(in), dimension(lowBeachHeightE) :: lowbruun !< lower beach 'Bruun' shape
    real(kind=double), intent(inout), dimension(nYsections) :: lowbeach_os !< Lower beach offset at each beach section
    real(kind=double), intent(inout), dimension(nce, nYsections) :: lowbeachy !< distance of surface of lower beach from baseline, active during storm conditions
    real(kind=double), intent(in), dimension(nce, nYsections) :: cliffy !< Discretised cliff and platform y value
    real(kind=double), intent(in), dimension(nce, nYsections) :: talusy !< Talus y value (horizontal distance of talus surface from baseline)
    real(kind=double), intent(inout), dimension(nce) :: beachdepth !< Depth of the beach
    real(kind=double), intent(in), dimension(nYsections) :: bvolume !< Beach volume
    real(kind=double), intent(in), dimension(nYsections) :: h !< Wave height
    real(kind=double), intent(in), dimension(nYsections) :: offshore_st !< Offshore sediment transport
    real(kind=double), intent(in) :: period !< Wave period
    real(kind=double), intent(in) :: tstep_secs !< Time step in seconds

    integer :: index
    real(kind=double) :: cbstrLookup
    
    if(stormflag(i) == 0)then ! this is the first timestep of 
        ! this sequence of storm waves
        stormflag(i) = 1
        ! Obtain values for the upper and lower beach volumes
        upbvolume(i) = bvolume(i)/2        !Upper beach volume
        lowbvolume(i) = upbvolume(i)    !lower beach volume
    else ! This is not the first timestep, and we have to 
    !account for the loss of material to the offshore, this 
    ! will have already been calculated in SEDIMENT_TRANSPORT, but the beach 
    ! volumes in the storm sequence 
    ! are disconnected from the rest of the model
        lowbvolume(i) = lowbvolume(i) + offshore_st(i)
        if(lowbvolume(i) < 0)then
            lowbvolume(i) = 0
        endif
    endif
    
    cbstrLookup = CROSSBEACH_STR_LOOKUP(h(i), period)
    crossbeach_st(i) = cbstrLookup * dx * tstep_secs / 3600 ! convert to units suitable for this timestep
    
    if(-1*crossbeach_st(i) > upbvolume(i))then
        crossbeach_st(i) = -1*upbvolume(i)
    endif
    upbvolume(i) = upbvolume(i) + crossbeach_st(i)
    lowbvolume(i) = lowbvolume(i) - crossbeach_st(i)

    
!*****************************************************    
!    Now find the positions of the UPPER beach
!*****************************************************
    
                call BEACH_LOCATION(i, nYsections, nce, upBeachHeightE, ubTop, ubBottom, &
                    upbruun, upbeach_os, upbeachy, &
                    cliffy, beachdepth, upbvolume, &
                    talusy, beach_upper)

!*****************************************************    
    ! Now find the positions of the LOWER beach
!*****************************************************
    
                call BEACH_LOCATION(i, nYsections, nce, lowBeachHeightE, lbTop, lbBottom, &
                    lowbruun, lowbeach_os, lowbeachy, &
                    cliffy, beachdepth, lowbvolume, &
                    upbeachy, beach_lower)
    
! End of crossbeach exchange calculations

end subroutine CROSSBEACH_EXCHANGE

end module platform_erosion_module
