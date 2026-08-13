!**********************************************
!*
!*       global_data.f90       
!*
!**********************************************
!
!    Purpose:
!>    Declares all the global data
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

    module global_data

    implicit none
    save
    
    integer, parameter :: double = kind(1D0)
    private double

    real(kind=double) :: tstep_secs !< Time step in seconds
    real(kind=double) :: pi, two_pi, half_pi
    real(kind=double) :: radcon
    real(kind=double) :: breakrat    !< Ratio of breaker height and water depth
    real(kind=double) :: rhow, g
    real(kind=double) :: msl_m !< Initial mean sea level from base of model (m)
!    real(kind=double) :: min_sect_height
    real(kind=double) :: msl_OD !< Mean sea level from Ordinance Datum (m)
    real(kind=double) :: baseline_OD !< Model baseline from Ordinance Datum (m)
    real(kind=double) :: sea_level_next_year !< Sea level at the end of the current year
    real(kind=double) :: slr_per_tide !< Sea level rise per tide (m)
    real(kind=double) :: block_size !< Largest depth that can be eroded in one tide
    real(kind=double) :: depthosc !< Depth of the offshore contour at high tide
    real(kind=double) :: a !< Bruun constant
    real(kind=double) :: a_up, a_low
    real(kind=double) :: dsf_erode !< The distance between points in the shape function This is as a fraction of the distance from the water line to the break point
    real(kind=double) :: dsf_drift !< The distance between points in the shape function This is as a fraction of the distance from the water line to the break point
    real(kind=double) :: sea_rise !< Offset (m) of curent sea lvel from height of msl_e platform element
    real(kind=double) :: tidecounter, mid_hourly_level
    real(kind=double), allocatable, dimension(:) :: tidal_amp
    real(kind=double) :: period
    real(kind=double) :: mean_dbreak,mean_half_tide,mean_setup,mean_setdown
    real(kind=double) :: totalbvolume
    real(kind=double) :: heightsurge !< Not used, left at zero
    real(kind=double) :: fines_volume

    

! MD 300304 for weighting of waves on the basis of distance
    real(kind=double), allocatable, dimension(:) :: hwp_temp,anglewp_temp
    real(kind=double), allocatable, dimension(:, :) :: distwp 

! MD 101003 for groyne bypassing

    real(kind=double), allocatable, dimension(:) :: beach_temp
    real(kind=double), allocatable, dimension(:) :: cliff_temp

    real(kind=double), allocatable, dimension(:,:) :: cliffy !< discretised cliff and platform y value
    real(kind=double), allocatable, dimension(:,:) :: beachy !< discretised platform y value (distance of beach surface from baseline)
    real(kind=double), allocatable, dimension(:,:) :: talusy !< talus y value (horizontal distance of talus surface from baseline)
    real(kind=double), allocatable, dimension(:,:) :: upbeachy !< distance of surface of upper beach from baseline, active during storm conditions
    real(kind=double), allocatable, dimension(:,:) :: lowbeachy !< distance of surface of lower beach from baseline, active during storm conditions
    real(kind=double), allocatable, dimension(:,:) :: beachdepth1
    real(kind=double), allocatable, dimension(:,:) :: startprofiles !< initial shape of cliffy (if more than one Y section)
    real(kind=double), allocatable, dimension(:,:) :: sf_drift_new,sf_erode_new
    real(kind=double), allocatable, dimension(:,:) :: tempslrc
    real(kind=double), allocatable, dimension(:) :: bruun !< 'Bruun' shape for beach (elevation,horizontal distance)
    real(kind=double), allocatable, dimension(:) :: upbruun !< upper beach 'Bruun' shape (elevation,horizontal distance)
    real(kind=double), allocatable, dimension(:) :: lowbruun !< lower beach 'Bruun' shape (elevation,horizontal distance)
    real(kind=double), allocatable, dimension(:,:) :: sect_drift
    real(kind=double), allocatable, dimension(:,:) :: sect_erode, sect_erode_new

    real(kind=double), allocatable, dimension(:) :: hwp,anglewp !< Wave height and angle at wave point
    real(kind=double), allocatable, dimension(:) :: seawall_os !< horizontal distance of the seawall (if 'active') from the baseline
    real(kind=double), allocatable, dimension(:) :: startprofile !< initial shape of cliffy (if only one Y section)
    real(kind=double), allocatable, dimension(:) :: cliffheights !< Heights of cliff elements
    real(kind=double), allocatable, dimension(:) :: berm
    real(kind=double), allocatable, dimension(:) :: beachos !< Beach offset at each beach section
    real(kind=double), allocatable, dimension(:) :: angleosc
    real(kind=double), allocatable, dimension(:) :: beachdepth !< depth of the beach, used to determine whether the platform below is protected
    real(kind=double), allocatable, dimension(:) :: upbeach_os !< Upper beach offset at each beach section
    real(kind=double), allocatable, dimension(:) :: lowbeach_os !< Lower beach offset at each beach section
    real(kind=double), allocatable, dimension(:) :: upbvolume
    real(kind=double), allocatable, dimension(:) :: lowbvolume
    real(kind=double), allocatable, dimension(:) :: sand_fraction
    real(kind=double), allocatable, dimension(:,:) :: resistance !< strength of the platform/cliff material, a calibration term
    real(kind=double), allocatable, dimension(:,:) :: sedimentContent !< cache sediment content in useVAgrid mode
    real(kind=double), allocatable, dimension(:) :: tvolume !< talus volume
    real(kind=double), allocatable, dimension(:) :: cliffTopElevations
    real(kind=double), allocatable, dimension(:) :: talus_offsets
    real(kind=double), allocatable, dimension(:) :: clifftoe_pos !< horizontal distance of the seawall (if 'active') from the baseline
    real(kind=double), allocatable, dimension(:) :: cliff_heights
    real(kind=double), allocatable, dimension(:) :: bvolume
    real(kind=double), allocatable, dimension(:) :: obvolume
    real(kind=double), allocatable, dimension(:) :: ubvolume
    real(kind=double), allocatable, dimension(:) :: offshore_st
    real(kind=double), allocatable, dimension(:) :: sect_height
    real(kind=double), allocatable, dimension(:) :: section_offsets !< initial distance from baseline to cliffy (gives the initial planshape of the coast)
    real(kind=double), allocatable, dimension(:) :: this_tide_drift
    real(kind=double), allocatable, dimension(:) :: cerck1
    real(kind=double), allocatable, dimension(:) :: groyneEffectAnnualSum !< Effect of groyne at a q point where zero means "no effect" and 1 means "completely blocks flow" 
    real(kind=double), allocatable, dimension(:) :: QAnn, QpAnn
    
    real(kind=double), allocatable, dimension(:) :: setup, setdown, h
    real(kind=double), allocatable, dimension(:) :: anglewob, cgroup, theta
    real(kind=double), allocatable, dimension(:) :: setup_qpoint, setdown_qpoint, h_qpoint
    real(kind=double), allocatable, dimension(:) :: anglewob_qpoint, cgroup_qpoint, theta_qpoint
    real(kind=double), allocatable, dimension(:) :: dbreak
    real(kind=double), allocatable, dimension(:) :: half_tide
        !< Half tide, this is a complicated one. Ultimately it is for the calculation of the tide-integrated cross-shore distributions of 
        !< drift and erosion. These use the static distribution of drift or erosion. These are shape functions (sf) sf_drift & sf_erode. 
        !< These have abscissa increments (dx) of 1/40th of the breaker depth. Half tide is the tidal amplitude measured in these increments.

    real(kind=double), allocatable, dimension(:) :: surface

    real(kind=double), allocatable, dimension(:) :: beach_addition
    real(kind=double), allocatable, dimension(:) :: total_beach_addition
    real(kind=double), allocatable, dimension(:) :: fines_addition !< Rate of fines addition (m3 / timestep)

    ! Modifications suggested by JT first implimented in 120208a
        !allocs pulled out of platform_erosion (JT221107a)
    real(kind=double), allocatable, dimension(:) :: this_erode !< Cross-shore distribution of erosion
    real(kind=double), allocatable, dimension(:) :: m !< slope of cliffy
    real(kind=double), allocatable, dimension(:) :: bdepth
    real(kind=double), allocatable, dimension(:) :: temp5 !< Related to attenuation factor
    real(kind=double), allocatable, dimension(:) :: battenuation !< factor representing the protection of the platform/cliff by sufficiently thick beaches
    real(kind=double), allocatable, dimension(:) :: crossbeach_st
    real(kind=double), allocatable, dimension(:) :: seawall_effect
    real(kind=double), allocatable, dimension(:) :: v_level_odn
    !end of platform allocs (JT221107a)

    ! More JT mods see 86(22)
    integer :: hwp_i, anglewp_i
    character :: T_c

    ! End of JT mods see 86(22)

    real(kind=double), allocatable, dimension(:, :) :: cglut, clut 
    real(kind=double), dimension(50, 2) :: sf_erode !< Shape function - vertical distribution of erosion with breaking waves ALL SHAPE FUNCTIONS MUST HAVE THE SAME dX
    real(kind=double), dimension(50, 2) :: sf_drift !< Shape function - distribution of alongshore drift across the surf zone ALL SHAPE FUNCTIONS MUST HAVE THE SAME dX
    real(kind=double), dimension(2) :: hour_section_limit !< hour_section_limits are the vertical limits
                !! in metres of the range of the tide in this hour, i.e. the tide moves from 
                !! hour_section_limit(1) to hour_section_limit(2) in one hour

    integer :: run_duration !< counts in tides at tidesPerYear tides per year
    integer :: beachHeightE !< height expressed as number of cliff elements
    integer :: upBeachHeightE !< height expressed as number of cliff elements
    integer :: lowBeachHeightE !< height expressed as number of cliff elements
    integer :: beachCrestLevelE  !< Storm Surge level, in cliff elements. 
      !! This is used to set the level of the berm, or top of beach
    integer :: msl_e !< Number of the platform element at msl
    integer :: ngroynes !< number of groynes
    integer :: nQpoints !< Number of Q points
    integer :: nYsections !< Number of profiles (number of beach sections)
    integer :: nce !< Number of cliff and platform vertical elements
    integer :: alloc_error
    integer :: groyne1,groyne2,groyne_flag_1
    integer :: np_sf_new_erode, np_sf_new_drift 
    integer :: waveleader1,waveleader2,waveleader3,waveleader4
    integer :: eventcounter
    integer :: ntide
!    integer :: tide_hour
!    integer :: tidal_period
    integer :: proceed_switch
    integer :: bottom, top
    logical :: outputYearFlag !< Output this year if true

    integer, allocatable, dimension(:) :: wpr !< Wave Point reference (wave point number for each Q point)
    integer, allocatable, dimension(:) :: top_drift,bot_drift
    integer, allocatable, dimension(:) :: beachelev !< Elevation of the berm as cliff element number
    integer, allocatable, dimension(:) :: bTop !< Elevation of the beach top as cliff element number
    integer, allocatable, dimension(:) :: bBottom !< Elevation of the beach bottom as cliff element number
    integer, allocatable, dimension(:) :: ubTop !< Elevation of the upper beach top as cliff element number
    integer, allocatable, dimension(:) :: ubBottom !< Elevation of the upper beach bottom as cliff element number
    integer, allocatable, dimension(:) :: lbTop !< Elevation of the lower beach top as cliff element number
    integer, allocatable, dimension(:) :: lbBottom !< Elevation of the lower beach bottom as cliff element number
    integer, allocatable, dimension(:) :: np_cliff
    integer, allocatable, dimension(:) :: top_erode, bot_erode 
    integer, allocatable, dimension(:) :: stormflag
    integer, allocatable, dimension(:) :: seawall_active !< existance of a seawall at this location
    integer, allocatable, dimension(:) :: seawall_construction
    integer, allocatable, dimension(:) :: seawall_removal
    integer, allocatable, dimension(:) :: groyne_construction
    integer, allocatable, dimension(:) :: groyne_removal
    integer, allocatable, dimension(:) :: v_ref
    
    real(kind=double), allocatable, dimension(:) :: profile_y !< Profile point Y values, as input
    real(kind=double), allocatable, dimension(:) :: profile_z !< Profile point Z values, as input
    integer, allocatable, dimension(:) :: base_profile_start !< Index in profile_y and profile_z where base geology profile for beach section starts
    integer, allocatable, dimension(:) :: base_profile_start2 !< Index in profile_y and profile_z of last point before YmaxH
    integer, allocatable, dimension(:) :: base_profile_end !< Index in profile_y and profile_z where base geology profile for beach section ends
    integer, allocatable, dimension(:, :) :: layer_profile_start !< Index in profile_y and profile_z where layer profile for layer within each beach section starts
    integer, allocatable, dimension(:, :) :: layer_profile_start2 !< Index in profile_y and profile_z of last point before YmaxH
    integer, allocatable, dimension(:, :) :: layer_profile_end !< Index in profile_y and profile_z where layer profile for layer within each beach section ends
    real(kind=double), allocatable, dimension(:, :) :: base_profile !< base geology profile values for each section at dy horizontal spacing
    real(kind=double), allocatable, dimension(:) :: base_strength !< Base geology strength value for each section
    real(kind=double), allocatable, dimension(:) :: base_sediment !< Base geology sediment content value for each section
    integer, allocatable, dimension(:) :: layer_count !< count of layer profiles for each section
    real(kind=double), allocatable, dimension(:, :, :) :: layer_profile !< layer profile depth values for each section at dy horizontal spacing
    real(kind=double), allocatable, dimension(:, :) :: layer_strength !< layer strength value for each section
    real(kind=double), allocatable, dimension(:, :) :: layer_sediment !< layer sediment content value for each section

    real(kind=double) :: nActiveLayers !< Number of sediment layers currently active.
        !< Initially set from the number of layer strength and sediment values present in the layerData.txt file.
        !< May change during the run as layers are added and exhausted.
    
    real(kind=double), allocatable, dimension(:) :: left_exclusion_angle !< angle anti-clockwise from shore normal, beyond which waves are excluded (radians).
    real(kind=double), allocatable, dimension(:) :: right_exclusion_angle !< angle clockwise from shore normal, beyond which waves are excluded (radians).
    
    integer, allocatable, dimension(:) :: vellinga_sections !< Section numbers of the beach sections subject to vellinga effects
        !< Only used if vellingaActive parameter is .true.
    end module global_data
