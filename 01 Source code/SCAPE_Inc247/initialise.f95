!**********************************************
!*       initialise.f95       
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


module initialise_module
    implicit none

contains

    !> FluidEarth entry point, also called from RUN_SCAPE
    !>
    !> Initialisation
    subroutine initialise(args, success_code)
    
        use setup_module
        use global_data
        use local_data
        use dump_module
        implicit none
        
		character(*), intent(in) :: args !< Argument not used
		integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

        pi = 4.0d0 * datan(1.0d0)
        two_pi = 2.0D0 * pi
        half_pi = 0.5D0 * pi
        radcon = pi / 180.0d0
    
        ! Control diagnostic dumps
        dump_id = 1
        dump_enable = .false.

		success_code = 0 
		
    end subroutine initialise
	

    !> Read seawall or groyne construction or removal data
    subroutine read_construction_data(values, nValues, unitNo, filename, success_code, moreValues)
        use local_data, only: FILE_PATH
        use utilities_module
        use exceptions
        implicit none
    
        integer, parameter :: double = kind(1D0)
    
        integer, intent(out), dimension(nValues) :: values !< Array of values read in
        integer, intent(in) :: nValues !< Number of values in the array
        integer, intent(in) :: unitNo !< FORTRAN IO unit number
        character(*), intent(in) :: filename !< Construction data file name
        integer, intent(inout) :: success_code !< zero if OK
        real(kind=double), intent(out), dimension(nValues), optional :: moreValues !< Array of additional values read in

        integer :: iostat
        integer :: j

        success_code = do_open_text_input_file(unitNo, FILE_PATH, filename)
        if (success_code < 0) return
            
        if (present(moreValues)) then
            read(unitNo,*, iostat=iostat) (values(j), moreValues(j), j = 1,nValues)
        else
            read(unitNo,*, iostat=iostat) values
        end if
        
        success_code = SUCCESS_CODE_FROM_IOSTAT(unitNo, iostat, 'Reading '//trim(filename))
        close(unitNo)

    end subroutine read_construction_data

    
    !> FluidEarth entry point, also called from RUN_SCAPE.
    
    !> Prepare for run, including:
    !>
    !> Initialise parameters.
    !>
    !> Log parameters and output metadata.
    !>
    !> Open, read and close input files.
    !>
    !> Open waves and tides file.
    !>
    !> Allocate and initialise arrays.
    !>
    !> Spin up model if required (only active if running as SCAPE ENGINE).
	subroutine prepare(success_code)
    
        use setup_module
        use exceptions
        use global_data
        use local_data
        use xshore_dist_t_data
        use dump_module
        use utilities_module
        use fromint_module
        use expand_sf_module
        use SLR_module
        use timestep_module
        use vertical_grid
        implicit none
        
        integer, allocatable, dimension(:) :: temp_lut
		integer, intent(out) :: success_code ! 0: ok, -ive fatal error, +ive warning

        integer :: iostat
        integer :: outputYearCount
        integer :: outputSectionCount
        real(double) :: talus_dy
        integer :: seawall_element ! index of cliffy element that defines the position of a seawall
        integer :: v_sect
        
        success_code = 0

        logFileUnit = get_file_unit()
        saveSedFluxLeftUnit = get_file_unit()
        savePotSedFluxLeftUnit = get_file_unit()
        saveSedFluxRightUnit = get_file_unit()
        savePotSedFluxRightUnit = get_file_unit()
        saveRockContourUnit = get_file_unit()
        saveShoreContourUnit = get_file_unit()
        saveAnnBvolUnit = get_file_unit()
        saveFinesVolUnit = get_file_unit()
        saveSedTransAnnUnit = get_file_unit()
        savePotTransAnnUnit = get_file_unit()
        saveBeachAddAnnUnit = get_file_unit()
        saveVellLevelsUnit = get_file_unit()
        saveSeawallActiveUnit = get_file_unit()
        saveGroyneEffectUnit = get_file_unit()
        saveRockProfilesUnit = get_file_unit()
        saveBeachProfilesUnit = get_file_unit()
        
        wavesAndTidesFileUnit = get_file_unit()
        tempInputFileUnit = get_file_unit()

    call INITIALISE_PARAMETERS(FILE_PATH, radcon, success_code)
    if (success_code < 0) return

    success_code = do_open_text_output_file(logFileUnit, FILE_PATH, 'log')
    if (success_code < 0) return
        
    ! To set the year at the start of the model run
    year = startYear
    nextOutputYear = startYear
    
    call GET_PERIOD_FACTOR(year, periodFactor, havePeriodFactor)
    call GET_HEIGHT_FACTOR(year, heightFactor, haveHeightFactor)
    call GET_ANGLE_CHANGE(year, angleChange, haveAngleChange)

    ! To determine the model run duration, in tides.
    run_duration =  tidesPerYear * (endYear - startYear + 1)
            
    logString = '*** SCAPE Version '//trim(SCAPE_VERSION)
    call SCAPE_LOG(.true.)
    logString = ''
    call SCAPE_LOG(.true.)
    
    logString = 'Timeframe:'
    call SCAPE_LOG(.true.)
    
    if (highResOutStartYear > year) then
        write(logString, fmt='(A,I6,A,I6,A)') 'From year', startYear, ' output every', lowResOutTimestep, ' years'
        call SCAPE_LOG(.true.)
        write(logString, fmt='(A,I6,A,I6,A)') 'From year', highResOutStartYear, ' output every', highResOutTimestep, ' years'
    else
        write(logString, fmt='(A,I6,A,I6,A)') 'From year', startYear, ' output every', highResOutTimestep, ' years'
    end if
        
    call SCAPE_LOG(.true.)
    write(logString, fmt='(A,I6)') 'Until year', endYear
    call SCAPE_LOG(.true.)
    write(logString, fmt='(A,I8,A)') 'Run duration:', run_duration, ' tides'
    call SCAPE_LOG(.true.)
    logString = ''
    call SCAPE_LOG(.true.)

    wave_count = waveTideFileStart

    !******************
    !     Set the grid
    !******************
    !*** Number of Q points, must be greater than one!
    nQpoints = nSections + 1
    nYsections = nSections

    !*** The number of cliff elements
    nce = nint(sectionHeight / dz)

    ! The height of the beach in cliff elements
    beachHeightE = nint(beachHeight / dz)

    logString = 'Grid Information:'
    call SCAPE_LOG()
    write(logString, fmt='(A,I6)') 'Numer of Q points:', nQpoints
    call SCAPE_LOG()
    write(logString, fmt='(A,I6)') 'Number of Y sections:', nYsections
    call SCAPE_LOG()
    write(logString, fmt='(A,G20.8)') 'Distance between Q points:', dx
    call SCAPE_LOG()
    logString = ''
    call SCAPE_LOG()
    
    call LOG_PARAMETERS()
    logString = ''
    call SCAPE_LOG()

    logString = 'Output Metadata:'
    call SCAPE_LOG()
    logString = ''
    call SCAPE_LOG()
    
    ! Write the count and list of output years to the log file
    outputYearCount = 0
    nextOutputYear = startYear

    do year = startYear, endYear
        if (year == nextOutputYear) then
            outputYearCount = outputYearCount + 1
            call SET_OUTPUT_YEAR(year)
        end if
    end do
    
    write(logString, fmt='(A,I6,A)') 'Outputs are for', outputYearCount, ' years, which are:'
    call SCAPE_LOG()
    nextOutputYear = startYear
    
    do year = startYear, endYear
        if (year == nextOutputYear) then
            write(logString, *) year
            call SCAPE_LOG()
            call SET_OUTPUT_YEAR(year)
        end if
    end do
    
    logString = ''
    call SCAPE_LOG()
    
    !*****************************************
    ! Open files to save data to
    success_code = do_open_text_output_file( &
        saveSedFluxLeftUnit, FILE_PATH, 'SED_FLUX_L_BOUND', tidesPerYear, 'm3 per tide')
    if (success_code < 0) return

	success_code = do_open_text_output_file( &
        savePotSedFluxLeftUnit, FILE_PATH, 'SED_FLUX_POT_L_BOUND', tidesPerYear, 'm3 per tide')
    if (success_code < 0) return

    success_code = do_open_text_output_file( &
        saveSedFluxRightUnit, FILE_PATH, 'SED_FLUX_R_BOUND', tidesPerYear, 'm3 per tide')
    if (success_code < 0) return

	success_code = do_open_text_output_file( &
        savePotSedFluxRightUnit, FILE_PATH, 'SED_FLUX_POT_R_BOUND', tidesPerYear, 'm3 per tide')
    if (success_code < 0) return

    success_code = do_open_text_output_file(saveRockContourUnit, FILE_PATH, 'CONTOUR_ROCK', nYsections, 'm')
    if (success_code < 0) return
    
    success_code = do_open_text_output_file(saveShoreContourUnit, FILE_PATH, 'CONTOUR_SHORE', nYsections, 'm')
    if (success_code < 0) return
    
    success_code = do_open_text_output_file(saveAnnBvolUnit, FILE_PATH, 'VOLUME_B_ANN_AVE', nYsections, 'm3')
    if (success_code < 0) return
    
    success_code = do_open_text_output_file(saveFinesVolUnit, FILE_PATH, 'LIBERATED_SED_FINES', 1, 'm3 per year')
    if (success_code < 0) return
    
    success_code = do_open_text_output_file(saveSedTransAnnUnit, FILE_PATH, 'SED_TRANS_ANN', nQpoints, 'm3 per year')
    if (success_code < 0) return
    
    success_code = do_open_text_output_file(savePotTransAnnUnit, FILE_PATH, 'SED_TRANS_POT_ANN', nQpoints, 'm3 per year')
    if (success_code < 0) return
    
    success_code = do_open_text_output_file(saveBeachAddAnnUnit, FILE_PATH, 'LIBERATED_SED_BEACH', nYsections, 'm3 per year')
    if (success_code < 0) return
    
    success_code = do_open_text_output_file(saveSeawallActiveUnit, FILE_PATH, 'SEAWALL_PRESENCE', nYsections)
    if (success_code < 0) return
    
    success_code = do_open_text_output_file(saveGroyneEffectUnit, FILE_PATH, 'GROYNE_FACTORS', nQpoints)
    if (success_code < 0) return
    
    logString = ''
    call SCAPE_LOG()

    if (firstYearProfileOutput /= -1000000) then
        ! Write the count and list of profile output years to the log file
        outputYearCount = lastYearProfileOutput - firstYearProfileOutput + 1
        outputSectionCount = lastSectionProfileOutput - firstSectionProfileOutput + 1
        write(logString, fmt='(A,I6,A,I6,A,I6)') 'Profile outputs are every', profileOutputTimestep, &
            ' year(s), from', firstYearProfileOutput, ' to',lastYearProfileOutput
        call SCAPE_LOG()
        
        logString = ''
        call SCAPE_LOG()
        write(logString, fmt='(A,I6,A,I6,A,I6)') 'Profile outputs are for', outputSectionCount, &
            ' sections, from ', firstSectionProfileOutput, ' to ',lastSectionProfileOutput
        call SCAPE_LOG()
        
        logString = ''
        call SCAPE_LOG()
        success_code = do_open_text_output_file(saveRockProfilesUnit, FILE_PATH, 'ROCK_PROFILES')
        if (success_code < 0) return
    
        success_code = do_open_text_output_file(saveBeachProfilesUnit, FILE_PATH, 'BEACH_PROFILES')
        if (success_code < 0) return
    
        call output_profiles_data(saveRockProfilesUnit, nce, success_code)
        call output_profiles_data(saveBeachProfilesUnit, beachHeightE, success_code)
        
        logString = ''
        call SCAPE_LOG()
    end if
    !*****************************************

    ! To set the year at the start of the model run
    year = startYear
    nextOutputYear = startYear
   
      tstep_secs = 44890    ! time step of one tide

!*** Some constants
    breakrat = 0.78D0        ! Ratio of breaker height and water depth 
    rhow = 1035.0D0            ! Density of water
    g = 9.81D0            

!*** For the call to WAVE_TRANSFORMATION
    ngroynes = 0

!***********************************
!
!     Set the number of wave points
!
!***********************************
    
    ! 'MWPM' is used as a marker meaning 'Multiple Wave Point Modifications' 
    ! it indicates where code has been changed to accommodate multiple wave points.
    ! The wave points are numbered, 1 to nWavePoints.
    ! The lowest number is the one closest to the lowest section number.

    ! Allocate arrays
    allocate(hwp(nWavePoints), hwp_temp(nQpoints), anglewp(nWavePoints), anglewp_temp(nQpoints), &
        wpr(nQpoints), ann_av_bvol(nYsections), tidal_amp(nYsections), &
        stat = alloc_error)  ! MD 300304
      if(alloc_error /= 0)then
        write(6,*)'*** Error : could not allocate space'
        write(6,*)'*** Error : for whp, anglewp, wpr'
        stop
      end if                    
    
    hwp = 0.0D0
    hwp_temp = 0.0D0        ! MD 300304
    anglewp_temp = 0.0D0    ! MD 300304
    anglewp = 0.0D0
    
    ! MWPM associate a wave point number to each profile section

    write(6,*)'Assigning wave point numbers '

    do wp = 1, nWavePoints
        if (wp < nWavePoints) then
            lastWPsection = fswp(wp + 1) - 1
        else
            lastWPsection = nQpoints
        end if
            
        do i = fswp(wp), lastWPsection
            wpr(i) = wp
        end do
    end do
    
!**********************************************************************
!
!    Set up the initial platform / cliff shape
!    and the offshore physical boundary conditions & sea-level rise
!
!**********************************************************************

!*** Allocate & initialise arrays

        allocate(cliffy(nce,nYsections),          &
                 beachy(nce,nYsections),          &
                 talusy(nce,nYsections),           &
                 upbeachy(nce,nYsections),           &
                 lowbeachy(nce,nYsections),           &
                 cliffheights(nce),           &
                 sect_erode(nce,nYsections),               &
                 sect_erode_new(nce,nYsections),         &
                 surface(nce),                         &
                 top_erode(nYsections),bot_erode(nYsections),   &
                 startprofile(nce),           &
                 startprofiles(nce, lastActiveSection),     &
                 sect_height(nYsections),np_cliff(nYsections),    &
                 ubvolume(nYsections), lowbvolume(nYsections),    &
                 stormflag(nYsections),angleosc(nQpoints),    &
                 initial_cliff_toe(nYsections), &
                 retreat(nYsections, startYear: endYear), &
                 stat = alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Error : could not allocate space'
          write(6,*)'*** Error : for cliffy, beachy, & cliffheights'
          stop
        end if

        allocate(setup(nYsections),setdown(nYsections),h(nYsections),       &
                 cgroup(nYsections),anglewob(nYsections),theta(nYsections),               &
                 setup_qpoint(nQpoints),setdown_qpoint(nQpoints),h_qpoint(nQpoints),       &
                 cgroup_qpoint(nQpoints),anglewob_qpoint(nQpoints),               &
                 theta_qpoint(nQpoints), &
                 top_drift(nYsections),bot_drift(nYsections),           &
                 sect_drift(nce, nYsections),              &
                 this_tide_drift(nce),                &
                 dbreak(nYsections),half_tide(nYsections),              &
                 beach_addition(nYsections),fines_addition(nYsections),total_beach_addition(nYsections), &
                 sediment_influx(nYsections), &
                 stat = alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Error : could not allocate space'
          write(6,*)'*** Error : for setup, setdown'
          stop
        end if                         

! MW to do: Produce a definition sketch of the zone influenced by wave action and the perched cliff.
! Define the model clearly in terms of the dynamically reshaping part, with the cliff being essentially outside the model.

    cliffy=0.0D0
    beachy=0.0D0
    talusy=0.0D0
    sect_erode=0.0D0
    sect_erode_new=0.0D0
    surface=0.0D0
    top_erode=0
    bot_erode=0
    startprofile=0.0D0
    startprofiles=0.0D0
    sect_height=0.0D0
    np_cliff=0
    setup=0.0D0
    setdown=0.0D0
    h=0.0D0
    cgroup=0.0D0
    anglewob=0.0D0
    theta=0.0D0
    setup_qpoint=0.0D0
    setdown_qpoint=0.0D0
    h_qpoint=0.0D0
    cgroup_qpoint=0.0D0
    anglewob_qpoint=0.0D0
    theta_qpoint=0.0D0
    top_drift=0
    bot_drift=0
    sect_drift=0.0D0
    this_tide_drift=0.0D0
    dbreak=0.0D0
    half_tide=0.0D0
    beach_addition=0.0D0
    total_beach_addition=0.0D0
    fines_addition=0.0D0
    fines_volume=0.0D0
    sediment_influx=0.0D0

!*** End array initialize
                             
    ! initial value for one year
    call SL_for_year(year, msl_OD, success_code)
    if (success_code < 0) return
    
    call SL_for_year(year + 1, sea_level_next_year, success_code)
    if (success_code < 0) return

    slr_per_tide = (sea_level_next_year - msl_OD) / tidesPerYear

!*** Produce an array containing all the heights of the elements
    do i=1, nce 
        cliffheights(i) = i * dz
    end do

! Each section is associated with a single offshore contour angle, which does not account for wave direction.
! This is very simplistic and should be improved.
!*** angle of offshore contour, converted to radians
    success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, oscAnglesFile)
    if (success_code < 0) return
        
    read(tempInputFileUnit, *, iostat=iostat)angleosc
    success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading angleosc')
    close(tempInputFileUnit)
    if (success_code < 0) return
        
    angleosc = (angleosc + 0.001) * radcon
    call dump_real("run_scape517", "angleosc", .true., angleosc)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!    Set the properties of the platform and beach
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !*** Allocate arrays
    allocate(resistance(nce, nYsections),stat = alloc_error)
     if(alloc_error /= 0)then
       write(6,*)'*** Error : could not allocate space'
       write(6,*)'*** Error : for Resistance'
       stop
     end if   
    
    if (useVAgrid) then
        resistance = -1.0D0 ! Flag resistance unknown
        call get_layers(success_code)
        if (success_code < 0) return
            
        msl_m = msl_OD - baseLine_OD
            
    else ! use rockStrength data
        msl_m = mslM            ! (m) above the base of the model
        baseLine_OD = msl_OD - msl_m
        
        !*** The 'resistance' of the clay as a vector
        success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, rockStrengthFile)
        if (success_code < 0) return
            
        do i=1,nce
            read(tempInputFileUnit, *, iostat=iostat) resistance(i, :)
            success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading resistance')
            if (success_code < 0) exit
        end do 
        close(tempInputFileUnit)
        if (success_code < 0) return
    
    end if   

    ! Establish the position of the mean sea level (msl) relative to the model 
    !    and to Ordinance Datum
    msl_e = int(msl_m / dz) + 1    ! Number of the platform element at msl

    if (trim(outputContourLevelUser) == '') then
        outputContourLevel_e = nce
    else
        outputContourLevel_e = msl_e + int(outputContourLevel / dz) + 1
    end if

!*** Characteristic block size
! This is an inelegant approximation.
! It represents the largest depth that can be eroded in one tide

    block_size = maxBlockSize !MW030304  !0.15    ! (m)
    
    if (nQpoints > 1)then    ! There is a beach so

        !*************** For the general beach *************** 
        allocate(cerck1(nQpoints), groyneEffectAnnualSum(nQpoints), seawall_os(nYsections), seawall_active(nYsections), &
            QAnn(nQpoints), QpAnn(nQpoints), left_exclusion_angle(nQpoints), right_exclusion_angle(nQpoints), &
            stat = alloc_error)
        success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'groyneEffectAnnualSum, seawall_os')
        if (success_code < 0) return

        if (includeGroynes) then
            allocate(groyne_construction(nQpoints), groyne_removal(nQpoints), stat = alloc_error)
            success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'groyne_construction, groyne_removal')
            if (success_code < 0) return
        end if

        if (includeSeawall) then
            allocate(seawall_construction(nYsections), seawall_removal(nYsections), stat = alloc_error)
            success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'seawall_construction, seawall_removal')
            if (success_code < 0) return
        end if

    QAnn = 0.0D0
    QpAnn = 0.0D0
    groyneEffectAnnualSum = 0.0D0

!**********************************************
!***    To describe the shape of the beach    ***
!**********************************************

      a = bruunConst        ! Bruun or Dean constant, sets the angle of the 
      ! curve of the beach surface.
      beachCrestLevelE = nint(beachCrestLevel / dz) ! beachCrestLevelE is the Storm Surge level, this is used to 
      ! set the level of the berm, or top of beach. The number here is in metres above 
      ! the very base of the model.  Returns the level in cliff elements

!*** During storms the beach is assumed to 'split' into an upper and lower parts
! The upper part is more steep than the average slope, the lower is more gentle.
! This helps the waves penetrate the beach to erode the platform. Material can 
! move from the upper to the lower beach.

!*** For upper beach 
      upBeachHeightE = beachHeightE ! The number is in metres, returns the 
      ! height in cliff elements
      a_up = a + 0.05D0    

!*** For lower beach
      lowBeachHeightE = beachHeightE 
      a_low = a - 0.05D0

    if (.not. useAltSedimentTransport) then
        !*** Input the CERC coefficients, constants for the sediment transport equation
        success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, cercConstantsFile)
        if (success_code < 0) return
            
        read(tempInputFileUnit, *, iostat=iostat) cerck1  ! NB, values needed on Q points
        success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading cerck1')
        close(tempInputFileUnit)
        if (success_code < 0) return
    end if

!*** Allocate arrays

          allocate(beachelev(nYsections), berm(nYsections), bTop(nYsections), bBottom(nYsections), &
              ubTop(nYsections), ubBottom(nYsections), lbTop(nYsections), lbBottom(nYsections), stat = alloc_error)
          if(alloc_error /= 0)then
            write(6,*)'*** Error : could not allocate space'
            write(6,*)'*** Error : for beachelev, berm'
            stop
          end if        

       beachelev = beachCrestLevelE
      berm = 0.0D0

    if (runMode == restart) then
        success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, bermWidthsFile)
        if (success_code < 0) return
            
        read(tempInputFileUnit, *, iostat=iostat) berm
        success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading berm')
        close(tempInputFileUnit)
        if (success_code < 0) return
        
    end if
     
    call dump_real("run_scape651", "rebermsistance", .false., berm)
      
!**************************************************************
!***    Input files to control engineering construction        ***
!**************************************************************

    seawall_active = 0
    seawall_os = 0D0
        
    if (includeSeawall) then
        !*** Input dates at which the seawalls are installed
        if (seawallMode == defineSeawallsByYear) then
            call read_construction_data(seawall_construction, nYsections, &
                tempInputFileUnit, seawallConstructionYearsFile, success_code)
        else ! seawallMode == defineSeawallsByYearAndPosition
            call read_construction_data(seawall_construction, nYsections, &
                tempInputFileUnit, seawallConstructionYearsFile, success_code, seawall_os)
        end if
        
        if (success_code < 0) return
    
        !*** Input dates at which the seawalls are removed
        call read_construction_data(seawall_removal, nYsections, tempInputFileUnit, seawallRemovalYearsFile, success_code)
        if (success_code < 0) return
    end if
            
    if (includeGroynes) then
        !*** Input dates at which the groynes are constructed
        call read_construction_data(groyne_construction, nQpoints, tempInputFileUnit, groyneConstructionYearsFile, success_code)
        if (success_code < 0) return
            
        !*** Input the file containing the dates at which the groynes are removed
        call read_construction_data(groyne_removal, nQpoints, tempInputFileUnit, groyneRemovalYearsFile, success_code)
        if (success_code < 0) return
    end if

    left_exclusion_angle = wave_exclusion_angle_L_radians
    right_exclusion_angle = wave_exclusion_angle_R_radians

    if (excludeWaves == 2) then
        success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, waveExclusionAnglesFile)
        if (success_code < 0) return
            
        read(tempInputFileUnit, *, iostat=iostat) (left_exclusion_angle(i), right_exclusion_angle(i), i = 1, nQpoints)
        success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading wave exclusion angles')
        close(tempInputFileUnit)
        if (success_code < 0) return
        
        left_exclusion_angle = radcon * left_exclusion_angle
        right_exclusion_angle = radcon * right_exclusion_angle

    end if
     
    write(6,*)'Structure control files loaded'
      
    end if    ! End of condition of (nQpoints > 1)

!*** Load in the look-up-tables for 
! wave celerity (clut), and
! wave group celerity (cglut).
!
! N.B. these might have to be increased if water becomes deeper or longer wave periods occur

    allocate(clut(periodValues, depthValues), cglut(periodValues, depthValues), stat = alloc_error)
    if(alloc_error /= 0)then
      write(6,*)'*** Error : could not allocate space'
      write(6,*)'*** Error : for clut, cglut'
      stop
    end if        

    clutUnit = get_file_unit() 
    cglutUnit = get_file_unit()

    if (celerityTablesAsTxtFiles) then
        success_code = do_open_text_input_file(clutUnit, FILE_PATH, 'clut.txt')
        if (success_code < 0) return
        success_code = do_open_text_input_file(cglutUnit, FILE_PATH, 'cglut.txt')
        if (success_code < 0) return
    
        do i = 1, depthValues
            read(clutUnit, *, iostat=iostat) clut(:, i)
            success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading clut')
            if (success_code < 0) return
    
            read(cglutUnit, *, iostat=iostat) cglut(:, i)
            success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading cglut')
            if (success_code < 0) return
        end do
    else
        allocate(temp_lut(periodValues), stat = alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Error : could not allocate space'
          write(6,*)'*** Error : for temp_lut'
          stop
        end if        
    
        success_code = do_open_binary_input_file(clutUnit, FILE_PATH, 'clut.bin')
        if (success_code < 0) return
        success_code = do_open_binary_input_file(cglutUnit, FILE_PATH, 'cglut.bin')
        if (success_code < 0) return
    
        do i = 1, depthValues
            read(clutUnit) temp_lut
    
            do j = 1, periodValues
                clut(j, i) = fromInt(temp_lut(j))
            end do
    
            read(cglutUnit)temp_lut
    
            do j = 1, periodValues
                cglut(j, i) = fromInt(temp_lut(j))
            end do
        end do

        deallocate(temp_lut, stat=alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Unexpected deallocation error'
          write(6,*)'*** for temp_lut'
        end if   
    end if

    close(clutUnit)
    close(cglutUnit)

    write(6,*)'look-up tables loaded'
      
!**********************************************
!    Set cliff y-positions & the offsets        ***
!**********************************************

    if (nQpoints == 1)then ! This might be the case when evolving from an initial plunging cliff 
        success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, startprofileFile) ! Only one section
        if (success_code < 0) return
            
        read(tempInputFileUnit, *, iostat=iostat) startprofile
        success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading startprofile')
        close(tempInputFileUnit)
        if (success_code < 0) return
        
    else if (runMode == restart .and. .not. useVAgrid) then
        ! MD START CONDITIONS
        success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, startProfilesFile)
        if (success_code < 0) return
            
        do i=1,nce
            read(tempInputFileUnit, *, iostat=iostat) startprofiles(i, :)
            success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading startprofiles')
            if (success_code < 0) exit
        end do 
        close(tempInputFileUnit)
        if (success_code < 0) return
    end if
    
    call dump_cliff_real("run_scape773", "startprofiles", startprofiles)
    
    if (nQpoints == 1)then
        cliffy(:,1) = startprofile
    else

!*** Allocate arrays

        allocate(section_offsets(nYsections), &
             beachos(nYsections),sand_fraction(nYsections),    &
             tvolume(nYsections), cliffTopElevations(nYsections), &
             cliff_heights(nYsections),clifftoe_pos(nYsections),    &
              talus_offsets(nYsections),upbeach_os(nYsections),lowbeach_os(nYsections),    &
              upbvolume(nYsections),                 stat = alloc_error)
        success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'section_offsets, beachos, sand_fraction')
        if (success_code < 0) return

        section_offsets=0.0d0
        beachos=0.0d0
     
        if (useVAgrid) then
            call get_y_values_from_profiles(nYsections, nce, cliffy)
        else
            ! START CONDITIONS
            success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, profOffsetsFile)
            if (success_code < 0) return
               
            read(tempInputFileUnit, *, iostat=iostat) section_offsets
            success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading section_offsets')
            close(tempInputFileUnit)
            if (success_code < 0) return
        
            call dump_real("run_sacpe806", "section_offsets", .false., section_offsets)
    
            do section = firstActiveSection, lastActiveSection
              cliffy(:,section) = startprofiles(:, section) + section_offsets(section) 
            end do
            do section = 1, firstActiveSection - 1
              cliffy(:,section) = startprofiles(:, firstActiveSection) + section_offsets(section)
            end do
            do section = lastActiveSection + 1, nYsections
              cliffy(:,section) = startprofiles(:, lastActiveSection) + section_offsets(section)
            end do
            do section = 1, nYsections
              beachos(section) = 0.0D0
            end do
        end if
    end if
    
    call dump_cliff_real("run_scape823", "cliffy", cliffy)
    call dump_real("run_sacpe824", "beachos", .false., beachos)

            ! The cliff toe position is saved once per year.
            write(saveRockContourUnit,"(F12.4)") cliffy(outputContourLevel_e, :)

            do section = 1, nYsections
                write(saveShoreContourUnit,"(F12.4)") &
                    max(cliffy(outputContourLevel_e, section), beachy(outputContourLevel_e, section))
            end do

!***************************************************
!     Parameters to describe the beach position
!***************************************************

!*** Allocate arrays
    
          allocate(bruun(beachHeightE),upbruun(upBeachHeightE), &
          lowbruun(lowBeachHeightE),bvolume(nYsections),obvolume(nYsections),    &
          beachdepth(nce),beachdepth1(nce,nYsections),        &
           offshore_st(nYsections),    &
                   stat = alloc_error)
          if(alloc_error /= 0)then
            write(6,*)'*** Error : could not allocate space'
            write(6,*)'*** Error : for bruun'
            stop
          end if   

      beachdepth1=0
      bvolume=0

!*************** Calc the Bruun beach shapes *************** 

    do i=1,beachHeightE
        bruun(i)=(dble(i)*dz/a)**1.5
    end do

    do i=1,upBeachHeightE
        upbruun(i)=(dble(i)*dz/a_up)**1.5
    end do

    do i=1,lowBeachHeightE
        lowbruun(i)=(dble(i)*dz/a_low)**1.5
    end do

    ! The range of the upper and lower beaches are
        bTop = beachelev
        bBottom = beachelev - beachHeightE + 1
        ubTop = beachelev
        ubBottom = beachelev - upBeachHeightE + 1
        lbTop = beachelev
        lbBottom = beachelev - lowBeachHeightE + 1

      do section = firstActiveSection, lastActiveSection

        beachy(bTop(section): bBottom(section): -1, section) =     &
            beachos(section) + berm(section) + bruun

        where (beachy(:,section) > cliffy(:,section))
            beachdepth1(:, section) = beachy(:,section) - cliffy(:,section)
        elsewhere
            beachdepth1(:, section) = 0
        end where

        bvolume(section) = sum(beachdepth1(:,section))*dz * dx
      end do


      if (runMode == restart) then
          ! Start Conditions, load the initial beach volume
          success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, beachVolumesFile)
          if (success_code < 0) return

          read(tempInputFileUnit, *, iostat=iostat) bvolume
         success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading bvolume')
         close(tempInputFileUnit)
         if (success_code < 0) return
        
          do section = 1, firstActiveSection - 1
            beachy(:,section) = beachy(:,firstActiveSection)+section_offsets(section)
            bvolume(section) = bvolume(firstActiveSection)
          end do
    
          do section = lastActiveSection + 1, nYsections
            beachy(:,section) = beachy(:,lastActiveSection)+section_offsets(section)
            bvolume(section) = bvolume(lastActiveSection)
          end do
      end if
      
     call dump_cliff_real("run_scape919", "beachy", beachy)
     call dump_cliff_real("run_scape920", "beachdepth1", beachdepth1)
     call dump_real("run_scape921", "bvolume", .false., bvolume)

     if (year < startYear + beachZeroDelay) then
         bvolume = 0.0D0
     end if
            
!********** Define the initial offsets of the upper  & lower beach ********** 
! Close to the cliff toe
    upbeach_os = beachos + 5
    upbvolume = 0
    lowbeach_os = beachos - 5
    lowbvolume = 0

!********** The initial shapes of the upper and lower beaches ********** 
    upbeachy = 0
    lowbeachy = 0
    
    do section = 1, nYsections
        upbeachy(ubTop(section): ubBottom(section): -1, section) = upbeach_os(section) + upbruun
        lowbeachy(lbTop(section): lbBottom(section): -1, section) = lowbeach_os(section) + lowbruun
    end do
    call dump_cliff_real("run_scape940", "upbeachy", upbeachy)
    call dump_cliff_real("run_scape941", "lowbeachy", lowbeachy)

!**************************************************************
!    Shape function for cross shore distribution of erosion   
!**************************************************************

    call read_sf(tempInputFileUnit, FILE_PATH, 'shapeFunctionErosion.txt', np_sf_erode, sf_erode, success_code)
    if (success_code < 0) return

    dsf_erode = 0.05D0    ! The distance between points in the shape function
    ! This is as a fraction of the distance from the water line to the break point

    ! ALL SHAPE FUNCTIONS MUST HAVE THE SAME dX!!!!!
    
    if(nQpoints > 1)then ! There is a beach 

!**************************************************************
!    Shape function for cross shore distribution of drift 
!**************************************************************

    call read_sf(tempInputFileUnit, FILE_PATH, 'shapeFunctionDrift.txt', np_sf_drift, sf_drift, success_code)
    if (success_code < 0) return

    dsf_drift = 0.05D0    ! Distance between points in the shape function
      
    ! ALL SHAPE FUNCTIONS MUST HAVE THE SAME dX!!!!!

    ! calc drift offset for cache
    driftOffset = (np_sf_erode + cacheCeil) * (cacheCeil - 1)

! Allocate memory for the xshorecache (JT241007)
! Memory storage optimised (JT221107)
    allocate(xshorecache(driftOffset+((np_sf_drift+cacheCeil)*(cacheCeil-1))),    &
        xshorecache_flag(2,cacheCeil), &
        stat = alloc_error)
    if(alloc_error /= 0)then
            write(6,*)'*** Error : could not allocate space'
            write(6,*)'*** Error : for xshorecache'
        stop
    end if

    xshorecache_flag = .false.

      !*** Allocate arrays pulled from platform erosion (JT221107a)
        allocate(this_erode(nce),                 &
         m(nce),    &
         temp5(nce),bdepth(nce),    &
         battenuation(nce),            &
         crossbeach_st(nYsections), &
         seawall_effect(nce), &          
                 stat = alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Error : could not allocate space'
          write(6,*)'*** Error : subroutine PLATFORM_EROSION'
          stop
        end if

!*** Allocate arrays pulled from xshore_dist_t (JT221107a)
        max_np_sf = max(np_sf_erode,np_sf_drift)

        allocate(distmarkers(nce), xshore(max_np_sf + 2*allocCeil),    &
         tide_duration(2*allocCeil+1), temp(max_np_sf + 2*allocCeil),    &
         start_fnctn(max_np_sf + 2*allocCeil), end_fnctn(max_np_sf + 2*allocCeil),    &
         locatesf(max_np_sf), duration_sub(max_np_sf), stat = alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Error : could not allocate space'
          write(6,*)'*** Error : subroutine xshore_dist_t'
          stop
        end if

!'''''' MD 101003 Allocate arrays for beach width / groyne effect stuff

        allocate(beach_temp(nce), cliff_temp(nce), &
         stat = alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Error : could not allocate space'
          write(6,*)'*** Error : subroutine SEDIMENT_TRANSPORT'
          stop
        end if

!**********************************
!***    Cliff and talus data    ***
!**********************************

    if (.not. useVAgrid) then
        !*** Load the sand fraction of the cliff sections 
        success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, beachContentFile)
        if (success_code < 0) return
            
        read(tempInputFileUnit, *, iostat=iostat) sand_fraction
        success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading sand_fraction')
        close(tempInputFileUnit)
        if (success_code < 0) return
    end if
            
    call dump_real("run_scape1179", "sand_fraction", .false., sand_fraction)

!*** Load the level of the cliff top 
    success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, cliffElevationsFile)  !MD080304 new boundaries
    if (success_code < 0) return
        
    read(tempInputFileUnit, *, iostat=iostat) cliffTopElevations
    success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading cliffTopElevations')
    close(tempInputFileUnit)
    if (success_code < 0) return
        
    call dump_real("run_scape1179", "cliffTopElevations", .false., cliffTopElevations)

!*** Work out the cliff heights 
! NB Only works if the cliff heights are not varied within the model, as is
! possible in the diffraction version, i.e. if sectHeight = min_sect_height.
    do section = 1, nYsections
        cliff_heights(section) = cliffTopElevations(section) - (sectionHeight - msl_m + msl_OD)
        
    ! '''''''''''''''''
    ! cliffTopElevations must not be below the highest point in cliffy
    
        if (cliff_heights(section) < 1.0D-6) then
            cliff_heights(section) = 0.0D0
        end if
            
        if (cliff_heights(section) < 0.0D0) then
            write(logString,*) 'cliffTopElevations below the highest point in cliffy at section ', section
            call RAISE_EXCEPTION(logString)
            success_code = -1
            return
        end if
    ! '''''''''''''''''
        
    end do
    call dump_real("run_scape1209", "cliff_heights", .false., cliff_heights)
    call dump_real("run_scape1210", "sand_fraction", .false., sand_fraction)

!*** Find the initial position of the cliff toe 
    clifftoe_pos = cliffy(nce,:)
    call dump_real("run_scape1216", "clifftoe_pos", .false., clifftoe_pos)

!*** Define the initial volumes of Talus 
    if (runMode == restart) then
        ! MD START CONDITIONS
        success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, talusVolumeFile)
        
        read(tempInputFileUnit, *, iostat=iostat) tvolume
        success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading tvolume')
        close(tempInputFileUnit)
        if (success_code < 0) return
        
    else
        tvolume = 0.0D0
    end if
    
    call dump_real("run_scape1225", "tvolume", .false., tvolume)

!*** Define the initial Talus offsets 
! A little way behind the cliff toe
    talus_offsets = clifftoe_pos - 2
    call dump_real("run_scape1232", "talus_offsets", .false., talus_offsets)

!*** Define the initial Talus vector
    talus_dy = dz * tan(talusSlope_radians)
    do section = 1, nYsections
        talusy(nce,section) = talus_offsets(section)
        do j = nce - 1, 1, -1
            talusy(j,section) = talusy(j+1,section) + talus_dy
        end do
    end do
    call dump_cliff_real("run_scape1242", "talus", talusy)


!********************************************************
! Information regarding the exchange of beach material 
! with the offshore bar
!********************************************************
    ! Offshore bar volume
    if (runMode == restart) then
        ! MD START CONDITIONS    
        success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, offshoreBarVolFile)
        if (success_code < 0) return
        
        read(tempInputFileUnit, *, iostat=iostat) obvolume
        success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading obvolume')
        close(tempInputFileUnit)
        if (success_code < 0) return
        
    else
        obvolume = 0.0D0
    end if
    
    call dump_real("run_scape1256", "obvolume", .false., obvolume)

    totalbvolume = sum(bvolume) !+ sum(obvolume)
    end if

!**************************************
!*** Preparations for the main loop ***
!**************************************

    waveleader1 = 0
    waveleader2 = 0
    waveleader3 = 0
    waveleader4 = 0
    tidecounter = 0
    eventcounter = 0
    sea_rise = mslOffset    
              
!*** Open tide and wave files
        
    success_code = do_open_text_input_file(wavesAndTidesFileUnit, FILE_PATH, wavesAndTidesFile)
    if (success_code < 0) return
          
    outputYearFlag = .false.

    if (runMode == restart) then ! Insert structures before the model starts running
        if (includeSeawall) then
            ! Seawalls and revetments on the Y sections
            if (seawallMode == defineSeawallsByYear) then
                seawall_element = seawallBaseLevel / dz
            end if
            
            where (seawall_construction > 0 .and. year >= seawall_construction .and. &
                (seawall_removal == 0 .or. year < seawall_removal))
                seawall_active = 1
            end where
                    
            if (seawallMode == defineSeawallsByYear) then
                where (seawall_active == 1)
                    seawall_os = cliffy(seawall_element, :)
                end where
            end if
        end if

        !Save info on where seawalls are
        write(saveSeawallActiveUnit,'(I1)') seawall_active
    end if ! runMode == restart

    ! Parameters for the VELLINGA profiles
    if (vellingaActive) then
        allocate(vellinga_sections(nSections), stat = alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Error : could not allocate space'
          write(6,*)'*** Error : for vellinga_sections'
          stop
        end if    
        success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, vellingaSectionsFile)
        nvsect = 0
        vellinga_sections = 0
        
        do
            read(tempInputFileUnit, *, iostat=iostat) v_sect
            if (iostat == -1) exit ! end of file
            success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading vellinga_sections')
            if (success_code < 0) exit

            if (v_sect > 0 .and. v_sect <= nSections) then
                nvsect = nvsect + 1
                vellinga_sections(nvsect) = v_sect
            else
                write(logString,*) 'Vellinga section number out of range', v_sect
                call RAISE_EXCEPTION(logString)
                success_code = -1
                return
            end if
        end do
        
        close(tempInputFileUnit)
        if (success_code < 0) return
        
        success_code = do_open_text_output_file(saveVellLevelsUnit, FILE_PATH, 'VELL_LEVELS', nvsect)
        if (success_code < 0) return
    
        !****   Allocate Vellinga arrays
        allocate(vxs(nce), vell_berm(nSections), vellinga(nce), v_ref(nvsect), v_level_odn(nvsect), stat = alloc_error)
              if(alloc_error /= 0)then
                write(6,*)'*** Error : could not allocate space'
                write(6,*)'*** Error : for Vellinga variables'
                stop
              end if    
        vell_berm = 0.0D0
        vellinga = 0.0D0
    end if ! vellingaActive


    ! capture the first cliff toe line
    initial_cliff_toe(:) = cliffy(nce, :);
    
    ann_av_bvol = 0.0D0

    time_start = to_mjd(julian_day(1, 1, startYear, 0, 0, 0))
    ntide = 1
    success_code = 0
        
    ! Take the model up to spinUpToYear.
    ! Important for FluidEarth/OpenMI that this is done during the prepare call.
    ! Makes no real difference for RUN_SCAPE.EXE runs.
    do while (year < spinUpToYear)
         call update(success_code)
         if (success_code < 0) return
    end do 

    end subroutine prepare
	

!> Output the output metadata required at the top of the 
!> ROCK_PROFILES.txt or BEACH_PROFILES.txt files.
subroutine output_profiles_data(unitNo, nElements, success_code)
    use setup_module, only: profileOutputTimestep, firstSectionProfileOutput, lastSectionProfileOutput, &
        firstYearProfileOutput, lastYearProfileOutput, dz, dx
    implicit none
    
    integer, parameter :: double = kind(1D0)
    
    integer, intent(in) :: unitNo !< FORTRAN unit number for already open output file
    integer, intent(in) :: nElements !< number of cliff elements to output
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

    write(unitNo,"(I6)") nElements
    write(unitNo,"(I6)") lastSectionProfileOutput - firstSectionProfileOutput + 1
    write(unitNo,"(I6)") lastYearProfileOutput - firstYearProfileOutput + 1
    write(unitNo,"(F12.4)") dz
    write(unitNo,"(F12.4)") dx
    write(unitNo,"(I6)") profileOutputTimestep
    write(unitNo,"(I6)") firstSectionProfileOutput
    write(unitNo,"(I6)") lastSectionProfileOutput
    write(unitNo,"(I6)") firstYearProfileOutput
    write(unitNo,"(I6)") lastYearProfileOutput

    success_code = 0

end subroutine  output_profiles_data

end module initialise_module