!**********************************************
!*       setup.f95       
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

module setup_module

implicit none
save
integer, parameter :: double = kind(1.0D0)

integer :: tidesPerYear = 706 !< Approximate number of tides per year
integer, parameter :: restart = 0 !< runMode value: restart run after reading initial conditions from files 
integer, parameter :: start = 1 !< runMode value: start run after setting default initial conditions
integer :: runMode = restart !< Controls whether run initial conditions are defaulted or read from files

integer :: startYear = -1000000 !< First (Gregorian) year simulated
integer :: endYear = -1000000 !< Last (Gregorian) year simulated
integer :: spinUpToYear = -1000000 !< In the FluidEarth/OpenMI context, the component time horizon starts at spinUpToYear.
    !< Default value -1000000 means use startYear.
    !< Has no effect on RUN_SCAPE.EXE runs.
integer :: beachZeroDelay = 0 !< Number of years that must elapse before the beach volume is allowed to become greater than zero.
integer :: lowResOutTimestep = 200 !< Number of years between model output before highResOutStartYear
    !< (intended for long simulations)
integer :: highResOutTimestep = 1 !< Number of years between model output after highResOutStartYear
    !< (intended for shorter simulations or the last part of long simulations)
integer :: highResOutStartYear = -1000000 !< Year (Gregorian) at which the resolution of the model output increases.
    !< Default value -1000000 means use startYear.
integer :: firstYearTransportOutput = -1000000 !< First year at which sediment transport is output.
    !< This may be used to avoid saving results from early in the simulation 
    !< Default value -1000000 means use startYear.
integer :: profileOutputTimestep = 1 !< Number of years between profile results outputs
integer :: firstYearProfileOutput = -1000000 !< First year at which profile results are output.
    !< Default value -1000000 means no profile output.
integer :: lastYearProfileOutput = -1000000 !< Last year at which profile results are output.
    !< Default value -1000000 means use endYear.
integer :: firstSectionProfileOutput = 1 !< First section for which profile results are output
integer :: lastSectionProfileOutput = -1 !< Last section for which profile results are output.
    !< Default value -1 means use nSections.
integer :: nSections = 10 !< Number of model cross sections 
real(kind=double) :: sectionWidth = 500.0D0 !< Model section width (distance between y points) (m)
real(kind=double) :: sectionHeight = 20.0 !< Model section height (m)
    !< sectionHeight is not the cliff height, but is the vertical extent of the material that may be 
    !< affected by wave action for a given sea level.
    !< This should be safely above the maximum level that wave action may reach.
real(kind=double) :: elementHeight = 0.05D0 !< Height of each model element (m)
    !< larger values of elementHeight will reduce run times, but will be less responsive to changes in sea level.
    
logical :: useVAgrid = .false. !< Use the vertical grid to determine strength and sediment content values
real(kind=double) :: elementWidth = 1.0D0 !< Width of each vertical grid element (m)
    !< larger values of elementWidth will reduce run times.
real(kind=double) :: YmaxH = -1.0D0 !< Width of model vertical grid (m)
integer :: nLayers = 10 !< Number of material layers over base goelogy
    !< nLayers limits the number of layers allowed during the run.
    !< The number of layers initially active is determined when reading the file layerData.txt.
    
real(kind=double) :: baselineAngle !< Model baseline angle (degrees)
logical, private :: baselineAngleRequired = .true.
integer :: firstActiveSection = 1 !< First active section
integer :: lastActiveSection = -1 !< Last active section. 
    !< Default value -1 mean use number of sections (nSections).
    !< firstActiveSection and lastActiveSection define the size of the model boundary regions.
    !< At model sections with numbers below firstActiveSection or above lastActiveSection erosion is extrapolated from further inside the model.

integer :: waveTideFileStart = 1 !< The point in the wave and tide input file at which data is read.
integer :: nWavePoints = -1 !< Number of wave points at which wave conditions are input.
    !< Default value -1 means use one wave point with a depth of 10.0 applied to all sections 
    !< The wave points will all be numbered, 1 to nWavePoints.
    !< The lowest number will be the one closest to the lowest section number.
integer, allocatable, dimension(:) :: fswp !< First section for each wave point
real(kind=double), allocatable, dimension(:) :: depthWPmsl !< Depth of each wave point below mean sea level (m)

integer, parameter :: defineSeawallsByYear = 0 !< seawallMode value: seawall position calculated from seawallBaseLevel
integer, parameter :: defineSeawallsByYearAndPosition = 1 !< seawallMode value: seawall position read from file
integer :: seawallMode = defineSeawallsByYear !< Controls whether seawall position is calculated or read from file
real(kind=double) :: seawallBaseLevel !< Height of the cliffy element used to define the position of the seawall
    !< in metres above mean sea level in the year that the seawall is installed (m)
    !< Required if includeSeawall is .true. and seawallMode is defineSeawallsByYear
logical, private :: seawallBaseLevelSupplied = .false.
real(kind=double) :: mslM = 15.0D0 !< Height of mean sea level above the base of the model (m)
    !< mslM is used only when useVAgrid is .false.
    !< When the vertically aligned grid is in use, the value is calculated from the highest seaward Z value in the base geology profiles
real(kind=double) :: mslOffset = 0.0D0 !< Model elements are shifted down(/up) to represent sea level rise(/fall).
    !< MSL is always at the element which happens to be msl_m / elementHeight from the base of the model.
    !< Every tide it rises up this element by a small amount.
    !< mslOffset records how far up the element it has risen.  (m)
real(kind=double) :: depthOSCMsl !< Depth of offshore contour below mean sea level (m)
logical, private :: depthOSCMslRequired = .true.
character(10) :: outputContourLevelUser = '' !< Height of the rock and shore contour relative to mean sea level (m).
    !< Default means set the level at the top element of cliffy
real(kind=double) :: minStormWaveHeight = 2.0D0 !< minimum storm wave height (m)
real(kind=double) :: minHs = 0.2D0 !< The wave height below which no erosion is assumed to occur (m)
real(kind=double) :: maxPlatformSlope = 70.0D0 !< maximum platform slope (degrees)
real(kind=double) :: bermSlope = -1.0D0 !< berm slope (1 in x).
    !< Default value -1 means no berm slope effects
real(kind=double) :: initialBermStep = 1.0D0 !< (m) Step size of the (horizontal) change in berm position.
    !< This will adapt in the model.
real(kind=double) :: minBermStep = 1.0D-6 !< (m) Minimum bermstep value.
real(kind=double) :: maxBermStep = 1.0D+6 !< (m) Maximum bermstep value.
real(kind=double) :: beachDisturbanceRatio = 0.2D0 !< beach disturbance (ratio)
    !< depth of the beach at which waves begin to erode / wave height
real(kind=double) :: beachReturnRatio = 0.01D0 !< beach return ratio
    !< beachReturnRatio is the proportion of the offshore bar that returns every tide, when conditions are calm.
    !< This is a simple behavioural rule, a candidate for improvement.
real(kind=double) :: talusStrengthRatio = 0.1D0 !< talus strength (ratio)
    !< strength of talus / strength of platform
real(kind=double) :: talusSlope = 45.0D0 !< angle of talus (degrees)
real(kind=double) :: maxBlockSize = 2.0D0 !< Maximum width of plarform material that could be removed in one tide (m)
    !< maxBlockSize applies across the model, but is rarely used by the model.
    !< When it is it tends to be at the foot of the cliff, where profile gradients are tight.
    !< Generally it represents the occasional loss of blocks of material from the toe of a cliff.
real(kind=double) :: beachHeight = 10.0D0 !< Height of the representation of the Bruun face of the beach (m)
    !< beachHeight extends down from the beachCrestLevel.
real(kind=double) :: bruunConst = 0.2D0 !< Bruun constant
real(kind=double) :: beachCrestLevel !< Elevation of the junction between the curving beach face and 
    !< the gently sloping berm on top of it (m). Relative to the base of the model.
logical, private :: beachCrestLevelRequired = .true.
real(kind=double) :: minBoundaryTransportRatio = 0.1D0 !< minimum boundary transport ratio
    !< minimum proportion of the potential alongshore sediment that may occur at the boundariues
real(kind=double) :: qpMaxBoundaryRightIn = 1.0D+9 !< Maximum sediment influx at the right model boundary
    !< in m^3 per tide
real(kind=double) :: qpMaxBoundaryLeftIn = 1.0D+9 !< Maximum sediment influx at the left model boundary
    !< in m^3 per tide
real(kind=double) :: qpMaxBoundaryRightOut = 1.0D+9 !< Maximum sediment outflux at the right model boundary
    !< in m^3 per tide
real(kind=double) :: qpMaxBoundaryLeftOut = 1.0D+9 !< Maximum sediment outflux at the left model boundary
    !< in m^3 per tide
integer :: excludeWaves = 0 !< Wave exclusion switch 
    !< Value 0 means not exclude
    !< Value 1 means exclude, same exclusion angles apply to all Q points
    !< Value 2 means exclude, exclusion angles for each Q point are read from waveExclusionAnglesFile
real(kind=double) :: waveExclusionAngleL = 22.5D0
    !< Angle anti-clockwise from shore normal, beyond which waves are excluded (degrees)
real(kind=double) :: waveExclusionAngleR = 22.5D0
    !< Angle clockwise from shore normal, beyond which waves are excluded (degrees)
integer :: switchTransportWhDiff = 0 !< Switch transport due to differential wave height 
    !< Value 1 means switch, 0 means do not switch
integer :: nCliffSimulations = 0 !< Number of cliff failure simulations 
integer :: slumpPeriod = 10 !< The number of erosion events that occurs before the cliff face slumps.
    !< The slumpPeriod parameter influences the frequency with which material is relased from the cliff.
    !< It does not influence the average rate at which material is released from the cliff.
    !< It is also disconnected from RegFalls3, the model used to calculate cliff top retreat.

logical :: useAltSedimentTransport = .false. !< Use sediment transport model described by French and Burningham.
    !< Taken from the paper: WAVE-DRIVEN SEDIMENT PATHWAYS ON A GRAVEL DOMINATED COAST
    !< SUBJECT TO A STRONGLY BI-MODAL WAVE CLIMATE, SUFFOLK, EASTERN UK.
    !< van Rijn, L.C. (2014). A simple general expression for longshore
    !< transport of sand, gravel and shingle. Coastal Engineering, 90, 23-39.
real(kind=double) :: medianGrainsize = 6.0D-3 !< Median sediment grain size (m).
    !< Used when useAltSedimentTransport is .true.

integer :: tidalAmpSection = 1 !< Y point (section) index at which the tidal_amp value read from the 
    !< Waves and Tides file is correct
real(kind=double) :: tidalAmpScaleLeft = 1.0D0 !< tidal_amp scale factor for left end of model looking offshore
real(kind=double) :: tidalAmpScaleRight = 1.0D0 !< tidal_amp scale factor for right end of model looking offshore

logical :: tidalRangeVariation = .true. !< Allow the tidal range to vary along the model
logical :: includeSeawall = .false. !< Include seawall effects
logical :: includeGroynes = .false. !< Include groyne effects
    
logical :: celerityTablesAsTxtFiles = .false. !< Read celerity tables from clut.txt and cglut.txt files
integer :: periodValues = 221 !< Celerity look up tables first dimension
    !< Wave period values start at 1.0 secs and increment 0.1 secs each column
integer :: depthValues = 5000 !< Celerity look up tables second dimension
    !< Water depth values start at 0.01m and increment 0.01m each row
    
logical :: vellingaActive = .false. !< Vellinga code active if .true.

real(kind=double), allocatable, dimension(:) :: cbstrWaveHeight !< Cross beach sediment transport lookup wave height
real(kind=double), allocatable, dimension(:,:) :: cbstrPeriod !< Cross beach sediment transport lookup period
real(kind=double), allocatable, dimension(:,:) :: cbstrValue !< Cross beach sediment transport lookup value
    
real(kind=double), allocatable, dimension(:) :: offstrWaveHeight !< Offshore sediment transport lookup wave height
real(kind=double), allocatable, dimension(:,:) :: offstrPeriod !< Offshore sediment transport lookup period
real(kind=double), allocatable, dimension(:,:) :: offstrValue !< Offshore sediment transport lookup value

real(kind=double), allocatable, dimension(:) :: groyneEffectBeachWidth !< Groyne effect lookup beach width
real(kind=double), allocatable, dimension(:) :: groyneEffectValue !< Groyne effect lookup value.
    !< Within the range zero (no effect) and 1 (groyne completely blocks drift)

integer, allocatable, dimension(:) :: modifyPeriodYear !< Year at which wave period factor applies.
    !< modifyPeriodYear values must increase monotonically
real(kind=double), allocatable, dimension(:) :: modifyPeriodFactor !< Wave period factor.
    !< modifyPeriodFactor values must be +ve

integer, allocatable, dimension(:) :: modifyHeightYear !< Year at which wave height factor applies.
    !< modifyHeightYear values must increase monotonically
real(kind=double), allocatable, dimension(:) :: modifyHeightFactor !< Wave height factor.
    !< modifyHeightFactor values must be +ve

integer, allocatable, dimension(:) :: modifyAngleYear !< Year at which wave angle factor applies.
    !< modifyAngleYear values must increase monotonically
real(kind=double), allocatable, dimension(:) :: modifyAngleAmount !< Wave angle change amount (degrees.
    !< modifyAngleAmount values may be +ve or -ve

! these values are derived
real(kind=double) :: maxPlatformSlope_radians !< maximum platform slope (radians) (derived value)
real(kind=double) :: wave_exclusion_angle_R_radians !< wave exclusion angle (radians) (derived value)
real(kind=double) :: wave_exclusion_angle_L_radians !< wave exclusion angle (radians) (derived value)
real(kind=double) :: baselineAngle_radians !< Model baseline angle (radians) (derived value)
real(kind=double) :: talusSlope_radians !< Angle of talus (radians) (derived value)
real(kind=double) :: outputContourLevel !< Height of the rock and shore contour relative to mean sea level (m) (derived value)
real(kind=double) :: dx !< sectionWidth paramaeter (derived value)
real(kind=double) :: dy !< elementWidth parameter (derived value)
real(kind=double) :: dz !< elementHeight parameter (derived value)

! Declare variables for input file names
character(64) :: oscAnglesFile = 'oscAngles.txt'
    !< File always used.
    !< Provides values for the angle of the offshore contour.
    !< Contents: A value for each beach section (1 : nSections).
character(64) :: cercConstantsFile = 'cercConstants.txt'
    !< Provides values for CERC coefficients, constants for the sediment transport equation.
    !< Contents: A value for each Q point (1 : nQpoints).
    !< Only used if the CERC equation is used (i.e. if useAltSedimentTransport is false)
character(64) :: cliffElevationsFile = 'cliffElevations.txt'
    !< File always used.
    !< Provides values for the level of the cliff top.
    !< Contents: A value for each beach section (1 : nSections).
character(64) :: wavesAndTidesFile = 'wavesAndTides.txt'
    !< File always used.
    !< Provides values for tidal amplitude, wave period, wave height and wave angle.
    !<
    !< The first set of values read from the file is at line waveTideFileStart.
    !< The file is kept open during the run and a set of values is read from the file each tide.
    !< When the end of file is reached, reading starts again at the beginning of the file.
    !<
    !< Contents (each line in the file): tidal amplitude value, wave period value, 
    !< then a wave height and wave angle value for each wavepoint (1 : nWavepoints).
character(64) :: seaLevelFile = 'seaLevel.txt'
    !< File always used.
    !< Provides values for the sea level, indexed by year.
    !< Contents (each line in the file): A year value and a sea level value.
    !< The file must include a sea level value for every year in the run.
character(64) :: bermWidthsFile = 'bermWidths.txt'
    !< File used if runMode is restart.
    !< Provides values for the initial width of the berm.
    !< Contents: A value for each beach section (1 : nSections).
character(64) :: beachVolumesFile = 'beachVolumes.txt'
    !< File used if runMode is restart.
    !< Provides values for initial beach volumes.
    !< Contents: A value for each beach section (1 : nSections).
character(64) :: talusVolumeFile = 'talusVolume.txt'
    !< File used if runMode is restart.
    !< Provides values for the initial talus volume.
    !< Contents: A value for each beach section (1 : nSections).
character(64) :: offshoreBarVolFile = 'offshoreBarVol.txt'
    !< File used if runMode is restart.
    !< Provides values for the initial offshore bar volume.
    !< Contents: A value for each beach section (1 : nSections).
character(64) :: rockStrengthFile = 'rockStrength.txt'
    !< File used if useVAgrid is .false.
    !< Provides values for the strength of the platform/cliff material.
    !< Contents: For each cliff element (1 : nce), a value for each beach section (1 : nSections).
character(64) :: beachContentFile = 'beachContent.txt'
    !< File used if useVAgrid is .false.
    !< Provides values for initial sand fraction (sediment content) ratio.
    !< Contents: A value (proportion in the range 0.0 to 1.0) for each beach section (1 : nSections).
character(64) :: startProfilesFile = 'startProfiles.txt'
    !< File used if runMode is restart and useVAgrid is .false.
    !< Provides values for the initial cliff profile.
    !< Contents: For each cliff element (1 : nce), a value (m) for each active beach section (1 : lastActiveSection).
character(64) :: profOffsetsFile = 'profOffsets.txt'
    !< File used if runMode is restart and useVAgrid is .false.
    !< Provides values for initial section offsets.
    !< Contents: A value (m) for each beach section (1 : nSections).
character(64) :: baseGeologyFile = 'baseGeology.txt'
    !< File used if useVAgrid is .true.
    !< Provides values for the base geology profiles.
    !< Contents: strength value and sediment content value (proportion in the range 0.0 to 1.0) , then a count of the number of profile points
    !< followed by pairs of Y (m) and Z (m from OD) values for the base geology profile. 
    !< This data repeated for each beach section in the model (1 : nSections).
character(64) :: layerDataFile = 'layerData.txt'
    !< File used if useVAgrid is .true.
    !< Provides values for the layer profiles.
    !< Contents: A count of the number of layers initially present, then for each layer, 
    !< a strength value and sediment content value (proportion in the range 0.0 to 1.0),
    !< followed by a count of the number of profile points, followed by pairs of Y (m) and depth (m) values. 
    !< This data repeated for each beach section in the model (1 : nSections).
character(64) :: cliffDataFile = 'cliffData.txt'
    !< File used if nCliffSimulations is greater than 0.
    !< Provides values for each cliff data section.
    !< Contents (each line in the file): cliff height, initial mean angle, initial angle stdeviation,
    !< mean falling angle, falling angle stdeviation, mean stable angle, stable angle stdeviation
    !< and beach section number for each cliff data section (1 : nSections).
character(64) :: seawallConstructionYearsFile = 'seawallConstructionYears.txt'
    !< File used if includeSeawall is .true.
    !< Provides dates at which the seawalls are installed and optionally provides seawall positions.
    !< Contents: A seawall construction year value or zero, followed by a seawall position value 
    !< if seawallMode is defineSeawallsByYearAndPosition, for each beach section (1 : nSections).
character(64) :: seawallRemovalYearsFile = 'seawallRemovalYears.txt'
    !< File used if includeSeawall is .true.
    !< Provides dates at which the seawalls are removed.
    !< Contents: A seawall removal year value or zero for each beach section (1 : nSections).
character(64) :: groyneConstructionYearsFile = 'groyneConstructionYears.txt'
    !< File used if includeGroynes is .true.
    !< Provides dates at which the groynes are installed.
    !< Contents: A groyne construction year value or zero for each Q point (1 : nQpoints).
character(64) :: groyneRemovalYearsFile = 'groyneRemovalYears.txt'
    !< File used if includeGroynes is .true.
    !< Provides dates at which the groynes are removed.
    !< Contents: A groyne removal year value or zero for each Q point (1 : nQpoints).
character(64) :: waveExclusionAnglesFile = 'waveExclusionAngles.txt'
    !< File used if excludeWaves is 2
    !< Provides angles from shore normal, beyond which waves are excluded (degrees).
    !< Contents: An exclusionAngleL and an exclusionAngleR for each Q point (1 : nQpoints).
    !< exclusionAngleL is an angle anti-clockwise from shore normal, beyond which waves are excluded (degrees).
    !< exclusionAngleR is an angle clockwise from shore normal, beyond which waves are excluded (degrees).
character(64) :: startprofileFile = 'startprofile.txt'
    !< File not used.
character(64) :: vellingaSectionsFile = 'vellingaSections.txt'
    !< File used if vellingaActive is .true.
    !< Identifies beach sections where vellinga effects can apply.
    !< Contents: A list of section numbers.

contains

!> Initialise parameters by reading from the setup.txt file.
!>
!> Apply defaults.
!>
!> Raise exceptions for missing or invalid parameters.
subroutine INITIALISE_PARAMETERS(path, radcon, success_code)

    use utilities_module
    use local_data
    use exceptions
    implicit none

    character(*), intent(in) :: path !< path to setup file
    real(kind=double), intent(in) :: radcon !< degrees to radians conversion factor
    integer, intent(inout) :: success_code !< 0: ok, -ive fatal error, +ive warning

    integer :: factorCount !< Number of year-factor pairs to read
    character(256) :: steering_file
    integer :: steering_file_unit
    character(32) :: keyword
    logical :: file_exists
    integer :: iostat
    integer :: alloc_error
    integer :: heights
    integer :: periods
    integer :: height
    integer :: widths

    steering_file = trim(path)//"setup.txt"
    feLogString = 'Opening '//steering_file
    call feLog()
    inquire(file=steering_file, exist=file_exists)
    
    if (file_exists) then
        feLogString = 'Opened '//steering_file
        call feLog()
        steering_file_unit = get_file_unit()
        open(steering_file_unit, file=steering_file, action='READ', iostat=iostat)
        
        if (iostat == 0) then
            do
                read(steering_file_unit, fmt=*, iostat = iostat) keyword

                if (iostat == -1) then
                    exit
                else
                    success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading keyword')
                    if (success_code < 0) return
                end if

                select case(UPPERCASE(trim(keyword)))
                    case('RESTARTMODE')
                        runMode = restart
                    case('STARTMODE')
                        runMode = start
                    case('STARTYEAR')
                        read(steering_file_unit, fmt=*, iostat= iostat) startYear
                    case('ENDYEAR')
                        read(steering_file_unit, fmt=*, iostat= iostat) endYear
                    case('SPINUPTOYEAR')
                        read(steering_file_unit, fmt=*, iostat= iostat) spinUpToYear
                    case('BEACHZERODELAY')
                        read(steering_file_unit, fmt=*, iostat= iostat) beachZeroDelay
                    case('LOWRESOUTTIMESTEP')
                        read(steering_file_unit, fmt=*, iostat= iostat) lowResOutTimestep
                    case('TIDESPERYEAR')
                        read(steering_file_unit, fmt=*, iostat= iostat) tidesPerYear
                    case('HIGHRESOUTSTARTYEAR')
                        read(steering_file_unit, fmt=*, iostat= iostat) highResOutStartYear
                    case('HIGHRESOUTTIMESTEP')
                        read(steering_file_unit, fmt=*, iostat= iostat) highResOutTimestep
                    case('FIRSTYEARTRANSPORTOUTPUT')
                        read(steering_file_unit, fmt=*, iostat= iostat) firstYearTransportOutput
                    case('PROFILEOUTPUTTIMESTEP')
                        read(steering_file_unit, fmt=*, iostat= iostat) profileOutputTimestep
                    case('FIRSTYEARPROFILEOUTPUT')
                        read(steering_file_unit, fmt=*, iostat= iostat) firstYearProfileOutput
                    case('LASTYEARPROFILEOUTPUT')
                        read(steering_file_unit, fmt=*, iostat= iostat) lastYearProfileOutput
                    case('FIRSTSECTIONPROFILEOUTPUT')
                        read(steering_file_unit, fmt=*, iostat= iostat) firstSectionProfileOutput
                    case('LASTSECTIONPROFILEOUTPUT')
                        read(steering_file_unit, fmt=*, iostat= iostat) lastSectionProfileOutput
                    case('NSECTIONS')
                        read(steering_file_unit, fmt=*, iostat= iostat) nSections
                    case('SECTIONWIDTH')
                        read(steering_file_unit, fmt=*, iostat= iostat) sectionWidth
                    case('SECTIONHEIGHT')
                        read(steering_file_unit, fmt=*, iostat= iostat) sectionHeight
                    case('ELEMENTHEIGHT')
                        read(steering_file_unit, fmt=*, iostat= iostat) elementHeight
                    case('USEVAGRID')
                        useVAgrid = .true.
                    case('ELEMENTWIDTH')
                        read(steering_file_unit, fmt=*, iostat= iostat) elementWidth
                    case('YMAXH')
                        read(steering_file_unit, fmt=*, iostat= iostat) YmaxH
                    case('BASELINEANGLE')
                        read(steering_file_unit, fmt=*, iostat= iostat) baselineAngle
                        baselineAngleRequired = .false.
                    case('FIRSTACTIVESECTION')
                        read(steering_file_unit, fmt=*, iostat= iostat) firstActiveSection
                    case('LASTACTIVESECTION')
                        read(steering_file_unit, fmt=*, iostat= iostat) lastActiveSection
                    case('WAVETIDEFILESTART')
                        read(steering_file_unit, fmt=*, iostat= iostat) waveTideFileStart
                    case('NWAVEPOINTS')
                        read(steering_file_unit, fmt=*, iostat= iostat) nWavePoints
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading nWavePoints')
                        if (success_code < 0) return

                        allocate(fswp(nWavePoints), depthWPmsl(nWavePoints), stat = alloc_error)
                        success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'fswp, depthWPmsl')
                        if (success_code < 0) return

                        read(steering_file_unit, fmt=*, iostat= iostat) fswp
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading fswp')
                        if (success_code < 0) return

                        read(steering_file_unit, fmt=*, iostat= iostat) depthWPmsl
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading depthWPmsl')
                        if (success_code < 0) return
                    case('DEFINESEAWALLSBYYEAR')
                        seawallMode = defineSeawallsByYear
                    case('DEFINESEAWALLSBYYEARANDPOSITION')
                        seawallMode = defineSeawallsByYearAndPosition
                    case('SEAWALLBASELEVEL')
                        read(steering_file_unit, fmt=*, iostat= iostat) seawallBaseLevel
                        seawallBaseLevelSupplied = .true.
                    case('MSLM')
                        read(steering_file_unit, fmt=*, iostat= iostat) mslM
                    case('MSLOFFSET')
                        read(steering_file_unit, fmt=*, iostat= iostat) mslOffset
                    case('DEPTHOSCMSL')
                        read(steering_file_unit, fmt=*, iostat= iostat) depthOSCMsl
                        depthOSCMslRequired = .false.
                    case('OUTPUTCONTOURLEVEL')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) outputContourLevelUser
                    case('MINSTORMWAVEHEIGHT')
                        read(steering_file_unit, fmt=*, iostat= iostat) minStormWaveHeight
                    case('MINHS')
                        read(steering_file_unit, fmt=*, iostat= iostat) minHs
                    case('MAXPLATFORMSLOPE')
                        read(steering_file_unit, fmt=*, iostat= iostat) maxPlatformSlope
                    case('BERMSLOPE')
                        read(steering_file_unit, fmt=*, iostat= iostat) bermSlope
                    case('INITIALBERMSTEP')
                        read(steering_file_unit, fmt=*, iostat= iostat) initialBermStep
                    case('MINBERMSTEP')
                        read(steering_file_unit, fmt=*, iostat= iostat) minBermStep
                    case('MAXBERMSTEP')
                        read(steering_file_unit, fmt=*, iostat= iostat) maxBermStep
                    case('BEACHDISTURBANCERATIO')
                        read(steering_file_unit, fmt=*, iostat= iostat) beachDisturbanceRatio
                    case('BEACHRETURNRATIO')
                        read(steering_file_unit, fmt=*, iostat= iostat) beachReturnRatio
                    case('TALUSSTRENGTHRATIO')
                        read(steering_file_unit, fmt=*, iostat= iostat) talusStrengthRatio
                    case('TALUSSLOPE')
                        read(steering_file_unit, fmt=*, iostat= iostat) talusSlope
                    case('MINBOUNDARYTRANSPORTRATIO')
                        read(steering_file_unit, fmt=*, iostat= iostat) minBoundaryTransportRatio
                    case('QPMAXBOUNDARYRIGHTIN')
                        read(steering_file_unit, fmt=*, iostat= iostat) qpMaxBoundaryRightIn
                    case('QPMAXBOUNDARYLEFTIN')
                        read(steering_file_unit, fmt=*, iostat= iostat) qpMaxBoundaryLeftIn
                    case('QPMAXBOUNDARYRIGHTOUT')
                        read(steering_file_unit, fmt=*, iostat= iostat) qpMaxBoundaryRightOut
                    case('QPMAXBOUNDARYLEFTOUT')
                        read(steering_file_unit, fmt=*, iostat= iostat) qpMaxBoundaryLeftOut
                    case('MAXBLOCKSIZE')
                        read(steering_file_unit, fmt=*, iostat= iostat) maxBlockSize
                    case('BEACHHEIGHT')
                        read(steering_file_unit, fmt=*, iostat= iostat) beachHeight
                    case('BRUUNCONST')
                        read(steering_file_unit, fmt=*, iostat= iostat) bruunConst
                    case('BEACHCRESTLEVEL')
                        read(steering_file_unit, fmt=*, iostat= iostat) beachCrestLevel
                        beachCrestLevelRequired = .false.
                    case('EXCLUDEWAVES')
                        read(steering_file_unit, fmt=*, iostat= iostat) excludeWaves
                    case('WAVEEXCLUSIONANGLES')
                        read(steering_file_unit, fmt=*, iostat= iostat) waveExclusionAngleL, waveExclusionAngleR
                    case('SWITCHTRANSPORTWHDIFF')
                        read(steering_file_unit, fmt=*, iostat= iostat) switchTransportWhDiff
                    case('NCLIFFSIMULATIONS')
                        read(steering_file_unit, fmt=*, iostat= iostat) nCliffSimulations
                    case('SLUMPPERIOD')
                        read(steering_file_unit, fmt=*, iostat= iostat) slumpPeriod
                    case('USEALTSEDIMENTTRANSPORT')
                        useAltSedimentTransport = .true.
                    case('MEDIANGRAINSIZE')
                        read(steering_file_unit, fmt=*, iostat= iostat) medianGrainsize
                    case('CROSSBEACH_STR')
                        read(steering_file_unit, fmt=*, iostat= iostat) heights, periods
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading heights, periods')
                        if (success_code < 0) return
  
                        allocate(cbstrWaveHeight(heights), &
                            cbstrPeriod(periods, heights), &
                            cbstrValue(periods, heights), &
                            stat = alloc_error)
                        success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'cbstrWaveHeight, cbstrPeriod, cbstrValue')
                        if (success_code < 0) return

                        read(steering_file_unit, fmt=*, iostat= iostat) cbstrWaveHeight
                            success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading cbstrWaveHeight')
                            if (success_code < 0) return
                        
                        do height = 1, heights
                            read(steering_file_unit, fmt=*, iostat= iostat) cbstrPeriod(:, height)
                            success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading cbstrPeriod')
                            if (success_code < 0) return

                            read(steering_file_unit, fmt=*, iostat= iostat) cbstrValue(:, height)
                            success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading cbstrValue')
                            if (success_code < 0) return
                        end do
                    case('OFFSHORE_STR')
                        read(steering_file_unit, fmt=*, iostat= iostat) heights, periods
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading heights, periods')
                        if (success_code < 0) return
  
                        allocate(offstrWaveHeight(heights), &
                            offstrPeriod(periods, heights), &
                            offstrValue(periods, heights), &
                            stat = alloc_error)
                        success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'offstrWaveHeight, offstrPeriod, offstrValue')
                        if (success_code < 0) return

                        read(steering_file_unit, fmt=*, iostat= iostat) offstrWaveHeight
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading offstrWaveHeight')
                        if (success_code < 0) return

                        do height = 1, heights
                            read(steering_file_unit, fmt=*, iostat= iostat) offstrPeriod(:, height)
                            success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading offstrPeriod')
                            if (success_code < 0) return
                            
                            read(steering_file_unit, fmt=*, iostat= iostat) offstrValue(:, height)
                            success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading offstrValue')
                            if (success_code < 0) return
                        end do
                    case('GROYNEEFFECT')
                        read(steering_file_unit, fmt=*, iostat= iostat) widths
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading widths')
                        if (success_code < 0) return
  
                        allocate(groyneEffectBeachWidth(widths), &
                            groyneEffectValue(widths), &
                            stat = alloc_error)
                        success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'groyneEffectBeachWidth, groyneEffectValue')
                        if (success_code < 0) return

                        read(steering_file_unit, fmt=*, iostat= iostat) groyneEffectBeachWidth
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading groyneEffectBeachWidth')
                        if (success_code < 0) return

                        read(steering_file_unit, fmt=*, iostat= iostat) groyneEffectValue
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading groyneEffectValue')
                        if (success_code < 0) return
                    case('OSCANGLES')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) oscAnglesFile
                    case('ROCKSTRENGTH')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) rockStrengthFile
                    case('CERCCONSTANTS')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) cercConstantsFile
                    case('BERMWIDTHS')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) bermWidthsFile
                    case('STARTPROFILE')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) startprofileFile
                    case('STARTPROFILES')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) startProfilesFile
                    case('PROFOFFSETS')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) profOffsetsFile
                    case('BEACHVOLUMES')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) beachVolumesFile
                    case('BEACHCONTENT')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) beachContentFile
                    case('CLIFFELEVATIONS')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) cliffElevationsFile
                    case('TALUSVOLUME')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) talusVolumeFile
                    case('OFFSHOREBARVOL')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) offshoreBarVolFile
                    case('WAVESANDTIDES')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) wavesAndTidesFile
                    case('SEALEVEL')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) seaLevelFile
                    case('CLIFFDATA')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) cliffDataFile
                    case('BASEGEOLOGY')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) baseGeologyFile
                    case('LAYERDATA')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) layerDataFile
                    case('SEAWALLCONSTRUCTIONYEARS')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) seawallConstructionYearsFile
                    case('SEAWALLREMOVALYEARS')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) seawallRemovalYearsFile
                    case('GROYNECONSTRUCTIONYEARS')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) groyneConstructionYearsFile
                    case('GROYNEREMOVALYEARS')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) groyneRemovalYearsFile
                    case('TIDALAMPSECTION')
                        read(steering_file_unit, fmt=*, iostat= iostat) tidalAmpSection
                    case('TIDALAMPSCALELEFT')
                        read(steering_file_unit, fmt=*, iostat= iostat) tidalAmpScaleLeft
                    case('TIDALAMPSCALERIGHT')
                        read(steering_file_unit, fmt=*, iostat= iostat) tidalAmpScaleRight
                    case('NOTIDALRANGEVARIATION')
                        tidalRangeVariation = .false.
                    case('INCLUDESEAWALL')
                        includeSeawall = .true.
                    case('INCLUDEGROYNES')
                        includeGroynes = .true.
                    case('CELERITYTABLESASTXTFILES')
                        celerityTablesAsTxtFiles = .true.
                        read(steering_file_unit, fmt=*, iostat= iostat) periodValues, depthValues
                    case('VELLINGAACTIVE')
                        vellingaActive = .true.
                    case('MODIFYWAVEPERIOD')
                        read(steering_file_unit, fmt=*, iostat= iostat) factorCount
                        
                        if (factorCount > 0) then
                            allocate(modifyPeriodYear(factorCount), modifyPeriodFactor(factorCount), stat = alloc_error)
                            success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'modifyPeriodYear, modifyPeriodFactor')
                            if (success_code < 0) return

                            call READ_MODIFY_YEAR_DATA(steering_file_unit, keyWord, factorCount, &
                                modifyPeriodYear, modifyPeriodFactor, success_code)
                            if (success_code < 0) return
                        end if
                    case('MODIFYWAVEHEIGHT')
                        read(steering_file_unit, fmt=*, iostat= iostat) factorCount
                        
                        if (factorCount > 0) then
                            allocate(modifyHeightYear(factorCount), modifyHeightFactor(factorCount), stat = alloc_error)
                            success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'modifyHeightYear, modifyHeightFactor')
                            if (success_code < 0) return

                            call READ_MODIFY_YEAR_DATA(steering_file_unit, keyWord, factorCount, &
                                modifyHeightYear, modifyHeightFactor, success_code)
                            if (success_code < 0) return
                        end if
                    case('MODIFYWAVEDIRECTION')
                        read(steering_file_unit, fmt=*, iostat= iostat) factorCount
                        
                        if (factorCount > 0) then
                            allocate(modifyAngleYear(factorCount), modifyAngleAmount(factorCount), stat = alloc_error)
                            success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'modifyAngleYear, modifyAngleAmount')
                            if (success_code < 0) return

                            call READ_MODIFY_YEAR_DATA(steering_file_unit, keyWord, factorCount, &
                                modifyAngleYear, modifyAngleAmount, success_code)
                            if (success_code < 0) return
                        end if

                    case default
                        call RAISE_EXCEPTION('Keyword not recognized: '//trim(keyWord))
                        success_code = -1
                        return
                end select
    
                success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'After keyword: '//keyword)
                if (success_code < 0) return

            end do

            close(steering_file_unit)
        end if
    end if
    
    if (.not. allocated(cbstrWaveHeight)) then
        ! Default CROSSBEACH_STR table
        allocate(cbstrWaveHeight(4), cbstrPeriod(3, 4), cbstrValue(3, 4), &
            stat = alloc_error)  ! MD 300304
        success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'cbstrWaveHeight, cbstrPeriod, cbstrValue')
        if (success_code < 0) return

        !**************************************************
        ! Find sed trans rate between upper and lower beach
        ! 130602, see WkBk p155(10)
        ! in metres cubed per m per hour
        ! These were all calculated for a 'typical' norfolk
        ! beach using COSMOS
        !**************************************************
        cbstrWaveHeight = (/ 3.0D0, 4.0D0, 5.0D0, -1.0D0 /)    
        cbstrPeriod(:, 1) = (/ 5.0D0, 6.0D0, -1.0D0 /)    
        cbstrValue(:, 1) = (/ -0.55D0, -0.74D0, -0.66D0 /)    
        cbstrPeriod(:, 2) = (/ 6.0D0, -1.0D0, 0.0D0 /)    
        cbstrValue(:, 2) = (/ -1.0D0, -1.27D0, 0.0D0 /)    
        cbstrPeriod(:, 3) = (/ 6.0D0, -1.0D0, 0.0D0 /)    
        cbstrValue(:, 3) = (/ -1.06D0, -1.46D0, 0.0D0 /)    
        cbstrPeriod(:, 4) = (/ 6.0D0, -1.0D0, 0.0D0 /)    
        cbstrValue(:, 4) = (/ -1.09D0, -1.55D0, 0.0D0 /)    
    end if

    if (.not. allocated(offstrWaveHeight)) then
        ! Default OFFSHORE_STR table
        allocate(offstrWaveHeight(4), offstrPeriod(3, 4), offstrValue(3, 4), &
            stat = alloc_error)  ! MD 300304
        success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'offstrWaveHeight, offstrPeriod, offstrValue')
        if (success_code < 0) return

        !**************************************************
        ! Find sed trans rate in metres cubed per m per hour
        ! These were all calculated for a 'typical' norfolk
        ! beach using COSMOS
        !**************************************************
        offstrWaveHeight = (/ 3.0D0, 4.0D0, 5.0D0, -1.0D0 /)    
        offstrPeriod(:, 1) = (/ 5.0D0, 6.0D0, -1.0D0 /)    
        offstrValue(:, 1) = (/ -0.03D0, -0.08D0, -0.09D0 /)    
        offstrPeriod(:, 2) = (/ 6.0D0, -1.0D0, 0.0D0 /)    
        offstrValue(:, 2) = (/ -0.16D0, -0.42D0, 0.0D0 /)    
        offstrPeriod(:, 3) = (/ 6.0D0, -1.0D0, 0.0D0 /)    
        offstrValue(:, 3) = (/ -0.17D0, -0.57D0, 0.0D0 /)    
        offstrPeriod(:, 4) = (/ 6.0D0, -1.0D0, 0.0D0 /)    
        offstrValue(:, 4) = (/ -0.18D0, -0.63D0, 0.0D0 /)    
    end if
    
    if (.not. allocated(groyneEffectBeachWidth)) then
        ! Default GROYNEEFFECT table
        allocate(groyneEffectBeachWidth(6), groyneEffectValue(6), &
            stat = alloc_error)  ! MD 300304
        success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'groyneEffectBeachWidth, groyneEffectValue')
        if (success_code < 0) return

        groyneEffectBeachWidth = (/ 50.0D0, 100.0D0, 150.0D0, 200.0D0, 250.0D0, -1.0D0 /)    
        groyneEffectValue = (/ 0.95D0, 0.90D0, 0.5D0, 0.4D0, 0.3D0, 0.2D0 /)    
    end if

    if (nWavePoints == -1) then
        nWavePoints = 1
        
        allocate(fswp(nWavePoints), depthWPmsl(nWavePoints), stat = alloc_error)
        success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'fswp, depthWPmsl')
        if (success_code < 0) return
        
        fswp(1) = 1
        depthWPmsl(1) = 10.0D0
    end if
    
    maxPlatformSlope_radians = radcon * maxPlatformSlope
    wave_exclusion_angle_L_radians = radcon * waveExclusionAngleL
    wave_exclusion_angle_R_radians = radcon * waveExclusionAngleR
    baselineAngle_radians = radcon * baselineAngle
    talusSlope_radians = radcon * talusSlope
    dx = sectionWidth
    dy = elementWidth
    dz = elementHeight

    if (lastActiveSection == -1) lastActiveSection = nSections

    if (baselineAngleRequired) then
        call RAISE_EXCEPTION('Please supply a value for baselineAngle')
        success_code = -1
        return
    end if

    if (startYear == -1000000) then
        call RAISE_EXCEPTION('Please supply a value for startYear')
        success_code = -1
        return
    end if

    if (endYear == -1000000) then
        call RAISE_EXCEPTION('Please supply a value for endYear')
        success_code = -1
        return
    end if

    if (spinUpToYear == -1000000) then
        spinUpToYear = startYear
    else if (spinUpToYear >= endYear .or. spinUpToYear < startYear) then
        write(logString,*) 'Invalid spinUpToYear: ', spinUpToYear, 'startYear: ', startYear, 'endYear: ', endYear
        call RAISE_EXCEPTION(logString)
        success_code = -1
        return
    end if

    if (highResOutStartYear == -1000000) highResOutStartYear = startYear
    
    if (firstYearTransportOutput == -1000000) firstYearTransportOutput = startYear

    if (lastYearProfileOutput == -1000000) lastYearProfileOutput = endYear
    
    if (lastSectionProfileOutput == -1) lastSectionProfileOutput = nSections

    if (firstYearProfileOutput /= -1000000) then
        if (firstYearProfileOutput < startYear .or. firstYearProfileOutput > endYear) then
           call RAISE_EXCEPTION('Invalid firstYearProfileOutput')
           success_code = -1
           return
        end if
        
        if (lastYearProfileOutput < firstYearProfileOutput .or. lastYearProfileOutput > endYear) then
           call RAISE_EXCEPTION('Invalid lastYearProfileOutput')
           success_code = -1
           return
        end if

        if (firstSectionProfileOutput < 1 .or. firstSectionProfileOutput > nSections) then
           call RAISE_EXCEPTION('Invalid firstSectionProfileOutput')
           success_code = -1
           return
        end if

        if (lastSectionProfileOutput < firstSectionProfileOutput .or. lastSectionProfileOutput > nSections) then
           call RAISE_EXCEPTION('Invalid lastSectionProfileOutput')
           success_code = -1
           return
        end if
    end if
    
    if (useVAgrid) then
        if (YmaxH == -1.0D0) then
            call RAISE_EXCEPTION('Please supply a value for YmaxH')
            success_code = -1
            return
        end if
    end if

    if (beachCrestLevelRequired) then
        call RAISE_EXCEPTION('Please supply a value for beachCrestLevel')
        success_code = -1
        return
    end if

    if (depthOSCMslRequired) then
        call RAISE_EXCEPTION('Please supply a value for depthOSCMsl')
        success_code = -1
        return
    end if

    if (beachCrestLevel < beachHeight) then
        write(logString,*) 'beachCrestLevel ', beachCrestLevel, ' must be at least beachHeight ', beachHeight
        call RAISE_EXCEPTION(logString)
        success_code = -1
        return
    end if

    if (includeSeawall .and. seawallMode == defineSeawallsByYear) then
        if (.not. seawallBaseLevelSupplied) then
            call RAISE_EXCEPTION('Please supply a value for seawallBaseLevel')
            success_code = -1
            return
        end if
    end if

    if (trim(outputContourLevelUser) == '') then
        outputContourLevel_e = nint(sectionHeight / elementHeight)
    else
        read(outputContourLevelUser, *) outputContourLevel

        if ((outputContourLevel + mslM) > sectionHeight .or. (outputContourLevel + mslM) < 0.0D0) then
           call RAISE_EXCEPTION('Invalid outputContourLevel')
           success_code = -1
           return
        end if

        outputContourLevel_e = int((outputContourLevel + mslM)/ elementHeight) + 1
    end if

end subroutine INITIALISE_PARAMETERS

!> Read data for modification of wave values by year.
!> Raise exceptions for missing or invalid parameters.
subroutine READ_MODIFY_YEAR_DATA(steering_file_unit, keyWord, valueCount, years, values, success_code)

    use exceptions
    implicit none

    integer, intent(in) :: steering_file_unit !< Steering file unit number
    character(*), intent(in) :: keyWord !< Keyword in setup.txt
    integer, intent(in) :: valueCount !< Number of year/value pairs expected
    integer, intent(out), dimension(valueCount) :: years !< Years at which values apply
    real(kind=double), intent(out), dimension(valueCount) :: values !< Value at given year
    integer, intent(inout) :: success_code !< 0: ok, -ive fatal error, +ive warning

    integer :: previous ! Previous year for validation
    integer :: index
    integer :: iostat

    read(steering_file_unit, fmt=*, iostat= iostat) years
    success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading '//trim(keyWord)//' years')
    if (success_code < 0) return

    read(steering_file_unit, fmt=*, iostat= iostat) values
    success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading '//trim(keyWord)//' values')
    if (success_code < 0) return

    if (valueCount > 1) then
        previous = years(1)

        do index = 2, valueCount
            if (years(index) <= previous) then
                call RAISE_EXCEPTION('Years must increase: '//trim(keyWord))
                success_code = -1
                return
            end if

            previous = years(index)
        end do
    end if

    if (trim(keyWord) /= 'MODIFYWAVEDIRECTION' .and. minval(values) < 0.0D0) then
        call RAISE_EXCEPTION('Negative values not permitted: '//trim(keyWord))
        success_code = -1
        return
    end if
    
end subroutine READ_MODIFY_YEAR_DATA


!> Write an array of double precision values to the log file.
subroutine LOG_DP_VALUES(values)
    use local_data, only: logString
    use utilities_module, only: SCAPE_LOG
    implicit none

    real(kind=double), intent(in), dimension(:) :: values !< array of values to be written

    integer :: fullLines
    integer :: tail
    integer :: index
    integer :: counter

    fullLines = size(values, 1) / 4
    tail = mod(size(values, 1), 4)
    counter = 1

    do index = 1, fullLines
        write(logString, fmt='(4G20.8)') values(counter : counter + 3)
        call SCAPE_LOG()
        counter = counter + 4
    end do

    if (tail > 0) then
        write(logString, fmt='(4G20.8)') values(counter : counter + tail - 1)
        call SCAPE_LOG()
    end if
    
end subroutine LOG_DP_VALUES


!> Write an array of integer values to the log file.
subroutine LOG_I_VALUES(values)
    use local_data, only: logString
    use utilities_module, only: SCAPE_LOG
    implicit none

    integer, intent(in), dimension(:) :: values !< array of values to be written

    integer :: fullLines
    integer :: tail
    integer :: index
    integer :: counter

    fullLines = size(values, 1) / 4
    tail = mod(size(values, 1), 4)
    counter = 1

    do index = 1, fullLines
        write(logString, fmt='(4I8)') values(counter : counter + 3)
        call SCAPE_LOG()
        counter = counter + 4
    end do

    if (tail > 0) then
        write(logString, fmt='(4I8)') values(counter : counter + tail - 1)
        call SCAPE_LOG()
    end if
    
end subroutine LOG_I_VALUES


!> Write all the parameters to the log file in a format which can be read by INITIALISE_PARAMETERS().
subroutine LOG_PARAMETERS()
    use local_data, only: logString
    use utilities_module, only: SCAPE_LOG
    implicit none

    integer :: index
    
    logString = 'Values of all parameters: '
    call SCAPE_LOG()
    logString = ''
    call SCAPE_LOG()

    if (runMode == restart) then
        logString = 'RESTARTMODE'
    else
        logString = 'STARTMODE'
    end if    
    call SCAPE_LOG()
    
    logString = 'STARTYEAR'
    call SCAPE_LOG()
    write(logString, *) startYear
    call SCAPE_LOG()

    logString = 'ENDYEAR'
    call SCAPE_LOG()
    write(logString, *) endYear
    call SCAPE_LOG()
    
    logString = 'SPINUPTOYEAR'
    call SCAPE_LOG()
    write(logString, *) spinUpToYear
    call SCAPE_LOG()
    
    logString = 'BEACHZERODELAY'
    call SCAPE_LOG()
    write(logString, *) beachZeroDelay
    call SCAPE_LOG()
    
    logString = 'LOWRESOUTTIMESTEP'
    call SCAPE_LOG()
    write(logString, *) lowResOutTimestep
    call SCAPE_LOG()
    
    logString = 'TIDESPERYEAR'
    call SCAPE_LOG()
    write(logString, *) tidesPerYear
    call SCAPE_LOG()
    
    logString = 'HIGHRESOUTSTARTYEAR'
    call SCAPE_LOG()
    write(logString, *) highResOutStartYear
    call SCAPE_LOG()
    
    logString = 'HIGHRESOUTTIMESTEP'
    call SCAPE_LOG()
    write(logString, *) highResOutTimestep
    call SCAPE_LOG()
    
    logString = 'FIRSTYEARTRANSPORTOUTPUT'
    call SCAPE_LOG()
    write(logString, *) firstYearTransportOutput
    call SCAPE_LOG()
    
    logString = 'PROFILEOUTPUTTIMESTEP'
    call SCAPE_LOG()
    write(logString, *) profileOutputTimestep
    call SCAPE_LOG()
    
    logString = 'FIRSTYEARPROFILEOUTPUT'
    call SCAPE_LOG()
    write(logString, *) firstYearProfileOutput
    call SCAPE_LOG()
    
    logString = 'LASTYEARPROFILEOUTPUT'
    call SCAPE_LOG()
    write(logString, *) lastYearProfileOutput
    call SCAPE_LOG()
    
    logString = 'FIRSTSECTIONPROFILEOUTPUT'
    call SCAPE_LOG()
    write(logString, *) firstSectionProfileOutput
    call SCAPE_LOG()
    
    logString = 'LASTSECTIONPROFILEOUTPUT'
    call SCAPE_LOG()
    write(logString, *) lastSectionProfileOutput
    call SCAPE_LOG()
    
    logString = 'NSECTIONS'
    call SCAPE_LOG()
    write(logString, *) nSections
    call SCAPE_LOG()
    
    logString = 'SECTIONWIDTH'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') sectionWidth
    call SCAPE_LOG()
    
    logString = 'SECTIONHEIGHT'
    call SCAPE_LOG()
    write(logString, *) sectionHeight
    call SCAPE_LOG()
    
    logString = 'ELEMENTHEIGHT'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') elementHeight
    call SCAPE_LOG()
    
    if (useVAgrid) then
        logString = 'USEVAGRID'
        call SCAPE_LOG()
    end if
    
    logString = 'ELEMENTWIDTH'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') elementWidth
    call SCAPE_LOG()
    
    if (YmaxH /= -1.0D0) then
        logString = 'YMAXH'
        call SCAPE_LOG()
        write(logString, fmt='(G20.8)') YmaxH
        call SCAPE_LOG()
    end if
    
    logString = 'BASELINEANGLE'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') baselineAngle
    call SCAPE_LOG()
    
    logString = 'FIRSTACTIVESECTION'
    call SCAPE_LOG()
    write(logString, *) firstActiveSection
    call SCAPE_LOG()
    
    logString = 'LASTACTIVESECTION'
    call SCAPE_LOG()
    write(logString, *) lastActiveSection
    call SCAPE_LOG()
    
    logString = 'WAVETIDEFILESTART'
    call SCAPE_LOG()
    write(logString, *) waveTideFileStart
    call SCAPE_LOG()
    
    logString = 'NWAVEPOINTS'
    call SCAPE_LOG()
    write(logString, *) nWavePoints
    call SCAPE_LOG()
    write(logString, *) fswp
    call SCAPE_LOG()
    call LOG_DP_VALUES(depthWPmsl)

    if (seawallMode == defineSeawallsByYear) then
        logString = 'DEFINESEAWALLSBYYEAR'
    else
        logString = 'DEFINESEAWALLSBYYEARANDPOSITION'
    end if    
    call SCAPE_LOG()
    
    if (seawallBaseLevelSupplied) then
        logString = 'SEAWALLBASELEVEL'
        call SCAPE_LOG()
        write(logString, fmt='(G20.8)') seawallBaseLevel
        call SCAPE_LOG()
    end if
    
    logString = 'MSLM'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') mslM
    call SCAPE_LOG()
    
    logString = 'MSLOFFSET'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') mslOffset
    call SCAPE_LOG()
    
    logString = 'DEPTHOSCMSL'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') depthOSCMsl
    call SCAPE_LOG()
    
    if (trim(outputContourLevelUser) /= '') then
        logString = 'OUTPUTCONTOURLEVEL'
        call SCAPE_LOG()
        write(logString, fmt='(A)') trim(outputContourLevelUser)
        call SCAPE_LOG()
    end if
    
    logString = 'MINSTORMWAVEHEIGHT'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') minStormWaveHeight
    call SCAPE_LOG()
    
    logString = 'MINHS'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') minHs
    call SCAPE_LOG()
    
    logString = 'MAXPLATFORMSLOPE'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') maxPlatformSlope
    call SCAPE_LOG()
    
    if (bermSlope > 0.0D0) then
        logString = 'BERMSLOPE'
        call SCAPE_LOG()
        write(logString, fmt='(G20.8)') bermSlope
        call SCAPE_LOG()
    end if
    
    logString = 'INITIALBERMSTEP'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') initialBermStep
    call SCAPE_LOG()
    
    logString = 'MINBERMSTEP'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') minBermStep
    call SCAPE_LOG()
    
    logString = 'MAXBERMSTEP'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') maxBermStep
    call SCAPE_LOG()
    
    logString = 'BEACHDISTURBANCERATIO'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') beachDisturbanceRatio
    call SCAPE_LOG()
    
    logString = 'BEACHRETURNRATIO'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') beachReturnRatio
    call SCAPE_LOG()
    
    logString = 'TALUSSTRENGTHRATIO'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') talusStrengthRatio
    call SCAPE_LOG()
    
    logString = 'TALUSSLOPE'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') talusSlope
    call SCAPE_LOG()
    
    logString = 'MINBOUNDARYTRANSPORTRATIO'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') minBoundaryTransportRatio
    call SCAPE_LOG()
    
    logString = 'QPMAXBOUNDARYRIGHTIN'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') qpMaxBoundaryRightIn
    call SCAPE_LOG()

    logString = 'QPMAXBOUNDARYLEFTIN'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') qpMaxBoundaryLeftIn
    call SCAPE_LOG()

    logString = 'QPMAXBOUNDARYRIGHTOUT'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') qpMaxBoundaryRightOut
    call SCAPE_LOG()

    logString = 'QPMAXBOUNDARYLEFTOUT'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') qpMaxBoundaryLeftOut
    call SCAPE_LOG()

    logString = 'MAXBLOCKSIZE'
    call SCAPE_LOG()
    write(logString, *) maxBlockSize
    call SCAPE_LOG()
    
    logString = 'BEACHHEIGHT'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') beachHeight
    call SCAPE_LOG()
    
    logString = 'BRUUNCONST'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') bruunConst
    call SCAPE_LOG()
    
    logString = 'BEACHCRESTLEVEL'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') beachCrestLevel
    call SCAPE_LOG()
    
    logString = 'EXCLUDEWAVES'
    call SCAPE_LOG()
    write(logString, *) excludeWaves
    call SCAPE_LOG()
    
    logString = 'WAVEEXCLUSIONANGLES'
    call SCAPE_LOG()
    write(logString, fmt='(2G20.8)') waveExclusionAngleL, waveExclusionAngleR
    call SCAPE_LOG()
    
    logString = 'SWITCHTRANSPORTWHDIFF'
    call SCAPE_LOG()
    write(logString, *) switchTransportWhDiff
    call SCAPE_LOG()
    
    logString = 'NCLIFFSIMULATIONS'
    call SCAPE_LOG()
    write(logString, *) nCliffSimulations
    call SCAPE_LOG()
    
    logString = 'SLUMPPERIOD'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') slumpPeriod
    call SCAPE_LOG()
    
    if (useAltSedimentTransport) then
        logString = 'USEALTSEDIMENTTRANSPORT'
        call SCAPE_LOG()
    end if
    
    logString = 'MEDIANGRAINSIZE'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') medianGrainsize
    call SCAPE_LOG()

    logString = 'CROSSBEACH_STR'
    call SCAPE_LOG()
    write(logString, *) size(cbstrWaveHeight, 1), size(cbstrPeriod, 1)
    call SCAPE_LOG()
    call LOG_DP_VALUES(cbstrWaveHeight)
        
    do index = 1, size(cbstrWaveHeight, 1)
        call LOG_DP_VALUES(cbstrPeriod(:, index))
        call LOG_DP_VALUES(cbstrValue(:, index))
    end do
    
    logString = 'OFFSHORE_STR'
    call SCAPE_LOG()
    write(logString, *) size(offstrWaveHeight, 1), size(offstrPeriod, 1)
    call SCAPE_LOG()
    call LOG_DP_VALUES(offstrWaveHeight)
    
    do index = 1, size(offstrWaveHeight, 1)
        call LOG_DP_VALUES(offstrPeriod(:, index))
        call LOG_DP_VALUES(offstrValue(:, index))
    end do
    
    logString = 'GROYNEEFFECT'
    call SCAPE_LOG()
    write(logString, *) size(groyneEffectBeachWidth, 1)
    call SCAPE_LOG()
    call LOG_DP_VALUES(groyneEffectBeachWidth)
    call LOG_DP_VALUES(groyneEffectValue)
    
    if (allocated(modifyPeriodYear)) then
        logString = 'MODIFYWAVEPERIOD'
        call SCAPE_LOG()
        write(logString, *) size(modifyPeriodYear, 1)
        call SCAPE_LOG()
        call LOG_I_VALUES(modifyPeriodYear)
        call LOG_DP_VALUES(modifyPeriodFactor)
    end if
    
    if (allocated(modifyHeightYear)) then
        logString = 'MODIFYWAVEHEIGHT'
        call SCAPE_LOG()
        write(logString, *) size(modifyHeightYear, 1)
        call SCAPE_LOG()
        call LOG_I_VALUES(modifyHeightYear)
        call LOG_DP_VALUES(modifyHeightFactor)
    end if
    
    if (allocated(modifyAngleYear)) then
        logString = 'MODIFYWAVEDIRECTION'
        call SCAPE_LOG()
        write(logString, *) size(modifyAngleYear, 1)
        call SCAPE_LOG()
        call LOG_I_VALUES(modifyAngleYear)
        call LOG_DP_VALUES(modifyAngleAmount)
    end if
    
    logString = 'OSCANGLES'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(oscAnglesFile)
    call SCAPE_LOG()
    
    logString = 'ROCKSTRENGTH'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(rockStrengthFile)
    call SCAPE_LOG()
    
    logString = 'CERCCONSTANTS'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(cercConstantsFile)
    call SCAPE_LOG()
    
    logString = 'BERMWIDTHS'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(bermWidthsFile)
    call SCAPE_LOG()
    
    logString = 'STARTPROFILE'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(startprofileFile)
    call SCAPE_LOG()
    
    logString = 'STARTPROFILES'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(startProfilesFile)
    call SCAPE_LOG()
    
    logString = 'PROFOFFSETS'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(profOffsetsFile)
    call SCAPE_LOG()
    
    logString = 'BEACHVOLUMES'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(beachVolumesFile)
    call SCAPE_LOG()
    
    logString = 'BEACHCONTENT'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(beachContentFile)
    call SCAPE_LOG()
    
    logString = 'CLIFFELEVATIONS'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(cliffElevationsFile)
    call SCAPE_LOG()
    
    logString = 'TALUSVOLUME'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(talusVolumeFile)
    call SCAPE_LOG()
    
    logString = 'OFFSHOREBARVOL'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(offshoreBarVolFile)
    call SCAPE_LOG()
    
    logString = 'WAVESANDTIDES'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(wavesAndTidesFile)
    call SCAPE_LOG()
    
    logString = 'SEALEVEL'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(seaLevelFile)
    call SCAPE_LOG()
    
    logString = 'CLIFFDATA'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(cliffDataFile)
    call SCAPE_LOG()
    
    logString = 'BASEGEOLOGY'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(baseGeologyFile)
    call SCAPE_LOG()
    
    logString = 'LAYERDATA'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(layerDataFile)
    call SCAPE_LOG()
    
    logString = 'SEAWALLCONSTRUCTIONYEARS'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(seawallConstructionYearsFile)
    call SCAPE_LOG()
    
    logString = 'SEAWALLREMOVALYEARS'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(seawallRemovalYearsFile)
    call SCAPE_LOG()
    
    logString = 'GROYNECONSTRUCTIONYEARS'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(groyneConstructionYearsFile)
    call SCAPE_LOG()
    
    logString = 'GROYNEREMOVALYEARS'
    call SCAPE_LOG()
    write(logString, fmt='(A)') trim(groyneRemovalYearsFile)
    call SCAPE_LOG()
    
    logString = 'TIDALAMPSECTION'
    call SCAPE_LOG()
    write(logString, *) tidalAmpSection
    call SCAPE_LOG()
    
    logString = 'TIDALAMPSCALELEFT'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') tidalAmpScaleLeft
    call SCAPE_LOG()
    
    logString = 'TIDALAMPSCALERIGHT'
    call SCAPE_LOG()
    write(logString, fmt='(G20.8)') tidalAmpScaleRight
    call SCAPE_LOG()

    if (.not. tidalRangeVariation) then
        logString = 'NOTIDALRANGEVARIATION'
        call SCAPE_LOG()
    end if

    if (includeSeawall) then
        logString = 'INCLUDESEAWALL'
        call SCAPE_LOG()
    end if

    if (includeGroynes) then
        logString = 'INCLUDEGROYNES'
        call SCAPE_LOG()
    end if
    
    if (celerityTablesAsTxtFiles) then
        logString = 'CELERITYTABLESASTXTFILES'
        call SCAPE_LOG()
        write(logString, *) periodValues, depthValues
        call SCAPE_LOG()
    end if
    
    if (vellingaActive) then
        logString = 'vellingaActive'
        call SCAPE_LOG()
    end if

end subroutine LOG_PARAMETERS


!> Returns the croosbeach sediment transport rate for a given wave height and period.
function CROSSBEACH_STR_LOOKUP(height, period) result(value)

    implicit none
    integer, parameter :: double = kind(1D0)

    real(kind=double), intent(in) :: height !< Wave height
    real(kind=double), intent(in) :: period !< Wave period
    real(kind=double) :: value !< Croosbeach sediment transport rate
    
    integer :: heightIndex
    integer :: periodIndex

    heightIndex = 1
    periodIndex = 1

    do while (cbstrWaveHeight(heightIndex) > 0.0D0 .and. height >= cbstrWaveHeight(heightIndex))
        heightIndex = heightIndex + 1
    end do

    do while (cbstrPeriod(periodIndex, heightIndex) > 0.0D0 .and. period >= cbstrPeriod(periodIndex, heightIndex))
        periodIndex = periodIndex + 1
    end do

    value = cbstrValue(periodIndex, heightIndex)
        
    return

end function CROSSBEACH_STR_LOOKUP

    
!> Returns the offshore sediment transport rate for a given wave height and period.
function OFFSHORE_STR_LOOKUP(height, period) result(value)

    implicit none
    integer, parameter :: double = kind(1D0)

    real(kind=double), intent(in) :: height !< Wave height
    real(kind=double), intent(in) :: period !< Wave period
    real(kind=double) :: value !< Offshore sediment transport rate
    
    integer :: heightIndex
    integer :: periodIndex

    heightIndex = 1
    periodIndex = 1

    do while (offstrWaveHeight(heightIndex) > 0.0D0 .and. height >= offstrWaveHeight(heightIndex))
        heightIndex = heightIndex + 1
    end do

    do while (offstrPeriod(periodIndex, heightIndex) > 0.0D0 .and. period >= offstrPeriod(periodIndex, heightIndex))
        periodIndex = periodIndex + 1
    end do

    value = offstrValue(periodIndex, heightIndex)
        
    return

end function OFFSHORE_STR_LOOKUP

    
!> Returns the groyne effect value for a given beach width.
!> Within the range zero (no effect) and 1 (groyne completely blocks drift)
function GROYNE_EFFECT_LOOKUP(width) result(value)

    implicit none
    integer, parameter :: double = kind(1D0)

    real(kind=double), intent(in) :: width !< Beach width
    real(kind=double) :: value !< groyne effect value
    
    integer :: widthIndex

    widthIndex = 1

    do while (groyneEffectBeachWidth(widthIndex) > 0.0D0 .and. width >= groyneEffectBeachWidth(widthIndex))
        widthIndex = widthIndex + 1
    end do

    value = groyneEffectValue(widthIndex)
        
    return

end function GROYNE_EFFECT_LOOKUP

    
!> Convert given string to upper case.
function UPPERCASE(string) result(value)
   character(*),intent(in) :: string
   character(len(string)) :: value 
   integer :: i,icode
   
   do i=1,len(string) 
      icode = iachar(string(i:i)) 
      if (icode >= 97 .and. icode <= 122) then
          value(i:i) = achar(icode-32)   
      else
          value(i:i) = achar(icode)
      end if
   end do
   
   return 
end function UPPERCASE 


!> Test if strings are equal ignoring case.
logical function EQUAL_NOCASE(string1,string2) 
   character(*),intent(in) :: string1,string2 
   equal_nocase = UPPERCASE(string1) == UPPERCASE(string2) 
   return 
end function EQUAL_NOCASE 


!> Deallocates arrays local to this module.
subroutine DEALLOCATE_PARAMETERS()

    integer :: alloc_error

    deallocate(fswp, depthWPmsl, cbstrWaveHeight, offstrWaveHeight, groyneEffectBeachWidth, &
           stat=alloc_error)
    if(alloc_error /= 0)then
      write(6,*)'*** Unexpected deallocation error'
      write(6,*)'*** for DEALLOCATE_PARAMETERS'
    endif

    return
    
end subroutine DEALLOCATE_PARAMETERS


!> Sets nextOutputYear to be the next year at which output is required.
subroutine SET_OUTPUT_YEAR(year)
    use local_data, only: nextOutputYear
    implicit none

    integer, intent(in) :: year !< Current year

    if (year >= highResOutStartYear) then
        nextOutputYear = year + highResOutTimestep
    else
        nextOutputYear = min(year + lowResOutTimestep, highResOutStartYear)
    end if

    nextOutputYear = min(nextOutputYear, endYear)

end subroutine SET_OUTPUT_YEAR
     
        
!> Takes a vertical element index in the range 1 : nce and returns element height relative to OD.
real(double) function relative_to_OD(element)
    use global_data, only: baseline_OD
    implicit none
    
    integer, parameter :: double = kind(1.0D0)
    
    integer, intent(in) :: element !< element index in range 1 to nce

    relative_to_OD = (element * elementHeight) + baseline_OD

end function relative_to_OD

    
!> Finds the factor to be applied to Wave Period values at the given year
subroutine GET_PERIOD_FACTOR(year, factor, haveFactor)
    implicit none

    integer, intent(in) :: year !< Current calendar year
    real(double), intent(out) :: factor !< Current factor
    logical, intent(out) :: haveFactor !< Current factor not 1.0

    integer :: index

    factor = 1.0D0
    haveFactor = .false.

    if (allocated(modifyPeriodYear)) then
        if (year >= modifyPeriodYear(1)) then
            do index = 1, ubound(modifyPeriodYear, 1)
              
                if (year < modifyPeriodYear(index)) then
                    exit
                end if
                
                if (year >= modifyPeriodYear(index)) then
                    factor = modifyPeriodFactor(index)
                    haveFactor = factor /= 1.0D0
                end if
            end do
        end if
    end if

end subroutine GET_PERIOD_FACTOR


!> Finds the factor to be applied to Wave Height values at the given year
subroutine GET_HEIGHT_FACTOR(year, factor, haveFactor)
    implicit none

    integer, intent(in) :: year !< Current calendar year
    real(double), intent(out) :: factor !< Current factor
    logical, intent(out) :: haveFactor !< Current factor not 1.0

    integer :: index

    factor = 1.0D0
    haveFactor = .false.

    if (allocated(modifyHeightYear)) then
        if (year >= modifyHeightYear(1)) then
            do index = 1, ubound(modifyHeightYear, 1)
              
                if (year < modifyHeightYear(index)) then
                    exit
                end if
                
                if (year >= modifyHeightYear(index)) then
                    factor = modifyHeightFactor(index)
                    haveFactor = factor /= 1.0D0
                end if
            end do
        end if
    end if

end subroutine GET_HEIGHT_FACTOR


!> Finds the change to be applied to Wave Angle values at the given year
subroutine GET_ANGLE_CHANGE(year, change, haveChange)
    implicit none

    integer, intent(in) :: year !< Current calendar year
    real(double), intent(out) :: change !< Current angle change
    logical, intent(out) :: haveChange !< Current angle change not 0.0

    integer :: index

    change = 0.0D0
    haveChange = .false.

    if (allocated(modifyAngleYear)) then
        if (year >= modifyAngleYear(1)) then
            do index = 1, ubound(modifyAngleYear, 1)
              
                if (year < modifyAngleYear(index)) then
                    exit
                end if
                
                if (year >= modifyAngleYear(index)) then
                    change = modifyAngleAmount(index)
                    haveChange = change /= 0.0D0
                end if
            end do
        end if
    end if

end subroutine GET_ANGLE_CHANGE


end module setup_module
