!**********************************************
!*       local_data.f95       
!**********************************************
!>    Declares all the local data
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


module local_data
    
    implicit none
    save
    
    integer, parameter :: double = kind(1D0)
    private double

!*** Local variables

    real(kind=double), allocatable, dimension(:) :: ann_av_bvol
    real(kind=double), allocatable, dimension(:) :: vell_berm
    real(kind=double), allocatable, dimension(:) :: vxs
    real(kind=double), allocatable, dimension(:) :: vellinga
    
    ! MW to do: Various variables (such as these) dealing with sea level rise will be changed at a later stage.
    real(kind=double), dimension(3) :: coeff_5, coeff_95
    
    ! MW to do: Various variables (such as these) used to save results will be changed at a later stage.
    real(kind=double), allocatable, dimension(:) :: initial_cliff_toe
    real(kind=double), allocatable, dimension(:, :) :: retreat
    
    real(kind=double) vell_y

    integer i,j
    integer :: wp    !< Wave point number
    integer :: lastWPsection
    integer year
    integer :: ierr
    integer wave_count !< Initialsed from the waveTideFileStart parameter.
        !< Determines the first record read from the wave file the first time through.
        !< Allows exploration of the effects of weather on simulations by varying waveTideFileStart.
    integer vell_top_e, vell_bott_e, nvsect
    integer month_counter
    character (20) :: arg
    logical :: use_random_seed !< For testing. If false, random number sequence is the same for all runs
    integer :: qPoint !< Q point index
    integer :: section !< Beach section index
    integer :: outputContourLevel_e !< Number of the platform element used for the rock contour and shore contour outputs

    logical :: feLogEnabled = .false. !< Enables logging used in SCAPE OpenMI engine
    logical :: feLogActive = .false. !< true if log file is open
    logical :: feLogAppend = .true. !< If true, append log message then close file
    character(256) :: feLogFile = '' !< FluidEarth log file name
    character(256) :: feLogString = '' !< Message buffer for FluidEarth logging call to feLog()
    integer, parameter :: feLogLimit = 50 !< number of messages that can be stored in feLogStrings
    integer :: feLogIndex = 0 !< index of next free message buffer within feLogStrings
    character(256), allocatable, dimension(:) :: feLogStrings
        !< buffers FluidEarth logging messages before logging file is opened

    character(256) :: logString = '' !< Message buffer for SCAPE logging call to SCAPE_LOG()
    
    integer :: logFileUnit !< FORTRAN IO unit number for SCAPE log file
    
    ! FORTRAN IO unit number for model output files
    integer :: saveSedFluxLeftUnit
    !< Sediment flux at left end of model
    integer :: savePotSedFluxLeftUnit
    !< Potential sediment flux at left end of model
    integer :: saveSedFluxRightUnit
    !< Sediment flux at right end of model
    integer :: savePotSedFluxRightUnit
    !< Potential sediment flux at right end of model
    integer :: saveRockContourUnit
    !< Rock contour
    integer :: saveShoreContourUnit
    !< Shore contour
    integer :: saveAnnBvolUnit
    !< Annual beach volume
    integer :: saveFinesVolUnit
    !< Fines volume
    integer :: saveSedTransAnnUnit
    !< Annual sediment transport
    integer :: savePotTransAnnUnit
    !< Potential annual sediment transport
    integer :: saveBeachAddAnnUnit
    !< Annual beach addition
    integer :: saveVellLevelsUnit
    !< Vellinga Levels
    integer :: saveSeawallActiveUnit
    !< Seawall active
    integer :: saveGroyneEffectUnit
    !< Groyne factors
    integer :: saveRockProfilesUnit
    !< Rock profiles
    integer :: saveBeachProfilesUnit
    !< Beach profiles
    
    ! FORTRAN IO unit number for input files
    integer :: wavesAndTidesFileUnit
    integer :: tempInputFileUnit
    integer :: clutUnit
    integer :: cglutUnit

    real(kind=double) :: time_start !< run start time in Julian Days

    ! Looking offshore, left end is at Q(nQpoints), right is at Q(1) and +ve Q is right to left

    real(kind=double) :: pot_left_sediment_flux !< Potential sediment flux at the left end of the model (m^3 / tide)
    real(kind=double) :: pot_right_sediment_flux !< Potential sediment flux at the right end of the model (m^3 / tide)
    real(kind=double) :: left_sediment_flux !< Sediment flux at the left end of the model (m^3 / tide)
    real(kind=double) :: right_sediment_flux !< Sediment flux at the right end of the model (m^3 / tide)
    logical :: have_left_sediment_input = .false.
    real(kind=double) :: left_sediment_input !< Sediment drift into the left end of the model (m^3 / tide)
    logical :: have_right_sediment_input = .false.
    real(kind=double) :: right_sediment_input !< Sediment drift into the right end of the model (m^3 / tide)

    real(kind=double), allocatable, dimension(:) :: sediment_influx !< Sediment influx at the model seaward boundary
        !< One value per section (m^3 / tide)
    
    character(256) :: FILE_PATH = "" !< Path used to prefix all file paths.
    !< Empty for RUN_SCAPE.EXE runs. May have a value for FluidEarth OpenMI runs.
    !< May be an absolute path or a relative path. If not empty, must end in a path separator character.
    !< Currennt working directory will be set by Pipistrelle (to the .chi file path) 
    !< for a SCAPE engine running in a composition.

    integer :: nextOutputYear !< next year at which output is required
    character(6) :: SCAPE_VERSION = '1.21'
    
    real(kind=double) :: periodFactor = 1.0D0 ! Wave period factor for current year
    logical :: havePeriodFactor = .false.
    real(kind=double) :: heightFactor = 1.0D0 ! Wave height factor for current year
    logical :: haveHeightFactor = .false.
    real(kind=double) :: angleChange = 0.0D0 ! Wave angle change for current year
    logical :: haveAngleChange = .false.

end module local_data