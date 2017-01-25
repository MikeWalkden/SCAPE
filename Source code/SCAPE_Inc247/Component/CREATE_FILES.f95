!**********************************************
!*       CREATE_FILES.f95       
!**********************************************
!
! Creates xml files describing the SCAPE OpenMI component.

! By default the SCAPE OpenMI component is captioned "SCAPE".
! The first command line argument, if given, is used to set the caption.

! This program generates the files <caption>.omi and <caption>.fcs which describe the SCAPE 
! OpenMI component to to the FluidEarth composition runner Pipistrelle.

! The output files are created from the files template.omi amd template.fcs by making
! substitutions for values that are particular to a given SCAPE model.
! Model information is taken from the file "setup.txt".
! Additional component information is taken from the file "component.txt".

! Model layout: Looking offshore, the origin of local coordinates is at the right-hand end of the 
! model baseline. The model baseline runs along the negative x-axis. Positive y-values increase 
! going offshore. Q points run from right to left with Q point 1 have an x-value of zero.
! Beach sections run from right to left between the Q points. The width of beach sections is 
! taken from the DX parameter in the setup.txt file.

! If translation of coordinates to "world coordinates" is required, give appropriate values to the 
! WORLDCOORDSID, LOCALORIGINX and LOCALORIGINY parameters in the component.txt file. 
! The angle of the model baseline is taken from the BASELINEANGLE parameter in the setup.txt file.

! The position of the model seaward boundary defaults to y-value of 1000 metres,
! but can be set using the SEAWARDBOUNDARY parameter in the component.txt file.
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

program main
    use setup_module, only: nSections, spinUpToYear, endYear, baselineAngle_radians, INITIALISE_PARAMETERS
    use component_parameters
    use exceptions
    use utilities_module
    use component_data
    use support

    implicit none

    character(len=256) :: cwd
#ifdef __GFORTRAN__
#else
    character(len=256) :: CURDIR@
#endif    
    character(keyLength), dimension(substitutions) :: keys ! keys to search for in the template files
    character(valueLength), dimension(substitutions) :: values ! values to substitute for the found keys
    real(kind=double) :: pi
    real(kind=double) :: radcon
    integer :: success_code
    logical :: file_exists
    integer :: iargc
    character(valueLength) :: arg = ''
    real(kind=double), dimension(3, 3) :: translate
    real(kind=double), dimension(3, 3) :: rotate
    real(kind=double), dimension(3, 3) :: transform

#ifdef __GFORTRAN__
    call getcwd(cwd)
#else
    cwd = CURDIR@()
#endif    

    print *, 'SCAPE .omi file creator'
    print *
    
    if (iargc() > 0) then
        call GETARG(1, arg)
    end if

    inquire(file='setup.txt', exist=file_exists)
    
    if (.not. file_exists) then
        call RAISE_EXCEPTION('File not found: setup.txt')
        stop
    end if
    
    pi = 4.0d0 * datan(1.0d0)
    radcon = pi / 180.0d0

    call INITIALISE_PARAMETERS('', radcon, success_code)
    if (success_code < 0) stop

    inquire(file='component.txt', exist=file_exists)
    
    if (.not. file_exists) then
        call RAISE_EXCEPTION('File not found: component.txt')
        stop
    end if
    
    call INITIALISE_COMPONENT_PARAMETERS(radcon, success_code)
    if (success_code < 0) stop

    ! Transformation matrix to translate world origin to local origin
    translate(1,1) = 1.0D0
    translate(2,1) = 0.0D0
    translate(3,1) = 0.0D0
    
    translate(1,2) = 0.0D0
    translate(2,2) = 1.0D0
    translate(3,2) = 0.0D0
    
    translate(1,3) = localOriginX
    translate(2,3) = localOriginY
    translate(3,3) = 1.0D0
    
    ! Transformation matrix to rotate world x-axis to local x-axis
    rotate(1,1) = sin(baselineAngle_radians)
    rotate(2,1) = cos(baselineAngle_radians)
    rotate(3,1) = 0.0D0
    
    rotate(1,2) = -cos(baselineAngle_radians)
    rotate(2,2) = sin(baselineAngle_radians)
    rotate(3,2) = 0.0D0
    
    rotate(1,3) = 0.0D0
    rotate(2,3) = 0.0D0
    rotate(3,3) = 1.0D0

    ! Transformation matric to convert local ccords to world coords
    transform = matmul(rotate, translate)
    
    keys(1) = '{{{START}}}' ! Time Horizon start as Modified Julian day
    write(values(1), fmt='(F10.0)') to_mjd(julian_day(1, 1, spinUpToYear, 0, 0, 0))
    
    keys(2) = '{{{END}}}' ! Time Horizon end as Modified Julian day
    write(values(2), fmt='(F10.0)') to_mjd(julian_day(1, 1, endYear, 0, 0, 0))

    keys(3) = '{{{YPOINTS_COUNT}}}'
    write(values(3), fmt='(I4)') nSections
    
    keys(4) = '{{{QPOINTS_COUNT}}}'
    write(values(4), fmt='(I4)') nSections + 1
    
    keys(5) = '{{{CAPTION}}}'
    values(5) = 'SCAPE'

    if (len_trim(arg) > 0) values(5) = arg
    
    keys(6) = '{{{DIRECTORY}}}'
    values(6) = cwd
    
    keys(7) = '{{{FCS_PATH}}}'
    values(7) = ''
    
    keys(8) = '{{{FLUIDEARTHPATH}}}'
    values(8) = fluidEarthPath
    
    keys(9) = '{{{ENGINEWRAPPERPATH}}}'
    values(9) = engineWrapperPath
    
    keys(10) = '{{{SEAWARD_BOUNDARY_COUNT}}}'
    write(values(10), fmt='(I4)') nSections
    
    keys(11) = '{{{BASELINE_SECTIONS_COUNT}}}'
    write(values(11), fmt='(I4)') nSections
    
    keys(12) = '{{{WORLD_COORDS_ID}}}'
    values(12) = worldCoordsId
    
    print *, 'Files generated '//trim(values(5))//'.omi and '//trim(values(5))//'.fcs'
    print *
    
    values(7) = trim(cwd)//'\'//trim(values(5))//'.fcs'

    call substituteFile(trim(templatesPath)//"template.omi", trim(values(5))//".omi", keys, values, transform)

    call substituteFile(trim(templatesPath)//"template.fcs", trim(values(5))//".fcs", keys, values, transform)

end program main
