!**********************************************
!*       vertical_grid.f95       
!**********************************************
!
!> Support functions added to support the SCAPE vertical grid enhancement.
!>
!> Base geology has a strength value and a sediment contant value and an extent defined by a Y - Z profile extending seaward to Z = 0.
!>
!> Layers above have a strength value and a sediment content value and an extent defined by a Y - Depth profile extending seaward to Y = YmaxH, 
!> where YmaxH is user-defined.
!>
!> Vertical grid constructed from the base geology and the layers extends horixontally from Y = 0 to Y = YmaxH.
!> The grid spacings are user-defined values, dy and dz.
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

module vertical_grid
use setup_module, only: double, dz, dy, nLayers

integer :: nHcells !< Number of horizontal cells

!< Variables used when reading profiles
real(kind=double), allocatable, dimension(:) :: temp_y_1 !< Temporary y value storage
real(kind=double), allocatable, dimension(:) :: temp_z_1 !< Temporary z value storage
real(kind=double), allocatable, dimension(:) :: temp_y_2 !< Temporary y value storage
real(kind=double), allocatable, dimension(:) :: temp_z_2 !< Temporary z value storage
integer :: active_temp_points !< Indicates which temp arrays are active (1 or 2)
integer :: index_limit !< Current dimension of temp arrays
integer :: next_index !< Next subscript in temp arrays

save

contains

!> Reads data giving layer properties.
!>
!> Reads, validates and interpolates the profiles for the base geology and for the layers above.
subroutine get_layers(success_code)
    use setup_module, only: layerDataFile, YmaxH
    use global_data, only: nce, base_strength, base_profile_start, base_profile_start2, base_profile_end, &
        layer_profile_start, layer_profile_start2, layer_profile_end, &
        base_sediment, base_profile, layer_count, layer_strength, layer_sediment, layer_profile, &
        nActiveLayers, nYsections, sedimentContent, baseline_OD
    use local_data, only: tempInputFileUnit, FILE_PATH
    use utilities_module
    use exceptions
    
    implicit none
    
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

    integer :: section ! Beach section index
    integer :: index ! Layer index
    integer :: iostat
    integer :: alloc_error

    baseline_OD = -500.0D0 ! Arbitrary large negative value to be updated as the base geology profiles are read in
    index = 1
    nHcells = nint(YmaxH / dy)

    allocate(base_profile_start(nYsections), base_profile_start2(nYsections), base_profile_end(nYsections), &
        layer_profile_start(nLayers, nYsections), layer_profile_start2(nLayers, nYsections), &
        layer_profile_end(nLayers, nYsections), &
        base_profile(0: nHcells, nYsections), base_strength(nYsections), base_sediment(nYsections), &
        layer_count(nYsections), layer_profile(0: nHcells, nLayers, nYsections), &
        layer_strength(nLayers, nYsections), layer_sediment(nLayers, nYsections), sedimentContent(nce, nYsections), &
        stat = alloc_error)
    success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'base_profile etc.')
    if (success_code < 0) return

    call init_profile_point_storage(success_code)
    if (success_code < 0) return
    
    call get_base_profile(success_code)
    if (success_code < 0) return

    success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, layerDataFile)
    if (success_code < 0) return

    do section = 1, nYsections
        read(tempInputFileUnit, fmt=*, iostat = iostat) layer_count(section)
        success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading '//trim(FILE_PATH)//trim(layerDataFile))
        if (success_code < 0) return

        do index = 1, layer_count(section)
            read(tempInputFileUnit, fmt=*, iostat = iostat) layer_strength(index, section), layer_sediment(index, section)
            success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading '//trim(FILE_PATH)//trim(layerDataFile))
            if (success_code < 0) return
    
            call read_profile(tempInputFileUnit, trim(FILE_PATH)//trim(layerDataFile), &
                layer_profile_start(index, section), layer_profile_start2(index, section), layer_profile_end(index, section), &
                .false., success_code)
            if (success_code < 0) return
        end do
    end do
    
    close(tempInputFileUnit)

    call save_profile_points(success_code)
    if (success_code < 0) return

    do section = 1, nYsections
        call interpolate_profile(base_profile_start(section), base_profile_end(section), &
            base_profile(:, section), success_code)

        do index = 1, layer_count(section)
            call interpolate_profile(layer_profile_start(index, section), layer_profile_end(index, section), &
                layer_profile(:, index, section), success_code)
            if (success_code < 0) return
        end do
    end do
        
end subroutine get_layers


!> Reads base geology or layer profile from a file of Y - Z or Y - Depth values. 
!>
!> Validates the profile and adds the profile points to temporary storage.
subroutine read_profile(unitNo, filename, profile_start, profile_start2, profile_end, base_geology, success_code)
    use setup_module, only: YmaxH
    use global_data, only: baseLine_OD
    use local_data, only: logString
    use utilities_module
    use exceptions
    
    implicit none
    
    integer, intent(in) :: unitNo !< File already open
    character(*), intent(in) :: filename !< input file name
    integer, intent(out) :: profile_start !< Start index within profile_y and _profile_z
    integer, intent(out) :: profile_start2 !< Index within profile_y and _profile_z of the last point before YmaxH
    integer, intent(out) :: profile_end !< End index within profile_y and _profile_z
    logical, intent(in) :: base_geology !< true if reading the base geology profile. false if reading a layer profile
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

    real(kind=double) :: yValue
    real(kind=double) :: zValue !< user profile value relative to Ordinance datum or layer depth value
    real(kind=double) :: lastyValue
    real(kind=double) :: lastzValue
    logical :: firstyValue
    integer :: pointsCount
    integer :: pointIndex
    integer :: iostat

    read(unitNo, fmt=*, iostat = iostat) pointsCount
    success_code = SUCCESS_CODE_FROM_IOSTAT(unitNo, iostat, 'Reading '//trim(filename))
    if (success_code < 0) return

    lastyValue = 0.0D0
    lastzValue = 50000.0D0
    firstyValue = .true.
    profile_start = next_index
    profile_start2 = profile_start

    do pointIndex = 1, pointsCount
        read(unitNo, fmt=*, iostat = iostat) yValue, zValue
        success_code = SUCCESS_CODE_FROM_IOSTAT(unitNo, iostat, 'Reading '//trim(filename))
        if (success_code < 0) return

        if (yValue < lastyValue) then
            write(logString, fmt='(A,F12.4,A)') 'Y value ', yValue, ' out of order in '//trim(filename)
            call RAISE_EXCEPTION(logString)
            success_code = -1
            return
        end if

        if (base_geology) then
            if (yValue > YmaxH .and. zValue > lastzValue) then
                write(logString, fmt='(A,F12.4,A)') 'Profile trough not allowed at Y value ', &
                    yValue, ' in '//trim(filename)
                call RAISE_EXCEPTION(logString)
                success_code = -1
                return
            end if
        else
            if (yValue > YmaxH) then
                write(logString, fmt='(A,F12.4,A)') 'Y value ', yValue, ' out of range in '//trim(filename)
                call RAISE_EXCEPTION(logString)
                success_code = -1
                return
            end if
    
            if (zValue < 0.0D0) then
                write(logString, fmt='(A,F12.4,A)') &
                    'Negative depth value not allowed ', zValue, ' in '//trim(filename)
                call RAISE_EXCEPTION(logString)
                success_code = -1
                return
            end if
        end if

        if (yValue <= YmaxH) then
            profile_start2 = next_index
        end if

        call add_profile_point(yValue, zValue, success_code)
        if (success_code < 0) return

        if (firstyValue) then
            if (yValue > 0.0D0) then
                call RAISE_EXCEPTION('Profile must start at Y = 0.0 in '//trim(filename))
                success_code = -1
                return
            end if
    
            firstyValue = .false.
        end if
        
        lastyValue = yValue
        lastzValue = zValue
    end do ! pointIndex = 1, pointsCount
    
    profile_end = next_index - 1
    
    if (lastyValue < YmaxH) then
        write(logString, fmt='(A,F12.4,A)') 'Profile must extend to ', YmaxH, ' in '//trim(filename)
        call RAISE_EXCEPTION(logString)
        success_code = -1
        return
    end if
    
    if (base_geology .and. lastzValue > baseLine_OD) then
        baseLine_OD = lastzValue
    end if

end subroutine read_profile


!> Interpolates between points in a base geology or layer profile to set values in the vertical grid.
subroutine interpolate_profile(profile_start, profile_end, profile, success_code)
    use setup_module, only: YmaxH
    use global_data, only: profile_y, profile_z
    use local_data, only: logString
    use utilities_module
    use exceptions
    
    implicit none
    
    integer, intent(in) :: profile_start !< Start index within profile_y and _profile_z
    integer, intent(in) :: profile_end !< End index within profile_y and _profile_z
    real(kind=double), dimension(0: ), intent(inout) :: profile !< profile values at dy horizontal spacing
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

    real(kind=double) :: lastyValue
    real(kind=double) :: lastzValue
    logical :: firstyValue
    integer :: pointIndex
    integer :: profileIndex
    integer :: nextProfileIndex
    integer :: profileIndexLimit
    
    success_code = 0

    lastyValue = 0.0D0
    lastzValue = 50000.0D0
    firstyValue = .true.
    profileIndex = 0
    profileIndexLimit = size(profile, 1)

    do pointIndex = profile_start, profile_end
        if (firstyValue) then
            profile(profileIndex) = profile_z(pointIndex)
            profileIndex = profileIndex + 1
            firstyValue = .false.
        else
            nextProfileIndex = int(profile_y(pointIndex) / dy)
            
!            if (nextProfileIndex >= profileIndex) then
                do while(profileIndex <= nextProfileIndex .and. profileIndex < profileIndexLimit)
                    call interpolate_profile_value(profile, profileIndex, lastyValue, lastzValue, &
                        profile_y(pointIndex), profile_z(pointIndex))
                    profileIndex = profileIndex + 1
                end do
!            end if
        end if
        
        lastyValue = profile_y(pointIndex)
        lastzValue = profile_z(pointIndex)
    end do ! pointIndex = profile_start, profile_end
    
end subroutine interpolate_profile


!> Interpolate values from a user base geology profile or layer profile
subroutine interpolate_profile_value(profile, index, y1, z1, y2, z2)
    implicit none

    real(kind=double), dimension(0:), intent(inout) :: profile !< profile values at dy horizontal spacing
    integer, intent(in) :: index !< profile vector index
    real(kind=double), intent(in) :: y1 !< first point y value
    real(kind=double), intent(in) :: z1 !< first point z value or depth
    real(kind=double), intent(in) :: y2 !< second point y value
    real(kind=double), intent(in) :: z2 !< second point z value or depth

    real(kind=double) :: yValue
    real(kind=double) :: factor

    yValue = dy * index
    factor = (z2 - z1) / (y2 - y1)
    profile(index) = z1 + factor * (yValue - y1)

end subroutine interpolate_profile_value

                   
!> Interpolate value y value along a y-z line
function interpolate_y_value(z, y1, z1, y2, z2) result(yValue)
    implicit none

    real(kind=double), intent(in) :: z !< z value at which interpolation is required
    real(kind=double), intent(in) :: y1 !< first point y value
    real(kind=double), intent(in) :: z1 !< first point z value or depth
    real(kind=double), intent(in) :: y2 !< second point y value
    real(kind=double), intent(in) :: z2 !< second point z value or depth
    real(kind=double) :: yValue !> interpolated yValue result

    real(kind=double) :: factor

    factor = (y2 - y1) / (z2 - z1)
    yValue = y1 + factor * (z - z1)

end function interpolate_y_value

                    
!> Interpolate value z value along a y-z line
function interpolate_z_value(y, y1, z1, y2, z2) result(zValue)
    implicit none

    real(kind=double), intent(in) :: y !< y value at which interpolation is required
    real(kind=double), intent(in) :: y1 !< first point y value
    real(kind=double), intent(in) :: z1 !< first point z value or depth
    real(kind=double), intent(in) :: y2 !< second point y value
    real(kind=double), intent(in) :: z2 !< second point z value or depth
    real(kind=double) :: zValue !> interpolated yValue result

    real(kind=double) :: factor

    factor = (z2 - z1) / (y2 - y1)
    zValue = z1 + factor * (y - y1)

end function interpolate_z_value

                    
!> Reads base geology data and profile from a file of Y - Z values. 
!>
!> Call read_profile to validates the profile and to get and store the profile points.
subroutine get_base_profile(success_code)
    use setup_module, only: baseGeologyFile, YmaxH
    use global_data, only: base_profile_start, base_profile_start2, base_profile_end, &
        base_strength, base_sediment, base_profile, nYsections
    use local_data, only: tempInputFileUnit, FILE_PATH, logString
    use utilities_module
    use exceptions

    implicit none
    
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

    integer :: section !< Beach section index
    integer :: index
    integer :: iostat

    success_code = do_open_text_input_file(tempInputFileUnit, FILE_PATH, baseGeologyFile)
    if (success_code < 0) return

    do section = 1, nYsections
        read(tempInputFileUnit, fmt=*, iostat = iostat) base_strength(section), base_sediment(section)
        success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Reading '//trim(FILE_PATH)//trim(baseGeologyFile))
        if (success_code < 0) return

        call read_profile(tempInputFileUnit, trim(FILE_PATH)//trim(baseGeologyFile), &
            base_profile_start(section), base_profile_start2(section), base_profile_end(section), .true., success_code)
        if (success_code < 0) return
    end do

    close(tempInputFileUnit, iostat = iostat)
    success_code = SUCCESS_CODE_FROM_IOSTAT(tempInputFileUnit, iostat, 'Closing base geology profile')
    if (success_code < 0) return

end subroutine get_base_profile


!> Initialises cliffy values from the base geology and layer profiles.
subroutine get_y_values_from_profiles(nYsections, nce, cliffy)
    use global_data, only: base_profile, base_profile_end, layer_count, layer_profile, profile_y, baseLine_OD
    use setup_module, only: nLayers, dy, dz
    use local_data, only: logString
    use utilities_module, only: SCAPE_LOG
    implicit none

    integer, intent(in) :: nYsections !< number of beach sections
    integer, intent(in) :: nce !< number of cliff elements
    real(kind=double), intent(inout), dimension(nce, nYsections) :: cliffy !< discretised cliff y value

    integer :: section !< Beach section index
    integer :: yIndex !< Horizontal index
    real(kind=double) :: leftY ! Y value at left of cell
    real(kind=double) :: rightY ! Y value at right of cell
    real(kind=double) :: leftZ ! Z value at left of cell at top of layer stack
    real(kind=double) :: rightZ ! Z value at right of cell at top of layer stack
    integer :: stackIndex ! Index of layer under consideration
    real(kind=double) :: last_base_profile_point_Y ! Y value of the last profile point for the beach section
    real(kind=double) :: change
    real(kind=double) :: leftZbase
    real(kind=double) :: zAbove
    real(kind=double) :: zBelow
    real(kind=double), parameter :: toler = 0.000001D0
    integer :: leftZindex
    
    cliffy = 0.0D0

    do section = 1, nYsections
        yIndex = 0
        last_base_profile_point_Y = profile_y(base_profile_end(section))
        
        do
            leftY = yIndex * dy
            rightY = leftY + dy

            if (leftY > last_base_profile_point_Y .or. rightY > last_base_profile_point_Y) then
                exit ! Profile done
            end if
        
            leftZ = base_profile_value(yIndex, section)
            rightZ = base_profile_value(yIndex + 1, section)
        
            do stackIndex = 1, layer_count(section)
                leftZ = leftZ + layer_profile_value(yIndex, stackIndex, section)
                rightZ = rightZ + layer_profile_value(yIndex + 1, stackIndex, section)
            end do
            
            change = rightZ - leftZ

            if (abs(change) > toler) then
                leftZbase = leftZ - baseLine_OD
                leftZindex = int((leftZbase + toler) / dz)

                if (leftZindex < 1) then
                    exit ! Gone below model baseline
                end if
            
                if (change > 0.0D0) then
                    ! Consolidated profile slopes up
                    leftZindex = leftZindex + 1
                    zAbove = leftZindex * dz + baseLine_OD

                    if (zAbove < leftZ) then
                        if ((leftZ - zAbove) < toler) then
                            zAbove = leftZ
                        else
                            write(6, *)'*** Error : zAbove out of range ', zAbove, leftY, leftZ, rightY, rightZ, section
                            stop
                        end if
                    end if

                    do while (zAbove <= rightZ)
                        cliffy(leftZindex, section) = interpolate_y_value(zAbove, leftY, leftZ, rightY, rightZ)
                        leftZindex = leftZindex + 1
                        zAbove = zAbove + dz
                    end do
                else
                    ! Consolidated profile slopes down
                    zBelow = leftZindex * dz + baseLine_OD

                    if (zBelow > leftZ) then
                        if ((zBelow - leftZ) < toler) then
                            zBelow = leftZ
                        else
                            write(6, *)'*** Error : zBelow out of range ', zBelow, leftY, leftZ, rightY, rightZ, section
                            stop
                        end if
                    end if

                    do while (zBelow >= rightZ)
                        cliffy(leftZindex, section) = interpolate_y_value(zBelow, leftY, leftZ, rightY, rightZ)
                        leftZindex = leftZindex - 1
                        zBelow = zBelow - dz
                    end do
                end if
            end if

            yIndex = yIndex + 1 ! Move to next cell further offshore
        end do
    end do ! section = 1, nYsections
    
end subroutine get_y_values_from_profiles


!> Cut the base geology and layer profiles to match current cliffy values.
!>
!> Only modifies the profiles in the range Y = 0 to YmaxH
subroutine cut_layer_profiles(nYsections, nce, cliffy)
    use global_data, only: base_profile, layer_count, layer_profile, msl_m, msl_OD
    use setup_module, only: nLayers, dy, dz
    implicit none

    integer, intent(in) :: nYsections !< number of beach sections
    integer, intent(in) :: nce !< number of cliff elements
    real(kind=double), intent(inout), dimension(nce, nYsections) :: cliffy !< discretised cliff y value

    real(kind=double), dimension(0: nHcells) :: profile ! Z value of profile from current cliffy values
    integer :: section !< Beach section index
    integer :: yIndex !< Horizontal index
    integer :: zIndex !< Vertical index
    real(kind=double) :: zValue ! Z value up the layer stack
    real(kind=double) :: midZ ! Z value midway up cell
    real(kind=double) :: previousY ! Y value at previous Y values point
    real(kind=double) :: previousZ ! Z valueat previous Y values point
    integer :: stackIndex !< Index of layer under consideration
    logical :: firstValueSet
    
    do section = 1, nYsections
        yIndex = 0
        firstValueSet = .false.
    
        ! Construct a profile to match the current Y values
        do zIndex = nce, 1, -1
            midZ = dz * (real(zIndex - 1, kind=double) + 0.5D0) - msl_m + msl_OD
    
            if (.not. firstValueSet .and. cliffy(zIndex, section) > 0.0D0) then
                profile(yIndex) = midZ
                firstValueSet = .true.
                yIndex = yIndex + 1
                previousY = 0.0D0
                previousZ = midZ
            end if
           
            if (firstValueSet) then
                do while ((yIndex * dy) <= cliffy(zIndex, section) .and. yIndex <= nHcells)
                    profile(yIndex) = interpolate_z_value(yIndex * dy, previousY, previousZ, cliffy(zIndex, section), midZ)
                    yIndex = yIndex + 1
                end do
            end if
                
            previousY = cliffy(zIndex, section)
            previousZ = midZ
        end do ! zIndex = nce, 1, -1

        if (yIndex < nHcells + 1) then ! extend the profile to YmaxH if needed
            profile(yIndex: nHcells) = previousZ
        end if
    
        do yIndex = 0, nHcells
            if (base_profile(yIndex, section) > profile(yIndex)) then
                base_profile(yIndex, section) = profile(yIndex)
            end if
    
            zValue = base_profile(yIndex, section)
    
            do stackIndex = 1, layer_count(section)
                if (zValue + layer_profile(yIndex, stackIndex, section) > profile(yIndex)) then
                    layer_profile(yIndex, stackIndex, section) = profile(yIndex) - zValue
                end if
    
                zValue = zValue + base_profile(yIndex, section)
            end do ! stackIndex = 1, layer_count(section)
        end do ! yIndex = 0, nHcells
    end do ! section = 1, nYsections
    
end subroutine cut_layer_profiles


!> Add a new layer and add material to it
subroutine add_new_layer(section, strength, sediment, volume, thickness, yPosition, success_code)
    use global_data, only: layer_count, layer_profile, layer_strength, layer_sediment
    use local_data, only: logString
    use setup_module, only: dx, dy
    use exceptions
    implicit none

    integer, intent(in) :: section !< Beach section index
    real(kind=double), intent(in) :: strength !< Strength of layer material
    real(kind=double), intent(in) :: sediment !< Sediment content of layer
    real(kind=double), intent(in) :: volume !< Volume of material to add to layer
    real(kind=double), intent(in) :: thickness !< Target maximum thickness of added material
    real(kind=double), intent(in) :: yPosition !< Target y position at middle of added material
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

    integer :: nLayer !< Layer index

    if (layer_count(section) == nLayers) then
        write(logString,*) 'Too many layers at beach section ', section
        call RAISE_EXCEPTION(logString)
        success_code = -1
        return
    end if

    nLayer = layer_count(section) + 1
    layer_count(section) = nLayer

    layer_strength(nLayer, section) = strength
    layer_sediment(nLayer, section) = sediment
    layer_profile(:, nLayer, section) = 0.0D0

    call add_material_to_layer(nLayer, section, volume, thickness, yPosition, success_code)

end subroutine add_new_layer


!> Add material to an existing layer.
!>
!> Calculates the number of cells required to achieve the required volume as a layer of the required thickness.
!> Centres the cells near the required yPosition. Adds a two cell wide taper to the layer at each end.
!>
!> The requested volume is always achieved exactly. The position and thickness arguments are treated as targets.
subroutine add_material_to_layer(nLayer, section, volume, thickness, yPosition, success_code)
    use global_data, only: layer_count, layer_profile
    use local_data, only: logString
    use setup_module, only: dx, dy
    use exceptions
    implicit none

    integer, intent(in) :: nLayer !< Layer index
    integer, intent(in) :: section !< Beach section index
    real(kind=double), intent(in) :: volume !< Volume of material to add to layer
    real(kind=double), intent(in) :: thickness !< Target maximum thickness of added material
    real(kind=double), intent(in) :: yPosition !< Target y position at middle of added material
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

    integer :: yIndex !< Y value index
    integer :: startIndex !< Y value index of first cell taking material
    integer :: endIndex !< Y value index of last cell taking material
    integer :: nCells !< Number of cells required to take the material, not including the startIndex and endIndex cells
    real(kind=double) :: area !< Area to be added to the profile for the layer
    real(kind=double) :: actualThickness !< Actual layer thickness, calculated to get the area right after the number of cells has been determined

    if (nLayer < 1 .or. nLayer > layer_count(section)) then
        write(logString,*) 'Layer index ', nLayer, ' out of range at beach section ', section
        call RAISE_EXCEPTION(logString)
        success_code = -1
        return
    end if

    if (volume <= 0.0D0) then ! No work to do
        return
    end if

    area = volume / dx
    nCells = int(area / (thickness * dy))

    if (nCells == 0) nCells = 1 ! Must have at least one cell
    
    yIndex = nint(yPosition / dy)
    startIndex = yIndex - (nCells / 2)

    if (mod(nCells, 2) == 1) then
        startIndex = startIndex - 1
    end if
    
    endIndex = yIndex + (nCells / 2) + 1

    do while (startIndex < 0) ! Make sure startIndex is within range
        startIndex = startIndex + 1

        if (endIndex < nHcells) then
            endIndex = endIndex + 1
        end if
    end do
    
    do while (endIndex > nHcells) ! Make sure endIndex is within range
        endIndex = endIndex - 1

        if (startIndex > 0) then
            startIndex = startIndex - 1
        end if
    end do

    nCells = endIndex - startIndex - 1 ! Equivalent number of full height cells 
    actualThickness = area / (nCells * dy)
    ! material added to startIndex cell tapers from 0 to half actuallThickness
    layer_profile(startIndex + 1, nLayer, section) = layer_profile(startIndex + 1, nLayer, section) + 0.5 * actualThickness 
    ! material added to endIndex cell tapers from half actuallThickness to 0
    layer_profile(endIndex, nLayer, section) = layer_profile(endIndex, nLayer, section) + 0.5 * actualThickness

    do yIndex = startIndex + 2, endIndex - 1
        ! material added to cells in the middle has depth actuallThickness
        layer_profile(yIndex, nLayer, section) = layer_profile(yIndex, nLayer, section) + actualThickness
    end do

end subroutine add_material_to_layer


!> Get the strength and sediment content values for a cell.
subroutine get_cell_values(zValue, yIndex, section, strength, sediment)
    use global_data, only: base_strength, base_sediment, layer_strength, layer_sediment
    use setup_module, only: nLayers
    use exceptions
    implicit none

    real(kind=double), intent(in) :: zValue !< Z value relatice to OD of a point within the cell
    integer, intent(in) :: yIndex !< Horizontal index to left of cell
    integer, intent(in) :: section !< Beach section index
    real(kind=double), intent(out) :: strength !< Strength value for cell
    real(kind=double), intent(out) :: sediment !< Secidment content value for cell

    integer :: dominantLayer

    if (yIndex < 0) then
        call RAISE_EXCEPTION('Negative Y Value index in HA grid')
    end if
    
    dominantLayer = get_cell_dominant_layer(zValue, yIndex, section)
    
    select case(dominantLayer)
        case(-1) ! Empty cell
            strength = 0.0D0
            sediment = 0.0D0
        case(0) ! Base geology
            strength = base_strength(section)
            sediment = base_sediment(section)
        case default ! Layer
            strength = layer_strength(dominantLayer, section)
            sediment = layer_sediment(dominantLayer, section)
    end select

end subroutine get_cell_values


!> Returns an int indicating whether the HA grid cell of interest is empty, or which of the base geology layer 
!> or one of the layers above covers most of the sell area.
function get_cell_dominant_layer(zValue, yIndex, section) result(layerSelected)
    use global_data, only: base_profile, layer_count, layer_profile, baseLine_OD
    use setup_module, only: nLayers
    implicit none

    real(kind=double), intent(in) :: zValue !< Z value relatice to OD of a point within the cell
    integer, intent(in) :: yIndex !< Horizontal index
    integer, intent(in) :: section !< Beach section index
    integer :: layerSelected !< Result: -1 indicates the cell is empty; 0 indicates the base geology layer predominates;
        !< Greater than 0 indicates which of the upper layers covers most of the cell area.

    real(kind=double) :: lowZ ! Z value at bottom of cell of interest
    real(kind=double) :: highZ ! Z value at top of cell of interest
    real(kind=double) :: leftY ! Y value at left of cell of interest
    real(kind=double) :: rightY ! Y value at right of cell of interest
    real(kind=double) :: cellArea ! Area of a cell
    real(kind=double) :: areaCovered ! Area of the cell of interest covered by layers found so far
    real(kind=double) :: areaBelow ! Area of the cell of interest below the layer currently under consideration
    real(kind=double) :: areaThisLayer ! Area of the cell of interest occupied by the layer under consideration
    real(kind=double) :: leftZ ! Z value at top of layer under consideration at left of cell of interest
    real(kind=double) :: rightZ ! Z value at top of layer under consideration at right of cell of interest
    real(kind=double) :: yInterpolated
    integer :: stackIndex !< Index of layer under consideration
    real(kind=double) :: selectedLayerArea !< Area of the cell covered by the selected layer

    leftY = yIndex * dy
    rightY = leftY + dy
    lowZ = int((zValue - baseLine_OD) / dz) * dz + baseLine_OD
    highZ = lowZ + dz

    cellArea = dy * dz
    stackIndex = 0
    areaCovered = 0.0D0
    layerSelected = -1
    selectedLayerArea = 0.0D0
    leftZ = base_profile_value(yIndex, section)
    rightZ = base_profile_value(yIndex + 1, section)

    do while (areaCovered < cellArea)
        ! Calculate the area of the cell of interest below the layer currently under consideration
        if (leftZ < lowZ) then
            if (rightZ < lowZ) then
                areaBelow = 0.0D0
            else if (rightZ >= lowZ .and. rightZ <= highZ) then
                areaBelow = 0.5D0 * (rightY - interpolate_y_value(lowZ, leftY, leftZ, rightY, rightZ)) * (rightZ - lowZ)
            else ! rightZ > highZ
                areaBelow = 0.5D0 * (interpolate_y_value(lowZ, leftY, leftZ, rightY, rightZ) + &
                    interpolate_y_value(highZ, leftY, leftZ, rightY, rightZ)) * dz
            end if
        else if (leftZ >= lowZ .and. leftZ <= highZ) then
            if (rightZ < lowZ) then
                areaBelow = 0.5D0 * interpolate_y_value(lowZ, leftY, leftZ, rightY, rightZ) * (leftZ - lowZ)
            else if (rightZ >= lowZ .and. rightZ <= highZ) then
                areaBelow = 0.5D0 * (leftZ - lowZ + rightZ - lowZ) * dy
            else ! rightZ > highZ
                yInterpolated = interpolate_y_value(highZ, leftY, leftZ, rightY, rightZ)
                areaBelow = 0.5D0 * (leftZ - lowZ + dz) * (yInterpolated - leftY) + dz * (rightY - yInterpolated)
            end if
        else ! leftZ > highZ
            if (rightZ < lowZ) then
                areaBelow = 0.5D0 * (interpolate_y_value(lowZ, leftY, leftZ, rightY, rightZ) + &
                    interpolate_y_value(highZ, leftY, leftZ, rightY, rightZ)) * dz
            else if (rightZ >= lowZ .and. rightZ <= highZ) then
                yInterpolated = interpolate_y_value(highZ, leftY, leftZ, rightY, rightZ)
                areaBelow = 0.5D0 * (rightZ - lowZ + dz) * (rightY - yInterpolated) + dz * (yInterpolated - leftY)
            else ! rightZ > highZ
                areaBelow = dy * dz
            end if
        end if

        areaThisLayer = areaBelow - areaCovered
        
        if (areaThisLayer > selectedLayerArea) then
            selectedLayerArea = areaThisLayer
            layerSelected = stackIndex

            if (selectedLayerArea >= 0.5D0 * cellArea) exit
        end if

        areaCovered = areaBelow
        
        stackIndex = stackIndex + 1

        if (stackIndex > layer_count(section) .or. yIndex + 1 > nHcells) exit
        
        leftZ = leftZ + layer_profile_value(yIndex, stackIndex, section)
        rightZ = rightZ + layer_profile_value(yIndex + 1, stackIndex, section)
    end do

end function get_cell_dominant_layer


!> Returns the base profile value at a cell boundary
real(kind=double) function base_profile_value(yIndex, section)
    use global_data, only: base_profile, base_profile_start2, base_profile_end, profile_y, profile_z
    use setup_module, only: dy
    use local_data, only: logString
    use exceptions
    implicit none

    integer, intent(in) :: yIndex !< Horizontal index
    integer, intent(in) :: section !< Beach section index

    real(kind=double) :: yValue
    integer :: pointIndex

    if (yIndex <= nHcells) then ! Within range set by YmaxH
        base_profile_value = base_profile(yIndex, section)
        return
    end if

    yValue = yIndex * dy

    do pointIndex = base_profile_start2(section), base_profile_end(section) - 1
        if (profile_y(pointIndex) <= yValue .and. profile_y(pointIndex + 1) >= yValue) then
            base_profile_value = interpolate_z_value(yValue, profile_y(pointIndex), profile_z(pointIndex), &
                profile_y(pointIndex + 1), profile_z(pointIndex + 1))
            return
        end if
    end do

    base_profile_value = 0.0D0

end function base_profile_value


!> Returns a layer profile value at a cell boundary
real(kind=double) function layer_profile_value(yIndex, layer, section)
    use global_data, only: layer_profile, layer_profile_start2, layer_profile_end, profile_y, profile_z
    use setup_module, only: dy
    use local_data, only: logString
    use exceptions
    implicit none

    integer, intent(in) :: yIndex !< Horizontal index
    integer, intent(in) :: layer !< Layer index
    integer, intent(in) :: section !< Beach section index

    real(kind=double) :: yValue
    integer :: pointIndex

    if (yIndex <= nHcells) then ! Within range set by YmaxH
        layer_profile_value = layer_profile(yIndex, layer, section)
        return
    end if

    yValue = yIndex * dy

    do pointIndex = layer_profile_start2(layer, section), layer_profile_end(layer, section) - 1
        if (profile_y(pointIndex) <= yValue .and. profile_y(pointIndex + 1) >= yValue) then
            layer_profile_value = interpolate_z_value(yValue, profile_y(pointIndex), profile_z(pointIndex), &
                profile_y(pointIndex + 1), profile_z(pointIndex + 1))
            return
        end if
    end do

    layer_profile_value = 0.0D0

end function layer_profile_value


!> Initialise the variables and arrays used when reading profiles.
subroutine init_profile_point_storage(success_code)
    use exceptions
    
    implicit none
    
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning
    
    integer :: alloc_error
    
    index_limit = 5000
    
    allocate(temp_y_1(index_limit), temp_z_1(index_limit), stat=alloc_error)
    success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'temp_y_1, temp_z_1')
    if (success_code < 0) return
        
    active_temp_points = 1
    next_index = 1
    
end subroutine init_profile_point_storage


!> Add the next point to the current base geology or layer profile.
subroutine add_profile_point(y, z, success_code)
    use exceptions
    
    implicit none
    
    real(kind=double), intent(in) :: y !< Y value
    real(kind=double), intent(in) :: z !< Z value
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning
    
    integer :: old_limit
    integer :: alloc_error
    
    success_code = 0
    
    if (next_index > index_limit) then ! Temporary storage is full. Allocate alternate storage twice the size and copy values across
        old_limit = index_limit
        index_limit = 2 * index_limit
        
        select case(active_temp_points)
        case(1)
            allocate(temp_y_2(index_limit), temp_z_2(index_limit), stat=alloc_error)
            success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'temp_y_2, temp_z_2')
            if (success_code < 0) return
                
            active_temp_points = 2
            temp_y_2(1: old_limit) = temp_y_1
            temp_z_2(1: old_limit) = temp_z_1
            deallocate(temp_y_1, temp_z_1)
        case(2)
            allocate(temp_y_1(index_limit), temp_z_1(index_limit), stat=alloc_error)
            success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'temp_y_1, temp_z_1')
            if (success_code < 0) return
                
            active_temp_points = 1
            temp_y_1(1: old_limit) = temp_y_2
            temp_z_1(1: old_limit) = temp_z_2
            deallocate(temp_y_2, temp_z_2)
        end select
    end if
    
    if (success_code < 0) return

    select case(active_temp_points)
    case(1)
        temp_y_1(next_index) = y
        temp_z_1(next_index) = z
    case(2)
        temp_y_2(next_index) = y
        temp_z_2(next_index) = z
    end select
    
    next_index = next_index + 1
    
end subroutine add_profile_point


!> Called after all the profile points have been read in.
!>
!> Allocates permanent storage for profile points and copies the points.
!>
!> Deallocates temporary storage.
subroutine save_profile_points(success_code)
    use global_data, only: profile_y, profile_z
    use exceptions

    implicit none
    
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning
    
    integer :: alloc_error
    
    allocate(profile_y(next_index - 1), profile_z(next_index - 1), stat=alloc_error)
    success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'profile_y, profile_z')
    if (success_code < 0) return
        
    select case(active_temp_points)
    case(1)
        profile_y = temp_y_1(1: next_index - 1)
        profile_z = temp_z_1(1: next_index - 1)
        deallocate(temp_y_1, temp_z_1)
    case(2)
        profile_y = temp_y_2(1: next_index - 1)
        profile_z = temp_z_2(1: next_index - 1)
        deallocate(temp_y_2, temp_z_2)
    end select
    
end subroutine save_profile_points


!> Deallocates vertical grid storage
subroutine finish_VA_grid()
    use global_data, only: base_strength, base_profile_start, base_profile_start2, base_profile_end, &
        layer_profile_start, layer_profile_start2, layer_profile_end, &
        base_sediment, base_profile, layer_count, layer_strength, layer_sediment, layer_profile, profile_y, profile_z, &
        sedimentContent
    use exceptions
    implicit none

    integer :: alloc_error
    integer :: success_code !< 0: ok, -ive fatal error, +ive warning
    
    deallocate(base_profile_start, base_profile_start2, base_profile_end, &
        layer_profile_start, layer_profile_start2, layer_profile_end, &
        base_profile, base_strength, base_sediment, &
        layer_count, layer_profile, layer_strength, layer_sediment, &
        profile_y, profile_z, sedimentContent, stat = alloc_error)
    success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'deallocate base_profile etc.')
    if (success_code < 0) return

end subroutine finish_VA_grid


end module vertical_grid
