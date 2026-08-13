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

module support

use component_data

implicit none

contains

!> Replace the first occurrence of key with value
!> Result is true if a replacement was made
logical function replace(line, key, value)
    implicit none
    character(*), intent(inout) :: line
    character(*), intent(in) :: key
    character(*), intent(in) :: value

    integer :: test

    test = index(line, key)
    replace = test > 0

    if (replace) then
        if (test == 1) then
            line = value//line(test + len(key):)
        else
            line = line(1 : test - 1)//value//line(test + len(key):)
        end if
    end if

end function replace


!> Substitute all keys in the given line
subroutine substituteLine(line, keys, values)
    implicit none

    character(*), intent(inout) :: line
    character(keyLength), dimension(substitutions) :: keys
    character(valueLength), dimension(substitutions) :: values

    integer :: sub

    outer: do

        inner: do sub = 1, size(keys, 1)
            if (replace(line, trim(keys(sub)), trim(values(sub)))) cycle outer
        end do inner
                
        exit
    end do outer

end subroutine substituteLine


subroutine ToWorld(localCoords, worldCoords, transform)
    use component_parameters, only: worldCoordsId

    implicit none

    real(kind=double), intent(in), dimension(3) :: localCoords
    real(kind=double), intent(out), dimension(3) :: worldCoords
    real(kind=double), intent(in), dimension(3, 3) :: transform

    if (trim(worldCoordsId) == '') then
        worldCoords = localCoords
    else
        worldCoords = matmul(transform, localCoords)
    end if
    
end subroutine ToWorld


!> Copy the template file to the destination file, making all required substitutions
subroutine substituteFile(sourceFile, destFile, keys, values, transform)
    use setup_module, only: dx, nSections
    use component_parameters, only: keyLength, valueLength, substitutions
    
    implicit none

    character(*), intent(in) :: sourceFile
    character(*), intent(in) :: destFile
    character(keyLength), intent(in), dimension(substitutions) :: keys
    character(valueLength), intent(in), dimension(substitutions) :: values
    real(kind=double), intent(in), dimension(3, 3) :: transform

    character(2000) :: line
    integer :: ios

    open(21, file=sourceFile, status='old')
    open(22, file=destFile)

    do
        read(21, fmt='(A)', iostat=ios) line

        if (ios == -1) exit

        if (index(line, '{{{YPOINTS}}}') > 0) then ! Completely replace line with y points xml
            call generatePoints(transform, -0.5D0 * dx, -dx, nSections)
        else if (index(line, '{{{QPOINTS}}}') > 0) then ! Completely replace line with q points xml
            call generatePoints(transform, 0.0D0, -dx, nSections + 1)
        else if (index(line, '{{{SEAWARD_BOUNDARY_IDS}}}') > 0) then ! Completely replace line with seaward boundary ids xml
            call generatePolylineIds(nSections)
        else if (index(line, '{{{SEAWARD_BOUNDARY_COORDS}}}') > 0) then ! Completely replace line with seaward boundary coords xml
            call generatePolylineCoords(transform, 0.0D0, -dx, seawardBoundary, nSections)
        else if (index(line, '{{{BASELINE_SECTIONS_IDS}}}') > 0) then ! Completely replace line with baseline sections ids xml
            call generatePolylineIds(nSections)
        else if (index(line, '{{{BASELINE_SECTIONS_COORDS}}}') > 0) then ! Completely replace line with baseline sections coords xml
            call generatePolylineCoords(transform, 0.0D0, -dx, 0.D0, nSections)
        else
            call substituteLine(line, keys, values) ! Make substituions within line
            write(22, fmt='(A)') line
        end if
    end do

    close(22)
    close(21)

end subroutine substituteFile


!> Generate an element set for y points or q points and write to the output file
subroutine generatePoints(transform, firstX, deltaX, count)
    implicit none
    integer, parameter :: double = kind(1.0D0)

    real(kind=double), intent(in), dimension(3, 3) :: transform
    real(kind=double), intent(in) :: firstX
    real(kind=double), intent(in) :: deltaX
    integer, intent(in) :: count

    integer :: point
    real(kind=double), dimension(3) :: localCoords
    real(kind=double), dimension(3) :: worldCoords

    do point = 1, count
        write(22, fmt='(A,I4,A,I4,A,I4,A)') &
            '<Identity id="', point, '"> <Describes caption="', point, '">', point, '</Describes> </Identity>'
    end do

    write(22, fmt='(A)') '<X><Values>'
    localCoords = (/ firstX, 0.0D0, 1.0D0 /)
    
    do point = 1, count - 1
        call ToWorld(localCoords, worldCoords, transform)
        write(22, fmt='(F12.2,A)') worldCoords(1), ','
        localCoords(1) = localCoords(1) + deltaX
    end do

    call ToWorld(localCoords, worldCoords, transform)
    write(22, fmt='(F12.2)') worldCoords(1)
    write(22, fmt='(A)') '</Values></X>'

    write(22, fmt='(A)') '<Y><Values>'
    localCoords = (/ firstX, 0.0D0, 1.0D0 /)
    
    do point = 1, count
        call ToWorld(localCoords, worldCoords, transform)
        write(22, fmt='(F12.2,A)') worldCoords(2), ','
        localCoords(1) = localCoords(1) + deltaX
    end do

    call ToWorld(localCoords, worldCoords, transform)
    write(22, fmt='(F12.2)') worldCoords(2)
    write(22, fmt='(A)') '</Values></Y>'

end subroutine generatePoints


!> Generate an element set for seaward boundary polyline coords and write to the output file
subroutine generatePolylineCoords(transform, firstX, deltaX, yValue, count)
    use component_parameters, only: worldCoordsId
    implicit none
    integer, parameter :: double = kind(1.0D0)

    real(kind=double), intent(in), dimension(3, 3) :: transform
    real(kind=double), intent(in) :: firstX
    real(kind=double), intent(in) :: deltaX
    real(kind=double), intent(in) :: yValue
    integer, intent(in) :: count

    integer :: point
    real(kind=double), dimension(3) :: localCoords1
    real(kind=double), dimension(3) :: worldCoords1
    real(kind=double), dimension(3) :: localCoords2
    real(kind=double), dimension(3) :: worldCoords2

    localCoords1 = (/ firstX, yValue, 1.0D0 /)
    localCoords2 = (/ firstX + deltaX, yValue, 1.0D0 /)

    do point = 1, count
        write(22, fmt='(A)') '<Coords>'
        write(22, fmt='(A)') '<X><Values>'
    
        call ToWorld(localCoords1, worldCoords1, transform)
        call ToWorld(localCoords2, worldCoords2, transform)
        
        write(22, fmt='(F12.2,A,F12.2)') worldCoords1(1), ',', worldCoords2(1)

        write(22, fmt='(A)') '</Values></X>'
        write(22, fmt='(A)') '<Y><Values>'
    
        write(22, fmt='(F12.2,A,F12.2)') worldCoords1(2), ',', worldCoords2(2)

        write(22, fmt='(A)') '</Values></Y>'
        write(22, fmt='(A)') '</Coords>'
        
        localCoords1(1) = localCoords1(1) + deltaX
        localCoords2(1) = localCoords2(1) + deltaX
    end do

end subroutine generatePolylineCoords


!> Generate an element set for seaward boundary polyline ids and write to the output file
subroutine generatePolylineIds(count)
    implicit none
    integer, parameter :: double = kind(1.0D0)

    integer, intent(in) :: count

    integer :: point

    do point = 1, count
        write(22, fmt='(A,I4,A,I4,A,I4,A)') &
            '<Identity id="', point, '"> <Describes caption="', point, '">', point, '</Describes> </Identity>'
    end do

end subroutine generatePolylineIds


end module support