!**********************************************
!*       dump.f95       
!**********************************************
!> Additional diagnostic dump support.

!> Dumps to txt files having names like dump<nnnnnn>.txt where <nnnnnn> increments for each call.
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

module dump_module

implicit none
save
integer, parameter :: double = kind(1.0D0)
private double
integer, private, parameter :: nce = 400
integer, private, parameter :: nYsections = 101
integer, private, parameter :: nQpoints = 102

integer :: dump_id !< Used in dump file name
logical :: dump_enable !< No output if this is false

contains

!> Dump q point or y point values to a text file
subroutine dump_real(location, item, qpoint_item, values)

    character(len=*), intent(in) :: location !< Unique string to link back to calling code
    character(len=*), intent(in) :: item !< Name of values to dump
    logical, intent(in) :: qpoint_item !< True if values are at Q points 
    real(kind=double), intent(in), dimension(:) :: values !< Array containing values

    integer :: index

    if (.not. dump_enable) return
    
    call open_file()
    
    call dump_values(location, item, qpoint_item, values)
    
    close(2014)

end subroutine dump_real


!> Dump y point values for entire cliff elevation to a text file
subroutine dump_cliff_real(location, item, values)

    character(len=*), intent(in) :: location !< Unique string to link back to calling code
    character(len=*), intent(in) :: item !< Name of values to dump
    real(kind=double), intent(in), dimension(:, :) :: values !< Array containing values

    integer :: index
    character(25) :: elevation

    if (.not. dump_enable) return
    
    call open_file()
    
    do index = 1, nce
        write(elevation,*) index
        call dump_values(location, item//elevation, .false., values(index, :))
    end do
    
    close(2014)

end subroutine dump_cliff_real


subroutine open_file()

    character(25) :: file_id
    character(25) :: filename

    write(file_id,*) dump_id
    dump_id = dump_id + 1

    filename = "dump_"//trim(file_id)//".txt"

    open(2014, file=filename)

end subroutine open_file


subroutine dump_values(location, item, qpoint_item, values)

    character(len=*), intent(in) :: location
    character(len=*), intent(in) :: item
    logical, intent(in) :: qpoint_item
    real(kind=double), intent(in), dimension(:) :: values

    integer :: index
    integer :: last

    write(2014,*) location
    write(2014,*) item

    if (qpoint_item) then
        last = nQpoints
    else
        last = nYsections
    end if

    do index = 1, last
        write(2014,*) index, values(index)
    end do
    
end subroutine dump_values

end module dump_module
