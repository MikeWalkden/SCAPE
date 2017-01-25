!**********************************************
!*       exceptions.f95       
!**********************************************
!> Support for exception handling when SCAPE+ is running freestanding.
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

module exceptions

implicit none

contains
 
    !> Test iostat to return success code.

    !> Log message if failure.
    integer function SUCCESS_CODE_FROM_IOSTAT(unit, iostat, message)
    integer, intent(in) :: unit !< FORTRAN I/O unit number
    integer, intent(in) :: iostat !< iostat value to test
    character(len=*), intent(in) :: message !< message to log on failure

    if (iostat /= 0)then
        write(6, *)'*** Error : unit: ', unit, ' unexpected iostat: ', iostat
        write(6, *)'*** Error : ', message
        SUCCESS_CODE_FROM_IOSTAT = -1
    else
        SUCCESS_CODE_FROM_IOSTAT = 0
    end if                    

    end function SUCCESS_CODE_FROM_IOSTAT
    

    !> Test allocation result to return success code.

    !> Log message if failure.
    integer function SUCCESS_CODE_FROM_ALLOCATION(alloc_error, message)
    integer, intent(in) :: alloc_error !< allocation result value to test
    character(len=*), intent(in) :: message !< message to log on failure

    if (alloc_error /= 0)then
        write(6, *)'*** Error : could not allocate space'
        write(6, *)'*** Error : for ', message
        SUCCESS_CODE_FROM_ALLOCATION = -2
    else
        SUCCESS_CODE_FROM_ALLOCATION = 0
    end if                    

    end function SUCCESS_CODE_FROM_ALLOCATION
 

    !> Log a failure message.
    subroutine RAISE_EXCEPTION(message)
    character(len=*), intent(in) :: message !< message to log

    write(6, *)'*** Error : ', message
    
    end subroutine RAISE_EXCEPTION

    
end module exceptions