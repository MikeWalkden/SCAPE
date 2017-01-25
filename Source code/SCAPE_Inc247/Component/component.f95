!**********************************************
!*       component.f95       
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

module component_parameters

use setup_module, only: UPPERCASE
use component_data

implicit none
save

contains

subroutine INITIALISE_COMPONENT_PARAMETERS(radcon, success_code)

    use utilities_module
    use exceptions
    implicit none

    real(kind=double), intent(in) :: radcon
    integer, intent(inout) :: success_code

    character(25) :: steering_file
    integer :: steering_file_unit
    character(25) :: keyword
    logical :: file_exists
    integer :: iostat

    steering_file = "component.txt"
    inquire(file=steering_file, exist=file_exists)
    
    if (file_exists) then
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
                    case('TEMPLATESPATH')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) templatesPath
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading templatesPath')
                        if (success_code < 0) return
                    case('ENGINEWRAPPERPATH')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) engineWrapperPath
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading engineWrapperPath')
                        if (success_code < 0) return
                    case('FLUIDEARTHPATH')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) fluidEarthPath
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading fluidEarthPath')
                        if (success_code < 0) return
                    case('SEAWARDBOUNDARY')
                        read(steering_file_unit, fmt=*, iostat= iostat) seawardBoundary
                    case('WORLDCOORDSID')
                        read(steering_file_unit, fmt='(A)', iostat= iostat) worldCoordsId
                        success_code = SUCCESS_CODE_FROM_IOSTAT(steering_file_unit, iostat, 'Reading worldCoordsId')
                        if (success_code < 0) return
                    case('LOCALORIGINX')
                        read(steering_file_unit, fmt=*, iostat= iostat) localOriginX
                    case('LOCALORIGINY')
                        read(steering_file_unit, fmt=*, iostat= iostat) localOriginY
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
    
end subroutine INITIALISE_COMPONENT_PARAMETERS


end module component_parameters
