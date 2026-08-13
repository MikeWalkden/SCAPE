!**********************************************
!*       SeaLevelRise.f95       
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

module SLR_module
    use utilities_module
    use local_data, only: FILE_PATH, feLogString
    use exceptions
    use setup_module, only: seaLevelFile
    
    implicit none

save

    integer, parameter :: double = kind(1.0D0)
    private double
    integer :: slFileUnit = 0
    integer :: slYear

contains


!> Get sea level for given year from Sea Level file.
subroutine SL_for_year(year, sl, success_code)
    integer, intent(in) :: year !< year for which sea level is required
    real(kind=double), intent(inout) :: sl !< sea level found
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

    integer :: ios
    
    success_code = 0

    if (slFileUnit == 0) then ! first time through
        slFileUnit = get_file_unit()
        success_code = do_open_text_input_file(slFileUnit, FILE_PATH, seaLevelFile)
        if (success_code < 0) return
    end if

    do
        read(slFileUnit, *, iostat = ios) slYear, sl

        if (ios == -1) then
            exit
        else
            success_code = SUCCESS_CODE_FROM_IOSTAT(slFileUnit, ios, 'Reading sea level data')
            if (success_code < 0) return
        end if

        if (slYear >= year) then
            exit
        end if
    end do
    
    if (slYear /= year) then
        write(feLogString, *) 'Could not find SL data for year ', year
        call RAISE_EXCEPTION(feLogString)
        success_code = -1
        return
    end if
    
end subroutine SL_for_year

end module SLR_module
