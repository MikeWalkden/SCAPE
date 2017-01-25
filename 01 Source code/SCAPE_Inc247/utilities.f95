!**********************************************
!*       utilities.f95       
!**********************************************
!
!> Utilities including file opening, date and time calculations and conversions, 
!> logging, array initialisation and avergare value calculation.
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

module utilities_module
    implicit none

save

integer, parameter :: double = kind(1D0)
private double

integer :: binaryOpenMode = 0 !< used by open_binary_file() to record successful binary file open
character(11) :: accessMode !< used by open_binary_file() to pass access parameter value to open
character(6) :: formMode = 'binary' !< used by open_binary_file() to pass form parameter value to open

contains

!> Opens a text file for output and return a success_code value.
!>
!> If valuesPerYear is present, that value is written to the output file.
integer function do_open_text_output_file(unitNo, path, filename, valuesPerYear, units)
    use exceptions
    use local_data
    implicit none
    
    integer, intent(in) :: unitNo !< FORTRAN IO unit number
    character(*), intent(in) :: path !< path to location of output files
    character(*), intent(in) :: filename !< output file name
    integer, intent(in), optional :: valuesPerYear !< number of values per year to be written to the output file
    character(*), intent(in), optional :: units !< units of values to be written to the output file
    
    character(256) :: message
    character(256) :: filePath
    integer :: ios
    integer :: success_code

    message = 'Opening '//trim(path)//trim(filename)//'.txt'
    open(unitNo, file = trim(path)//trim(filename)//'.txt', iostat=ios)
    success_code = SUCCESS_CODE_FROM_IOSTAT(unitNo, ios, trim(message))

    if (success_code >= 0 .and. trim(filename) /= 'log') then
        filePath = trim(path)//trim(filename)//'.txt'
        write(logString, fmt='(A)') filePath
        call SCAPE_LOG()
        
        if (present(valuesPerYear)) then
            write(unitNo, *) valuesPerYear

            if (present(units) .and. units /= '') then
                write(logString, fmt='(A,I6,A)') '    Values per Year:', valuesPerYear, ' Units: '//trim(units)
            else
                write(logString, fmt='(A,I6)') '    Values per Year:', valuesPerYear
            end if

            call SCAPE_LOG()
        end if
    end if

    do_open_text_output_file = success_code

end function do_open_text_output_file


!> Opens a text file for input and return a success_code value.
integer function do_open_text_input_file(unitNo, path, filename)
    use exceptions
    use local_data
    implicit none
    
    integer, intent(in) :: unitNo !< FORTRAN IO unit number
    character(*), intent(in) :: path !< path to location of input files
    character(*), intent(in) :: filename !< input file name
    
    character(256) :: message
    integer :: ios

    message = 'Opening '//trim(path)//trim(filename)
    feLogString = message
    call feLog()
    open(unitNo, file = trim(path)//trim(filename), status='old', action='read', iostat=ios)
    do_open_text_input_file = SUCCESS_CODE_FROM_IOSTAT(unitNo, ios, trim(message))

end function do_open_text_input_file


!> Opens a binary file for input and return a success_code value.
integer function do_open_binary_input_file(unitNo, path, filename)
    use exceptions
    use local_data
    implicit none
    
    integer, intent(in) :: unitNo !< FORTRAN IO unit number
    character(*), intent(in) :: path !< path to location of input files
    character(*), intent(in) :: filename !< input file name
    
    character(256) :: message
    integer :: ios
    logical :: file_exists
    
    message = 'Opening '//trim(path)//trim(filename)
    feLogString = message
    call feLog()
    inquire(file=trim(path)//trim(filename), exist=file_exists)

    if (.not. file_exists) then
        call RAISE_EXCEPTION('*** Error : '//trim(path)//trim(filename)//' File does not exist')
        do_open_binary_input_file = -1
        return
    end if

    call open_binary_file(unitNo, trim(path)//trim(filename), ios)
    do_open_binary_input_file = SUCCESS_CODE_FROM_IOSTAT(unitNo, ios, trim(message))

end function do_open_binary_input_file


!> Opens a binary file for input, allowing for different FORTRAN compilers.
!>
!> The file is assumed to exist.
subroutine open_binary_file(unitNo, filename, ios)
    implicit none
    
    integer, intent(in) :: unitNo !< FORTRAN IO unit number
    character(*), intent(in) :: filename !< input file name
    integer, intent(inout) :: ios !< iostat returned by open
    
    if (binaryOpenMode == 0) then ! This works for Silverfrost FTN95
        accessMode = 'transparent'
        open(unitNo, file=trim(filename), form='unformatted', access=accessMode, status='old', iostat=ios)
    
        if (ios == 0) return
            
        binaryOpenMode = binaryOpenMode + 1
    end if
    
    if (binaryOpenMode == 1) then ! This works for gFORTRAN
        accessMode = 'stream'
        open(unitNo, file=trim(filename), form='unformatted', access=accessMode, status='old', iostat=ios)
    
        if (ios == 0) return
            
        binaryOpenMode = binaryOpenMode + 1
    end if
    
    if (binaryOpenMode == 2) then ! Original from MW
        open(unitNo, file=trim(filename), form=formMode, status='old', iostat=ios)

        if (ios == 0) return
            
        binaryOpenMode = 0
    end if

end subroutine open_binary_file


!> Computes elapsed time in days, hours, minutes and seconds
!> for a time interval defined by start and stop times in Julian Days.
subroutine elapsed_time(start_day, stop_day, days, hours, minutes, seconds)
    implicit none
    
    real(kind=double), intent(in) :: start_day !< interval start time in Julian days
    real(kind=double), intent(in) :: stop_day !< interval stop time in Julian days
    integer, intent(out) :: days !< elapsed days
    integer, intent(out) :: hours !< elapsed hours
    integer, intent(out) :: minutes !< elapsed minutes
    integer, intent(out) :: seconds !< elapsed seconds
    
    real(kind=double) :: time

    ! compute elapsed time in days, hours, minutes and seconds
    time = (stop_day - start_day)
    days = int(time)
    time = time - days
    time = time * 24
    hours = int(time)
    time = time - hours
    time = time * 60
    minutes = int(time)
    time = time - minutes
    time = time * 60
    seconds = nint(time)

    ! handle rollover
    if (seconds .eq. 60) then
          seconds = 0
          minutes = minutes + 1
          end if

    if (minutes .eq. 60) then
          minutes = 0
          hours = hours + 1
          end if

    if (hours .eq. 24) then
          hours = 0
          days = days + 1
          end if

end subroutine elapsed_time
 
     
!> Calculates and returns Julian day from Gregorian date.
function julian_day(month, day, year, hour, minute, second) result (jd)
    implicit none
    
    integer, intent(in) :: month !< Gregorian month
    integer, intent(in) :: day !< Gregorian day
    integer, intent(in) :: year !< Gregorian year
    integer, intent(in) :: hour !< Gregorian hour
    integer, intent(in) :: minute !< Gregorian minute
    integer, intent(in) :: second !< Gregorian second
    real(kind=double) :: jd !< Julian day

    real(kind=double) :: day_fraction
    integer :: y
    integer :: m
    integer :: a
    integer :: b
    
    day_fraction = dble(day) + hour / 24.0D0 + minute / 1440.0D0 + second / 86400.0D0

    if (month .gt. 2) then
        y = year
        m = month
    else
        y = year - 1
        m = month + 12
    end if

    a = y / 100
    b = 2 - a + a / 4

    jd = int(365.25D0 * (y + 4716)) + int(30.6001D0 * (m + 1)) + day_fraction + b - 1524.5D0

    return

end function 
 
 
!> Calculates Gregorian date from Julian day.
subroutine gregorian_date(julian_day, month, day, year)
    implicit none
    
    real(kind=double), intent(in) :: julian_day !< Julian day
    integer, intent(out) :: year !< Gregorian year
    integer, intent(out) :: month !< Gregorian month
    real(kind=double), intent(out) :: day !< Gregorian day

    real(kind=double) :: jd
    integer :: z
    real(kind=double) :: f
    integer :: a
    integer :: alpha
    integer :: b
    integer :: c
    integer :: d
    integer :: e
    integer :: y
    
    jd = julian_day + 0.5D0 
    
    z = int(jd)
    f = jd - z

    if (z .lt. 2299161) then
       a = z
    else
       alpha = int((z - 1867216.25d0) / 36524.25d0)
       a = z + 1 + alpha - alpha / 4
    end if

    b = a + 1524
    c = int((b - 122.1d0) / 365.25d0)
    d = int(365.25d0 * c)
    e = int((b - d) / 30.6001d0)

    day = b - d - int(30.6001d0 * e) + f

    if (e .lt. 14) then
       month = e - 1
    else
       month = e - 13
    end if

    if (month .gt. 2) then
       y = c - 4716
    else
       y = c - 4715
    end if

    year = y

end subroutine gregorian_date
     
 
!> Calculates and returns Modified Julian day from Julian day.
function to_mjd(jd) result(mjd)
    implicit none
    
    real(kind=double), intent(in) :: jd !< Julian day
    real(kind=double) :: mjd !< Modified Julian day
    
    mjd = jd - 2400000.5D0

end function to_mjd
  
 
!> Calculates and returns Julian day from Modified Julian day.
function to_jd(mjd) result(jd)
    implicit none
    
    real(kind=double), intent(in) :: mjd !< Modified Julian day
    real(kind=double) :: jd !< Julian day
    
    jd = mjd + 2400000.5D0

end function to_jd
  
 
!> Returns a unit number that is not in use.
integer function get_file_unit()
    integer, save :: lu = 9 !< Avoid pre-assigned units
    integer :: iostat
    logical :: opened

    do while (lu < 100)
        lu = lu + 1
        inquire (unit=lu, opened=opened, iostat=iostat)
        if (iostat.ne.0) cycle
        if (.not.opened) exit
    end do

    get_file_unit = lu
    return
    
end function get_file_unit


!> Write logString to the log file and optionally to the screen.
subroutine SCAPE_LOG(toScreen)
    use local_data, only: logFileUnit, logString

    implicit none

    logical, optional :: toScreen !< write logString to screen is present and true

    logical :: opened
    integer :: ios

    inquire(unit=logFileUnit, opened=opened)
    
    if (opened) then
        write(logFileUnit, fmt='(A)', iostat=ios) trim(logString)
    end if

    if (present(toScreen) .and. toScreen) then
        print*, trim(logString)
    end if

end subroutine SCAPE_LOG


!> Support for FluidEarth (OpenMI) logging.
subroutine feLog()
    use local_data, only: feLogEnabled, feLogFile, feLogString, feLogStrings, feLogLimit, &
        feLogIndex, feLogAppend, feLogActive, FILE_PATH
    
    implicit none

    logical :: opened
    integer :: ios

    if (feLogEnabled) then
        if (feLogActive) then
            if (feLogAppend) then
                open(unit=98, file=trim(feLogFile), position='append', iostat=ios)
                    
                inquire(unit=98, opened=opened)
                
                if (opened) then
                    write(98, fmt='(A)', iostat=ios) trim(feLogString)
                    close(98)
                end if
            else ! feLogAppend
                inquire(unit=98, opened=opened)
                
                if (opened) then
                    write(98, fmt='(A)', iostat=ios) trim(feLogString)
                end if
            end if ! feLogAppend
        else if (feLogIndex < feLogLimit) then
            feLogIndex = feLogIndex + 1
            feLogStrings(feLogIndex) = feLogString
        end if ! feLogActive
    end if ! feLogEnabled

end subroutine feLog


!> Initialise a range of values in an array by interpolation.
subroutine INTERPOLATE_VALUES(dest, value1, value2)
    implicit none

    real(kind=double), dimension(:), intent(inout) :: dest !< the array to be initialised
    real(kind=double) :: value1 !< value at first index
    real(kind=double) :: value2 !< value at last index

    real(kind=double) :: offset
    integer :: dest_size
    integer :: index

    dest_size = size(dest, 1)

    if (dest_size > 0) then
        dest(1) = value1
    
        if (dest_size > 2) then
            offset = (value2 - value1) / (dest_size - 1)

            do index = 2, dest_size - 1
                dest(index) = dest(index - 1) + offset
            end do
        end if
    
        dest(dest_size) = value2
    end if
    
end subroutine INTERPOLATE_VALUES


!> Average the values in an array.
function AVERAGE_VALUE(values) result(avg)
    implicit none

    real(kind=double), intent(in), dimension(:) :: values !< the array containing values to be averaged
    real(kind=double) :: avg !< calculated average value

    if (size(values, 1) == 0) then
        avg = 0.0d0
    else
        avg = SUM(values) / size(values, 1)
    end if

end function AVERAGE_VALUE


!> Average the values in a masked array array.
function MASKED_AVERAGE_VALUE(values, mask) result(avg)
    implicit none

    real(kind=double), intent(in), dimension(:) :: values !< the array containing values to be averaged
    logical, intent(in), dimension(:) :: mask !< The mask
    real(kind=double) :: avg !< calculated average value
    
    integer :: number
    
    number = COUNT(mask)

    if (number == 0) then
        avg = 0.0d0
    else
        avg = SUM(values, mask) / number
    end if

end function MASKED_AVERAGE_VALUE



end module utilities_module
 