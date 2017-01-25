!**********************************************
!*
!*       RUN_SCAPE.f90       
!*
!**********************************************
!
!> Main program to run SCAPE freestanding.
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

    use global_data
    use utilities_module
    use dump_module
    use setup_module
    use local_data
    use timestep_module
    use initialise_module
    use finish_module

    implicit none
    
    character (LEN = 12) REAL_CLOCK(3)
    integer date_time1(8)
    integer date_time2(8)
    real :: cpu_time1
    real :: cpu_time2
    real(kind=double) start_time
    real(kind=double) stop_time
    integer :: duration_days
    integer :: duration_hours
    integer :: duration_minutes
    integer :: duration_seconds
    integer :: iargc
    character(1) :: args
    integer :: success_code

    ! To calculate the total run time
    call cpu_time(cpu_time1)
    call date_and_time (REAL_CLOCK (1), REAL_CLOCK (2), REAL_CLOCK (3), date_time1)

    ! Use time to seed random number generator
    use_random_seed = .true.

    if (iargc() > 0) then
        call GETARG(1, arg)
        if (arg .eq. "NOSEED") use_random_seed = .false.
    end if

	call initialise(args, success_code)

	call prepare(success_code)
    if (success_code < 0) stop

    do while (ntide <= run_duration)
        ! Run duration is counted in tides
        call update(success_code)
        if (success_code < 0) stop
    end do 

    write(6,*)'*** Program finishing'

    call date_and_time (REAL_CLOCK (1), REAL_CLOCK (2), REAL_CLOCK (3), date_time2)
    call cpu_time(cpu_time2)
    start_time = julian_day(date_time1(2), date_time1(3), date_time1(1), date_time1(5), date_time1(6), date_time1(7))
    stop_time = julian_day(date_time2(2), date_time2(3), date_time2(1), date_time2(5), date_time2(6), date_time2(7))
    call elapsed_time(start_time, stop_time, duration_days, duration_hours, duration_minutes, duration_seconds)
    print*,  'run duration ', duration_seconds, '  -seconds', duration_minutes, '  -minutes',&
        duration_hours, '  -hours', duration_days, '  -days'
    print*, 'CPU time is ', (cpu_time2 - cpu_time1) / 60.0, ' minutes'
    
    print*, 'Final mslOffset: ', sea_rise

    call finish(success_code)

    stop
    end
