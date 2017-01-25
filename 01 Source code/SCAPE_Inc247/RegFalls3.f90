!****************************************************************************
!
! RegFalls3.F90
!
!****************************************************************************
!> This is a probabalistic model of landsliding that is separate from SCAPE,
!> runs after it. It takes the projections of cliff toe recession from SCAPE
!> and uses them to predict cliff top location. It does this many times, 
!> randomly selecting sequences of landslide period and size. In this way
!> an ensemble set of cliff top recession estimates is created.
!> See Dawson, R., Dickson, M., Nichols, R., Hall, J., Walkden, M., Stansby, P.,
!> Mokrech, M., Richards, J., Zhou, J., Milligan, J., Jordan, A., Pearson, S.,
!> Rees, S., Bates, P., Koukoulas, S., Watkinson, A., (2009)
!> Integrated analysis of risks of coastal flooding and cliff erosion under
!> scenarios of long term change. Climatic Change, Climatic Change (2009) 
!> 95:249-288. DOI 10.1007/s10584-008-9532-8. Springer.
!****************************************************************************
!
!   Copyright (C) James Hall 2014
!   This work is the intellectual property of the author.
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


module regfalls3_module
    implicit none

contains

    !> Estimate cliff top recession.
    subroutine RegFalls3(nYsections, nyears, recession, filename, use_random_seed)
        
    use normal_func
    use setup_module
    use utilities_module
    use local_data, only: FILE_PATH
    use exceptions

    implicit none

    integer, intent(in) :: nYsections !< number of beach sections
    integer, intent(in) :: nyears  !< number of years in run
    real(kind=double), intent(inout), dimension(nYsections, nyears) :: recession  !< horizontal retreat (of cliffy) during one year
    character (*), intent(in) :: filename !< name of cliff data file
    logical, intent(in) :: use_random_seed !< For testing. Always uses same sequence of random numbers if false

    integer :: i, j, m
    integer :: date_time (8)
    integer :: cpsect
     
    integer, dimension(:), allocatable ::  seed
    integer :: seed_dimension


    real(kind=double) :: meaninia, sdinia, meanfalla, sdfalla, meanstaba, sdstaba
    
    real(kind=double), dimension(nCliffSimulations) :: rand_ns

    real(kind=double), dimension(8) :: cliffdata 
    real(kind=double), dimension(nyears) :: toeret
    real(kind=double), dimension(nCliffSimulations) :: iniasim
    real(kind=double), dimension(nCliffSimulations) :: cht
    real(kind=double), dimension(nCliffSimulations) :: angle
        
    real(kind=double), dimension(nCliffSimulations) :: fallasim, stabasim
    logical :: haveFallInfo
    real(kind=double), dimension(nyears, nCliffSimulations) :: ctop

    integer :: cliffDataUnit
    integer :: saveCliffTopUnit
    integer :: iostat
    integer :: success_code
    
    cliffDataUnit = get_file_unit()
    saveCliffTopUnit = get_file_unit()
    
    success_code = do_open_text_input_file(cliffDataUnit, FILE_PATH, cliffDataFile)
    if (success_code < 0) return

    success_code = do_open_text_output_file(saveCliffTopUnit, FILE_PATH, filename)
    if (success_code < 0) return

    call date_and_time (values = date_time)

    call random_seed(size = seed_dimension)
    allocate(seed(1:seed_dimension))
    call random_seed(get = seed)

    if (use_random_seed) seed(1) = date_time(8)*date_time(7)

    call random_seed(put=seed)

! Read in data describing the cliff
! cliffdata.txt  (8 by nCliffSections) is [cliff height; current angle mean; current angle sd; pre-fall angle mean;
! sd; post-fall angle mean; sd; predsect]
    
    do
        ctop = 0
        
        read(cliffDataUnit, *, iostat=iostat) cliffdata
        if (iostat == -1) exit ! end of file
        success_code = SUCCESS_CODE_FROM_IOSTAT(cliffDataUnit, iostat, 'Reading cliffData')
        if (success_code < 0) exit
        
        cht = cliffdata(1)            ! cliff height
        meaninia = cliffdata(2)        ! Initial mean angle
        sdinia = sqrt(cliffdata(3))    ! Initial angle stdeviation
        meanfalla = cliffdata(4)        ! Mean falling angle
        sdfalla = sqrt(cliffdata(5))    ! Stdeviation of falling angle 
        meanstaba = cliffdata(6)        ! Mean stable angle, post fall
        sdstaba = sqrt(cliffdata(7))    ! Stdeviation stable angle
        cpsect = cliffdata(8)            ! SCAPE section number 

        toeret = recession(cpsect, :) !Toe retreat

        !simulate cliff angles before and after fall
        call Normal(rand_ns, nCliffSimulations)

        iniasim = (rand_ns * sdinia) + meaninia
    
        do i = 1, nCliffSimulations
            ctop(:, i) = cht(i) / tan(iniasim(i) * 3.142 / 180) !  Distance behind the local cliff toe
        enddo

        haveFallInfo = .false.

        do i = 1,nyears
            !for each simulation calculate new angle to cliff top based on how much the toe has retreated
            angle = atan(cht / (ctop(i, :) + toeret(i))) * 180 / 3.142   ! toe retreat is negative

            ! If a failure occurs then recalculate the cliff top position
            do j = 1, nCliffSimulations
                if (.not. haveFallInfo) then
                    call Normal(rand_ns, nCliffSimulations)
                    fallasim = (rand_ns * sdfalla) + meanfalla
                    
                    ! To simulate cliff stable angle stabasim 
                    call Normal(rand_ns, nCliffSimulations)
                    stabasim = (rand_ns * sdstaba) + meanstaba ! cliff stable angle

                    haveFallInfo = .true.
                end if
                    
                if ((angle(j) .gt. fallasim(j) .and. stabasim(j) < angle(j)) .or. angle(j) < 0) then
                    ! Finds the difference in horizontal position associated with the two angles
                    ctop(i,j) = ctop(i,j) + (cht(j) / (tan(stabasim(j) * 3.142 / 180)) - cht(j) / tan(angle(j) * 3.142 / 180))
                    
                    do m = i + 1, nyears
                        ctop(m, j) = ctop(i, j)
                    enddo
                                        
                    !ctop(i:nyears,j) = ctop(i,j) + [cht(j) / (tan(stabasim(j)*3.142/180)) - cht(j)/tan(angle(j)*3.142/180)]

                    haveFallInfo = .false.
                endif
            enddo

        enddo


        do i = 1,nyears ! years
            do j = 1, nCliffSimulations ! simulations
                if (i == 1) then ! year 1 so save initial data
                    write(saveCliffTopUnit,"(I3)") cpsect    ! seciton
                    write(saveCliffTopUnit,"(I3)") j ! simulation
                    write(saveCliffTopUnit,"(I3)") i    ! year
                    write(saveCliffTopUnit,"(F5.1)") ctop(i,j) ! offset
                elseif (ctop(i,j) /= ctop(i-1,j)) then
                    write(saveCliffTopUnit,"(I3)") cpsect    
                    write(saveCliffTopUnit,"(I3)") j ! simulation
                    write(saveCliffTopUnit,"(I3)") i    ! year
                    write(saveCliffTopUnit,"(F5.1)") ctop(i,j) ! offset
                endif
            enddo
        enddo

    enddo

    close(cliffDataUnit)
    
    close(saveCliffTopUnit)
    
    END SUBROUTINE RegFalls3

end module regfalls3_module
