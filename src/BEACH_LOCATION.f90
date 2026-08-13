!**********************************************
!*       BEACH_LOCATION.f90       
!**********************************************
!
!>    Purpose:
!>    To adjust the discretised description of 
!>    the beach (beachy(,)) at each section until 
!>    it approximately matches the calculated beach 
!>    volume (bvolume(i)).
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

module beach_location_module
    implicit none

    !> Permitted values for the beach_type argument.
    !> Changes the behaviour of subroutine BEACH_LOCATION according to the context in
    !> which it is called.
    integer, parameter :: beach_regular = 0
    integer, parameter :: beach_upper = 1
    integer, parameter :: beach_lower = 2
    integer, parameter :: beach_vellinga = 3

contains
    
    subroutine BEACH_LOCATION(iProfile, nYsections, nce, beachHeightE, top, bottom, &
        bruun, beachos, beachy, cliffy, beachdepth, bvolume, &
        surface, beach_type, beachelev, berm, &
        seawall_os, vellinga, vxs, vell_top_e)

    use setup_module
    implicit none
    
    integer, intent(in) :: iProfile !< Beach section index or -1 to loop over all beach sections
    integer, intent(in) :: nYsections !< Number of beach sections
    integer, intent(in) :: nce !< Number of cliff elements
    integer, intent(in) :: beachHeightE !< height expressed as number of cliff elements
    integer, intent(in), dimension(nYsections) :: top !< Top active cliff element at each beach section
    integer, intent(in), dimension(nYsections) :: bottom !< Bottom active cliff element at each beach section
    real(kind=double), intent(in), dimension(beachHeightE) :: bruun !< 'Bruun' shape for beach
    real(kind=double), intent(inout), dimension(nYsections) :: beachos !< Beach offset at each beach section
    real(kind=double), intent(inout), dimension(nce, nYsections) :: beachy !< discretised platform y value (distance of beach surface from baseline)
    real(kind=double), intent(in), dimension(nce, nYsections) :: cliffy !< discretised cliff and platform y value
    real(kind=double), intent(inout), dimension(nce) :: beachdepth !< depth of the beach
    real(kind=double), intent(in), dimension(nYsections) :: bvolume !< Beach volume at each beach section
    real(kind=double), intent(in), dimension(nce, nYsections) :: surface !< Effective beach surface
    
    integer, intent(in) :: beach_type !< takes values 0 - 3, as defined above
    
    !< optional arguemnts, controlled by beach_type argument
    integer, intent(in), optional, dimension(nYsections) :: beachelev !< only used if beach_regular
    real(kind=double), intent(inout), optional, dimension(nYsections) :: berm !< only used if beach_regular or beach_vellinga
    real(kind=double), intent(in), optional, dimension(nYsections) :: seawall_os !< only used if beach_vellinga
    real(kind=double), dimension(:), intent(in), optional :: vxs !< only used if beach_vellinga
    real(kind=double), dimension(:), intent(inout), optional :: vellinga !< only used if beach_vellinga
    integer, intent(in), optional :: vell_top_e !< only used if beach_vellinga
    
!*** Local variables

    real(kind=double), dimension(nce) :: exposed_surface
    real(kind=double) :: bVolumeApprox
    real(kind=double) :: bermstep

    integer :: bermflag
    integer :: bermflag1
    integer :: index
    integer :: first
    integer :: last
    integer :: section

    if (iProfile == -1) then ! loop over all nYsections
        first = 1
        last = nYsections
    else
        first = iProfile
        last = iProfile
    end if            
    
    do section = first, last

        if (beach_type == beach_vellinga) then
            
        	! Find where the beach overlies the platform
            do index = vell_top_e + 1, nce
                vellinga(index) = 0
            enddo
            
            do index = 1, vell_top_e
                vellinga(index) = seawall_os(section) + berm(section) + vxs(index)
            enddo
    
            ! SCAPE describes the beach volume twice, once with precision
            ! in the sediment transport / sediment budget calculations using the name:
            !        bvolume(i)
            ! and a second approximate value:
            !        bVolumeApprox 
            ! which is the volume calculated by summing the difference between the 
            ! discretised VELLINGA PROFILE and the discretised platform (cliffy(,)). 
            ! bVolumeApprox is changed by moving the discretised beach
            ! horizontally, until it is similar to bvolume(i).
        
            bVolumeApprox = 0.0d0
            
            do index = 1, nce
                if ((vellinga(index) - cliffy(index, section)) .gt. 0)then
                    bVolumeApprox = bVolumeApprox + (vellinga(index) - cliffy(index, section)) * dz * dx
                endif
            enddo
            
        else ! regular beach, upper beach or lower beach
            exposed_surface(bottom(section): top(section)) = &
                max(cliffy(bottom(section): top(section), section), surface(bottom(section): top(section), section))

            beachdepth = 0

            ! Find where the beach overlies the platform
            do index = bottom(section), top(section)
                if (beachy(index, section) > exposed_surface(index)) then
                    beachdepth(index) = beachy(index, section) - exposed_surface(index)
                end if
            end do
          
            ! SCAPE describes the beach volume twice, once with precision
            ! in the sediment transport / sediment budget calculations using the name:
            !        bvolume(i)
            ! and a second approximate value:
            !        bVolumeApprox 
            ! which is the volume calculated by summing the difference between the 
            ! discretised beach (beachy(,)) and the discretised platform (cliffy(,)). 
            ! bVolumeApprox is changed by moving the discretised beach
            ! horizontally, until it is similar to bvolume(i).
            
            bVolumeApprox = sum(beachdepth(bottom(section): top(section))) * dz * dx
        
            if (beach_type == beach_regular .and. bermSlope > 0.0D0) then
                if (beachy(beachelev(section) - 1, section) > cliffy(beachelev(section) - 1, section)) then
                    bVolumeApprox = bVolumeApprox + ((0.5D0 * dx * (beachy(beachelev(section) - 1, section) - &
                    cliffy(beachelev(section) - 1, section))**2) / bermSlope)
                    ! This has been added 80(13) to account for the volume of material assumed by 
                    ! PLATFORM_EROSION to be perched on top of any berm
                    ! This is only a rough estimate, it should be improved laer, in particular to 
                    ! take account of the presence of any talus.
                    ! NB the same calculation is made again below
                end if
            endif
        endif
        
        bermflag = 0 ! flag for direction of change
        bermflag1 = 0 ! flag for change in change direction
        bermstep = initialBermStep

          do while(bVolumeApprox > 1.05 * bvolume(section) .or. bVolumeApprox < 0.95 * bvolume(section))
            if(bVolumeApprox < bvolume(section))then ! The estimate of volume is too small
              if(bermflag == -1)then ! Last time the app beach was too large
                call HALVE_BERM_STEP(bermstep) ! Modify step size
                bermflag1 = 1        ! there has been a change in direction
              else if(bermflag == 1 .and. bermflag1 == 0)then
                call DOUBLE_BERM_STEP(bermstep) ! Modify step size
              end if

              ! Move the approximate beach out
              if (beach_type == beach_upper .or. beach_type == beach_lower) then
                  beachos(section) = beachos(section) + bermstep 
              else
                  berm(section) = berm(section) + bermstep
              end if
                  
              bermflag = 1
            else if(bVolumeApprox > bvolume(section))then ! The approximate beach is too large
              if(bermflag == 1)then ! last time the app beach was too small
                call HALVE_BERM_STEP(bermstep) ! Modify step size
                bermflag1 = 1 ! there has been a change in direction, 
              else if(bermflag == -1 .and. bermflag1 == 0)then
                call DOUBLE_BERM_STEP(bermstep) ! Modify step size
              end if

              ! Move the approximate beach in
              if (beach_type == beach_upper .or. beach_type == beach_lower) then
                  beachos(section) = beachos(section) - bermstep 
              else
                  berm(section) = berm(section) - bermstep
              end if
                  
              bermflag = -1
            end if

        if (beach_type == beach_vellinga) then
            
        	! Find where the beach overlies the platform
            do index = vell_top_e + 1, nce
                vellinga(index) = 0
            enddo
            
            do index = 1, vell_top_e
                vellinga(index) = seawall_os(section) + berm(section) + vxs(index)
            enddo
    
            bVolumeApprox = 0.0d0
            
            do index = 1, nce
                if ((vellinga(index) - cliffy(index, section)) .gt. 0)then
                    bVolumeApprox = bVolumeApprox + (vellinga(index) - cliffy(index, section)) * dz * dx
                endif
            enddo
            
        else ! regular beach, upper beach or lower beach
            ! Now the approximate beach has moved, check the approximate beach volume
            ! First 'redraw' the approximate beach
              if (beach_type == beach_regular) then
                  beachy(top(section) + 1: bottom(section) + 1: -1, section) = beachos(section) + berm(section) + bruun
              else
                  beachy(top(section): bottom(section): -1, section) = beachos(section) + bruun
              end if

            ! earlier verrsion in WkBk 9(97)
            
            beachdepth = 0

            ! Find where the beach overlies the platform
            do index = bottom(section), top(section)
                if (beachy(index, section) > exposed_surface(index)) then
                    beachdepth(index) = beachy(index, section) - exposed_surface(index)
                end if
            end do
            
            ! bVolumeApprox is the aparent beach volume (see above), it has to be 
            ! adjusted to (approximately) match the actual beach volume
            bVolumeApprox = sum(beachdepth(bottom(section): top(section))) * dz * dx
        
            if (beach_type == beach_regular .and. bermSlope > 0.0D0) then
                if (beachy(beachelev(section)-1,section) > cliffy(beachelev(section)-1, section))then
                    bVolumeApprox = bVolumeApprox + ((0.5D0 * dx*(beachy(beachelev(section) - 1,section) - &
                    cliffy(beachelev(section) - 1, section))**2) / bermSlope)
                    ! The extra volume has been added 80(13) to account for the volume of material assumed by 
                    ! PLATFORM_EROSION to be perched on top of any berm
                    ! This is only a rough estimate, it should be imporved laer, in particular to 
                    ! take account of the presence of any talus
                    ! NB the same calculation is made again above
                end if
            end if
        end if
        
            ! Rounding check to prevent over calculation
            if(nint(bVolumeApprox) == nint(bvolume(section)))then
              bVolumeApprox = bvolume(section)
            end if
          end do
      end do

    return
    end subroutine BEACH_LOCATION


    subroutine HALVE_BERM_STEP(bermstep)
    use setup_module
    implicit none

    real(kind=double), intent(inout) :: bermstep

    bermstep = bermstep / 2.0D0

    if (bermstep < minBermStep) then
        write(6,*) 'bermstep too small: ', bermstep
        stop
    end if
    
    end subroutine HALVE_BERM_STEP


    subroutine DOUBLE_BERM_STEP(bermstep)
    use setup_module
    implicit none

    real(kind=double), intent(inout) :: bermstep

    bermstep = bermstep * 2.0D0

    if (bermstep > maxBermStep) then
        write(6,*) 'bermstep too large: ', bermstep
        stop
    end if
    
    end subroutine DOUBLE_BERM_STEP

end module beach_location_module
