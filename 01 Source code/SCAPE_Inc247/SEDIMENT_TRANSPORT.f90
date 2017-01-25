!**********************************************
!*
!*       SEDIMENT_TRANSPORT.f90       
!*
!**********************************************
!
!>    This module of SCAPE contains a one-line beach model, with some 
!>    additional algorithms. 
!>
!>    A 'bar' is assumed to exist offshore. This does nothing except store sediment
!>    moved offshore from the beach. Sediment transport rates offshore during 
!>    storms are determined by COSMOS output. A simple function has been assumed 
!>    that moves material back during calm periods.
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

module sediment_transport_module
    implicit none

    contains


    !> Implements the one-line beach model.
    subroutine SEDIMENT_TRANSPORT(nce, nYsections, nQpoints, year, &
        offshore_st, bvolume, sediment_influx, obvolume, top_drift, bot_drift, h, &
        sect_drift, beachy, cliffy, QAnn, QPAnn, &
        period, tstep_secs, h_qpoint, cgroup_qpoint, theta_qpoint, cerck1)

    use setup_module
    use global_data, only: g, beachCrestLevelE, beachHeightE, &
        beach_temp, cliff_temp, groyne_construction, groyne_removal, groyneEffectAnnualSum, &
        left_exclusion_angle, right_exclusion_angle
    use local_data, only: pot_left_sediment_flux, pot_right_sediment_flux, logString, &
        left_sediment_flux, right_sediment_flux, &
        have_left_sediment_input, left_sediment_input, have_right_sediment_input, right_sediment_input
    use local_data, only: feLogString
    use utilities_module

    implicit none

    integer, intent(in) :: nce !< number of cliff elementss
    integer, intent(in) :: nYsections !< number of beach sections
    integer, intent(in) :: nQpoints !< number of Q points
    integer, intent(in) :: year !< Current year
!    integer, intent(in) :: ngroynes !< number of groynes
    real(kind=double), intent(out), dimension(nYsections) :: offshore_st !< Offshore sediment transport
    real(kind=double), intent(inout), dimension(nYsections) :: bvolume !< beach volume
    real(kind=double), intent(in), dimension(nYsections) :: sediment_influx !< Sediment influx at the model seaward boundary
    real(kind=double), intent(inout), dimension(nYsections) :: obvolume !< Offshore bar volume
    integer, intent(in), dimension(nYsections) :: top_drift !< Upper vertical limit of the drift distribution
    integer, intent(in), dimension(nYsections) :: bot_drift !< Lower vertical limit of the drift distribution
    real(kind=double), intent(in), dimension(nYsections) :: h !< Wave height
    real(kind=double), intent(inout), dimension(nce, nYsections) :: sect_drift !< Cross-shore distribution of longshore drift
!    integer, intent(in), dimension(ngroynes) :: groynesection
!    real(kind=double), intent(in), dimension(nce, ngroynes) :: groyney
    real(kind=double), intent(in), dimension(nce, nYsections) :: beachy !< Discretised platform y value (distance of beach surface from baseline)
    real(kind=double), intent(in), dimension(nce, nYsections) :: cliffy !< Discretised cliff and platform y value
    real(kind=double), intent(inout), dimension(nQpoints) :: QAnn !< Annual sediment transport at each Q point
    real(kind=double), intent(inout), dimension(nQpoints) :: QPAnn !< Potential annual sediment transport at each Q point
    real(kind=double), intent(in) :: period !< Wave period
    real(kind=double), intent(in) :: tstep_secs !< Time step in seconds
    real(kind=double), intent(in), dimension(nQpoints) :: h_qpoint !< Wave height at Q points
    real(kind=double), intent(in), dimension(nQpoints) :: cgroup_qpoint !< Group celerity at Q points
    real(kind=double), intent(in), dimension(nQpoints) :: theta_qpoint !< Wave angle at Q points
    real(kind=double), intent(in), dimension(nQpoints) :: cerck1 !< CERC coefficients at Q points
    
!*** Local variables

    real(kind=double), allocatable, dimension(:) :: Q, Qp
    real(kind=double), allocatable, dimension(:) :: drift_fraction
    real(kind=double), allocatable, dimension(:) :: bvolnew
    real(kind=double), allocatable, dimension(:) :: tanBeachAlt ! Beach tangent estimated from beachy values for active part of beach

    real(kind=double) :: volcheck
    real(kind=double) :: factor ! Used to take invariant part of alternate longshore sediment transport equation out of the qPoint loop
    real(kind=double) :: tanBeachQpoint ! Beach tangent at Qpoint estimated from tangent at adjacent beach sections 
    integer :: jmax ! index of top vertical element for which sect_drift is not zero
    integer :: jmin ! index of bottom vertical element for which sect_drift is not zero
 
    integer section, qPoint, j
    integer alloc_error
    real(kind=double) :: beachWidth
    real(kind=double) :: bermWidth
    real(kind=double) :: groyneEffect
 
!*** Allocate arrays

    allocate(Q(nQpoints), Qp(nQpoints), drift_fraction(nQpoints), bvolnew(nYsections), stat = alloc_error)
    if(alloc_error /= 0)then
      write(6,*)'*** Error : could not allocate space'
      write(6,*)'*** Error : subroutine SEDIMENT_TRANSPORT'
      stop
    end if                               

    if (useAltSedimentTransport) then
        allocate(tanBeachAlt(nYsections), stat = alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Error : could not allocate space'
          write(6,*)'*** Error : subroutine SEDIMENT_TRANSPORT'
          stop
        end if                               
    end if                               

    Q=0
    Qp=0
    drift_fraction=0
    bvolnew=0

    !*************************************************
    !    Account for elements where there is no beach
    !*************************************************

    do section = 1, nYsections
        do j = bot_drift(section), top_drift(section)    ! Work on the elements within the range of drift
            if (beachy(j,section) - cliffy(j,section) <= 0) sect_drift(j, section) = 0
        end do
    end do


    if (useAltSedimentTransport) then
        factor = 0.00018D0 * (g**0.5) / (medianGrainsize**0.6)
    
        do section = 1, nYsections ! Calculate beach tangents at beach sections
            jmax = 0
            jmin = 0

            do j = 1, nce
                if (sect_drift(j, section) > 0.0D0) jmax = j
                if (jmin == 0 .and. sect_drift(j, section) > 0.0D0) jmin = j
            end do

            if (jmin == 0 .or. jmax == jmin) then
                ! No range of drift or only one beach elements within the range of drift.
                ! Use the whole beach beach tangent, up to 3m beach height.
                jmax = beachCrestLevelE
                jmin = jmax - min(beachHeightE, nint(3.0D0 / dz))
            end if

            tanBeachAlt(section) = ((jmax - jmin) * dz) / (beachy(jmin, section) - beachy(jmax, section))
        end do
    end if                               

    call EXCHANGE_MATERIAL(nYsections, offshore_st, bvolume, obvolume, h, period, tstep_secs)

    ! This code excludes waves from inappropriate directions

    do qPoint = 1, nQpoints
        if ((theta_qpoint(qPoint) < 998.0D0) .and. &
              ((excludeWaves == 0) .or. ((theta_qpoint(qPoint) <= right_exclusion_angle(qPoint)) .and. &
              (theta_qpoint(qPoint) >= -left_exclusion_angle(qPoint))))) then
              
            if (useAltSedimentTransport) then ! Use Alternate sediment transport calculation
                if (qPoint == 1) then
                    tanBeachQpoint = tanBeachAlt(1)
                else if (qPoint == nQpoints) then
                    tanBeachQpoint = tanBeachAlt(nYsections)
                else
                    tanBeachQpoint = (tanBeachAlt(qPoint - 1) + tanBeachAlt(qPoint)) * 0.5D0
                end if
            
                ! Alternate longshore sediment transport equation.
                Qp(qPoint) = factor * (h_qpoint(qPoint)**3.1) * (tanBeachQpoint**0.4) * sin(2.0D0 * theta_qpoint(qPoint))
                
            else ! Use the CERC equation of longshore sediment transport.
                Qp(qPoint) = h_qpoint(qPoint) * h_qpoint(qPoint) * cgroup_qpoint(qPoint) * &
                    sin(2.0D0 * theta_qpoint(qPoint)) * cerck1(qPoint) / 37.8D0
            end if
        end if
    end do

    !*************************************************
    !    update effect of groynes
    !*************************************************

    if (includeGroynes) then
        ! Groynes on the Q sections
        do qPoint = 1, nQpoints
            ! Installation of groynes,
                     
            ! ------------- MC: 101003, Groyne Bypassing
                   
            if (year >= groyne_construction(qPoint) .and. (groyne_construction(qPoint) > 0) &
                .and. year < groyne_removal(qPoint))then
  
                ! MC: 10 Oct 2003
                ! For sections in which groynes have been built, this bit of code
                ! checks what the beach width at that section is and adjusts the  
                ! sediment transport rate through the groyne effect value
                    
                ! beach width was calculated on the Y sections
                ! but that data needed at the Q sections

                ! Find updrift section
                if (Qp(qPoint) > 0.0D0) then
                    section = qPoint - 1
                else
                    section = qPoint
                end if

                if (section /= 0 .and. section /= nQpoints) then
                    !Updrift section is inside model boundary
                    beach_temp = 0.0D0
                    cliff_temp = 0.0D0
                    where (beachy(:, section) >= cliffy(:, section)) 
                        beach_temp = beachy(:, section)
                        cliff_temp = cliffy(:, section) 
                    end where
                    
                    beachWidth = maxval(beach_temp) - minval(beach_temp, mask = beach_temp > 0.0D0)
                    bermWidth = minval(beach_temp, mask = beach_temp > 0.0D0) - minval(cliff_temp, mask = cliff_temp > 0.0D0)
                    beachWidth = beachWidth + bermWidth
        
                    groyneEffect = GROYNE_EFFECT_LOOKUP(beachWidth)
                    
                    if (groyneEffect /= 0.0D0) then
                        Qp(qPoint) = Qp(qPoint) * (1.0D0 - groyneEffect)
                        groyneEffectAnnualSum(qPoint) = groyneEffectAnnualSum(qPoint) + groyneEffect
                    end if
                end if
            end if
        end do ! qPoint = 1, nQpoints
    end if

    ! Extrapolate the boundary conditions from the two sections just inside the model
    !
        Qp(1) = Qp(2) - ((Qp(3)-Qp(2)))
        Qp(nQpoints) = Qp(nQpoints-1) - ((Qp(nQpoints-2)-Qp(nQpoints-1)))
    !

    Qp = Qp * tstep_secs

    if (Qp(1) > 0.0D0) then
        ! influx at right boundary
		if (Qp(1) > qpMaxBoundaryRightIn)then
        	Qp(1) = qpMaxBoundaryRightIn
        end if
    else
        ! outflux at right boundary
		if (Qp(1) < -qpMaxBoundaryRightOut)then
        	Qp(1) = -qpMaxBoundaryRightOut
        end if
    end if
        
    if (Qp(nQpoints) > 0.0D0) then
        ! outflux at left boundary
		if (Qp(nQpoints) > qpMaxBoundaryLeftOut)then
        	Qp(nQpoints) = qpMaxBoundaryLeftOut
        end if
    else
        ! influx at left boundary
		if (Qp(nQpoints) < -qpMaxBoundaryLeftIn)then
        	Qp(nQpoints) = -qpMaxBoundaryLeftIn
        end if
    end if
        

!****************************************************************
!    Add up what is left of the distribution of longshore drift
!****************************************************************

    do section = 1, nYsections 
      drift_fraction(section) = sum(sect_drift(:,section))
    end do


! The drift fraction has been calculated at the Y points, but must be 
! applied to transport at the Q points

    do qPoint = 2, nQpoints - 1    ! Not able to deal with the boundaries here
      if(Qp(qPoint) < 0)then
        Q(qPoint) = Qp(qPoint) * drift_fraction(qPoint)
      else
        Q(qPoint) = Qp(qPoint) * drift_fraction(qPoint - 1)
      endif
    end do

!*****************************************************************************
! The boundaries are constricted by the drift fraction just inside the model
!*****************************************************************************
!
! First limit the smallness of the drift fraction for the boundary edges.
! When the beach is empty at the edges some sand must be allowed to come in.

    if(drift_fraction(1) < minBoundaryTransportRatio) drift_fraction(1) = minBoundaryTransportRatio
    if(drift_fraction(nYsections) < minBoundaryTransportRatio) drift_fraction(nYsections) = minBoundaryTransportRatio

    Q(nQpoints) = Qp(nQpoints) * drift_fraction(nYsections)
    Q(1) = Qp(1) * drift_fraction(1)

    ! Get values from OpenMI.
    ! Looking offshore, left end is at Q(nQpoints), right is at Q(1) and +ve Q is right to left
    if (have_left_sediment_input) then
        Q(nQpoints) = -left_sediment_input
    end if

    if (have_right_sediment_input) then
        Q(1) = right_sediment_input
    end if

! Update the beach volumes
    do section = 1, nYsections
      bvolnew(section) = bvolume(section) + Q(section) - Q(section + 1)
    end do

    volcheck = minval(bvolnew)
! The sediment balancing algorithms have been checked out on page 40 WkBk 6
    do while(volcheck < 0) ! if this is less than zero
      do section = 1, nYsections ! Find out the (first) empty section
        if(bvolnew(section) < 0)then ! It is section profile
        ! For the volume to have dropped below zero, there must have been flow away.
         ! There are three possibilities, flow away downchainage, upchainage and in both directions
         
          if(Q(section) < 0)then ! Material is flowing away downchain
            if(Q(section+1) > 0)then ! Material is flowing away upchain as well
              Q(section) = -1*floor((-1)*  bvolume(section) * ((Q(section)/Q(section+1)) / (1 - (Q(section)/Q(section+1))))  )
              Q(section+1) = floor(bvolume(section) + Q(section)) 
            else ! Only flow away downchainage
              Q(section) =  -1*floor((-1)* (Q(section+1) - bvolume(section)))
            end if
          else ! only flow away upchainage
            Q(section+1) = floor(bvolume(section) + Q(section))
          end if
          exit     ! so that it only corrects one beach volume at a time
        end if
      end do
      do section = 1, nYsections
        bvolnew(section) = nint(bvolume(section) + Q(section) - Q(section + 1))
      end do
      volcheck = minval(bvolnew)
    end do

    ! Save values for OpenMI.
    ! Looking offshore, left end is at Q(nQpoints), right is at Q(1) and +ve Q is right to left
    pot_left_sediment_flux = Qp(nQpoints)
    left_sediment_flux = Q(nQpoints)
    
    pot_right_sediment_flux = Qp(1)
    right_sediment_flux = Q(1)
    
    ! Finish updating the beach volumes
    bvolume = bvolnew

    ! Inject sediment_influx (if any) from OpenMI
    bvolume = bvolume + sediment_influx

    ! code to record the average sediment transport each year
      QAnn = QAnn + Q
      QpAnn = QpAnn + Qp

!*** Deallocate arrays
    if (useAltSedimentTransport) then
        deallocate(tanBeachAlt, stat = alloc_error)
    end if

       deallocate(Q, Qp, drift_fraction, bvolnew, stat=alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Unexpected deallocation error'
          write(6,*)'*** for subroutine SEDIMENT_TRANSPORT'
        endif

        return
        end subroutine SEDIMENT_TRANSPORT

!**************************************************
!> Exchange material between the bar and the beach
!> New, 130602, see p151(10)
!**************************************************
subroutine EXCHANGE_MATERIAL(nYsections, offshore_st, bvolume, obvolume, &
    h, period, tstep_secs)

    use setup_module
    
    implicit none

    integer, intent(in) :: nYsections !< number of beach sections
    real(kind=double), intent(out), dimension(nYsections) :: offshore_st !< Offshore sediment transport
    real(kind=double), intent(inout), dimension(nYsections) :: bvolume !< beach volume
    real(kind=double), intent(inout), dimension(nYsections) :: obvolume !< Offshore bar volume
    real(kind=double), intent(in), dimension(nYsections) :: h !< Wave height
    real(kind=double), intent(in) :: period !< Wave period
    real(kind=double), intent(in) :: tstep_secs !< Time step in seconds
    
    real(kind=double) :: offstrLookup
    integer :: section

    offshore_st = 0
    
    do section = 1, nYsections
        if (h(section) > minStormWaveHeight) then
            offstrLookup = OFFSHORE_STR_LOOKUP(h(section), period)
            offshore_st(section) = dx * offstrLookup * tstep_secs / 3600 ! convert 
            ! to units suitable for this timestep
            if(-1*offshore_st(section) > bvolume(section))then
                offshore_st(section) = -1*bvolume(section)
            endif
        else ! Not a storm condition so the material moves onshore 
            offshore_st(section) = obvolume(section) * beachReturnRatio 
        endif

        bvolume(section) = bvolume(section) + offshore_st(section)
        obvolume(section) = obvolume(section) - offshore_st(section)

        if (obvolume(section) < 1.0D-6) obvolume(section) = 0.0D0

    enddo


end subroutine EXCHANGE_MATERIAL

end module sediment_transport_module
