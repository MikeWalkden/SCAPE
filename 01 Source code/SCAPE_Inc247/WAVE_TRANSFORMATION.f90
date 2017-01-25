!**********************************************
!*
!*       WAVE_TRANSFORMATION.f90       
!*
!**********************************************
!
!> The input wave conditions within the model are assumed to be constant throughout a tidal timestep. 
!> Wave transformation processes including; refraction and shoaling are also included within the model 
!> and are represented using linear wave theory (Kamphuis, 1992, Kamphuis, 2000). 
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

module wave_transformation_module
    implicit none

contains

    !> Transforms waves from the wave points to the shore.
    subroutine WAVE_TRANSFORMATION(nQpoints, nYsections, nce, &
            hwp, t,anglewp, depthwp_msl, tidal_amp, &  ! MD 300304
            heightsurge, angleosc, periodValues, depthValues, clut,    &
            cglut,breakrat, beachy,cliffy,        &
            msl_e,setup,setdown,h,anglewob,cgroup, theta,      &
            setup_qpoint,setdown_qpoint,h_qpoint,anglewob_qpoint,cgroup_qpoint, &
            theta_qpoint, wpr)
        
    use setup_module, only: double, nWavePoints, dx, baseLineAngle_radians, depthOscMsl, fswp
    use utilities_module, only: AVERAGE_VALUE
    use global_data, only: pi, two_pi, half_pi
    use local_data, only: feLogString
    use utilities_module, only: feLog

    implicit none

    integer, intent(in) :: nQpoints !< number of Q points
    integer, intent(in) :: nYsections !< number of beach sections
    integer, intent(in) :: nce !< number of cliff elements
!    integer, intent(in) :: ngroynes !< number of groynes
    real(kind=double), intent(inout), dimension(nce, nYsections) :: cliffy !< Discretised cliff y value
    real(kind=double), intent(inout), dimension(nce, nYsections) :: beachy !< Discretised platform y value (distance of beach surface from baseline)
    integer, intent(in) :: periodValues !< Celerity look up tables first dimension
    integer, intent(in) :: depthValues !< Celerity look up tables second dimension
    real(kind=double), dimension(periodValues, depthValues), intent(in) :: clut !< wave celerity
    real(kind=double), dimension(periodValues, depthValues), intent(in) :: cglut !< wave group celerity
    real(kind=double), intent(inout), dimension(nYsections) :: setup
    real(kind=double), intent(inout), dimension(nYsections) :: setdown
    real(kind=double), intent(inout), dimension(nYsections) :: h !< Wave height at beach sections
    real(kind=double), intent(inout), dimension(nYsections) :: anglewob
    real(kind=double), intent(inout), dimension(nYsections) :: cgroup !< Group celerity at beach sections
    real(kind=double), intent(inout), dimension(nYsections) :: theta !< Wave angle at beach sections
    real(kind=double), intent(inout), dimension(nQpoints) :: setup_qpoint
    real(kind=double), intent(inout), dimension(nQpoints) :: setdown_qpoint
    real(kind=double), intent(inout), dimension(nQpoints) :: h_qpoint !< Wave height at Q points
    real(kind=double), intent(inout), dimension(nQpoints) :: anglewob_qpoint
    real(kind=double), intent(inout), dimension(nQpoints) :: cgroup_qpoint !< Group celerity at Q points
    real(kind=double), intent(inout), dimension(nQpoints) :: theta_qpoint !< Wave angle at Q points
    real(kind=double), intent(in), dimension(nWavePoints) :: hwp !< Wave height at wave point
    real(kind=double), intent(in), dimension(nWavePoints) :: anglewp !< Wave angle at wave point
!    real(kind=double), intent(in), dimension(ngroynes) :: groyneoffset
!    real(kind=double), intent(in), dimension(ngroynes) :: groynelength
    real(kind=double), intent(in), dimension(nQpoints) :: angleosc !< Offshore contour angle
    real(kind=double), intent(in), dimension(nWavePoints) :: depthwp_msl !< Depth of each wave point below mean sea level (m)
      
    real(kind=double), intent(in) :: t !< Wave period
    real(kind=double), intent(in), dimension(nYsections) :: tidal_amp !< Tidal amplitude
!    real(kind=double), intent(in) :: msl_m,
    real(kind=double), intent(in) ::heightsurge !< Not used, left at zero
    real(kind=double), intent(in) :: breakrat !< Ratio of breaker height and water depth
!    real(kind=double), dimension(2), intent(in) :: hour_section_limit !< hour_section_limits are the vertical limits

    integer, intent(in), dimension(nQpoints) :: wpr !< MWPM wpr is Wave Point Reference (wave point number for each Q point)
!    real(kind=double), intent(in) :: radcon

    integer, intent(in) :: msl_e !< Number of the platform element at msl

!*** Local variables

    real(kind=double), allocatable, dimension(:) :: angbeachbp !< Beach angle at q point
    real(kind=double), allocatable, dimension(:) :: beach_breaker_os !< Beach breaker offshore distance at profile
    real(kind=double), allocatable, dimension(:) :: avg_beach_breaker_os
    real(kind=double), allocatable, dimension(:) :: diffraction_theta
    real(kind=double), allocatable, dimension(:) :: hiosc
    real(kind=double), allocatable, dimension(:) :: angleoscbp
    real(kind=double), allocatable, dimension(:) :: dbreak_wt
    real(kind=double), allocatable, dimension(:) :: depthwp !< Depth at wave point

    real(kind=double) :: k,ks,kr
    real(kind=double) :: kd
    real(kind=double) :: alphaquad
    real(kind=double) :: alphaoosc, alphaiosc
    real(kind=double) :: alphaobp
    real(kind=double) :: directionfactor
    real(kind=double) :: c, cg, hbp, alphaibp, dbreak_new
    real(kind=double) :: depthosc

    integer, allocatable, dimension(:) :: active_groynes

    integer :: i
    integer :: idepth,idepthwp,idepthosc,it
    integer :: quadswitch,quadrant
    integer :: counter1
    integer :: repeatswitch,repeatcounter
    integer :: alloc_error
    integer :: qPoint
    integer :: section
    integer :: lastWPsection
    integer :: wp

!*** Allocate temp arrays
        allocate(angbeachbp(nQpoints), beach_breaker_os(nYsections), avg_beach_breaker_os(nQpoints),    &
         diffraction_theta(nYsections), active_groynes(nQpoints), hiosc(nQpoints),    &
         angleoscbp(nQpoints),dbreak_wt(nQpoints),depthwp(nWavePoints),                    &
         stat = alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Error : could not allocate space'
          write(6,*)'*** Error : subroutine WAVE_TRANSFORMATION'
          stop
        end if                          

    angbeachbp=0
    beach_breaker_os=0
    diffraction_theta=0
    active_groynes=0

      do wp = 1, nWavePoints
        if (wp < nWavePoints) then
            lastWPsection = fswp(wp + 1) - 1
        else
            lastWPsection = nQpoints
        end if
            
          ! find the depth of the wave point, at high tide
          depthwp(wp) = depthwp_msl(wp) + &
              AVERAGE_VALUE(tidal_amp(min(fswp(wp), nYsections): min(lastWPsection, nYsections))) + heightsurge
      end do

      ! find the depth of the offshore contour at high tide
      depthosc = depthOSCMsl + AVERAGE_VALUE(tidal_amp) + heightsurge
      
      if(nQpoints == 1)then
        i=1
        angbeachbp(i) = baselineAngle_radians
      else
        do section = 1, nYsections
          ! Estimate the offset of the breaker point
          beach_breaker_os(section) = max(beachy(msl_e,section),cliffy(msl_e,section))
        end do
        do qPoint = 2, nQpoints - 1
          ! Estimate the angle of the beach at q points
          angbeachbp(qPoint) = baselineAngle_radians - atan((beach_breaker_os(qPoint - 1)-beach_breaker_os(qPoint))/ dx)
        end do
        angbeachbp(1) = angbeachbp(2)
        angbeachbp(nQpoints) = angbeachbp(nQpoints - 1)
      end if

    ! For the look-up tables
    idepthosc = nint(depthosc * 100.0D0)
    it = nint(((nint(t*10.0D0)/10.0D0)-0.9D0)*10.0D0)
    
    do counter1 = 1, nQpoints    ! Loop on Q sections
      ! For the look-up table
      idepthwp = nint(depthwp(wpr(counter1)) * 100.0D0)    ! MWPM
      ! Shoaling coefficient
      ks = sqrt(cglut(it, idepthwp) / cglut(it, idepthosc))
      ! To find the quadrant that the wave approaches form
      alphaquad = anglewp(wpr(counter1)) - (angleosc(counter1) - pi) ! MWPM
      !alphaquad = anglewp_temp(counter1) - (angleosc(counter1) - pi) ! MD 300304
      if(alphaquad > two_pi)then
        alphaquad  = alphaquad - two_pi
      endif

      quadrant = 1
      quadswitch = 1

      if(alphaquad > half_pi)then
        quadrant = 2
        quadswitch = -1;
      end if
      
      if(sin(alphaquad) < 0.000001)then ! Quadrant 3 or 4, wave conditions assumed unchanged. 
        hiosc(counter1) = hwp(wpr(counter1)) ! MWPM
        !hiosc(counter1) = hwp_temp(counter1) ! MD 300304
         angleoscbp(counter1) = anglewp(wpr(counter1))+pi 
         !angleoscbp(counter1) = anglewp_temp(counter1)+3.1416 ! MD 300304
         ! Adds PI to convert 
        ! angle from 'direction coming from' to 'direction travelling to' ! MWPM
      else
        alphaoosc = quadswitch * (half_pi - alphaquad)
        alphaiosc = asin(sin(alphaoosc)*clut(it, idepthosc)/clut(it, idepthwp))
        ! Refraction coefficient
        kr = sqrt(cos(alphaoosc)/cos(alphaiosc))
        hiosc(counter1) = hwp(wpr(counter1)) * ks * kr    ! MWPM
        !hiosc(counter1) = hwp_temp(counter1) * ks * kr    ! MD 300304
        if(quadrant == 1)then
          angleoscbp(counter1) = angleosc(counter1) - alphaiosc + half_pi
        else
          angleoscbp(counter1) = angleosc(counter1) + alphaiosc + half_pi
        end if
        if(hiosc(counter1) > depthosc * breakrat)then
          hiosc(counter1) = hwp(wpr(counter1))    ! MWPM
          !hiosc(counter1) = hwp_temp(counter1)    ! MD 300304
          angleoscbp(counter1) = anglewp(wpr(counter1))+pi    ! MWPM
          !angleoscbp(counter1) = anglewp_temp(counter1)+pi    ! MD 300304
                    
          write(6,*)'Breaking beyond the offshore contour and seaward limit of the beach' 
        end if
      end if

      dbreak_wt(counter1) = hiosc(counter1)/breakrat
      dbreak_wt(counter1) = (nint(dbreak_wt(counter1) * 1000.0d0))/1000.0d0
    end do
! End of the new loop to turn the osc into a vector

    do counter1 = 1, nQpoints
!      repeatswitch_kd = 1
!      kd_old = 1
      kd = 1
      alphaquad = angleoscbp(counter1) - angbeachbp(counter1)
      quadrant = 1
      quadswitch = 1
      if(alphaquad > half_pi)then
        quadrant = 2
        quadswitch = -1
      end if

      if(sin(alphaquad) < 0.00001D0)then
        h_qpoint(counter1) = 0.01D0
        cgroup_qpoint(counter1) = 0.001D0
        anglewob_qpoint(counter1) = 999.0D0
        setup_qpoint(counter1) = 0.001D0
        setdown_qpoint(counter1) = 0.001D0
        theta_qpoint(counter1) = 999.0D0
      else
        alphaobp = quadswitch * (half_pi - alphaquad)
        repeatswitch = 1
        repeatcounter = 0
        do while (repeatswitch > 0)
          repeatcounter = repeatcounter + 1

          dbreak_new = hiosc(counter1)*kd/breakrat
          idepth = nint(100.0D0 * dbreak_new)
          
          if (idepth > idepthosc)then
            idepth = idepthosc
            repeatswitch = 0
          endif
          
          
          if(idepth == 0)idepth = 1
          c = clut(it, idepth)
          cg = cglut(it, idepth)
          ks = sqrt(cglut(it, idepthosc) / cg)
          alphaibp = asin(sin(alphaobp) * c / clut(it, idepthosc))
          kr = sqrt(cos(alphaobp)/cos(alphaibp))
          hbp = hiosc(counter1) * kd * ks * kr

          dbreak_new = hbp/breakrat
          if(abs(dbreak_wt(counter1) - dbreak_new) < 0.01D0)then
            repeatswitch = 0
          else
             dbreak_wt(counter1) = (nint(dbreak_new * 1000.0D0))/1000.0D0
          end if
        end do 

        h_qpoint(counter1) = hbp
        cgroup_qpoint(counter1) = cg
        if(quadrant == 1)then
          anglewob_qpoint(counter1) = angbeachbp(counter1) - alphaibp + half_pi
          theta_qpoint(counter1) = -1.0D0 * alphaibp
        else
          anglewob_qpoint(counter1) = angbeachbp(counter1) + alphaibp + half_pi
          theta_qpoint(counter1) = alphaibp
        end if

        k = two_pi/(c*t)

        if( h_qpoint(counter1) /= 0.0D0)then
          directionfactor = cos(alphaibp)
          setdown_qpoint(counter1) = -0.125D0*k* h_qpoint(counter1)* h_qpoint(counter1)/sinh(2.0D0*k* h_qpoint(counter1)/breakrat)
          setup_qpoint(counter1) = directionfactor*(setdown_qpoint(counter1) + (0.375D0* h_qpoint(counter1)*breakrat))
        else
          setdown_qpoint(counter1) = 0
          setup_qpoint(counter1) = 0
        end if
      end if
    end do





    if(nQpoints == 1)then
      i=1
    else
      do section = 1, nYsections
        anglewob(section) = (anglewob_qpoint(section)+anglewob_qpoint(section+1))/2.0d0
        h(section) = (h_qpoint(section)+h_qpoint(section+1))/2.0d0
        cgroup(section) = (cgroup_qpoint(section)+cgroup_qpoint(section+1))/2.0d0
        setup(section) = (setup_qpoint(section)+setup_qpoint(section+1))/2.0d0
        setdown(section) = (setdown_qpoint(section)+setdown_qpoint(section+1))/2.0d0
        theta(section) = (theta_qpoint(section)+theta_qpoint(section+1))/2.0d0
      end do
    end if
          
        
!*** Deallocate temp arrays

        deallocate(angbeachbp,beach_breaker_os,diffraction_theta,        &
            active_groynes,hiosc,angleoscbp,dbreak_wt,depthwp, stat=alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Unexpected deallocation error', alloc_error
          write(6,*)'*** for subroutine WAVE_TRANSFORMATION'
        endif                   

        return
        end subroutine WAVE_TRANSFORMATION             

end module wave_transformation_module
