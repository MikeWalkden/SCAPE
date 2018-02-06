!**********************************************
!*
!*       xshore_dist_t.f90       
!*
!**********************************************
!**********************************************
!
!>    This module of SCAPE calculates cross-shore 
!>    distributions of:
!>    a) Longshore drift
!>    b) Erosion
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

module xshore_dist_t_module
    implicit none

contains

!>    This routine does a numerical integration, stepping the 
!>    shape function through the tidal cycle and summing it
!>    each step to produce a distribution representing the 
!>    whole tidal cycle.
!>
!>    This process has to represent the different amount of 
!>    time that the tide spends at the different levels, more 
!>    at the extremes thena in the middle.
    subroutine xshore_dist_t(half_tide, sf, nce, dbreak,    &
        dsf, setup, setdown, heightsurge, np_sf, lowlim, uplim,    &
        cliffheights, msl_m, this_dist, bottom, top, sfid, success_code)
    
    use xshore_dist_t_data ! (JT241007)
    use setup_module
    use exceptions
    use local_data, only: logString
    use local_data, only: feLogString
    use utilities_module, only: feLog
    implicit none

    integer, intent(in) :: half_tide !< Tidal amplitude measured in shape function increments
    real(kind=double), dimension(50, 2), intent(in) :: sf !< The shape function (of drift or erosion,
        !< under breaking waves with a static water level
    integer, intent(in) :: nce !< Number of cliff elements
    real(kind=double), intent(in) :: dbreak !< Depth of breaking
    real(kind=double), intent(in) :: dsf !< Dx of the shape function (see elsewhere)
    real(kind=double), intent(in) :: setup
    real(kind=double), intent(in) :: setdown
    real(kind=double), intent(in) :: heightsurge !< Not used, left at zero
    integer, intent(in) :: np_sf !<  Number of points in the shape function
    integer, intent(in) :: lowlim !< The lower limit of the model (element 2)
    integer, intent(in) :: uplim !< The upper limit of the model (nce)
    real(kind=double), dimension(nce), intent(in) :: cliffheights !< Heights of cliff elements
    real(kind=double), intent(in) :: msl_m !< Mean sea level in metres
    real(kind=double), dimension(nce), intent(inout) :: this_dist !< The derived distribution
    integer, intent(inout) :: bottom !< Lower limit of the derived distribution
    integer, intent(inout) :: top !< Upper limit of the derived distribution
    !JT080108 optimisation for caching on half_tide        
    integer, intent(in) :: sfid !< Shape function identifier, 1=sf_drift & 2=sf_erode (JT241007)
    integer, intent(out) :: success_code !< 0: ok, -ive fatal error, +ive warning

!*** Local variables

    real(kind=double) :: first,second
    real(kind=double) :: jj
    real(kind=double) :: distlow,disthigh
    real(kind=double) :: dzdist

    integer :: j, k
    integer :: npp
    integer :: counter1, counter2
    integer :: alloc_error
    integer :: xShoreCacheOffset
    logical :: alreadyCached

    success_code = 0
    
    ! npp the number of points in the distribution, with dx = the 
    ! dx of the shape function. 
    npp = np_sf + 2 * half_tide

    !*** Allocate arrays
    
    !*** Allocate arrays if they not big enough (JT221107a)
    if (half_tide > allocCeil)then
        deallocate(xshore, tide_duration, start_fnctn,        &
                   end_fnctn, temp, stat=alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Unexpected deallocation error'
          write(6,*)'*** for subroutine xshore_dist_t'
        endif

        allocate(xshore(npp),    &
         tide_duration(2 * half_tide + 1), temp(npp),    &
         start_fnctn(npp), end_fnctn(npp),    &
         stat = alloc_error)
        success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'subroutine xshore_dist_t')
        if (success_code < 0) return
     
    endif
    !*** end of alloc  (JT221107a)               

    ! calculate the offset where the cache will be stored (JT221107)
    if (sfid == 1)then
        xShoreCacheOffset = driftOffset+((np_sf+half_tide)*(half_tide-1))
    else
        xShoreCacheOffset = (np_sf+half_tide)*(half_tide-1)
    endif

    ! retrieve xshore if already cached (JT241007) mem optimised (JT221107)
    alreadyCached = half_tide < cacheCeil
    if (alreadyCached) alreadyCached = xshorecache_flag(sfid, half_tide)
    
    if (alreadyCached) then
        xshore(1:npp)=xshorecache((xShoreCacheOffset+1):(xShoreCacheOffset+npp))
    else
        ! otherwise calculate xshore (JT241007)
    
        ! Account for the sin variation of the tide
        ! The tide is symmetrical so only have to deal with a part of it
        ! Calculate the duration of the time spent at each elevation
    
        tide_duration = 0.0D0
        counter1 = half_tide
        counter2 = half_tide+1
        tide_duration(1) = 1
    
        do j=1,half_tide
            jj = j
            counter1 = counter1+1
            counter2 = counter2-1
            first = (jj-1)/half_tide
            second = jj/half_tide
            tide_duration(counter1) = 2.0d0*(asin(second)-asin(first))
            tide_duration(counter2) = tide_duration(counter1)
        end do
    
        tide_duration = tide_duration / 6.283D0
    
        ! The start and end functions define start and end positions in the shape function.
        start_fnctn(1:npp) = 1
        do j=1,np_sf+1
            start_fnctn(npp-np_sf-1+j) = j
        end do
        do j=1,npp
            if(j>=np_sf)then
                end_fnctn(j) = np_sf
            else
                end_fnctn(j) = j
            end if
        end do
    
        xshore(1:npp) = 0.0d0
     
        ! start of corrected version see 105 wkbk 9
        counter1 = 0
        do j = 1, npp - 1
            counter1 = counter1+1
            if(counter1 > 2*half_tide) then
                counter1 = 2*half_tide
            endif
            !    check for unreal situation of no tide
            if(counter1 == 0)counter1=1
            
            locatesf(1:(end_fnctn(j)-start_fnctn(j)+1))=sf(start_fnctn(j):end_fnctn(j),1)
    
            duration_sub(1:(end_fnctn(j)-start_fnctn(j)+1)) = tide_duration(counter1:counter1-(end_fnctn(j)-start_fnctn(j)):-1)
    
            do k=1,(end_fnctn(j)-start_fnctn(j)+1)
                xshore(j) = xshore(j) + (locatesf(k)*duration_sub(k));
            enddo
        enddo
        ! end of corrected version see 105 wkbk 9
     
        ! store xshore if under cache ceiling (JT241007) mem optimised (JT221107)
        if (half_tide < cacheCeil) then
            xshorecache((xShoreCacheOffset + 1): (xShoreCacheOffset + npp)) = xshore(1:npp)
            xshorecache_flag(sfid, half_tide) = .true.
        endif
    
    ! end of caching brackets (JT241007)
    endif

    do j = 1, npp
        temp(j) = xshore(npp + 1 - j)
    end do
    xshore(1: npp) = temp(1: npp)

    this_dist=0.0d0

    dzdist = dbreak * dsf
    distlow = msl_m + heightsurge + setdown - ((half_tide + np_sf)*dzdist)
    disthigh = msl_m + heightsurge + setup + (half_tide * dzdist)
    top = nint(disthigh / dz)
    bottom = nint(distlow / dz)
    
    if ((top < lowlim) .or. (bottom > uplim - 1)) then
        write(logString, fmt='(A,I4,A,I4,A,I4,A,I4,A)') &
            ' xshore_dist_t : distribution (', bottom, ':', top, ') model (', lowlim, ':', uplim, ')'
        call RAISE_EXCEPTION(logString)
        success_code = -1
        return
    end if
    
    if (bottom < lowlim) bottom = lowlim
    if (top > uplim-1) top = uplim - 1

    distmarkers(bottom:top) = nint((cliffheights(bottom:top) - distlow) / dzdist)
    do j = bottom, top ! Remove possible rounding errors
        if(distmarkers(j) < 1) distmarkers(j)=1
        if(distmarkers(j) > npp-1) distmarkers(j)=npp-1
    end do

    this_dist(bottom:top) = xshore(distmarkers(bottom:top))/        &
        sum(xshore(distmarkers(bottom:top)))

    ! If have alloced larger than allocCeil reset to allocCeil (JT221107a)
    if (half_tide > allocCeil)then
        deallocate(xshore, tide_duration, start_fnctn, end_fnctn, temp, stat=alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Unexpected deallocation error'
          write(6,*)'*** for subroutine xshore_dist_t'
        endif

        npp = max_np_sf + 2 * allocCeil
        allocate(xshore(npp),    &
         tide_duration(2 * allocCeil + 1), temp(npp),    &
         start_fnctn(npp), end_fnctn(npp),    &
         stat = alloc_error)
        success_code = SUCCESS_CODE_FROM_ALLOCATION(alloc_error, 'subroutine xshore_dist_t end')
        if (success_code < 0) return
    endif
    ! end of dealloc

    return
    end subroutine xshore_dist_t

end module xshore_dist_t_module
