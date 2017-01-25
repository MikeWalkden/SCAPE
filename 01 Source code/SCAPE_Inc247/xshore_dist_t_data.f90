!********************************************************
!*
!*       xshore_dist_t_data.f90     
!*
!********************************************************
!
!>    Global data shared by xshore_dist_t and RUN_cliffSCAPE
!>    to enable caching in xshore_dist_t.
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


    module xshore_dist_t_data
    
    implicit none
    save

    integer, parameter :: double = kind(1D0)
    private double

    ! Maximum half_tide to cache, higher values result in generally faster runtime
    ! but cache memory storage requirement increases proportionally to the square of cacheCeil.
    integer :: cacheCeil = 5000

    ! Record of which xshores have been cached for a given shapeFunctionID & half_tide
    logical, allocatable, dimension(:,:) :: xshorecache_flag

    ! Cache of xshore
    ! Originally (JT241007), three dimensions shapeFunctionID, half_tide & index within xshore.
    ! which is very wasteful each cached xshore stored in a space big enough for the largest
    ! as xshore is determined the number of points in the shape function + 2* half_tide this means
    ! each record was assigned (np_sf_drift+2*cacheCeil) space i.e. 10041 when the smallest needs only 42!
    ! New approach (JT221107) stores each element in the minimum space needed.
    ! This array will be of length (np_sf_erode+cacheCeil)*(cacheCeil-1) + ((np_sf_drift+cacheCeil)*(cacheCeil-1))
    ! erosion storage as (np_sf_erode+half_tide)*(half_tide-1)
    ! drift storage will be driftOffset+((np_sf_drift+half_tide)*(half_tide-1))
    real(kind=double), allocatable, dimension(:) :: xshorecache
    
    ! Index offset after which the drift data can be foumd in the cache (JT221107)
    ! i.e. (np_sf_erode+cacheCeil)*(cacheCeil-1)
    integer :: driftOffset

    ! Number of points in each shape function saved, to allow better cache storage (JT221107)
    integer :: np_sf_erode, np_sf_drift

    ! for use in allocating below variables (JT221107a)
    integer :: max_np_sf !largest np_sf i.e. 41
    integer :: allocCeil = 5000 ! if half tide is below this no allocation is needed,
    ! otherwise the below variables will need to be allocated

    ! local allocs from xshore_dist_t (JT221107a)
    real(kind=double), allocatable, dimension(:) :: xshore, tide_duration
    integer, allocatable, dimension(:) :: start_fnctn,end_fnctn
    real(kind=double), allocatable, dimension(:) :: temp
    
    real(kind=double), allocatable, dimension(:) :: locatesf
    real(kind=double), allocatable, dimension(:) :: duration_sub

    integer, allocatable, dimension(:) :: distmarkers  

    end module xshore_dist_t_data
