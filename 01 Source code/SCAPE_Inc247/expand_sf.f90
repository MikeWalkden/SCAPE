!**********************************************
!*
!*       expand_sf.f90       
!*
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

module expand_sf_module

implicit none

contains

    !> Widens an input shape function, which may describe the distribution of erosion or alongshore transport under a breaking wave field.
    subroutine expand_sf(sf, np_sf, sf_new, dsf, np_sf_new)

    integer, parameter :: double = kind(1D0)

    real(kind=double), dimension(:,:), intent(inout) :: sf !< original shape function
    integer, intent(in) :: np_sf !< Number of points in the shape function
    real(kind=double), dimension(:,:), intent(out) :: sf_new !< new shape function
    real(kind=double), intent(in) :: dsf !< Dx of the shape function
    integer, intent(in) :: np_sf_new !< Number of points in the new shape function

!*** Local variables
                       
    real(kind=double) :: ctod,ab,ac,be
    real(kind=double) :: abscissa,toler 

    integer :: i
    integer :: repeatswitch,counter

    toler=1.0d-6  

    ! Turn the shape function into a two dimensional array containing both coordinates
    do i=1,np_sf
      sf(i,2) = sf(i,1)
      sf(i,1) = (i-1)*dsf
    end do

    sf_new = 0

    do i = 1, np_sf_new  ! Number of points in the new Shape Function
      sf_new(i,1) = sf(1,1) + ((i-1)*(sf(np_sf,1) - sf(1,1)))/(np_sf_new-1)
    end do
    
    do i=1,np_sf_new
      abscissa = sf_new(i,1)
      repeatswitch = 1
      counter = 1
      do while(repeatswitch > 0 )
        if (abs(abscissa-sf(counter,1)) < toler)then
          sf_new(i,2) = sf(counter,2)
          repeatswitch = 0
        else if(sf(counter,1) > abscissa)then
          ctod = sf(counter,2) - sf(counter-1,2)
          ab = abscissa - sf(counter-1,1)
          ac = sf(counter,1) - sf(counter-1,1)
          be = ctod*ab/ac
          repeatswitch = 0
          sf_new(i,2) = be + sf(counter-1,2)
        end if
        counter = counter+1
      end do
    end do

    return
    end subroutine expand_sf


    !> Read shape function and set number of points from the number of values in the file.
    
    !> Number of points is limited by the first dimension of the sf array.
    subroutine read_sf(fileUnit, path, fileName, np_sf, sf, success_code)
    use exceptions
    use utilities_module
    implicit none

    integer, parameter :: double = kind(1D0)

    integer, intent(in) :: fileUnit !< FORTRAN file unit number
    character(*), intent(in) :: path !< Path to file
    character(*), intent(in) :: fileName !< name of the file containing the shape function definition
    integer, intent(out) :: np_sf !< Number of points in the shape function
    real(kind=double), dimension(:,:), intent(inout) :: sf !< shape function
    integer, intent(inout) :: success_code !< 0 indicates success

    integer :: i
    integer :: ierr
    
    np_sf = size(sf, 1) ! Assume there are enough values in the file to fill the array

    success_code = do_open_text_input_file(fileUnit, path, fileName)
        
    do i = 1, np_sf
        read(fileUnit, *, iostat = ierr) sf(i, 1)
        sf(i, 2) = 0.0d0

        if (ierr == -1) then ! Set number of points from number of values in the file
            np_sf = i - 1
            exit
         else
            success_code = SUCCESS_CODE_FROM_IOSTAT(fileUnit, ierr, 'Reading shape function from '//fileName)
            if (success_code < 0) return
         end if
    end do

    close(fileUnit)

    end subroutine read_sf

end module expand_sf_module
