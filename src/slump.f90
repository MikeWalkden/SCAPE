!**********************************************
!*
!*       slump.f90       
!*
!**********************************************
!
!>    This subroutine of SCAPE slumps undercut cliff
!>
!>    The cliff is represented in two parts, as the upper 
!>    region of the platform, represented by elements, and
!>    as a simple block of material on top of the upper 
!>    element.
!>
!>    Material slumps from the cliff and adds to the talus
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


module slump_module
    implicit none

contains

    !> Slump the given beach section
    subroutine slump(section)

    use global_data
    use setup_module
    implicit none

    integer, intent(in) :: section !< Beach section index
    
!*** Local variables

    real(kind=double) :: overhang

    integer :: j

    cliffy(1, section) = cliffy(2, section) + 1

    ! Slump the cliff toe (upper platform)
    do j=1,nce-1
        overhang = cliffy(j+1, section) -  cliffy(j, section)

        if(overhang > 0)then
            cliffy(j+1, section) = cliffy(j, section)
            tvolume( section) = tvolume( section) + overhang * dx * dz
        end if
        
    end do
    
    ! Slump the main body of the cliff
    tvolume(section) = tvolume(section) + ((clifftoe_pos(section) - cliffy(nce, section)) &
        * dx * cliff_heights(section))

    clifftoe_pos(section) = cliffy(nce, section)

    return
    end subroutine slump
      
end module slump_module
