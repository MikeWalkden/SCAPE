!**********************************************
!*       component_data.f95       
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

module component_data

implicit none
save

integer, parameter :: double = kind(1.0D0)

integer, parameter :: substitutions = 12
integer, parameter :: keyLength = 30
integer, parameter :: valueLength = 256

! Path to the folder containing template.omi and templates.fcs files
character(valueLength) :: templatesPath = &
    'C:\\FortranProjects\\MESO_i\\Component\\'

! Path to the dll for the C# wrapper for the MESO_i engine
character(valueLength) :: engineWrapperPath = &
    'C:\\FortranProjects\\SCAPE\\Engine_Wrapper\\Engine_Wrapper\\bin\\x86\\Debug\\Engine_Wrapper.dll'

! Path to the FluidEarth SDK installed file FluidEarth2_Sdk.exe
character(valueLength) :: fluidEarthPath = &
    'C:\\Source\\FluidEarth2_Sdk\\bin\\x86\\Debug\\FluidEarth2_Sdk.exe'

! y-value for the model seaward boundary
real(double) :: seawardBoundary = 1000.0D0

! world coordinates identified
character(valueLength) :: worldCoordsId = ''

! local origin x-value in world coordinates
real(double) :: localOriginX = 0.0D0

! local origin y-value in world coordinates
real(double) :: localOriginY = 0.0D0
    

end module component_data
