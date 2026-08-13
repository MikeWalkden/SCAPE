!**********************************************
!*       normal.f90
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



module normal_func
    implicit none

contains
        
        !> This is a SUBROUTINE that generates random numbers from a Normal distribution with
        !> mean zero and variance one.

        ! Stats of 1000 numbers checked by MW 061205
        
        SUBROUTINE Normal(random, nsim) ! code from Andrew Jeffrey
        
        implicit none
        
        integer :: nsim !< number of elements in the random vector
        integer, parameter :: double = kind(1D0)

        real(kind=double), dimension(nsim), intent(out) ::    random
        !< Each element of this output vector will contain a N(0,1) independent random number.
        
        real(kind=double), dimension(min(size(random),nsim)) ::    W
        real(kind=double), dimension(2*size(W)) ::                V
        logical, dimension(size(W)) ::                    flag

        integer N,i,M,number,dim1,dim2,dim3,index

        dim1=size(W)
        dim2=dim1+1
        dim3=dim1+dim1
        number=size(random)

        i=0
        do
            call random_number(V)
            V=2.0d0*V-1.0d0
            W=V(1:dim1)*V(1:dim1) + V(dim2:dim3)*V(dim2:dim3)
            flag=.false.

            do index = 1, dim1
                if (W(index) <= 1.0d0 .and. W(index) > 0.0d0) then
                    W(index) = log(W(index)) / W(index)
                    W(index) = sqrt(-W(index) - W(index))
                    V(index) = V(index) * W(index)
                    V(dim1 + index) = V(dim1 + index) * W(index)
                    flag(index) = .true. !true means it is part of N(0,1) series
                end if
            end do
            
            !where (W <= 1.0d0 .and. W > 0.0d0)
            !    W=log(W)/W
            !    W=sqrt(-W-W)
            !    V(1:dim1)=V(1:dim1)*W
            !    V(dim2:dim3)=V(dim2:dim3)*W
            !    flag=.true. !true means it is part of N(0,1) series
            !end where
            
            N=count(flag)
            V(1:N)=pack(V(1:dim1),flag)
            V(N+1:N+N)=pack(V(dim2:dim3),flag)
            M=min(N+N,number-i)
            random(i+1:i+M)=V(1:M)
            i=i+M
            if(i.ge.number) exit
        end do
        
        return

        END SUBROUTINE Normal

end module normal_func
