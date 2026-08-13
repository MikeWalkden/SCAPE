!**********************************************
!*       finish.f95       
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

module finish_module
    implicit none

integer, parameter :: double = kind(1D0)
private double

contains

    !> FluidEarth entry point, also called from RUN_SCAPE.
    
    !> Find cliff top positions.
    !> Close all files.
    !> Deallocate memory.
    subroutine finish(success_code)
    
        use setup_module
        use global_data
        use local_data
        use xshore_dist_t_data
        use regfalls3_module
        use SLR_module
        use vertical_grid
        implicit none
        
        integer, intent(out) :: success_code ! 0: ok, -ive fatal error, +ive warning

        logical :: opened
        
        write(6, fmt='(A,A,A,I6)') ' SCAPE+ ', trim(FILE_PATH), ' Finishing at year:', year ! To screen

! Now find cliff top positions
    if (nCliffSimulations > 0) then
        call RegFalls3(nYsections, endYear - startYear, retreat, 'CLIFF_TOP_POS', use_random_seed)
    end if

    !******************************************
    !
    ! close files 
    !
    !******************************************

    
    close(saveRockContourUnit)
    close(saveShoreContourUnit)
    close(saveAnnBvolUnit)
    close(saveFinesVolUnit)
    close(saveSedTransAnnUnit)
    close(savePotTransAnnUnit)
    close(saveBeachAddAnnUnit)
    
    if (vellingaActive) then
        close(saveVellLevelsUnit)
    end if
    
    close(saveSeawallActiveUnit)
    close(saveGroyneEffectUnit)
    
    if (firstYearProfileOutput /= -1000000) then
        close(saveRockProfilesUnit)
        close(saveBeachProfilesUnit)
    end if

    close(wavesAndTidesFileUnit)
    close(slFileUnit)
    close(saveSedFluxLeftUnit)
    close(savePotSedFluxLeftUnit)
    close(saveSedFluxRightUnit)
    close(savePotSedFluxRightUnit)

!*** Deallocate arrays

    if(nQpoints > 1)then

          deallocate(bruun, upbruun, lowbruun, bvolume, tvolume, sand_fraction, &
          cliffTopElevations, cliff_heights, clifftoe_pos,        &
          talus_offsets, offshore_st, resistance,        &
          cerck1, groyneEffectAnnualSum, seawall_os, &
          seawall_active, &
          QAnn,QpAnn,        &
          stat=alloc_error)
          if(alloc_error /= 0)then
            write(6,*)'*** Unexpected deallocation error'
            write(6,*)'*** for bruun'
          end if          

          if (includeGroynes) then
              deallocate(groyne_construction, groyne_removal, stat=alloc_error)
              if(alloc_error /= 0)then
                write(6,*)'*** Unexpected deallocation error'
                write(6,*)'*** for deallocategroyne_construction'
              end if          
          end if
    
          if (includeSeawall) then
              deallocate(seawall_construction, seawall_removal, stat=alloc_error)
              if(alloc_error /= 0)then
                write(6,*)'*** Unexpected deallocation error'
                write(6,*)'*** for seawall_construction'
              end if          
          end if

          deallocate(beachelev, berm, beachdepth, beachdepth1, &
              bTop, bBottom, ubTop, ubBottom, lbTop, lbBottom, stat=alloc_error)
          if(alloc_error /= 0)then
            write(6,*)'*** Unexpected deallocation error'
            write(6,*)'*** for beachelev, berm'
          end if

          deallocate(section_offsets, beachos, upbeach_os,        &
          lowbeach_os,upbvolume,obvolume, stat=alloc_error) ! MW030304
          if(alloc_error /= 0)then
            write(6,*)'*** Unexpected deallocation error'
            write(6,*)'*** for section_offsets, beachos'   
          end if          

    end if

    deallocate(clut, cglut, stat=alloc_error)
    if(alloc_error /= 0)then
      write(6,*)'*** Unexpected deallocation error'
      write(6,*)'*** for clut, cglut'   
    end if          

    if (useVAgrid) call finish_VA_grid()
        
    deallocate(setup,setdown,theta,cgroup,h,anglewob, &
               setup_qpoint,setdown_qpoint,theta_qpoint,cgroup_qpoint,h_qpoint,anglewob_qpoint,     &
               top_drift,bot_drift,sect_drift,dbreak,half_tide,     &
               beach_addition,fines_addition,total_beach_addition,this_tide_drift, sediment_influx, &
               stat=alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Unexpected deallocation error'
          write(6,*)'*** for setup'
        end if                 

    deallocate(cliffy,beachy,talusy,upbeachy,lowbeachy,cliffheights,     &
           sect_erode_new,sect_erode,surface,        &
           top_erode,bot_erode,np_cliff,    &
           sect_height,startprofile,startprofiles,angleosc,    &
           ubvolume,lowbvolume,stormflag, initial_cliff_toe, retreat, &
           stat=alloc_error)
    if(alloc_error /= 0)then
      write(6,*)'*** Unexpected deallocation error'
      write(6,*)'*** for cliffy, beachy & cliffheights'
    end if

    deallocate(wpr,hwp,hwp_temp,anglewp,anglewp_temp,ann_av_bvol,        &   ! MD 300304
           stat=alloc_error)
    if(alloc_error /= 0)then
      write(6,*)'*** Unexpected deallocation error'
      write(6,*)'*** for wpr,hwp,hwp_temp,anglewp,anglewp_temp,ann_av_bvol'
    end if

    deallocate(beach_temp, cliff_temp, &
           stat=alloc_error)
    if(alloc_error /= 0)then
      write(6,*)'*** Unexpected deallocation error'
      write(6,*)'*** for beach_temp, cliff_temp'
    end if

    if (vellingaActive) then
        deallocate(vellinga_sections, vxs, vell_berm, vellinga, v_ref, v_level_odn, &
               stat=alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Unexpected deallocation error'
          write(6,*)'*** for vxs,vell_berm,vellinga,v_ref,v_level_odn'
        end if
    end if

    call DEALLOCATE_PARAMETERS()


! JT221007
    deallocate(xshorecache, &
        xshorecache_flag, &
        stat=alloc_error)
    if(alloc_error /= 0)then
      write(6,*)'*** Unexpected deallocation error'
      write(6,*)'*** for xshorecache'
    end if

    !*** Deallocate arrays from platform erosion (JT221107a)

       deallocate( this_erode, m,            &
                   temp5,bdepth,battenuation,        &
                   crossbeach_st, &
                   seawall_effect, &
                   stat=alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Unexpected deallocation error'
          write(6,*)'*** subroutine PLATFORM_EROSION'
        end if

    !*** Deallocate arrays from xshore_dist_t (JT221107a)
       deallocate(distmarkers, xshore, tide_duration, start_fnctn,        &
                   end_fnctn, locatesf, duration_sub, temp, stat=alloc_error)
        if(alloc_error /= 0)then
          write(6,*)'*** Unexpected deallocation error'
          write(6,*)'*** for subroutine xshore_dist_t'
        end if

        if (feLogActive) then
            inquire(unit=98, opened=opened)
            
            if (opened) then
                close(98)
            end if
        end if

        close(logFileUnit)
        success_code = 0
        
    end subroutine finish

end module finish_module