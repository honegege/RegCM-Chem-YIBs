      module disturbance

      use yibs_types

      implicit none


      contains
      !*********************************************************************

      subroutine fire_frequency_cell(dtsec,time, ycell)
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: time  !Greenwich Mean Time
      type(ycelltype) :: ycell

      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------


      end subroutine fire_frequency_cell

      !*********************************************************************
      subroutine calc_cell_disturbance_rates(dtsec,time,ycell)
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: time  !Greenwich Mean Time
      type(ycelltype) :: ycell
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------

      end subroutine calc_cell_disturbance_rates


      !*********************************************************************

      end module disturbance
