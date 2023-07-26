      module yibs_make_struct
!@sum yibs_make_struct  Off-line module for generating an YIBS vegetation
!@+   data structure, given ascii input file of ycell-patch-cohort structure.
      use yibs_types
      use yibs_const
      use cohorts
      use patches
      use ycells

      public yibs_struct_readcsv
      
      contains

!************************************************************************
      subroutine skipstar (iu_ystruct)
      !* Skip a comment line in an YIBS structure csv file.
      integer :: iu_ystruct
      !---
      character :: check

      read(iu_ystruct,*) check
      !write(*,*) 'ss',check
      do while (check.eq.'*') 
        read(iu_ystruct,*) check
      end do
      backspace(iu_ystruct)
      end subroutine skipstar

!************************************************************************


      subroutine read_ycell_struct ( ecp, iu_ystruct )
      type(ycelltype) :: ecp
      integer :: iu_ystruct
      !---
      !real*8 :: stext1,stext2,stext3,stext4,stext5 !soil_texture
      
      !read(iu_ystruct,*) stext1,stext2,stext3,stext4,stext5 !soil_texture
      !ecp%soil_texture(1) = stext1
      !ecp%soil_texture(2) = stext2
      !ecp%soil_texture(3) = stext3
      !ecp%soil_texture(4) = stext4
      !ecp%soil_texture(5) = stext5
      call skipstar(iu_ystruct)
      read(iu_ystruct,*) ecp%soil_texture 
      end subroutine read_ycell_struct

!************************************************************************
      

      subroutine read_patch_struct (pp,iu_ystruct )
      type(patch) :: pp
      integer :: iu_ystruct
      !---
      integer :: layer
!      real*8 :: age, area, soil_type
!     real*8 :: May also add soil texture
!      real*8 :: tc1,tc2,tc3,tc4,tc5,tc6,tc7,tc8,tc9,tc10,tc11,tc12 !TpoolC
!      real*8 :: tn1,tn2,tn3,tn4,tn5,tn6,tn7,tn8,tn9,tn10,tn11,tn12 !TpoolN

      call skipstar(iu_ystruct)
      read(iu_ystruct,*) pp%age, pp%area, pp%soil_type
      do layer = 1,N_CASA_LAYERS
        call skipstar(iu_ystruct)
        read(iu_ystruct,*) pp%Tpool(CARBON,:,layer)
        read(iu_ystruct,*) pp%Tpool(NITROGEN,:,layer)
      end do
      end subroutine read_patch_struct

!************************************************************************

      subroutine read_cohort_struct ( cop, iu_ystruct )
      type(cohort) :: cop
      integer :: iu_ystruct
      !---

      call skipstar(iu_ystruct)
      read(iu_ystruct,*) cop%pft
      call skipstar(iu_ystruct)
      read(iu_ystruct,*) cop%n,cop%h,cop%crown_dx,cop%crown_dy,
     &     cop%dbh,cop%root_d,cop%clump
      end subroutine read_cohort_struct


      end module yibs_make_struct
