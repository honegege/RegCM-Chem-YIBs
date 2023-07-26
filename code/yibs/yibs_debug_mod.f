      module yibs_debug_mod
      implicit none

      integer, parameter :: NPFT=16
      integer, parameter :: SIZE_YIBS_DEBUG=NPFT*(16+11 +1)

      type yibs_debug
        sequence
        real*8 :: vf(NPFT)      ! 1
        real*8 :: Anet(NPFT)    ! 2
        real*8 :: Atot(NPFT)    ! 3
        real*8 :: Rd(NPFT)      ! 4
        real*8 :: GCANOPY(NPFT) ! 5
        real*8 :: TRANS_SW(NPFT)! 6
        real*8 :: LAI(NPFT)     ! 7
        real*8 :: Resp_fol(NPFT)! 8
        real*8 :: Resp_sw(NPFT) ! 9
        real*8 :: Resp_lab(NPFT) ! 10 
        real*8 :: Resp_root(NPFT)! 11
        real*8 :: Resp_maint(NPFT) ! 12
        real*8 :: Resp_growth_1(NPFT) ! 13
        real*8 :: Resp_growth(NPFT) ! 14
        real*8 :: GPP(NPFT)     ! 15
        real*8 :: R_auto(NPFT)  !16

        real*8 :: total(NPFT)  ! 17
        real*8 :: Resp_soil(NPFT)  ! 25
        real*8 :: phenofactor(NPFT) ! 26
        real*8 :: betad(NPFT) ! 27

        ! cell vars................. 28
        real*8 :: soiltemp_10d ! 1
        real*8 :: airtemp_10d  ! 2
        real*8 :: par_10d      ! 3
        real*8 :: gdd          ! 4
        real*8 :: ncd          ! 5
        real*8 :: sgdd         ! 6
        real*8 :: ld           ! 7
        real*8 :: ld0          ! 8
        real*8 :: foo(NPFT-8)

       end type yibs_debug

      type(yibs_debug), target :: yibs_d
      !real*8 :: yibs_dl(SIZE_YIBS_DEBUG)
      !equivalence(yibs_d,yibs_dl)

      contains

      subroutine get_yibs_debug_ptr(ptr)
      use ISO_C_BINDING
      real*8, pointer :: ptr(:)
      !---
      type (C_PTR) :: cptr

      cptr = c_loc(yibs_d)
      call c_f_pointer( cptr, ptr, (/SIZE_YIBS_DEBUG/) )

      end subroutine get_yibs_debug_ptr

      end module yibs_debug_mod
