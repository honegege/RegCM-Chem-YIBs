#include "rundeck_opts.h"

      module yibs_prescribed_drv

      !*********************************************************************
!@sum !*    SUBROUTINES TO READ IN prescribed VEGETATION DATA SETS 
!@+   !+      for initialization and updating prescribed cover like crops.
!@+   !*    Array data only, no ycells or patches info.
!@+   !*    Interfaces with yibs_prescr_veg for Matthews pft-level calculations
!@+   !*      or with yibs_prescribed_drv_geo for geographic initialization.
      !*********************************************************************

      use yibs_const
      use yibs_pfts
      use yibs_prescr_veg
      use mod_dynparam
      use mod_regcm_types
	  
      implicit none
      private
      save

      public 
     &     yibs_init_params,
     &     init_canopy_physical,
     &     prescr_vegdata,
     &     prescr_veg_albedodata

      public init_yibs_laidata, init_yibs_hdata
     &     ,prescr_get_yibs_plant  ,prescr_get_soilpools

      public prescr_get_laidata,
     &     prescr_get_cropdata,
     &     prescr_get_soil_C_total
      public prescr_get_hdata
      public prescr_calc_canopy_geometry, prescr_get_carbonplant
      public prescr_get_pft_vars


#ifdef MIXED_CANOPY
      public yibs_struct_get_phys
#endif

      contains

!***************************************************************************
      subroutine yibs_init_params()
!@sum Initialize some YIBS parameters.
      use yibs_prescr_veg, only : init_params

      !call prescr_calcconst !renamed init_params
      call init_params()

      end subroutine yibs_init_params
!***************************************************************************
      subroutine init_canopy_physical(
     & I0,I1,J0,J1,Ci_ini, CNC_ini, Tcan_ini, Qf_ini)
!@sum For old Friend & Kiang (2005) biophysics. Initialize LSM outputs.
      integer,intent(in) :: I0,I1,J0,J1
      real*8, DIMENSION(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini

      Ci_ini(:,:) = 0.0127d0
      CNC_ini(:,:) = 0.d0
      Tcan_ini(:,:) = 0.d0        !Should be a forcing from land surface model.
      Qf_ini(:,:) = 0.d0          !Should be a forcing from land surface model.

      end subroutine init_canopy_physical
      
!***************************************************************************
      subroutine prescr_get_soilpools(IM,JM,I0,I1,J0,J1,
     &     soil_C_total, Tpool_ini, lms)
!@sum Prescribe initial partitioning of soil carbon pools given total.
      !* For global runs, this routine gets soil_C_total from subroutine 
      !get_soil_C_total,which reads in total soil pool amounts (measured). 
      !Then individual soil pool fractions (modeled pft-dependent values 
      !from spinup runs by PK),are used to prescribe individual amounts.
      !**all carbon amounts should be in g/m2** -PK 12/07, NK 7/11
      !
      !* For site runs, this routine reads in all soil carbon fractions.

      use FILEMANAGER, only : openunit,closeunit
      integer,intent(in) :: IM,JM,I0,I1,J0,J1
      type(lm_state) , intent(in) :: lms
      real*8,intent(in) ::
     &     soil_C_total(N_CASA_LAYERS,I0:I1,J0:J1)
      real*8,intent(out) :: 
     &      Tpool_ini(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,  
     &                I0:I1,J0:J1)!prescribed soil pools, g/m2
      !-----Local------
!      first 3 for eventually reading in globally gridded dataset, e.g. ISRIC-WISE
      integer :: iu_SOILCARB, iunit
      integer :: n,p,nn,k
      real*8, dimension(N_PFT,NPOOLS-NLIVE,N_CASA_LAYERS) :: Cpool_fracs
      real*4 :: cpool1(IM,JM)
      real*8 :: cpool_pft(IM,JM,N_PFT)
      character*80 :: head,fname

      Tpool_ini(:,:,:,:,:,:) = 0.d0  !initialize all pools to zero (g-C/m^2)

#ifdef PFT_MODEL_YIBS
!YK - temp. values, modified from 8 GISS pfts below
!NK - later these arrays should be moved to yibs_pfts
        Cpool_fracs(1,:,1) = (/ !ever_ES_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(2,:,1) = (/ !ever_LS_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(3,:,1) = (/ !ever_ES_needle
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(4,:,1) = (/ !ever_LS_needle
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(5,:,1) = (/ !cold_ES_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(6,:,1) = (/ !cold_LS_broad
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(7,:,1) = (/ !drought_broad
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989 /)
        Cpool_fracs(8,:,1) = (/ !decid_needle
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(9,:,1) = (/ !shrub_cold
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(10,:,1) = (/ !shrub_arid
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989 /)
        Cpool_fracs(11,:,1) = (/ !c3grass
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(12,:,1) = (/ !c4grass
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(13,:,1) = (/ !c3grass_ann
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(14,:,1) = (/ !c3grass_arctic
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(15,:,1) = (/ !cropsc4
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(16,:,1) = (/ !cropstree
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
#else
!***  for now define for 8 GISS pfts, one-layer only -PK 1/23/08***
        Cpool_fracs(1,:,1) = (/ !tundra (for now=C3 grass)
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(2,:,1) = (/ !C3 grass (Vaira)
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
        Cpool_fracs(3,:,1) = (/ !shrub (for now=savanna)
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989 /)
        Cpool_fracs(4,:,1) = (/ !savanna (Tonzi)
     &       0.0026084,0.077080104,0.001512116,0.059312743,0.064831966,
     &       0.007522776,0.022795146,0.57388576,0.190450989 /)
        Cpool_fracs(5,:,1) = (/ !decid broadl (MMSF)
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(6,:,1) = (/ !evergr needl (for now=decid broadl)
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(7,:,1) = (/ !trop rainf (for now=decid broadl)
     &       0.005387981,0.062119495,0.004371574,0.050484497,0.280263607
     &       ,0.007488613,0.032152921,0.406417963,0.151313349 /)
        Cpool_fracs(8,:,1) = (/ !crops (for now=C3 grass)
     &       0.001891469,0.078919906,0.000456561,0.016762327,0.,
     &       0.009848682,0.014675189,0.692995043,0.184450822 /)
#endif

#ifdef RESTART_CPOOLS
      fname=trim(dirter)//pthsep//'SOILCARB_restart'
      call openunit(trim(fname),iunit,.true.,.true.)
      do k=1,N_PFT
        read(iunit) head, cpool1
        cpool_pft(:,:,k) = cpool1
      enddo
      call closeunit(iunit)
#endif

      !assign Tpool_ini values (pft-specific)
      do p=1,N_PFT
        do n=1,N_CASA_LAYERS 
          do nn=NLIVE+1,NPOOLS
#ifdef RESTART_CPOOLS
            Tpool_ini(p,CARBON,nn-NLIVE,n,I0:I1,J0:J1) =
     &           Cpool_fracs(p,nn-NLIVE,n)
     &           * cpool_pft(I0:I1,J0:j1,p)*1.d3
#else
            Tpool_ini(p,CARBON,nn-NLIVE,n,I0:I1,J0:J1) =
     &           Cpool_fracs(p,nn-NLIVE,n)
c     &           * soil_C_total(n,I0:I1,J0:J1)*1.d3  
     &           * lms%cpool(I0:I1,J0:J1,p)*1.d3   !!!!! restart cpool
#endif
          end do
        end do
      end do
ccc#endif

      end subroutine prescr_get_soilpools


      subroutine prescr_get_soil_C_total(IM,JM,I0,I1,J0,J1,
     &     soil_C_total,lm,lms)
!@sum (gC/m2) Read map of total soil carbon from file.
      use FILEMANAGER, only : openunit,closeunit
      integer,intent(in) :: IM,JM,I0,I1,J0,J1
      type(lm_exchange) , intent(in) :: lm
      type(lm_state) , intent(in) :: lms
      real*8,intent(out) ::
     &     soil_C_total(N_CASA_LAYERS,I0:I1,J0:J1)
      !---
      real*4 :: buf(N_CASA_LAYERS,144,90)
      character*80 :: title,fname
      integer :: iu_SOILCARB
      real*8 :: lons(144), lats(90)
      integer k,i,j, ii, jj

      !fname=trim(dirter)//pthsep//'SOILCARB_global'
      !call openunit(trim(fname),iu_SOILCARB,.true.,.true.)
      !read (iu_SOILCARB) title, buf !data in kg/m2 (converted to g/m2 below)
      !call closeunit(iu_SOILCARB)

      If (IM .ne. 360 .and. IM .ne. 270) stop "Wrong IM for soil carbon"

      do i = 1, 144
         lons(i) = -178.75d0 + dble(i-1)*2.5d0
      enddo
      do j = 1, 90
         lats(j) = -89.0d0 + dble(j-1)*2.0d0
      enddo

      do j = j0, j1
      do i = i0, i1
         do k=1,N_CASA_LAYERS
            soil_C_total(k,i,j)=lms%soil_c(i,j)
         enddo
      enddo
      enddo

      end subroutine prescr_get_soil_C_total


      subroutine read_soilcarbon_site(I0,I1,J0,J1,Tpool_ini)
!@sum (gC/m2 by soil layer) Read site values for soil carbon pools from file.
!External file should be named as below and should be organized as follows:
!(1) there should be 1 or 2 columns (corresponding to each soil bgc layer);
!(2) first non-header row should have total site-measured pool (in g/m2);
!(3) 9 subsequent rows correspond to modeled 9 soil pool fractions
      use FILEMANAGER, only : openunit,closeunit,nameunit

      integer,intent(in) :: I0,I1,J0,J1
      real*8,intent(out) :: Tpool_ini(N_PFT,PTRACE,NPOOLS-NLIVE
     &     ,N_CASA_LAYERS, I0:I1,J0:J1)!prescribed soil pools, g/m2
      !---Local------------------------
      !Site total measured soil C_org:
      integer :: iu_SOILCARB  !File ID
      real*8, dimension(N_CASA_LAYERS) :: total_Cpool  !g-C/m^2
      real*8, dimension(NPOOLS-NLIVE,N_CASA_LAYERS) :: Cpool_fracs_in !fraction

      !Variables for calculating soil carbon pools
      real*8, dimension(N_PFT,NPOOLS-NLIVE,N_CASA_LAYERS) :: Cpool_fracs !fraction
      integer :: nn,p,n
      character*40 fname

      fname=trim(dirter)//pthsep//'SOILCARB_site'
      call openunit(trim(fname),iu_SOILCARB,.false.,.true.)  !csv dataset
      read(iu_SOILCARB,*)  !skip optional header row(s)
      read(iu_SOILCARB,*) total_Cpool(:)
      do nn=1,NPOOLS-NLIVE
        read(iu_SOILCARB,*) Cpool_fracs_in(nn,:)
      end do

      do p=1,N_PFT      
       do n=1,N_CASA_LAYERS 
        do nn=NLIVE+1,NPOOLS
         Cpool_fracs(p,nn-NLIVE,n) = Cpool_fracs_in(nn-NLIVE,n)
         Tpool_ini(p,CARBON,nn-NLIVE,n,I0:I1,J0:J1) =
     &          Cpool_fracs(p,nn-NLIVE,n)*total_Cpool(n)
        end do
       end do
      end do
      end subroutine read_soilcarbon_site

!***************************************************************************
      ! This module is not used for GCM runs, but only for standalone runs.
      subroutine prescr_vegdata(lm,lms,jday, year, IM,JM,I0,I1,J0,J1,
     &     vegdata,albedodata,laidata,hdata,nmdata,popdata,dbhdata,
     &     craddata,cpooldata,rootprofdata,soil_color,soil_texture,
     &     Tpooldata, vegdata2, cropdata2, elev, do_clmpft,
     &     do_soilinit,do_phenology_activegrowth,do_init_geo
     &     ,do_read_from_files)
!@sum prescr_vegdata - File reading and Matthews prescribed calculations
!@+   of vegetation structure. Off-line (standalone) runs only.
!     This is a general driver routine that calls routines
!     in this module only to read or calculate vegetation and soil structure
!     and carbon pools; the module subroutines call routines in yibs_prescr_veg,
!     which do the explicit calculations and are not available to outside
!     drivers.  
!     This routine may be imitated by drivers used for off-line runs
!     or coupled runs to GCMs.
      use yibs_prescribed_drv_geo, only : init_yibsvegdata_geo
     &     ,init_yibs_laimax_geo
      implicit none
      integer,intent(in) :: jday, year
      integer,intent(in) :: IM,JM,I0,I1,J0,J1 !long/lat grid number range
      real*8,intent(in)  :: vegdata2(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(in)  :: cropdata2(I0:I1,J0:J1)
      real*8,intent(in)  :: elev(I0:I1,J0:J1)
      type(lm_exchange) , intent(in) :: lm
      type(lm_state) , intent(in) :: lms
      real*8,intent(out) :: vegdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: albedodata(N_BANDS,N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: laidata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: hdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: nmdata(N_COVERTYPES)
      real*8,intent(out) :: popdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: dbhdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: craddata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: cpooldata(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1)
      real*8,intent(out) :: rootprofdata(N_COVERTYPES,N_DEPTH)
      integer,intent(out) :: soil_color(N_COVERTYPES)
      real*8,intent(out) :: soil_texture(N_SOIL_TEXTURES,I0:I1,J0:J1)
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &     I0:I1,J0:J1):: Tpooldata !in g/m2 -PK
      logical,intent(in) :: do_clmpft
      logical,intent(in) :: do_soilinit
      logical,intent(in) :: do_phenology_activegrowth, do_init_geo
      logical,intent(in) :: do_read_from_files

      !-----Local------
      integer :: i,j,k, jeq, p
      integer hemi(I0:I1,J0:J1)
      REAL*8 :: soil_C_total(N_CASA_LAYERS,I0:I1,J0:J1)
      real*8 :: laimaxdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8 :: laimaxdata1(N_COVERTYPES,I0:I1,J0:J1)
      real*8 :: laidata1(N_PFT,I0:I1,J0:J1)
      real*8 :: hdata1(N_COVERTYPES,I0:I1,J0:J1)

      jeq = JM/2
      do j=J0,J1
        hemi(:,j) = 1

        if (j <= jeq) hemi(:,j) = -1
      enddo

      call init_vfraction(IM,JM,I0,I1,J0,J1,vegdata,vegdata2)   !veg fractions
         do j=J0,J1
         do i=I0,I1
         do k = 1,18
         if (vegdata(k,i,j) .gt. 1.) then
         print*,'here after init_vfraction: ',k,i,j,vegdata(k,i,j)
         endif
         enddo
         enddo
         enddo
      if (.not.do_clmpft)
     &  call prescr_update_vegcrops(year,IM,JM,I0,I1,J0,J1
     &  ,vegdata,cropdata2,elev)
         do j=J0,J1
         do i=I0,I1
         do k = 1,18
         if (vegdata(k,i,j) .gt. 1.) then
         print*,'apfter prescr_update_vegcrops: ',k,i,j,vegdata(k,i,j)
         endif
         enddo
         enddo
         enddo

      if ((.not.do_phenology_activegrowth).and.
     &     (.not.do_read_from_files)) then
         call prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laidata) !lai
         call prescr_get_hdata(I0,I1,J0,J1,hdata) !height
         call prescr_calc_canopy_geometry(I0,I1,J0,J1
     i        ,hdata
     o        ,dbhdata,popdata,craddata)
         call prescr_get_carbonplant(I0,I1,J0,J1,
     &           laidata,hdata,dbhdata,popdata,cpooldata)
!---------------------------------------------------
      else  !if do_phenology_activegrowth=true or do_read_from_files
         laimaxdata(:,:,:) = 0.d0
         hdata(:,:,:) = 0.d0
         laidata(:,:,:) = 0.d0
         laidata1(:,:,:) = 0.d0
         laimaxdata1(:,:,:) = 0.d0
         hdata1(:,:,:) = 0.d0
         print*,"here 111"
         call getlaix_m16(laimaxdata1,hdata1,I0,
     &       I1,J0,J1,do_clmpft,lm,lms)
         print*,"here 222"
         call getlai_m16(jday, laidata1,I0,I1,J0,J1, do_clmpft,lm,lms)
         print*,"here 333"
         call prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laidata) !lai
         call prescr_get_laimaxdata(I0,I1,J0,J1,laimaxdata)
         call prescr_get_hdata(I0,I1,J0,J1,hdata)
         print*,"here 444"
         do j=J0,J1
           do i=I0,I1
             do k=1,16
               if (laidata1(k,i,j) .gt. 10.d0) then
                  print*,"laidata1 error: ", laidata1(k,i,j)
               endif
               if (laidata1(k,i,j) .gt. 0.d0) then
                  laidata(k,i,j) = laidata1(k,i,j)
               endif
               if (hdata1(k,i,j) .gt. 0.d0) then
                  hdata(k,i,j) = hdata1(k,i,j)
               endif
               if (laimaxdata1(k,i,j) .gt. 0.d0) then
                  laimaxdata(k,i,j) = laimaxdata1(k,i,j)
               endif
             enddo
           enddo
         enddo
         !update diameter, population density, carbon plant &  crown rad
         if (.not.do_init_geo) then
          print *,'Calling prescr_get_yibs_plant in yibs_prescribed_drv'
          call prescr_get_yibs_plant(I0,I1,J0,J1, 
     i           vegdata, laidata,hdata,laimaxdata,
     o           dbhdata,popdata,craddata,cpooldata)
         do j=J0,J1
         do i=I0,I1
         do k = 1,18
         if (vegdata(k,i,j) .gt. 1.) then
         print*,'here prescr_get_yibs_plant: ',k,i,j,vegdata(k,i,j)
         endif
         enddo
         enddo
         enddo

         else
          print *, 'Initializing geographic veg data.'
          call init_yibs_laidata(IM,JM,I0,I1,J0,J1,laidata) !lai
          call init_yibsvegdata_geo( IM,JM,I0,I1,J0,J1
     i           ,laidata,hdata,laimaxdata
     o           ,popdata ,dbhdata,craddata,cpooldata) 
         endif
      endif

      call prescr_veg_albedodata(jday,hemi,I0,I1,J0,J1,albedodata)

      call prescr_get_pft_vars(nmdata,rootprofdata,soil_color)

      if ( do_read_from_files ) then
      do j = j0, j1
      do i = i0, i1
         do k=1,N_SOIL_TEXTURES
            soil_texture(k,i,j) = lms%stext(i,j,k)
            if(isnan(soil_texture(k,i,j))) soil_texture(k,i,j)=0
            if((soil_texture(k,i,j)).gt.100) soil_texture(k,i,j)=0
            if((soil_texture(k,i,j)).lt.-100) soil_texture(k,i,j)=0
         enddo
       enddo
       enddo
      endif
C        !&     call prescr_get_soiltexture(IM,JM,I0,I1,J0,J1,
C        !&     soil_texture,lms)
#ifdef SET_SOILCARBON_GLOBAL_TO_ZERO
      Tpooldata(:,:,:,:,:,:) = 0.d0
#else
      if ( do_soilinit ) then
#ifdef SOILCARB_SITE
         !Site soil carbon pools
           print *,"Getting site soil carbon"
           call read_soilcarbon_site(I0,I1,J0,J1,Tpooldata)
#else
         !Global soil carbon pools
          print *,'Reading global soil carbon data'
          call prescr_get_soil_C_total(IM,JM,I0,
     &      I1,J0,J1,soil_C_total,lm,lms)
          call prescr_get_soilpools(IM,JM,I0,I1,J0,J1,
     &                              soil_C_total,Tpooldata,lms)
#endif
      else
         Tpooldata(:,:,:,:,:,:) = 0.d0
      endif
#endif

      end subroutine prescr_vegdata


!***************************************************************************
      subroutine init_vfraction(im,jm,I0f,I1f,J0f,J1f,vfraction,vfrac)
!@sum Read in vegetation and soil cover from file (mosaicked fractions).
      use FILEMANAGER, only : openunit,closeunit,nameunit
      integer, intent(in) :: im,jm,I0f,I1f,J0f,J1f
      real*8, intent(in)  :: vfrac(N_COVERTYPES,I0f:I1f,J0f:J1f) 
      real*8, intent(out) :: vfraction(N_COVERTYPES,I0f:I1f,J0f:J1f) 
      !------Local---------------------
      !1    2    3    4    5    6    7    8    9   10   11    12
      !BSAND TNDRA GRASS SHRUB TREES DECID EVERG RAINF CROPS BDIRT ALGAE GRAC4
      character*80 :: title
      real*4 :: buf(im,jm)
      integer :: iu_VEG
      integer :: k,i,j
      real*8 :: s

      ! Make sure that unused fractions are set to 0
      vfraction(:,:,:) = 0.d0
      do i=I0f,I1f
         do j=J0f,J1f
            vfraction(:,i,j) = vfrac(:,i,j)
         enddo
      enddo

      !NK - This looks like it could be setting water to SAND??
      ! make sure that veg fractions are reasonable
      do j=J0f,J1f
         do i=I0f,I1f
          do k=1,N_COVERTYPES
            ! get rid of unreasonably small fractions
            if ( vfraction(k,i,j) < 1.d-4 ) vfraction(k,i,j) = 0.d0
          enddo
          s = sum( vfraction(:,i,j) )
          if ( s > .99d0 ) then !Keep if at least 90% specified, scale to 1.
             vfraction(:,i,j) = vfraction(:,i,j)/s
          else if ( s < .1d0 ) then
             print *, "missing veg data at ",i,j,"assume bare soil",s
             vfraction(:,i,j) = 0.d0
             vfraction(COVER_SAND,i,j) = 1.d0
          else
c             stop "Incorrect data in VEG file"
          endif
        enddo
      enddo

      print *,'End of init_vfraction'
      end subroutine init_vfraction

!**************************************************************************
      subroutine prescr_get_cropdata(year,IM,JM,I0,I1,J0,J1,cropdata)
!@sum Read in crop cover fractions from file.
!@+   Calculates crop fraction for given year.
      use FILEMANAGER, only : openunit,closeunit,nameunit
      integer, intent(in) :: year
      integer, intent(in) :: IM, JM, I0, I1, J0, J1
      real*8, intent(out) :: cropdata(I0:I1,J0:J1)
      integer i
      !----------
      integer :: iu_CROPS
      integer :: year1, year2
      real*4 crop4(im,jm)
      real*8 wt, crop1(I0:I1,J0:J1), crop2(I0:I1,J0:J1)
      character*80 title,fname

      !* Calculate fraction for given gcmtime:  interpolate between years*/

      year1 = -32768 ; crop1(:,:) = 0.d0
      year2 = -32767 ; crop2(:,:) = 0.d0
      wt = 1.d0

      fname=trim(dirter)//pthsep//'CROPS'
      call openunit(trim(fname),iu_CROPS,.true.,.true.)
      do while( year2 < year )
        year1 = year2
        crop1(:,:) = crop2(:,:)
        read (iu_CROPS,end=10) title , crop4
        read(title,*) year2 !Read year integer out of character array title
        crop2(I0:I1,J0:J1) = crop4(I0:I1,J0:J1)
      enddo
      wt = (year-year1)/(real(year2-year1,kind=8))
 10   continue
      call closeunit(iu_CROPS)

      cropdata(:,:) = max(0.d0, crop1(:,:)
     &     + wt * (crop2(:,:) - crop1(:,:)))  !Set min to zero, since no land mask yet -nyk 1/22/08

      end subroutine prescr_get_cropdata

!**************************************************************************

      subroutine prescr_update_vegcrops(year,IM,JM,I0,I1,J0,J1,
     &     vegdata, vcrop, elev)
!@sum Rescales natural vegetation cover fractions given new crop cover.
      integer,intent(in) :: year
      integer, intent(in) :: IM,JM,I0,I1,J0,J1
      real*8, intent(in) :: vcrop(I0:I1,J0:J1)
      real*8, intent(in) :: elev(I0:I1,J0:J1)
      real*8, intent(inout) :: vegdata(N_COVERTYPES,I0:I1,J0:J1)
      !--------
      real*8,ALLOCATABLE,dimension(:,:) :: cropdata !grid array
      integer :: i,j,k
      real*8 crops_old

      ALLOCATE(cropdata(I0:I1,J0:J1))

      !* Loop *!
      print *,"Getting crop cover."
cxyue  call prescr_get_cropdata(year,IM,JM,I0,I1,J0,J1,cropdata) !crop fraction
      cropdata(I0:I1,J0:J1) = vcrop(I0:I1,J0:J1)

      !* If cropdata was prepared somewhere else, then cover is as simple as
      !* modifying the vegetation fractions.  Need to update cohort and
      !* patch summary variables.
      do j=J0,J1
        do i=I0,I1
          if (elev(i,j) .lt. 0.) then
            cropdata(i,j) = 0.d0
          endif
          if ( cropdata(i,j) == 1.d0 ) then
            vegdata(:,i,j) = 0.d0
            vegdata(CROPS+COVEROFFSET,i,j) = 1.d0
          else
            crops_old = vegdata(CROPS+COVEROFFSET,i,j)
            if ( crops_old .gt. 0.99d0 ) then
              crops_old=0.99d0
             ! call stop_model("incompatible crops: old=1, new<1",255)
            endif
         do k = 1,18
         if (vegdata(k,i,j) .gt. 1.) then
         print*,'before change crop: ',k,i,j,vegdata(k,i,j)
         endif
         enddo
            vegdata(:,i,j) = vegdata(:,i,j)
     $           * (1.d0-cropdata(i,j))/(1.d0-crops_old)
            vegdata(CROPS+COVEROFFSET,i,j) = cropdata(i,j)
         do k = 1,18
         if (vegdata(k,i,j) .gt. 1.) then
         print*,'after change crops: ',k,i,j,vegdata(k,i,j)
         endif
         enddo
          endif
        end do
        !write(*,*) vegdata(CROPS+COVEROFFSET,:,j)
      end do

      DEALLOCATE(cropdata)
 
      end subroutine prescr_update_vegcrops

!**************************************************************************

      subroutine prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laidata)
!@sum Returns prescr GCM leaf area index for entire grid and given jday.
!@+   Based on seasonal sinusoidally varying LAI by veg type from R&A 1997.
      use yibs_const,only : N_COVERTYPES
      integer,intent(in) :: jday
      integer, intent(in) :: I0,I1,J0,J1
      real*8 :: laidata(N_COVERTYPES,I0:I1,J0:J1) 
      integer :: hemi(I0:I1,J0:J1) !@var hemi =1 in N. hemisphere, =-1 in South
      !----------
      integer :: n !@var cover type
      integer i,j,jeq

      !jeq = JM/2

      do j=J0,J1
        !hemi = 1
        !if (j <= jeq) hemi = -1
        do i=I0,I1
          do n=1,N_COVERTYPES
            laidata(n,i,j) = prescr_calc_lai(n,jday,hemi(i,j))
          enddo
        enddo
      enddo

      !* Return lai for each vegetation type.
      end subroutine prescr_get_laidata

!**************************************************************************

      subroutine init_yibs_laidata(IM,JM,I0,I1,J0,J1,laidata)
!@sum Read the initial LAI for prognostic vegetation 
!@+   (i.e., prognostic phenology/growth with YIBS PFTs)

      use FILEMANAGER, only : openunit,closeunit,nameunit !for VEG_PROGNOSTIC
      use yibs_const,only : N_COVERTYPES
      integer, intent(in) :: IM,JM,I0,I1,J0,J1
      real*8 :: laidata(N_COVERTYPES,I0:I1,J0:J1) 
      !----------
      character*80 :: title,fname
      real*4 :: buf(IM,JM)
      integer :: iu_LAI
      integer :: k !@var cover type

      fname=trim(dirter)//pthsep//'LAIinit'
      call openunit(trim(fname),iu_LAI,.true.,.true.)

      do k=1,N_COVERTYPES
        read(iu_LAI) title , buf
        laidata(k,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
      end do
      !print *,"laidata", laidata(:,I0:I1,J0:J1) !#DEBUG

      call closeunit(iu_LAI)

      !* Return lai for each vegetation type.
      end subroutine init_yibs_laidata

!**************************************************************************
      subroutine prescr_veg_albedodata(jday,hemi,I0,I1,J0,J1,albedodata)
!@sum Assign seasonal prescribed albedo by cover type (Matthews, 1983).
      integer,intent(in) :: jday
      integer, intent(in) :: I0,I1,J0,J1
      real*8 :: albedodata(N_BANDS,N_COVERTYPES,I0:I1,J0:J1)
      integer :: hemi(I0:I1,J0:J1)
      !----------
      !integer :: pft
      integer :: ncov
      integer i,j,jeq
      
      !jeq = JM/2

      do j=J0,J1
        !hemi = 1
        !if (j <= jeq) hemi = -1
        do i=I0,I1
          do ncov = 1, N_COVERTYPES
            call prescr_veg_albedo(hemi(i,j),ncov,jday,
     &           albedodata(:,ncov,i,j))
          end do
        enddo
      enddo

      end subroutine prescr_veg_albedodata

!**************************************************************************

      subroutine init_yibs_hdata(IM,JM,I0,I1,J0,J1,hdata3d)
!@sum Read the initial height for prognostic vegetation. 
!@+   (i.e., prognostic phenology/growth with YIBS PFTs)
      use FILEMANAGER, only : openunit,closeunit,nameunit !for VEG_PROGNOSTIC
      use yibs_const,only : N_COVERTYPES
      integer, intent(in) :: IM,JM,I0,I1,J0,J1
      real*8 :: hdata3d(N_COVERTYPES,I0:I1,J0:J1) 
      !----------
      character*80 :: title,fname
      real*4 :: buf(IM,JM)
      integer :: iu_HITE
      integer :: k !@var cover type
     
      print *, 'Reading height data'
      fname=trim(dirter)//pthsep//'HITE'
      call openunit(trim(fname),iu_HITE,.true.,.true.)

      do k=1,N_COVERTYPES
        read(iu_HITE) title , buf
        hdata3d(k,I0:I1,J0:J1) = buf(I0:I1,J0:J1)
      end do
      !print *,"hdata", hdata(:,I0:I1,J0:J1) !#DEBUG

      call closeunit(iu_HITE)

      end subroutine init_yibs_hdata

!**************************************************************************

      subroutine prescr_get_laimaxdata(I0,I1,J0,J1
     &     ,laimaxdata)!,laimindata)
!@sum Make arrays of prescribed LAI: alamax, alamin from Matthews (1983).
      !* This should access alamax and alamin via yibs_prescr_veg, but
      !* so many levels is annoying for a trivial assignment. ...NK
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*8,intent(out) :: laimaxdata(N_COVERTYPES,I0:I1,J0:J1) 
      !real*8,intent(out) :: laimindata(N_COVERTYPES,I0:I1,J0:J1) 
      !--- Local ----
      integer :: i,j
      do i=I0,I1
         do j=J0,J1
            laimaxdata(:,i,j) = alamax(:)
            !laimindata(:,i,j) = alamin(:)
         enddo
      enddo
      end subroutine prescr_get_laimaxdata

!**************************************************************************
      subroutine prescr_get_hdata(I0,I1,J0,J1,hdata) !height
!@sum Return array of prescribed canopy heights (m) by vegetation type
!@+   (Matthews, 1983).
      use yibs_const,only : N_COVERTYPES
      use yibs_prescr_veg, only : prescr_calc_hdata
      integer, intent(in) :: I0,I1,J0,J1
      real*8,intent(out) :: hdata(N_COVERTYPES,I0:I1,J0:J1)
      !---Local----
      integer :: i,j

      !* Read file of height data - TBA

      !* Calculate prescribed model tree heights
      do j=J0,J1
         do i=I0,I1
            call prescr_calc_hdata(hdata(:,i,j))
         enddo
      enddo
      
      end subroutine prescr_get_hdata

!**************************************************************************
      subroutine prescr_calc_canopy_geometry(I0,I1,J0,J1
     i     ,hdata,dbhdata,popdata,craddata)
!@sum Calculate canopy structure variables that would not be provided 
!@+   by data files but are determined through allometric relations to
!@+   height.  
!@+   Subject to change according to data.
      use yibs_prescr_veg, only: prescr_calc_woodydiameter
     &     ,prescr_get_pop,prescr_get_crownrad
      integer, intent(in) :: I0,I1,J0,J1
      real*8,intent(in) :: hdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: dbhdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: popdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: craddata(N_COVERTYPES,I0:I1,J0:J1)
      !---Local---
      integer :: i, j
      real*8 :: laimaxdata(N_COVERTYPES,I0:I1,J0:J1)
      
      call prescr_get_laimaxdata(I0,I1,J0,J1,laimaxdata)

      do j=J0,J1
         do i=I0,I1
            call prescr_calc_woodydiameter(
     &           hdata(:,i,j), dbhdata(:,i,j))
            call prescr_get_pop(
     &           dbhdata(:,i,j), laimaxdata(:,i,j), popdata(:,i,j))
            !## Need to re-do prescr_get_crownrad fot non-closed canopy.
            call prescr_get_crownrad(popdata(:,i,j), craddata(:,i,j))
         enddo
      enddo
      end subroutine prescr_calc_canopy_geometry
      
!**************************************************************************
!**************************************************************************
      subroutine prescr_get_carbonplant(I0,I1,J0,J1,
     &     laidata, hdata, dbhdata,popdata
     &     , cpooldata)
!@sum Calculate per plant carbon pools (g-C/plant) from geometry and 
!@+   alamax for Matthews (1983) vegetation structure arrays.
      use allometryfn, only : allom_plant_cpools,init_Clab
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*8,intent(in) :: laidata(N_COVERTYPES,I0:I1,J0:J1) 
      !real*8,intent(in) :: laimaxdata(N_COVERTYPES,I0:I1,J0:J1) 
      real*8,intent(in) :: hdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(in) :: dbhdata(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(in) :: popdata(N_COVERTYPES,I0:I1,J0:J1) 
      real*8,intent(out) :: cpooldata(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1)
      !Array: 1-foliage, 2-sapwood, 3-hardwood, 4-labile,5-fine root, 6-coarse root
      !----Local----
      integer :: p,pft !@var pft vegetation type
      integer i,j,jeq

      cpooldata(:,:,:,:) = 0.d0  !Zero initialize
      do j=J0,J1
        do i=I0,I1
          do pft=1,N_PFT
            p = pft + COVEROFFSET
            call allom_plant_cpools(pft,laidata(p,i,j)
     &           , hdata(p,i,j)
     &           ,dbhdata(p,i,j), popdata(p,i,j),cpooldata(p,:,i,j))
            call init_Clab(pft,dbhdata(p,i,j)
     &           ,hdata(p,i,j) ,cpooldata(p,LABILE,i,j))
          enddo
        enddo
      enddo
      
      end subroutine prescr_get_carbonplant
!**************************************************************************
      subroutine prescr_get_yibs_plant(I0,I1,J0,J1
     i     ,vegdata, laidata, hdata3d, laimax
     o     ,dbhdata3d, popdata3d, craddata3d,cpooldata)
!@sum Calculate woody diameter, population density, crown geometry,
!@+   carbon pools for the YIBS prognostic vegetation and geographic variation.
!@+   Density is obtained from max LAI and geometric allometry.
      use allometryfn, only: allom_plant_cpools, height2dbh,nplant
     &     ,Crown_rad_max_from_density,Crown_rad_allom, init_Clab
      implicit none
      integer,intent(in) :: I0,I1,J0,J1
      real*8,intent(in) :: vegdata(N_COVERTYPES,I0:I1,J0:J1) 
      real*8,intent(in) :: laidata(N_COVERTYPES,I0:I1,J0:J1) 
      real*8,intent(in) :: hdata3d(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(in) :: laimax(N_COVERTYPES,I0:I1,J0:J1) 
      real*8,intent(out) :: dbhdata3d(N_COVERTYPES,I0:I1,J0:J1)
      real*8,intent(out) :: popdata3d(N_COVERTYPES,I0:I1,J0:J1) 
      real*8,intent(out) :: craddata3d(N_COVERTYPES,I0:I1,J0:J1) 
      real*8,intent(out) :: cpooldata(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1)
      !Array: 1-foliage, 2-sapwood, 3-hardwood, 4-labile,5-fine root, 6-coarse root
      !----Local----
      integer :: pft !@var pft vegetation type
      integer :: i,j, n

!Zero initialize.
      dbhdata3d(:,:,:) = 0.0 
      popdata3d(:,:,:) = 0.0 
      craddata3d(:,:,:) = 0.0
      cpooldata(:,:,:,:) = 0.d0 
      do j=J0,J1
        do i=I0,I1
          do pft = 1,N_PFT
            n = pft + COVEROFFSET
            if (pfpar(pft)%woody)  dbhdata3d(n,i,j) 
     &                             =height2dbh(pft,hdata3d(n,i,j))  
            !update population denisty
            popdata3d(n,i,j) = nplant(pft,dbhdata3d(n,i,j), 
     &                       hdata3d(n,i,j), laimax(n,i,j))
            !update crown rad
            craddata3d(n,i,j) = min(
     &           Crown_rad_max_from_density(popdata3d(n,i,j))
     &           ,Crown_rad_allom(pft,hdata3d(n,i,j)))
            !update cabon pools
            call allom_plant_cpools(pft,laidata(n,i,j)
     &           ,hdata3d(n,i,j),dbhdata3d(n,i,j), popdata3d(n,i,j)
     &           ,cpooldata(n,:,i,j))

            call init_Clab(pft,dbhdata3d(n,i,j),
     &           hdata3d(pft,i,j),cpooldata(n,LABILE,i,j))
          enddo
        enddo
      enddo
      
      end subroutine prescr_get_yibs_plant

!*************************************************************************
      subroutine prescr_get_pft_vars(nmdata,rootprofdata,soil_color)
!@sum Routine to assign pft-dependent misc. arrays that are not gridded.
      use yibs_prescr_veg, only : prescr_calc_initnm,
     &     prescr_calc_rootprof_all,prescr_calc_soilcolor

      real*8,intent(out) :: nmdata(N_COVERTYPES)
      real*8,intent(out) :: rootprofdata(N_COVERTYPES,N_DEPTH)
      integer,intent(out) :: soil_color(N_COVERTYPES)

      call prescr_calc_initnm(nmdata) !nm ! mean canopy nitrogen
      call prescr_calc_rootprof_all(rootprofdata)
      call prescr_calc_soilcolor(soil_color)
      
      end subroutine prescr_get_pft_vars
!*************************************************************************
      subroutine prescr_get_soiltexture(IM,JM,I0,I1,J0,J1,
     &     soil_texture,lms)
!@sum Return arrays of GISS soil color and texture.
      use FILEMANAGER, only : openunit,closeunit
      integer, intent(in) :: IM,JM,I0,I1,J0,J1
      type(lm_state) , intent(in) :: lms
      real*8, intent(out) :: soil_texture(N_SOIL_TEXTURES,I0:I1,J0:J1)
      !------
      real*8 :: buf(144,90,N_SOIL_TEXTURES)
      real*8 :: long(IM), latg(JM), lons(144), lats(90)
      real*8 :: dlon, dlat, lon0
      integer :: iu_SOIL
      integer k,i,j, ii, jj
      character*50 :: fname

      fname=trim(dirter)//pthsep//'soiltextures'
      print *, 'here reading soil textures',trim(fname) !#DEBUG
      call openunit(trim(fname),iu_SOIL,.true.,.true.)
      print *,IM,JM,N_COVERTYPES !#DEBUG
      read(iu_SOIL) buf
      call closeunit(iu_SOIL)
      print *,"soil fractions:",buf(I0,J0,:)!#DEBUG
      If (IM .ne. 360 .and. IM .ne. 270) stop "Wrong IM for soil carbon"
      If (IM .eq. 360) then
        do i = 1, IM
           long(i) = -180.0d0 + dble(i-1)*1.d0 + 1.d0/2.d0
        enddo
      else
        do i = 1, IM
           long(i) = -180.0d0 + dble(i-1)*4.d0/3.d0 + 1.d0/3.d0
        enddo
      Endif

      do j = 1, JM
         latg(j) = -90.0d0 + dble(j-1)*1.0d0
      enddo
      do i = 1, 144
         lons(i) = -178.75d0 + dble(i-1)*2.5d0
      enddo
      do j = 1, 90
         lats(j) = -89.0d0 + dble(j-1)*2.0d0
      enddo

      dlon = lons(10) - lons(9)
      dlat = lats(10) - lats(9)
      do j = j0, j1
      do i = i0, i1
         lon0 = long(i)
         If (lon0 .gt. 180.0d0) lon0 = lon0 - 360.0d0
         ii  = nint((lon0-lons(1))/dlon)+1
         jj  = nint((latg(j)-lats(1))/dlat)+1
         ii  = max( 1, min( ii, 144) )
         jj  = max( 1, min( jj, 90) )
         do k=1,N_SOIL_TEXTURES
            soil_texture(k,i,j) = buf(ii,jj,k)
         enddo
      enddo
      enddo

      end subroutine prescr_get_soiltexture
!*************************************************************************
#ifdef MIXED_CANOPY
      subroutine yibs_struct_get_phys(lms,IM,JM,I0, I1, J0, J1,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini,
     &     soil_texture,soil_C_total,Tpooldata,
     &     do_soilinit,do_read_from_files)
!@sum Initialize canopy and soil physics for mixed canopies.
      integer,intent(in) :: IM,JM,I0,I1,J0,J1
      type(lm_state) , intent(in) :: lms
      real*8, dimension(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini
      real*8, dimension(N_CASA_LAYERS,I0:I1,J0:J1) :: soil_C_total
!      integer, dimension(N_COVERTYPES) :: soil_color !soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &                  I0:I1,J0:J1):: Tpooldata  !g/m2
      logical,intent(in) :: do_soilinit
      logical,intent(in) :: do_read_from_files
      !------
      real*8 :: soil_C_total(N_CASA_LAYERS,I0:I1,J0:J1)
      call init_canopy_physical(I0, I1, J0, J1,
     &     Ci_ini, CNC_ini, Tcan_ini, Qf_ini)
      
      if ( do_read_from_files )
     &     call prescr_get_soiltexture(IM,JM,I0,I1,J0,J1,
     &     soil_texture)
      !soil_color !Don't need to do here. Gets read in structure file.

#ifdef SET_SOILCARBON_GLOBAL_TO_ZERO
      Tpooldata(:,:,:,:,:,:) = 0.d0
#else
      if ( do_soilinit ) then
        call prescr_get_soil_C_total(IM,JM,I0,I1,J0,J1,soil_C_total,lm,lms)
        call prescr_get_soilpools(IM,JM,I0,I1,J0,J1,
     &                            soil_C_total,Tpooldata,lms)
      else
        Tpooldata = 0.d0
      endif
#endif
      end subroutine yibs_struct_get_phys

#endif
!#MIXED_CANOPY
!*************************************************************************
      end module yibs_prescribed_drv

