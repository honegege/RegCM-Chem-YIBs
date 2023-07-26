      module run_yibs

      use yibs_mod
      !use yibs_prescrveg
      use filemanager
      use yibs_const 
      use vegfor_com, only: mvegid
      use geos_com
      use yibs_prescribed_drv, only: prescr_get_laidata
      use flux_com
      use mod_domain
      use mod_date
      use mod_runparams
      use mod_regcm_types
      use mod_dynparam
      
      
      implicit none
      public
      integer IM, JM, I0, I1, J0, J1, i0f, i1f, j0f, j1f
      integer jday, year, jday2, year2
      integer :: lyear , lmonth , lday , lhour      
      real*8 yibs_dt !seconds
      real*8  lon1, lat1, lon2, lat2
      real*8  lon, lat, c3_alpha
      logical :: force_VEG
      logical :: do_soilinit,do_soilresp
      logical :: do_phenology_activegrowth, do_structuralgrowth
      logical :: do_frost_hardiness
      logical :: do_patchdynamics, do_init_geo, do_spinup
      logical :: do_clmpft, do_islscp, c3c4_mix, pft_cover
      integer :: skip !#HACK to skip records at beginning of forcing file
      integer :: id_veg

      !---Local----
      integer :: jdaycount  !Only needed for prognostic phenology
      integer :: nh, nd, nm, ny, ny_old, totdays, nm_old
      integer, parameter :: mdays(12) = (/
     &    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      
!      real*8 :: max_time !In seconds
      real*8,parameter :: max_days = 365 !In days
      real*8, parameter :: save_interval=86400.d0*10d0 !save every 10 days

      real*8 :: long(360), latg(181), dlon, dlat
! GCM forcings
      real*8, pointer, dimension(:,:) :: o3s
! WFDEI forcings
      real*8, pointer, dimension(:,:) :: tam
      real*8, pointer, dimension(:,:) :: qfm
      real*8, pointer, dimension(:,:) :: pmm
      real*8, pointer, dimension(:,:) :: wdm
      real*8, pointer, dimension(:,:) :: prm
      real*8, pointer, dimension(:,:,:) :: stm
      real*8, pointer, dimension(:,:,:) :: smm

      real*8 :: tdy, cosz1, cosz2, cosz
      real*8 :: rtmp, co2p 
      real*8, pointer, dimension(:,:) ::
     &     lat2d, lon2d                !2-d latitude - for phenology
     &     ,TairC               !Air temperature (Celsius) !KIM - for phenology
     &     ,TcanopyC            !Canopy temperature (Celsius)
     &     ,Qf                  !Foliage surface specif humidity (kg vapor/ kg air)
     &     ,P_mbar              !Atmospheric pressure (mb)
     &     ,Ca                  !@Atmos CO2 conc at surface height (mol/m3).
     &     ,Ch                  !Ground to surface heat transfer coefficient 
     &     ,U                   !Surface layer wind speed (m s-1)
     &     ,IPARdif             !Incident diffuse PAR (vis) (W m-2)
     &     ,IPARdir             !Incident direct PAR (vis) (W m-2)
     &     ,CosZen              !cos of solar zenith angle
     &     ,O3
     &     ,O3c
     &     ,Cosz2d
      real*8, pointer, dimension(:,:,:) ::  !Changed to GCM layering -NK
     &      Soiltemp           !soil temp 
     &     ,Soilmoist          !soil volum moist 
     &     ,Soilmp             !Soil matric potential
     &     ,fice               !Fraction of soil water that is ice.
      real*8, dimension(:,:,:), pointer ::
     &     LAI                  ! prescribed LAI (for all pft's in the cell)
     &     ,height              ! prescribed height (for all pft's in the cell)
      

      !Coupling and diagnostic variables specific to YIBS (output)
      real*8, pointer, dimension(:,:,:) :: laidata
      real*8, pointer, dimension(:,:,:) :: laitmp1
      real*8, pointer, dimension(:,:,:)        :: laitmp
      real*8, pointer, dimension(:,:,:)        :: laidatam
      real*8, pointer, dimension(:,:,:) :: fv0
      real*8, pointer, dimension(:,:,:) :: lai_pfts
      real*8, pointer, dimension(:,:,:) :: phen_pfts
      real*8, pointer, dimension(:,:,:) :: ht_pfts
      real*8, pointer, dimension(:,:,:) :: gpp_pfts
      real*8, pointer, dimension(:,:,:) :: ipp_pfts
      real*8, pointer, dimension(:,:,:) :: mtp_pfts
      real*8, pointer, dimension(:,:,:) :: Clive_pfts
      real*8, pointer, dimension(:,:,:) :: Cdead_pfts
      real*8, pointer, dimension(:,:) ::
     &     GCANOPY              !Canopy conductance of water vapor (m s-1). 
     &     ,Ci                  !Internal foliage CO2 concentration (mol/m3)
     &     ,NPP                 !Net primary productivity (kg[C]/m2/s).
     &     ,GPP                 !Gross primary productivity (kg[C]/m2/s).
     &     ,GPPD                !Gross primary productivity (kg[C]/m2/s).
     &     ,IPP                 !Isoprene Emission (kg[C]/m2/s).
     &     ,MTP                 !Monoterpene Emission (kg[C]/m2/s).
     &     ,FO3                 !ozone flux to stomata (nmol O3 m-2 s-1).
     &     ,dFO3                !excess ozone flux to stomata (nmol O3 m-2 s-1).
     &     ,TRANS_SW            !Transmittance of shortwave through canopy to soil 
     &     ,z0, CO2flux         !Variables from YIBS to GCM
     &     ,R_auto              !Variables from YIBS to GCM
     &     ,lai0                !Variables from YIBS to GCM
     &     ,vegfrac
     &     ,lai_p               !Prognostic LAI
     &     ,ht_p                !Prognostic HT
     &     ,phenf               !Prognostic phenology
     &     ,landf               !Land fraction
     &     ,resp_r              !Root respiration
     &     ,resp_l              !Leaf respiration
     &     ,resp_w              !Wood respiration
     &     ,carb_tot            !Total Land carbon storage
     &     ,carb_soil           !Total Soil carbon storage
                     !&     ,CO2c  
      real*8, pointer, dimension(:,:,:) ::  !Changed to GCM layering -NK
     &      soilt              !soil temp 
     &     ,soilm              !soil volum moist 
     &     ,soilp              !Soil matric potential
     &     ,soili              !Fraction of soil water that is ice.
      real*8, pointer, dimension(:,:,:) ::
     &      laipfts, htpfts, gpppfts, ipppfts, mtppfts, 
     &      fvpfts, phenpfts, clivepfts, cdeadpfts

#ifdef OUTPUT_DAILY
      real*8, pointer, dimension(:,:) :: agcanopy, aci, anpp, agpp, 
     &    agppd, aco2flux, ar_auto, acoszen, ao3c, alai0, avegfrac,
     &    afo3, adfo3, alai_p, aht_p, aphenf, aipp, amtp, 
     &    aresp_r, aresp_l, aresp_w, acarb_tot, acarb_soil
      real*8, pointer, dimension(:,:,:) ::
     &    alaipfts, ahtpfts, agpppfts, afvpfts, aphenpfts,
     &    aipppfts, amtppfts, aclivepfts, acdeadpfts
#ifdef OUTPUT_FORCING
      real*8, pointer, dimension(:,:) :: atairc, atcanopyc, aqf,
     &    ap_mbar, aca, ach, au, aipardir, aipardif
      real*8, pointer, dimension(:,:,:) ::  
     &    asoilt, asoilm, asoilp, asoili
#endif
      real*8 nhrs
#endif

      real*8, pointer, dimension(:,:,:) ::
     &     betadl
      real*8, pointer, dimension(:,:,:) ::
     &     albedo               !Variables from YIBS to GCM

      !output for debug by xiaoxie
      real*8, pointer, dimension(:,:) :: in_lat !cohort
      real*8, pointer, dimension(:,:,:) :: in_veg !cohort
      real*8, pointer, dimension(:,:,:) :: in_pop  !cohort
      real*8, pointer, dimension(:,:,:) :: in_h
      real*8, pointer, dimension(:,:,:) :: in_dbh
      real*8, pointer, dimension(:,:,:) :: in_crad !cohort
      real*8, pointer, dimension(:,:,:) :: in_cpool !cohort
      real*8, pointer, dimension(:,:,:) :: in_albedo !cohort
      real*8, pointer, dimension(:,:,:) :: in_texture !cohort
      real*8, pointer, dimension(:,:,:) :: in_tpool  !cohort
      real*8, pointer, dimension(:,:) :: in_plant
      real*8, pointer, dimension(:,:) :: in_harvest      
      real*8, pointer, dimension(:,:,:) :: in_lai
      !---------------------------------------------------------------------

      
      type(ycelltype_public), pointer, dimension(:,:) :: cells
      !integer iu_forcings
      real*8 yibs_time, time_since_last_save
      integer, pointer, dimension(:,:):: hemi !hemisphere flags = 1 for N., =-1 for S.
      logical :: update_day    !For prescribed phenology, litter
      real*8, pointer, dimension(:,:):: fw ! fraction of wet canopy
      real*8 fv_tot      
      real(rk8) :: extime, timeend
      type(rcm_time_interval) :: tdif
  
      contains
      
      subroutine yibs_initialize(lm,lms)
      use mod_mppparam
      use mod_memutil
      use mod_outvars
      implicit none
      type(lm_exchange) , intent(in) :: lm
      type(lm_state) , intent(in) :: lms
      integer :: i,j,k,iflag      
      !* Default configuration
      yibs_dt = dtsrf              !second. Default value may be over-ridden by yibs_input
      force_VEG = .false.
      do_soilinit = .true.
      do_soilresp = .true.
      do_phenology_activegrowth = .false.
      do_structuralgrowth = .false.
      do_frost_hardiness = .true.
      do_patchdynamics = .false.
      do_init_geo = .false.
      do_spinup = .false.
      do_islscp = .false.
      do_clmpft = .true.
      c3c4_mix  = .true.
      pft_cover = .false.
      skip = 0
      id_veg = 4
      c3_alpha = 0.5 

      !* Set world bounds (should correspond to format of input files)
      IM = 360; JM = 181

      !* Default, entire grid.
      i0 = 1; i1 = IM
      j0 = 1; j1 = JM

      !* dims to default forcings file
      i0f = 1; i1f = IM
      j0f = 1; j1f = JM

      !* Default date start for GISS GCM 10-day forcings
      jday = 152                !June 1
      year = 1980
      jday2 = jday + 10
      year2 = -1 !Initialize

      write(*,*) 'Got here before read_input_parameters'
      call split_idate(idate1,lyear,lmonth,lday,lhour)      
      year = lyear
      jday = lday + sum(mdays(1:lmonth-1))

      if (year2.eq.-1) year2=year  !If year2 was not specified, default same as year.
      
      do i = 1, 360
         long(i) = -180.0d0 + dble(i-1)*1.d0 + 1.d0/2.d0
      enddo
      do j = 1, 181
         latg(j) = -90.0d0 + dble(j-1)*1.0d0
      enddo

      dlon = long(10) - long(9)
      dlat = latg(10) - latg(9)      

      I0   = jci1
      J0   = ici1
      I1   = jci2
      J1   = ici2

      !if ( myid == italk ) then
      print*,"jce1,ice1,jce2,ice2:",jce1,ice1,jce2,ice2,
     & "jci1,ici1,jci2,ici2:",jci1,ici1,jci2,ici2
      !endif

      I0f  = I0
      I1f  = I1
      J0f  = J0
      J1f  = J1

      !* Read elevation
      call get_elev(IM, JM, I0, I1, J0, J1, long, latg, lms)

      !* Read in [CO2]
      call read_co2
      
! GCM forcings call getmem2d(ncid%i4buf,jout1,jout2,iout1,iout2,'clm_createfile')
      call getmem2d(o3s,I0,I1,J0,J1,'yibs')
! WFDEI forcings
      call getmem2d(tam,I0,I1,J0,J1,'yibs') 
      call getmem2d(qfm,I0,I1,J0,J1,'yibs') 
      call getmem2d(pmm,I0,I1,J0,J1,'yibs') 
      call getmem2d(wdm,I0,I1,J0,J1,'yibs') 
      call getmem2d(prm,I0,I1,J0,J1,'yibs') 
      call getmem3d(stm,1,N_DEPTH,I0,I1,J0,J1,'yibs') 
      call getmem3d(smm,1,N_DEPTH,I0,I1,J0,J1,'yibs') 

      call getmem2d(lat2d,I0,I1,J0,J1,'yibs')
      call getmem2d(lon2d,I0,I1,J0,J1,'yibs')
      call getmem2d(TairC,I0,I1,J0,J1,'yibs')
      call getmem2d(TcanopyC,I0,I1,J0,J1,'yibs')
      call getmem2d(Qf,I0,I1,J0,J1,'yibs')
      call getmem2d(P_mbar,I0,I1,J0,J1,'yibs')
      call getmem2d(Ca,I0,I1,J0,J1,'yibs')
      call getmem2d(Ch,I0,I1,J0,J1,'yibs')
      call getmem2d(U,I0,I1,J0,J1,'yibs')
      call getmem2d(IPARdif,I0,I1,J0,J1,'yibs')
      call getmem2d(IPARdir,I0,I1,J0,J1,'yibs')
      call getmem2d(CosZen,I0,I1,J0,J1,'yibs')
      call getmem2d(O3,I0,I1,J0,J1,'yibs')
      call getmem2d(O3c,I0,I1,J0,J1,'yibs')
      call getmem2d(Cosz2d,I0,I1,J0,J1,'yibs')
      
      ! output for debug by xiaoxie
      call getmem2d(in_lat,I0,I1,J0,J1,'yibs')
      call getmem2d(in_plant,I0,I1,J0,J1,'yibs')
      call getmem2d(in_harvest,I0,I1,J0,J1,'yibs')
      call getmem3d(in_veg,1,N_COVERTYPES,I0,I1,J0,J1,'yibs')
      call getmem3d(in_pop,1,N_COVERTYPES,I0,I1,J0,J1,'yibs')
      call getmem3d(in_h,1,N_COVERTYPES,I0,I1,J0,J1,'yibs')
      call getmem3d(in_dbh,1,N_COVERTYPES,I0,I1,J0,J1,'yibs')
      call getmem3d(in_crad,1,N_COVERTYPES,I0,I1,J0,J1,'yibs')
      call getmem3d(in_texture,1,N_SOIL_TEXTURES,I0,I1,J0,J1,'yibs')
      call getmem3d(in_lai,1,N_COVERTYPES,I0,I1,J0,J1,'yibs')
      call getmem3d(in_cpool,1,N_COVERTYPES,I0,I1,J0,J1,'yibs')
      call getmem3d(in_albedo,1,N_COVERTYPES,I0,I1,J0,J1,'yibs')
      call getmem3d(in_tpool,1,N_PFT,I0,I1,J0,J1,'yibs')
      
      call getmem3d(Soiltemp,1,N_DEPTH,I0,I1,J0,J1,'yibs')
      call getmem3d(Soilmoist,1,N_DEPTH,I0,I1,J0,J1,'yibs')
      call getmem3d(Soilmp,1,N_DEPTH,I0,I1,J0,J1,'yibs')
      call getmem3d(fice,1,N_DEPTH,I0,I1,J0,J1,'yibs')
      
      !Coupling and diagnostic variables specific to YIBS (output)
      call getmem3d(laidata,1,N_COVERTYPES,I0,I1,J0,J1,'yibs') 
      call getmem3d(laitmp1,1,N_COVERTYPES,I0,I1,J0,J1,'yibs') 
      call getmem3d(laitmp,1,N_PFT,I0,I1,J0,J1,'yibs') 
      call getmem3d(laidatam,1,N_PFT,I0,I1,J0,J1,'yibs') 
      call getmem3d(fv0,1,N_COVERTYPES,I0,I1,J0,J1,'yibs') 
      call getmem3d(lai_pfts,1,N_COVERTYPES,I0,I1,J0,J1,'yibs') 
      call getmem3d(phen_pfts,1,N_COVERTYPES,I0,I1,J0,J1,'yibs') 
      call getmem3d(ht_pfts,1,N_COVERTYPES,I0,I1,J0,J1,'yibs') 
      call getmem3d(gpp_pfts,1,N_COVERTYPES,I0,I1,J0,J1,'yibs') 
      call getmem3d(ipp_pfts,1,N_COVERTYPES,I0,I1,J0,J1,'yibs') 
      call getmem3d(mtp_pfts,1,N_COVERTYPES,I0,I1,J0,J1,'yibs') 
      call getmem3d(Clive_pfts,1,N_COVERTYPES,I0,I1,J0,J1,'yibs') 
      call getmem3d(Cdead_pfts,1,N_COVERTYPES,I0,I1,J0,J1,'yibs') 
      
      call getmem2d(GCANOPY,I0,I1,J0,J1,'yibs')
      call getmem2d(Ci,I0,I1,J0,J1,'yibs')
      call getmem2d(NPP,I0,I1,J0,J1,'yibs')
      call getmem2d(GPP,I0,I1,J0,J1,'yibs')
      call getmem2d(GPPD,I0,I1,J0,J1,'yibs')
      call getmem2d(IPP,I0,I1,J0,J1,'yibs')
      call getmem2d(MTP,I0,I1,J0,J1,'yibs')
      call getmem2d(FO3,I0,I1,J0,J1,'yibs')
      call getmem2d(dFO3,I0,I1,J0,J1,'yibs')
      call getmem2d(TRANS_SW,I0,I1,J0,J1,'yibs')
      call getmem2d(z0,I0,I1,J0,J1,'yibs')
      call getmem2d(CO2flux,I0,I1,J0,J1,'yibs')
      call getmem2d(R_auto,I0,I1,J0,J1,'yibs')
      call getmem2d(lai0,I0,I1,J0,J1,'yibs')
      call getmem2d(vegfrac,I0,I1,J0,J1,'yibs')
      call getmem2d(lai_p,I0,I1,J0,J1,'yibs')
      call getmem2d(ht_p,I0,I1,J0,J1,'yibs')
      call getmem2d(phenf,I0,I1,J0,J1,'yibs')
      call getmem2d(landf,I0,I1,J0,J1,'yibs')
      call getmem2d(resp_r,I0,I1,J0,J1,'yibs')
      call getmem2d(resp_l,I0,I1,J0,J1,'yibs')
      call getmem2d(resp_w,I0,I1,J0,J1,'yibs')
      call getmem2d(carb_tot,I0,I1,J0,J1,'yibs')
      call getmem2d(carb_soil,I0,I1,J0,J1,'yibs')
         !call getmem2d(CO2c,I0,I1,J0,J1,'yibs')

      call getmem3d(soilt,I0,I1,J0,J1,1,N_DEPTH,'yibs')
      call getmem3d(soilm,I0,I1,J0,J1,1,N_DEPTH,'yibs')
      call getmem3d(soilp,I0,I1,J0,J1,1,N_DEPTH,'yibs')
      call getmem3d(soili,I0,I1,J0,J1,1,N_DEPTH,'yibs')

      call getmem3d(laipfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(htpfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(gpppfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(ipppfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(mtppfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(fvpfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(phenpfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(clivepfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(cdeadpfts,I0,I1,J0,J1,1,N_PFT,'yibs')

#ifdef OUTPUT_DAILY
      call getmem2d(agcanopy,I0,I1,J0,J1,'yibs')
      call getmem2d(aci,I0,I1,J0,J1,'yibs') 
      call getmem2d(anpp,I0,I1,J0,J1,'yibs')
      call getmem2d(agpp,I0,I1,J0,J1,'yibs') 
      call getmem2d(agppd,I0,I1,J0,J1,'yibs')
      call getmem2d(aco2flux,I0,I1,J0,J1,'yibs')
      call getmem2d(ar_auto,I0,I1,J0,J1,'yibs')
      call getmem2d(acoszen,I0,I1,J0,J1,'yibs')
      call getmem2d(ao3c,I0,I1,J0,J1,'yibs')
      call getmem2d(alai0,I0,I1,J0,J1,'yibs')
      call getmem2d(avegfrac,I0,I1,J0,J1,'yibs')
      call getmem2d(afo3,I0,I1,J0,J1,'yibs')
      call getmem2d(adfo3,I0,I1,J0,J1,'yibs')
      call getmem2d(alai_p,I0,I1,J0,J1,'yibs')
      call getmem2d(aht_p,I0,I1,J0,J1,'yibs')
      call getmem2d(aphenf,I0,I1,J0,J1,'yibs')
      call getmem2d(aipp,I0,I1,J0,J1,'yibs')
      call getmem2d(amtp,I0,I1,J0,J1,'yibs')
      call getmem2d(acarb_tot,I0,I1,J0,J1,'yibs')
      call getmem2d(acarb_soil,I0,I1,J0,J1,'yibs')
      call getmem2d(aresp_r,I0,I1,J0,J1,'yibs')
      call getmem2d(aresp_l,I0,I1,J0,J1,'yibs')
      call getmem2d(aresp_w,I0,I1,J0,J1,'yibs')
      
      call getmem3d(alaipfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(ahtpfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(agpppfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(afvpfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(aphenpfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(aipppfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(amtppfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(aclivepfts,I0,I1,J0,J1,1,N_PFT,'yibs')
      call getmem3d(acdeadpfts,I0,I1,J0,J1,1,N_PFT,'yibs')
#ifdef OUTPUT_FORCING
      call getmem2d(atairc,I0,I1,J0,J1,'yibs')
      call getmem2d(atcanopyc,I0,I1,J0,J1,'yibs')
      call getmem2d(aqf,I0,I1,J0,J1,'yibs')
      call getmem2d(ap_mbar,I0,I1,J0,J1,'yibs')
      call getmem2d(aca,I0,I1,J0,J1,'yibs')
      call getmem2d(ach,I0,I1,J0,J1,'yibs')
      call getmem2d(au,I0,I1,J0,J1,'yibs')
      call getmem2d(aipardir,I0,I1,J0,J1,'yibs')
      call getmem2d(aipardif,I0,I1,J0,J1,'yibs')
      call getmem3d(asoilt,I0,I1,J0,J1,1,N_DEPTH,'yibs')
      call getmem3d(asoilm,I0,I1,J0,J1,1,N_DEPTH,'yibs')
      call getmem3d(asoilp,I0,I1,J0,J1,1,N_DEPTH,'yibs')
      call getmem3d(asoili,I0,I1,J0,J1,1,N_DEPTH,'yibs')
#endif
#endif

      call getmem3d(betadl,1,N_DEPTH,I0,I1,J0,J1,'yibs')
      call getmem3d(albedo,1,N_BANDS,I0,I1,J0,J1,'yibs')
      !---------------------------------------------------------------------
      allocate(cells(I0:I1,J0:J1))
      call getmem2d(hemi,I0,I1,J0,J1,'yibs') !hemisphere flags = 1 for N., =-1 for S.
      call getmem2d(fw,I0,I1,J0,J1,'yibs') ! fraction of wet canopy
      
      
      if (.not.do_phenology_activegrowth .and.
     &   do_structuralgrowth) then
         print*,"impossible combinations of input parameters"
         stop
      endif

      print *,"yibs_input: "
     &     , jday, year, jday2, year2, yibs_dt
     &     , im, jm, i0f, i1f, j0f, j1f, force_VEG
     &     , do_soilinit, do_soilresp, do_phenology_activegrowth
     &     , do_structuralgrowth, do_frost_hardiness
     &     , do_patchdynamics, do_init_geo,do_spinup
     &     , do_islscp, do_clmpft
     &     , c3c4_mix, pft_cover, c3_alpha !, mixed_VEG


125   format(1x,a8, 6x, i8, i8, i8, f13.4, f13.4)
      
      fw(:,:) = 0.d0   ! all canopy is dry

      dlon = 360.d0/dble(IM)
      dlat = 180.d0/dble(JM-1)
      do i = 1,IM
         long(i) = -180.0d0 + dble(i-1)*dlon + 1.d0/2.d0
      enddo
      do j = 1, JM
         latg(j) = -90.0d0 + dble(j-1)*dlat
      enddo

      !* Now allocate optional arrays
      if ( force_VEG ) then
        call getmem3d( LAI,1,N_PFT,I0,I1,J0,J1,'yibs' )
        call getmem3d( height,1,N_PFT,I0,I1,J0,J1,'yibs' )
        !call getmem2d( heightm,1,N_PFT,I0,I1,J0,J1,'yibs' )
      else
        nullify( LAI )
        nullify( height )
        !nullify(heightm)
      endif

      !* Now YIBS should be initialized before any other calls to yibs_*
      call yibs_init_config(
     &     do_soilresp=do_soilresp
     &     ,do_phenology_activegrowth=do_phenology_activegrowth
     &     ,do_structuralgrowth=do_structuralgrowth
     &     ,do_frost_hardiness=do_frost_hardiness
     &     ,do_patchdynamics=do_patchdynamics)
!     &     ,mixed_VEG=mixed_VEG)

      !* Set hemisphere flags.
      !if ( J0<=JM/2 )   hemi(:,J0:min(JM/2,J1))   = -1    ! S.
      !if ( J1>=JM/2+1 ) hemi(:,max(JM/2+1,J0):J1) =  1    ! N.
      hemi = 1 ! all N.
      print *,"set hemi ok"
      
      !* Initialize yibs cells (makes head ycell).
      call yibs_cell_construct( cells )
      print *,"yibs_cell_construct passed ok" 

      !* Initialize vegetation cover.  
      !call readlai_m16(year,long,latg,IM,JM)
      !if ( associated(xlon_out) ) then
      !lat2d = xlat_out
      !lon2d = xlon_out
      lat2d = lm%xlat(I0:I1,J0:J1)
      lon2d = lm%xlon(I0:I1,J0:J1)
      !end if
      !call assignpnt(lm%xlat,lat2d)
      !call assignpnt(lm%xlon,lon2d)

      print*,"lm xlat 1,2:",size(lm%xlat,1),size(lm%xlat,2),
     & size(lm%xlon,1),size(lm%xlon,2)
      print*, IM, JM, I0, I1, J0, J1, jday, year
      call yibs_init_vegstruct( cells, IM, JM, I0, I1, J0, J1,
     &        jday, year, lat2d, lon2d, lm, lms,
c     &        in_lat,in_veg,in_pop,in_h,in_dbh,
c     &        in_crad,in_cpool,in_albedo,in_texture,
c     &        in_tpool,in_plant,in_harvest,in_lai,
     &        do_islscp, id_veg,
     &        do_clmpft, c3c4_mix, c3_alpha, pft_cover,
     &        do_soilinit,do_phenology_activegrowth,
     &        reinitialize=.true.)!,mixed_VEG)
      print *,"yibs_init_vegstruct passed ok"
      
!      max_time = 86400.*(365*(year2-year) + jday2)
!      print *,"max_time: ",max_time
      !time = 0.d0
      yibs_time = 86400.*(jday-1) !Allows starting run in middle of year at jday1.
      time_since_last_save = 0.d0
      jdaycount = jday
      ny     = year
      ny_old = -1
      nm_old = -1
      update_day = .true. !Initialize

#ifdef OUTPUT_DAILY
      nhrs = 86400./yibs_dt
#endif

      tdif = idate2 - idate1
      timeend = tohours(tdif) * secph
      extime = 0
      end subroutine yibs_initialize
      
      subroutine yibs(lm,lms,sfs,cpsb,nyear,nmm,ndd)
c      subroutine yibs(lm,lms,sfs,chia,atms,cpsb,sdelt)
      use flux_com
      use mod_dynparam
      use mod_mppparam
      use mod_memutil
      use mod_outvars
      
      implicit none
      type(lm_exchange) , intent(in) :: lm
      type(lm_state) , intent(in) :: lms
      type(surfstate) , intent(in) :: sfs
c      real(rk8) , pointer , intent(in) :: chia(:,:,:,:) 
      real(rk8) , pointer , intent(in) :: cpsb(:,:)
      integer, intent(in) :: nyear,nmm,ndd
      integer :: i,j,k,iflag
      integer :: ico2=38 , io3=7  !!! CBMZ
c        print *,"---------------------------------------------------"
c        print *,"started time step, time=", time, "yibs_dt=", yibs_dt
        !totdays = 0
        !nm  = 1
        !do while (totdays < jdaycount) 
        !   totdays = totdays + mdays(nm)
        !   nm = nm + 1
        !end do
        !nm = nm - 1
        !totdays = totdays - mdays(nm)
        !nd = jdaycount - totdays
        nh = mod(yibs_time/3600,24.)
        ny = nyear
        nm = nmm
        nd = ndd 
c        if ( myid == italk ) then
c        print *,"yt,dt,hour,day,month,year =",
c     &         nh,nd,nm,ny,nyear
c        end if

        ! get CO2 concentrations
        If (ny .ne. ny_old) then 
            call get_co2(ny, co2p)
            ny_old = ny
        Endif



#ifdef DEBUG        
      call assign_dummyvals(I0,I1,J0,J1,N_DEPTH,
     &  TairC,TcanopyC,Qf,P_mbar,Ca,Ch,U,
     &  IPARdif,IPARdir,CosZen, Soilmoist, Soilmp,fice)
#endif        


        tdy = jdaycount+(dble(nh))/24.0d0        
        call calc_solarzen(tdy,lat,cosz1)
        tdy = jdaycount+(dble(nh)+0.5)/24.0d0        
        call calc_solarzen(tdy,lat,cosz2)
        cosz = (cosz1+cosz2)/2.0d0
        Do i = I0, I1
        Do j = J0, J1
          tdy = jdaycount + dble(nh) + long(i)/15.0d0
          call calc_solarzen(tdy,latg(j),Cosz2d(i,j))
        Enddo
        Enddo

!       Get daily soil fraction
        call get_landf(jdaycount,landf,I0,I1,J0,J1,lm,lms)

!        call read_gcm(IM,JM,I0,I1,J0,J1,nm,nd,nh,
!     &      tag, tcang, qfg, pmg, cag, chg, wdg, ptg, prg, coszg, 
!     &      fwg, o3s, soiltg, soilmg, soilpg, soilig)

       ! print *,"yibs_time =",yibs_time
       if(mod(yibs_time,3600.).eq.0.and.yibs_time.ge.0) then
           P_mbar = (sfs%psa+5)*10 !!
           !P_mbar = ps_out/100. !!
           !print *,"yibs_time pmbar",yibs_time,maxVal(P_mbar)
           CosZen  = srf_cosz_out(:,:,1)
           Do k = 1,n_depth
           Do i = I0,I1
           DO j = J0,J1
           soiltemp(k,i,j)  = srf_stm_out(i,j,k)
           soilmoist(k,i,j) = srf_sq_out(i,j,k)
           soilmp(k,i,j)    = srf_smp_out(i,j,k)
           fice(k,i,j)      = srf_sif_out(i,j,k)
           if(soiltemp(k,i,j).gt.100) soiltemp(k,i,j) = -30
           if(soilmoist(k,i,j).gt.100) soilmoist(k,i,j) = 0
           if(fice(k,i,j).gt.100) fice(k,i,j) = 0
           enddo
           enddo
           enddo
           !print *,"yibs_time soiltemp",yibs_time,maxVal(soiltemp)

           Do i = I0,I1
           DO j = J0,J1
           Ch(i,j) = srf_ch_out(i,j,1)*1.e-6 !! gb------ umol/m2/s --> mol/m2/s
           if(Ch(i,j).gt.100) Ch(i,j) = 1.e-6
           if(Ch(i,j).lt.1.e-6) Ch(i,j) = 1.e-6
           TcanopyC(i,j) = srf_tcn_out(i,j,1)
           if(TcanopyC(i,j).gt.100) TcanopyC(i,j) = -30
           Qf(i,j) = srf_qcn_out(i,j,1)
           TairC(i,j) = srf_t2m_out(i,j,1)-273.15           
           U(i,j)  = sqrt(srf_u10m_out(i,j,1)**2 + 
     &         srf_v10m_out(i,j,1)**2)
        IPARdir(i,j)=2.3d0*max(srf_pard_out(i,j,1),zero)/4.55
        IPARdif(i,j)=2.3d0*max(srf_parf_out(i,j,1),zero)/4.55
           !IPARdir(i,j) = srf_pard_out(i,j,1)
           !IPardif(i,j) = srf_parf_out(i,j,1)
           fw(i,j) = srf_fwt_out(i,j,1)
c           Ca(i,j) = chia(i,j,17,ico2) / 
c     &          cpsb(i,j)*1.0e2*P_mbar(i,j)/8.314/(TairC(i,j)+tfrz)
         rtmp = 1.0d3*P_mbar(i,j)*100.0d0/287.05/(TcanopyC(i,j)+tfrz)
c         O3c(i,j) = chia(i,j,17,io3) / cpsb(i,j)
c         o3s(i,j) = O3c(i,j)*1.0e11*P_mbar(i,j)/8.314/(TairC(i,j)+tfrz)
         Ca(i,j) = CO2p*rtmp/28.9*1.0e-6
         enddo
         enddo
c         O3 = 0. !!!O3s
         O3 = O3s
c           print*,"here in yibs",maxVal(soiltemp),maxVal(srf_stm_out)     !maxVal(Qf),maxVal(TairC),maxVal(fw)
c     & ,maxVal(CosZen),maxVal(IPardif),maxVal(Ch),maxVal(P_mbar)
c     & ,maxVal(soiltemp),maxVal(soilmoist),maxVal(fice),maxVal(TcanopyC)
c     & ,maxVal(Ipardir),maxVal(ca),
c     &     maxVal(o3s),o3s(5,5)

c        !if ( myid == italk ) then
C       print *,"size",size(TairC,1),size(TairC,2),
C     &   "air_temperature",TairC(69,32),srf_t2m_out(69,32,1),
C     &   "canopy_temperature=",TcanopyC(69,32),srf_tcn_out(69,32,1),
C     &       "canopy_air_humidity=",Qf(69,32),srf_qcn_out(69,32,1),
C     &       "surf_pressure=",P_mbar(69,32),ps_out(69,32),
C     &       "surf_CO2=",Ca(69,32),
C     &       "surf_O3=",O3(69,32),
C     &       "heat_transfer_coef=",Ch(69,32),
C     &       "wind_speed=",U(69,32),
C     &       "dif_visible_rad=",IPARdif(69,32),
C     &       "direct_visible_rad=",IPARdir(69,32),
C     &       "cos_solar_zenith_angle=",CosZen(69,32),
C     &       "canopy_wet_fraction=",fw(69,32),
C     &       "soil_temp=",Soiltemp(1,69,32),  
C     &       "soil_moist=",Soilmoist(1,69,32),
C     &       "soil_matric_pot=",Soilmp(1,69,32),
C     &       "soil_ice_fraction=",fice(1,69,32)
        !end if


        call yibs_set_forcings( cells,
     &       air_temperature=TairC,
     &       canopy_temperature=TcanopyC,
     &       canopy_air_humidity=Qf,
     &       surf_pressure=P_mbar,
     &       surf_CO2=Ca,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     &       surf_O3=O3,
#endif
     &       heat_transfer_coef=Ch,
     &       wind_speed=U,
     &       total_visible_rad=IPARdir+IPARdif,
     &       direct_visible_rad=IPARdir,
     &       cos_solar_zenith_angle=CosZen,
     &       canopy_wet_fraction=fw,
     &       soil_temp=Soiltemp,  
     &       soil_moist=Soilmoist,
     &       soil_matric_pot=Soilmp,
     &       soil_ice_fraction=fice
     &       )
      !  print *, 'call yibs_set_forcings ok'  
       endif
        !print *,"yibs_time after call set_forcings",yibs_time
      !* NEW STREAMLINED CONTROL *!
      if (update_day) then
       !   print *, 'Got here, updating veg.'
        if (force_VEG) then
       !   print *, 'Got here, updating veg.1111111'
          call yibs_prescribe_vegupdate(cells,hemi,jdaycount,year,
     &         do_giss_phenology=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         do_giss_albedo=.true.,
     &         do_giss_lai=
     &         (.not.do_phenology_activegrowth).and.(.not.force_VEG),
     &         update_crops=.false.,
     &         laidata=LAI,     ! pft LAI
     &         hdata=height, init=.false. )   ! pft height
        !  print *,'yibs_prescribe_vegupdate with LAI and height'
        else !Don't include laidata and hdata in parameter list.
        !  print *, 'Got here, updating veg.2222222222'
          call getlai_m16(jdaycount,laitmp,I0,I1,J0,J1,do_clmpft,lm,lms)
        !  print *, 'Got here, getlai_m16'
          call prescr_get_laidata(jday,hemi,I0,I1,J0,J1,laitmp1)
        !  print *, 'Got here, prescr_get_laidata'
          do j=J0,J1
           do i=I0,I1
             do k=1,N_PFT
               if (laitmp(k,i,j) .le. 0.d0) then
                  laitmp(k,i,j) = laitmp1(k,i,j)
               endif
             enddo
           enddo
          enddo
          If (do_islscp) then
             do j = J0,J1
             do i = I0,I1
             laidatam(:,i,j) = 0.d0
             If (lms%mvegid(i,j) .gt. 0) 
     &          laidatam(lms%mvegid(i,j),i,j) = 3.0d0
             if (c3c4_mix .and. lms%mvegid(i,j) .eq. 15) then
                laidatam(13,i,j) = 3.0d0
                laidatam(15,i,j) = 3.0d0
             endif
             enddo 
             enddo
             laitmp = laidatam
          Endif
          !!!print*,'max lai: ',maxval(laitmp(13,15,:))
          call yibs_prescribe_vegupdate(cells,hemi,jdaycount,year,
     &         do_giss_phenology=.false.,
     &         do_giss_albedo=.true.,
     &         do_giss_lai=.false.,
     &         laidata=laitmp,
     &         update_crops=.false.,init=.false.)
        endif
      endif                     !update_day

      call yibs_run(cells,yibs_dt,update_day) !Change name to yibs_processes, calls yibs_integrate
 
#ifdef DEBUG
       ! print *,"Got here before yibs_get_exports."
#endif
        !* Extract results from YIBS structure for GCM or diagnostics.
        call yibs_get_exports( cells,
     &       canopy_conductance=GCANOPY,
     &       beta_soil_layers=betadl,
     &       shortwave_transmit=TRANS_SW,
     &       leafinternal_CO2=Ci,
     &       foliage_humidity=Qf,
     &       canopy_npp=NPP,
     &       canopy_gpp=GPP,
     &       canopy_gpp0=GPPD,
     &       canopy_resp_r=Resp_r,
     &       canopy_resp_l=Resp_l,
     &       canopy_resp_w=Resp_w,
     &       canopy_carb_tot=Carb_tot,
     &       canopy_carb_soil=Carb_soil,
     &       canopy_gpp_pfts=GPP_pfts,
     &       canopy_ipp_pfts=IPP_pfts,
     &       canopy_mtp_pfts=MTP_pfts,
     &       canopy_clive_pfts=Clive_pfts,
     &       canopy_cdead_pfts=Cdead_pfts,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     &       ozone_flux=FO3,
     &       excess_ozone_flux=dFO3,
#endif
#ifdef ACTIVE_GROWTH
     &       LAI_prognostic=lai_p,
     &       height_prognostic=ht_p,
     &       Phenf_pfts=phen_pfts,
#endif
     &       LAI_pfts=lai_pfts,
     &       height_pfts=ht_pfts,
     &       roughness_length=z0,
     &       flux_CO2=CO2flux,
     &       R_auto=R_auto,
     &       albedo=albedo,
     &       leaf_area_index=lai0,
     &       vegetation_fractions=fv0,
     &       canopy_ipp=IPP,
     &       canopy_mtp=MTP
     &     )
      
        do j = j0,j1
        do i = i0,i1
           npp(i,j)     = npp(i,j)*landf(i,j) !*6.31/9.98
           gpp(i,j)     = gpp(i,j)*landf(i,j) !*6.31/9.98   !!!! xiaoxie 2018-11-02 for year 2013
           gppd(i,j)    = gppd(i,j)*landf(i,j)
           ipp(i,j)     = ipp(i,j)*landf(i,j)
           mtp(i,j)     = mtp(i,j)*landf(i,j)
           resp_r(i,j)  = resp_r(i,j)*landf(i,j)
           resp_l(i,j)  = resp_l(i,j)*landf(i,j)
           resp_w(i,j)  = resp_w(i,j)*landf(i,j)
           r_auto(i,j)  = r_auto(i,j)*landf(i,j)
           co2flux(i,j) = co2flux(i,j)*landf(i,j)  !*1.4/2.08 !!!! xiaoxie 2018-11-02 for year 2013
           carb_tot(i,j)= carb_tot(i,j)*landf(i,j)
           carb_soil(i,j)= carb_soil(i,j)*landf(i,j)
           vegfrac(i,j) = sum(fv0(:,i,j))
           do k = 1,n_depth
              soilt(i,j,k) = soiltemp(k,i,j)
              soilm(i,j,k) = soilmoist(k,i,j)
              soilp(i,j,k) = soilmp(k,i,j)
              soili(i,j,k) = fice(k,i,j)
           enddo
           do k = 1,n_pft
              laipfts(i,j,k)   = lai_pfts(k,i,j)
              in_pop(k,i,j)   = lai_pfts(k,i,j)  !!!! xiaoxie 2018-6-22
              htpfts(i,j,k)    = ht_pfts(k,i,j)
              in_h(k,i,j)    = ht_pfts(k,i,j)  !!! xie
              clivepfts(i,j,k) = clive_pfts(k,i,j)
              cdeadpfts(i,j,k) = cdead_pfts(k,i,j)
              in_dbh(k,i,j) = cdead_pfts(k,i,j) !!! xie
              gpppfts(i,j,k)   = gpp_pfts(k,i,j)*landf(i,j)
              in_crad(k,i,j)   = gpp_pfts(k,i,j)*landf(i,j) !!!! xiaoxie 2018-6-22
              ipppfts(i,j,k)   = ipp_pfts(k,i,j)*landf(i,j)
              mtppfts(i,j,k)   = mtp_pfts(k,i,j)*landf(i,j)
              fvpfts(i,j,k)    = fv0(k,i,j)
              in_veg(k,i,j)    = fv0(k,i,j)
#ifdef ACTIVE_GROWTH
              phenpfts(i,j,k) = phen_pfts(k,i,j)
#endif
           enddo
        enddo
        enddo



        !* Write yibs state to a restart file.
        time_since_last_save = time_since_last_save + yibs_dt
        if ( time_since_last_save > save_interval ) then
          call yibs_write_state( cells )
          time_since_last_save = 0.d0
        endif

        yibs_time = yibs_time + yibs_dt  !Moved yibs_time update to before check of update_day-NK
        update_day=(dmod(yibs_time,86400.d0) .eq. 0.d0)
     &  .and.(yibs_time.ne.0.d0)

        if (update_day) then
          jdaycount = jdaycount + 1
          if (jdaycount>365) then !numdays=365 for 1-year data, etc.
!YK - to accomodate the case in which the starting jday is not 1.
!          if (jdaycount-jday+1 > 365) then
            jdaycount = 1       !reset
            ny = ny + 1
          endif
          !jdaycount = mod(jdaycount,365) !mod(365) goes from 0 to 364 for days 365 and 1-364.  Use this if jday is not read in with forcings.
        endif

      extime = extime + yibs_dt  
      end subroutine yibs

!************************************************************************
      subroutine yibs_read_state( cells )
!@sum read yibs state from the file
      type(ycelltype_public), intent(out) :: cells(:,:)
      !---
      integer, parameter :: MAX_BUFFER=100000 ! need realistic estimate
      real*8 buffer(MAX_BUFFER)
      integer iu_yibsstate
      integer ic, jc, i, j

      ic = size(cells,1)
      jc = size(cells,2)

      call openunit('yibs_state',iu_yibsstate,.true.,.true.)
      do j=1,jc
        do i=1,ic
          read(iu_yibsstate) buffer
          ! check length of buffer : if( buffer(1) > MAX_BUFFER ) ??
          call yibs_cell_unpack(buffer, cells(i,j))
        enddo
      enddo
      call closeunit(iu_yibsstate)

      end subroutine yibs_read_state

!************************************************************************
      subroutine yibs_write_state( cells )
!@sum write yibs state to the file
      type(ycelltype_public), intent(in) :: cells(:,:)
      !---
      real*8, pointer :: buffer(:)
      integer iu_yibsstate
      integer ic, jc, i, j

      ic = size(cells,1)
      jc = size(cells,2)

      call openunit('yibs_state_new',iu_yibsstate,.true.,.false.)
      do j=1,jc
        do i=1,ic
          call yibs_cell_pack(buffer, cells(i,j))
          write(iu_yibsstate) buffer
          deallocate(buffer)
        enddo
      enddo
      call closeunit(iu_yibsstate)

      end subroutine yibs_write_state


!************************************************************************

      subroutine yibs_init_vegstruct( cells,
     &     IM, JM, I0, I1, J0, J1, jday, year, lat, lon,
     &     lm, lms,
c     &     in_lat,in_veg,in_pop,in_h,in_dbh,
c     &     in_crad,in_cpool,in_albedo,in_texture,
c     &     in_tpool,in_plant,in_harvest,in_lai,
     &     do_islscp, id_veg, 
     &     do_clmpft, c3c4_mix,  c3_alpha, pft_cover, 
     &     do_soilinit, do_phenology_activegrowth, reinitialize)!,mixed_VEG)
      use yibs_prescribed_drv, only: init_canopy_physical,prescr_vegdata
     &     ,yibs_init_params
      !use yibs_prescr_veg, only : prescr_calcconst
#ifdef YIBS_1D_DIAG
      use yibs_prescr_veg, only : print_yibs_pfts
      use yibs_pfts, only : alamin, alamax
#endif
      use yibs_const
      use vegfor_com, only: mfrac, mcrop, mvegid, mcrop_cal, 
     $                      mfrac_clm, mfrac_isl
      use geos_com,   only: elev

      implicit none
      type(ycelltype_public), intent(inout) :: cells(I0:I1,J0:J1)
      integer, intent(in) :: IM, JM, I0, I1, J0, J1, jday
      integer, intent(in) :: year
      type(lm_exchange) , intent(in) :: lm
      type(lm_state) , intent(in) :: lms
      real*8, dimension(I0:I1,J0:J1) :: lat, lon
      integer, intent(in) :: id_veg
      logical, intent(in) :: do_islscp, do_clmpft
      logical, intent(in) :: c3c4_mix, pft_cover, reinitialize
      logical, intent(in) :: do_soilinit, do_phenology_activegrowth
      real*8, intent(in)  :: c3_alpha
! output for debug by xiaoxie
      !real*8, intent(out) :: in_lat(I0:I1,J0:J1) !cohort
      !real*8, intent(out) :: in_veg(N_COVERTYPES,I0:I1,J0:J1) !cohort
      !real*8, intent(out) :: in_pop(N_COVERTYPES,I0:I1,J0:J1)  !cohort
      !real*8, intent(out) :: in_h(N_COVERTYPES,I0:I1,J0:J1)
      !real*8, intent(out) :: in_dbh(N_COVERTYPES,I0:I1,J0:J1)
      !real*8, intent(out) :: in_crad(N_COVERTYPES,I0:I1,J0:J1) !cohort
      !real*8, intent(out) :: in_cpool(N_COVERTYPES,I0:I1,J0:J1) !cohort
      !real*8, intent(out) :: in_albedo(N_COVERTYPES,I0:I1,J0:J1) !cohort
      !real*8, intent(out) :: in_texture(N_SOIL_TEXTURES,I0:I1,J0:J1) !cohort
      !real*8, intent(out) :: in_tpool(N_PFT,I0:I1,J0:J1)  !cohort
      !real*8, intent(out) :: in_plant(I0:I1,J0:J1)
      !real*8, intent(out) :: in_harvest(I0:I1,J0:J1)      
      !real*8, intent(out) :: in_lai(N_COVERTYPES,I0:I1,J0:J1)
!      logical, intent(in) :: mixed_VEG
      !---Local variables-----
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: vegdata !cohort
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: vegdata2 !cohort
      real*8, dimension(N_BANDS,N_COVERTYPES,I0:I1,J0:J1) :: albedodata !patch, NOTE:snow
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: laidata  !cohort
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: hdata    !cohort
      real*8, dimension(N_COVERTYPES) :: nmdata    !cohort
      real*8, dimension(N_COVERTYPES,N_DEPTH) :: rootprofdata !Root fraction of veg type.
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: popdata !Dummy population density:  0-bare soil, 1-vegetated
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: dbhdata !Diameter at breast height for woody veg.(cm)
      real*8, dimension(N_COVERTYPES,I0:I1,J0:J1) :: craddata !Crown radius (m)
      real*8, dimension(N_COVERTYPES,N_BPOOLS,I0:I1,J0:J1) :: cpooldata !Carbon pools in individuals
      integer, dimension(N_COVERTYPES) :: soildata ! soil types 1-bright 2-dark
      real*8, dimension(N_SOIL_TEXTURES,I0:I1,J0:J1) :: soil_texture
      real*8, dimension(I0:I1,J0:J1) :: Ci_ini,CNC_ini,Tcan_ini,Qf_ini
      real*8, dimension(I0:I1,J0:J1) :: plant_date,harvest_date
      real*8, dimension(I0:I1,J0:J1) :: cropdata2
      real*8, dimension(N_PFT,PTRACE,NPOOLS-NLIVE,N_CASA_LAYERS,
     &                  I0:I1,J0:J1):: Tpooldata  !in g/m2 
      !------
      integer :: iu,i,j, k,ii,jj
      real*8  :: s,dlon,dlat
      
      call yibs_init_params()

#ifdef YIBS_1D_DIAG
      call print_yibs_pfts()
      write(*,*) "alamin: ", alamin
      write(*,*) "alamax: ", alamax
#endif
 
#ifdef PFT_MODEL_YIBS
      write(*,*) "PFT_MODEL_YIBS defined"
#endif

      if (c3c4_mix) then
         if (c3_alpha .gt. 1. .or. c3_alpha .lt. 0.) then
            print*,"C3_alpha should be between 0-1"
            stop
         endif
      endif


         !* Mosaicked subgrid fractions. *!
         !* NOTE:  This initializes with default LAI.
         dlon = long(10) - long(9)
         dlat = latg(10) - latg(9)
         print*,'dlon,dlat: ',dlon,dlat
         vegdata2 = 0.0d0
         do j=J0,J1
         do i=I0,I1
         cropdata2(i,j)=lms%vcrop(i,j)
         plant_date(i,j)=lms%crop_cal(i,j,1)
         harvest_date(i,j)=lms%crop_cal(i,j,2)
         if(plant_date(i,j).gt.365) plant_date(i,j)=0
         if(harvest_date(i,j).gt.365) harvest_date(i,j)=0
         do k=1,18
         If (do_clmpft) then
         vegdata2(k,i,j)=lms%clmfrac(i,j,k)
         if(vegdata2(k,i,j).gt.100)vegdata2(k,i,j)=0
         else
         vegdata2(k,i,j)=lms%vfrac(i,j,k)
         Endif
         enddo
         if (cropdata2(i,j) .gt. 1.) then
         print*,'here before prescr_vegdata: ',k,i,j,cropdata2(i,j)
         endif
         enddo
         enddo

         print*,"input vegdata: ",I0,I1,J0,J1,
     &  vegdata(3,15,15),vegdata2(3,15,15),lms%clmfrac(15,15,3)

         call prescr_vegdata(lm,lms,jday, year, 
     &     IM,JM,I0,I1,J0,J1,vegdata,albedodata,laidata,hdata,nmdata,
     &     popdata,dbhdata,craddata,cpooldata,rootprofdata,
     &     soildata,soil_texture,Tpooldata, 
     &     vegdata2, cropdata2, elev, do_clmpft,
     &     do_soilinit,do_phenology_activegrowth,do_init_geo=.false.,
     &     do_read_from_files=.true.)

         print*,"input vegdata: ",I0,I1,J0,J1,
     &  vegdata(3,15,15),vegdata2(3,15,15)
         do j=J0,J1
         do i=I0,I1
         if (laidata(16,i,j).gt.10.or.laidata(16,i,j).lt.0) then
         print*,"lai data error: ",i,j,laidata(16,i,j)
         endif
         enddo
         enddo
         
         If (do_islscp) then
         do j=J0,J1
         do i=I0,I1
           vegdata(:,i,j) = 0.d0
           laidata(:,i,j) = 0.d0
           If (pft_cover) then
              do k=1,16
                vegdata(k,i,j)=lms%frac_isl(i,j,k)
              enddo
           else
              If (lms%mvegid(i,j) .gt. 0) 
     &     vegdata(lms%mvegid(i,j),i,j) = 1.0d0
           Endif
           do k = 1,16
              If (vegdata(k,i,j) .gt. 0.0d0) laidata(k,i,j) = 2.0d0
           enddo
         enddo
         enddo
         Endif


         do j=J0,J1
         do i=I0,I1
           if (c3c4_mix .and. vegdata(15,i,j) .gt. 0.) then
              vegdata(13,i,j) = vegdata(13,i,j)
     &                        + c3_alpha*vegdata(15,i,j)
              vegdata(15,i,j) = (1-c3_alpha)*vegdata(15,i,j)
           endif
           do k = 1,18
              if (vegdata(k,i,j) .gt. 1.) then
               print*,'here before sum: ',k,i,j,vegdata(k,i,j)
              endif
           enddo
           s = sum( vegdata(:,i,j) )
           If (s .gt. 1.0d0) vegdata(:,i,j) = vegdata(:,i,j)/s
           do k = 1,18
              if (vegdata(k,i,j) .gt. 1.) then 
               print*,'here before yibs_cell set: ',k,i,j,vegdata(k,i,j)
              endif
           enddo
         enddo
         enddo
                 
         !Translate gridded data to YIBSdata structure
         if (reinitialize) 
     &   call init_canopy_physical(I0, I1, J0, J1,
     &        Ci_ini, CNC_ini, Tcan_ini, Qf_ini)

         call yibs_cell_set(cells, jday, lat, vegdata, popdata,laidata,
     &        hdata, dbhdata, craddata, cpooldata, nmdata, rootprofdata,
     &        soildata, albedodata, soil_texture,
     &        Ci_ini, CNC_ini, Tcan_ini, Qf_ini, Tpooldata,
     &        plant_date, harvest_date,
     &        reinitialize)  

         do j=J0,J1
         do i=I0,I1
         if (laidata(16,i,j).gt.10.or.laidata(16,i,j).lt.0) then
         print*,"lai data error2222: ",i,j,laidata(16,i,j)
         endif
         enddo
         enddo     
      !in_lat = lat
      !in_veg = vegdata
      !print*,"input vegdata: ",I0,I1,J0,J1,vegdata(3,15,15)
      !in_pop = popdata
      !in_h = hdata
      !in_dbh = dbhdata
      !in_crad = craddata
      !in_cpool = cpooldata(:,1,:,:)
      !in_nm = nmdata
      !in_root = rootprofdata
      !in_soil = soildata
      !in_albedo = albedodata(1,:,:,:)
      !in_texture = soil_texture
      !in_tpool = tpooldata(:,1,1,1,:,:)
      !in_plant = plant_date
      !in_harvest =  harvest_date
      !in_lai = laidata
      
      end subroutine yibs_init_vegstruct

      subroutine calc_solarzen(td,latdegrees,sbeta1)
      !* Calculate solar zenith angle **in radians**
      !* From Spitters, C. J. T. (1986), AgForMet 38: 231-242.
      implicit none
      real*8,intent(in) :: td             ! day(to minute fraction)
      real*8,intent(in) :: latdegrees     ! latitude in degrees
      real*8,parameter :: pi = 3.1415926535897932d0 !@param pi    pi
      real*8,parameter :: rad = pi/180.d0 ! Conversion from degrees to radians.
      real*8 :: hour,latrad
      real*8 :: delta                     ! declination angle
      real*8 :: td0
      real*8,intent(out) :: sbeta1        ! sbeta1=cos(zen angle)=sin(elev angle)
!      real*8,intent(out) :: solarelev    ! solar elevation angle (rad)
!      real*8,intent(out) :: solarzen     ! solar zenith angle (rad)
      
      td0 = td
      If (td0 .lt. 0.d0) td0 = td0 + 365.0d0
      If (td0 .gt. 365.0d0) td0 = td0 - 365.0d0
      hour = (td0-floor(td0))*24.d0
      latrad = latdegrees*rad
      delta = asin(-sin(rad*23.45d0)*cos(2.d0*pi*(td0+10.d0)/365.d0))
      sbeta1 = sin(latrad)*sin(delta)+
     &     cos(latrad)*cos(delta)*cos(rad* 15.d0*(hour-12.d0))
c      if (sbeta1 < 0.d0) sbeta1 = 0.d0  !**GCM does this too** 
 
      end subroutine calc_solarzen
      
      end module run_yibs
