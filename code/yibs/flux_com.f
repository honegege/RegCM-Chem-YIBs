      module flux_com

      implicit none
      save

      integer ncyid, start(3), ncount(3)
      integer lonyid, latyid, timeyid, vlonyid, vlatyid, vtimyid
      integer gcyid, ciyid, coszyid, nppyid, gppyid, gpp0yid, royid
      integer fo3yid, dfo3yid, ippyid, mtpyid, cflxyid
      integer laisyid, htsyid, gppsyid, fvsyid, ippsyid, mtpsyid
      integer clsyid, cdsyid
      integer respryid, resplyid, respwyid, ctotyid, csoilyid
      integer start3(4), ncut3(4), pftyid, vpftyid, npt(4)
#ifdef ACTIVE_GROWTH
      integer laipyid, htpyid, phenyid
#endif
#ifdef OUTPUT_FORCING
      integer tayid, tcyid, qcyid, psyid, cayid, chyid, wdyid
      integer pdiryid, pdifyid
      integer styid, smyid, spyid, siyid
      Integer start2(4), ncut2(4), levyid, vlevyid, nzt(4)
#endif
      integer laiyid, fvyid, o3yid
      integer, parameter :: nrecmax = 12000
      integer nhour(nrecmax), nday(31)
      integer :: ny_old2 = -1
      logical :: ifirst = .true.

      public :: open_output_netcdf, close_output_netcdf
      public ::  write_output_netcdf
	  contains


      subroutine open_output_netcdf(I0,I1,J0,J1,lon2d,lat2d,ny,nm)


      implicit none
      include 'netcdf.inc'

      integer, intent(in)  ::  I0, I1, J0, J1, ny, nm
      real*8, intent(in)   :: lon2d(I0:I1,J0:J1), lat2d(I0:I1,J0:J1)
      integer vardim(3)
      integer nstat
      integer lenc
      character*120 fsite
      character*4 cyear
      character*2 cmon

      write(cyear,'(I4.4)') ny
      write(cmon,'(I2.2)') nm
      fsite = 'flux.'//cyear//cmon//'.nc'

      If (nm .ne. ny_old2 .or. ifirst) then
      
      If (.not.ifirst) call close_output_netcdf
      ny_old2 = nm
      ifirst = .false.

      nstat = NF_CREATE(fsite(1:lenc(fsite)),NF_CLOBBER,ncyid)
      IF(nstat .NE. NF_NOERR) Call handle_err(nstat)

         !print*, 'here start creating  output : ',ncyid
CCCCCC   DEFINE DIMENSIONS    CCCCCCCC

      nstat = NF_DEF_DIM(ncyid, 'lon', I1-I0+1, lonyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_DEF_DIM(ncyid, 'lat', J1-J0+1, latyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_DEF_DIM(ncyid, 'pfts', 16, pftyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_DEF_DIM(ncyid, 'time',NF_UNLIMITED, timeyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      vardim(1) = lonyid
      vardim(2) = latyid
      vardim(3) = timeyid
      start(1)  = 1
      start(2)  = 1
      start(3)  = 1
      ncount(1) = I1-I0+1
      ncount(2) = J1-J0+1
      ncount(3) = 1

      npt(1) = lonyid
      npt(2) = latyid
      npt(3) = pftyid
      npt(4) = timeyid
      start3(1) = 1
      start3(2) = 1
      start3(3) = 1
      start3(4) = 1
      ncut3(1)  = I1-I0+1
      ncut3(2)  = J1-J0+1
      ncut3(3)  = 16
      ncut3(4)  = 1

#ifdef OUTPUT_FORCING
      nstat = NF_DEF_DIM(ncyid, 'depth', 10, levyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nzt(1) = lonyid
      nzt(2) = latyid
      nzt(3) = levyid
      nzt(4) = timeyid
      start2(1) = 1
      start2(2) = 1
      start2(3) = 1
      start2(4) = 1
      ncut2(1)  = I1-I0+1
      ncut2(2)  = J1-J0+1
      ncut2(3)  = 10
      ncut2(4)  = 1
#endif

      nstat = NF_PUT_ATT_TEXT(ncyid,NF_GLOBAL,'title',34,
     &       'Gryidded standalone vegetation flux')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_DEF_VAR(ncyid,'lon',NF_DOUBLE,2
     &  ,(/lonyid,latyid/),vlonyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_DEF_VAR(ncyid,'lat',NF_DOUBLE,2
     &  ,(/lonyid,latyid/),vlatyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_DEF_VAR(ncyid,'time',NF_INT,1,timeyid,vtimyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
#ifdef OUTPUT_FORCING
      nstat = NF_DEF_VAR(ncyid,'depth',NF_INT,1,levyid,vlevyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
#endif

      nstat = NF_PUT_ATT_TEXT(ncyid,vlonyid,'long_name',9,
     &       'Longitude')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid,vlonyid,'units',9,
     &       'degrees_e')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid,vlatyid,'long_name',8,
     &       'Latitude')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid,vlatyid,'units',9,
     &       'degrees_n')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid,vtimyid,'long_name',4,
     &       'Time')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
#ifndef OUTPUT_DAILY
      nstat = NF_PUT_ATT_TEXT(ncyid,vtimyid,'units',30,
     &       'hours since '//cyear//'-1-1 00:00:0.0')
#else 
      nstat = NF_PUT_ATT_TEXT(ncyid,vtimyid,'units',29,
     &       'days since '//cyear//'-1-1 00:00:0.0')
#endif
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid,vtimyid,'calendar',7,
     &       '365_day')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

#ifdef OUTPUT_FORCING
      nstat = NF_PUT_ATT_TEXT(ncyid,vlevyid,'long_name',5,
     &       'Depth')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
#endif

      nstat = NF_DEF_VAR(ncyid,'gc',NF_FLOAT,3,vardim,gcyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, gcyid,'long_name',33,
     &       'Canopy conductance of water vapor')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, gcyid,'units',5,
     &       'm s-1')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'ci',NF_FLOAT,3,vardim,ciyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, ciyid,'long_name',34,
     &       'Internal foliage CO2 concentration')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, ciyid,'units',7,
     &       'mol m-3')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         !print*, 'here start creating  output ci : ',ciyid
		 
      nstat = NF_DEF_VAR(ncyid,'cosz',NF_FLOAT,3,vardim,coszyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, coszyid,'long_name',28,
     &       'cosine of solar zenith angle')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, coszyid,'units',4,
     &       'none')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'NPP',NF_FLOAT,3,vardim, nppyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, nppyid,'long_name',24,
     &       'Net primary productivity')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, nppyid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'GPP',NF_FLOAT,3,vardim,gppyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, gppyid,'long_name',26,
     &       'Gross primary productivity')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, gppyid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'GPP0',NF_FLOAT,3,vardim,gpp0yid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, gpp0yid,'long_name',34,
     &       'Gross primary productivity offline')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, gpp0yid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

#ifdef PS_BVOC

      nstat = NF_DEF_VAR(ncyid,'Isoprene',NF_FLOAT,3,vardim,ippyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, ippyid,'long_name',17,
     &       'Isoprene emission')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, ippyid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Monoterpene',NF_FLOAT,3,vardim,mtpyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, mtpyid,'long_name',20,
     &       'Monoterpene emission')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, mtpyid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

#endif

      nstat = NF_DEF_VAR(ncyid,'R_auto',NF_FLOAT,3,vardim,royid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, royid,'long_name',23,
     &       'Autotrophic respiration')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, royid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'CO2flux',NF_FLOAT,3,vardim,cflxyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, cflxyid,'long_name',15,
     &       'Net CO2 flux up')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, cflxyid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'resp_r',NF_FLOAT,3,vardim, respryid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, respryid,'long_name',16,
     &       'Root respiration')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, respryid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'resp_l',NF_FLOAT,3,vardim, resplyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, resplyid,'long_name',16,
     &       'Leaf respiration')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, resplyid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'resp_w',NF_FLOAT,3,vardim, respwyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, respwyid,'long_name',16,
     &       'Wood respiration')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, respwyid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Carbon_total',NF_FLOAT,
     &              3,vardim, ctotyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, ctotyid,'long_name',25,
     &       'Total Land Carbon Storage')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, ctotyid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Carbon_soil',NF_FLOAT,
     &      3,vardim, csoilyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, csoilyid,'long_name',25,
     &       'Total Soil Carbon Storage')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, csoilyid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'lai_pfts',NF_FLOAT,4,npt,laisyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, laisyid,'long_name',17,
     &       'PFT-specified LAI')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, laisyid,'units',5,
     &       'm2/m2')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'ht_pfts',NF_FLOAT,4,npt,htsyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, htsyid,'long_name',16,
     &       'PFT-specified HT')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, htsyid,'units',1,
     &       'm')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'clive_pfts',NF_FLOAT,4,npt,clsyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, clsyid,'long_name',31,
     &       'PFT-specified Total Live Carbon')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, clsyid,'units',5,
     &       'Kg/m2')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'cdead_pfts',NF_FLOAT,4,npt,cdsyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, cdsyid,'long_name',31,
     &       'PFT-specified Total Dead Carbon')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, cdsyid,'units',5,
     &       'Kg/m2')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'gpp_pfts',NF_FLOAT,4,npt,gppsyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, gppsyid,'long_name',17,
     &       'PFT-specified GPP')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, gppsyid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

#ifdef PS_BVOC

      nstat = NF_DEF_VAR(ncyid,'ipp_pfts',NF_FLOAT,4,npt,ippsyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, ippsyid,'long_name',22,
     &       'PFT-specified Isoprene')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, ippsyid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'mtp_pfts',NF_FLOAT,4,npt,mtpsyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, mtpsyid,'long_name',22,
     &       'PFT-specified Monoterpene')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, mtpsyid,'units',10,
     &       'kg[C]/m2/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

#endif

      nstat = NF_DEF_VAR(ncyid,'fv_pfts',NF_FLOAT,4,npt,fvsyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, fvsyid,'long_name',16,
     &       'PFT-specified HT')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, fvsyid,'units',8,
     &       'fraction')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)

      nstat = NF_DEF_VAR(ncyid,'FO3',NF_FLOAT,3,vardim,fo3yid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, fo3yid,'long_name',21,
     &       'ozone flux to stomata')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, fo3yid,'units',15,
     &       'nmol O3 m-2 s-1')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'dFO3',NF_FLOAT,3,vardim,dfo3yid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, dfo3yid,'long_name',28,
     &       'excess ozone flux to stomata')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, dfo3yid,'units',15,
     &       'nmol O3 m-2 s-1')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

#endif

#ifdef ACTIVE_GROWTH

      nstat = NF_DEF_VAR(ncyid,'lai_p',NF_FLOAT,3,vardim,laipyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, laipyid,'long_name',14,
     &       'prognostic LAI')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, laipyid,'units',6,
     &       'm2 m-2')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'ht_p',NF_FLOAT,3,vardim,htpyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, htpyid,'long_name',17,
     &       'prognostic height')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, htpyid,'units',1,
     &       'm')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'phen_pfts',NF_FLOAT,4,npt,phenyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, phenyid,'long_name',25,
     &       'PFT-specified phenofactor')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, phenyid,'units',8,
     &       'fraction')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

#endif

#ifdef OUTPUT_FORCING
      nstat = NF_DEF_VAR(ncyid,'Tair',NF_FLOAT,3,vardim,tayid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, tayid,'long_name',23,
     &       'Surface air temperature')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, tayid,'units',8,
     &       'degree_C')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Tcanopy',NF_FLOAT,3,vardim,tcyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, tcyid,'long_name',22,
     &       'Canopy air temperature')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, tcyid,'units',8,
     &       'degree_C')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Qcanopy',NF_FLOAT,3,vardim,qcyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, qcyid,'long_name',19,
     &       'Canopy air humyidity')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, qcyid,'units',4,
     &       'g/kg')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Psurf',NF_FLOAT,3,vardim,psyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, psyid,'long_name',16,
     &       'Surface pressure')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, psyid,'units',3,
     &       'hPa')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Ca',NF_FLOAT,3,vardim,cayid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, cayid,'long_name',13,
     &       'Surface [CO2]')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, cayid,'units',7,
     &       'umol/m3')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Ch',NF_FLOAT,3,vardim,chyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, chyid,'long_name',24,
     &       'Heat trasfer coefficient')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Wind',NF_FLOAT,3,vardim,wdyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, wdyid,'long_name',18,
     &       'Surface wind speed')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, wdyid,'units',3,
     &       'm/s')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Pardir',NF_FLOAT,3,vardim,pdiryid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, pdiryid,'long_name',10,
     &       'Direct PAR')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, pdiryid,'units',4,
     &       'W/m2')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Pardif',NF_FLOAT,3,vardim,pdifyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, pdifyid,'long_name',11,
     &       'Diffuse PAR')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, pdifyid,'units',4,
     &       'W/m2')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Soiltemp',NF_FLOAT,4,nzt,styid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, styid,'long_name',16,
     &       'Soil temperature')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, styid,'units',8,
     &       'degree_C')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Soilmoist',NF_FLOAT,4,nzt,smyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, smyid,'long_name',13,
     &       'Soil moisture')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, smyid,'units',10,
     &       'percentage')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'Soilmp',NF_FLOAT,4,nzt,spyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, spyid,'long_name',21,
     &       'Soil matrix potential')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'fice',NF_FLOAT,4,nzt,siyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, siyid,'long_name',17,
     &       'Soil ice fraction')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, siyid,'units',10,
     &       'percentage')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

#endif

      nstat = NF_DEF_VAR(ncyid,'O3conc',NF_FLOAT,3,vardim,o3yid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, o3yid,'long_name',16,
     &       'O3 concentration')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, o3yid,'units',3,
     &       'ppb')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'lai',NF_FLOAT,3,vardim,laiyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, laiyid,'long_name',15,
     &       'Leaf area index')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, laiyid,'units',6,
     &       'm2 m-2')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_DEF_VAR(ncyid,'vfrac',NF_FLOAT,3,vardim,fvyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, fvyid,'long_name',18,
     &       'Total veg fraction')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_ATT_TEXT(ncyid, fvyid,'units',8,
     &       'fraction')
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_ENDDEF(ncyid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      nstat = NF_PUT_VAR_DOUBLE(ncyid,vlonyid,lon2d(I0:I1,J0:J1))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VAR_DOUBLE(ncyid,vlatyid,lat2d(I0:I1,J0:J1))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      Endif

      end subroutine open_output_netcdf


      subroutine close_output_netcdf


      implicit none
      include 'netcdf.inc'

      integer nstat

      nstat=NF_CLOSE(ncyid)
      IF(nstat .NE. NF_NOERR) Call handle_err(nstat)

      end subroutine close_output_netcdf


      subroutine write_output_netcdf(I0,I1,J0,J1,nm,nd,nh, 
     &       GCANOPY, Ci, NPP, GPP, GPP0, IPP, MTP, CO2flux, 
     &       resp_r, resp_l, resp_w, c_tot, c_soil,
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
     &       FO3, dFO3,
#endif
#ifdef ACTIVE_GROWTH
     &       lai_p, ht_p, phenpfts,
#endif
#ifdef OUTPUT_FORCING
     &       TairC, TcanopyC, Qf, P_mbar, Ca, Ch, U, IPARdir, IPARdif,
     &       soilt, soilm, soilp, soili,
#endif
     &       laipfts, htpfts, gpppfts, ipppfts, mtppfts, fvpfts,
     &       clivepfts, cdeadpfts, R_auto, coszen, o3s, lai0, fv0)


      implicit none
      include 'netcdf.inc'

      integer, intent(in) :: I0, I1, J0, J1, nm, nd, nh
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: GCANOPY
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: Ci
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: NPP
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: GPP
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: GPP0
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: CO2flux
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: IPP
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: MTP
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: resp_r
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: resp_l
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: resp_w
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: c_tot
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: c_soil
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: FO3
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: dFO3
#endif
#ifdef ACTIVE_GROWTH
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: lai_p
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: ht_p
      real*8, dimension(I0:I1,J0:J1,16), intent(in):: phenpfts
#endif
#ifdef OUTPUT_FORCING
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: TairC
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: TcanopyC
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: Qf
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: P_mbar
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: Ca
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: Ch
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: U
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: IPARdir
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: IPARdif
      real*8, dimension(I0:I1,J0:J1,10), intent(in) :: soilt
      real*8, dimension(I0:I1,J0:J1,10), intent(in) :: soilm
      real*8, dimension(I0:I1,J0:J1,10), intent(in) :: soilp
      real*8, dimension(I0:I1,J0:J1,10), intent(in) :: soili
#endif

      real*8, dimension(I0:I1,J0:J1,16), intent(in):: laipfts
      real*8, dimension(I0:I1,J0:J1,16), intent(in):: htpfts
      real*8, dimension(I0:I1,J0:J1,16), intent(in):: gpppfts
      real*8, dimension(I0:I1,J0:J1,16), intent(in):: ipppfts
      real*8, dimension(I0:I1,J0:J1,16), intent(in):: mtppfts
      real*8, dimension(I0:I1,J0:J1,16), intent(in):: fvpfts
      real*8, dimension(I0:I1,J0:J1,16), intent(in):: clivepfts
      real*8, dimension(I0:I1,J0:J1,16), intent(in):: cdeadpfts
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: R_auto
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: coszen
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: o3s
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: lai0
      real*8, dimension(I0:I1,J0:J1), intent(in)   :: fv0

      integer, parameter :: mdays(12) = (/
     &    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      integer totdays, n, nstat

      If (start(3) .gt. nrecmax) then
         print*, 'Record length exceed limit, increase nrecmax!'
         stop
      Endif
      ! print*,'here  nm, nd, nh :', nm, nd, nh
      totdays = 0
      do n = 1, nm
         totdays = totdays + mdays(n)
      enddo
      totdays = totdays - mdays(nm) + nd
#ifndef OUTPUT_DAILY
      nhour(start(3)) = (totdays-1)*24 + nh
      !print*,'here start(3), nhour',start(3),nhour(start(3))
#else
      nday(start(3))  = totdays
#endif

       !  print*, 'here start writing output : ',ncyid
      nstat = NF_PUT_VARA_REAL(ncyid,gcyid,start,
     &  ncount,real(gcanopy))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
        ! print*, 'here start writing  output ci : ',start,ncount,ciyid
      nstat = NF_PUT_VARA_REAL(ncyid,ciyid,start,ncount,real(ci))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,nppyid,start,ncount,real(NPP))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,gppyid,start,ncount,real(GPP))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,gpp0yid,start,ncount,real(GPP0))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,cflxyid,start,ncount,real(CO2flux))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,o3yid,start,ncount,real(o3s))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,royid,start,ncount,real(R_auto))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,respryid,start,ncount,real(resp_r))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,resplyid,start,ncount,real(resp_l))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,respwyid,start,ncount,real(resp_w))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,ctotyid,start,ncount,real(c_tot))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,csoilyid,start,ncount,real(c_soil))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,coszyid,start,ncount,real(coszen))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,laiyid,start,ncount,real(lai0))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,fvyid,start,ncount,real(fv0))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
#ifdef PS_BVOC
      nstat = NF_PUT_VARA_REAL(ncyid,ippyid,start,ncount,real(IPP))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,mtpyid,start,ncount,real(MTP))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,ippsyid,start3,ncut3,real(ipppfts))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,mtpsyid,start3,ncut3,real(mtppfts))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
#endif
#if (defined O3DEP_UPTAKE) || (defined O3DEP_UPTAKE_OFFLINE)
      nstat = NF_PUT_VARA_REAL(ncyid,fo3yid,start,ncount,real(FO3))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,dfo3yid,start,ncount,real(dFO3))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
#endif
#ifdef ACTIVE_GROWTH
      nstat = NF_PUT_VARA_REAL(ncyid,laipyid,start,ncount,real(lai_p))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,htpyid,start,ncount,real(ht_p))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,phenyid,start3,ncut3,
     &   real(phenpfts))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
#endif
      nstat = NF_PUT_VARA_REAL(ncyid,laisyid,start3,ncut3,
     &   real(laipfts))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,htsyid,start3,ncut3,real(htpfts))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,clsyid,start3,ncut3,
     &         real(clivepfts))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,cdsyid,start3,ncut3,
     &           real(cdeadpfts))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,gppsyid,start3,ncut3,real(gpppfts))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,fvsyid,start3,ncut3,real(fvpfts))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      start3(4) = start3(4) + 1
#ifdef OUTPUT_FORCING
      nstat = NF_PUT_VARA_REAL(ncyid,tayid,start,ncount,real(TairC))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,tcyid,start,ncount,real(TcanopyC))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,qcyid,start,ncount,real(Qf*1000.))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,psyid,start,ncount,real(P_mbar))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,cayid,start,ncount,real(Ca))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,chyid,start,ncount,real(Ch))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,wdyid,start,ncount,real(U))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,pdiryid,start,ncount,real(IPARdir))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,pdifyid,start,ncount,real(IPARdif))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,styid,start2,ncut2,real(soilt))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,smyid,start2,ncut2,real(soilm))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,spyid,start2,ncut2,real(soilp))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat = NF_PUT_VARA_REAL(ncyid,siyid,start2,ncut2,real(soili))
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      start2(4) = start2(4) + 1
#endif

#ifndef OUTPUT_DAILY
      !print*,'output nhour:',nhour(start(3))
      nstat = NF_PUT_VAR_INT(ncyid,vtimyid,nhour(1:start(3)))
#else
      nstat = NF_PUT_VAR_INT(ncyid,vtimyid,nday(1:start(3)))
#endif
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      start(3) = start(3) + 1

      end subroutine write_output_netcdf

      end module flux_com
