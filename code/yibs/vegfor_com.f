      
      module vegfor_com
      use mod_dynparam 
      use mod_regcm_types  
      implicit none
      save

      integer, parameter :: nlon=144
      integer, parameter :: nlat=90
      integer, parameter :: nim =360
      integer, parameter :: njm =181
      integer, parameter :: mdays(12) = (/
     &    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      integer, PARAMETER :: JMPERY = 12
      integer :: JDmidOfM(0:JMPERY+1) = (
     *     /-15,16,45,75,106,136,167,197,228,259,289,320,350,381/)
      integer, Parameter :: id2id(18) = (/
     &    -1,4,2,6,6,4,10,10,10,12,13,10,15,-1,15,-1,-1,-1/) 
      real*8, parameter :: canopy_ht(16)  = (/
     &    19.0, 19.0, 16.5, 16.5, 19.0, 19.0, 19.0, 16.5, 
     &    1.0, 1.0, 0.8, 1.26, 0.8, 0.8, 0.8, 1.26/)
      integer, Parameter :: clmid(11) = (/
     &    9,13,10,-1,6,4,2,15,12,17,18/) 

      real*8, dimension(nlon,nlat,16,12)  :: vlai
      real*8, dimension(nlon,nlat,18)     :: vfrac
      real*8, dimension(nlon,nlat,10)     :: vfrac0
      real*8, dimension(nlon,nlat,16)     :: vheight
      real*8, dimension(nlon,nlat,16)     :: vlaix
      real*8, dimension(nlon,nlat)        :: vcrop
      real*8, dimension(nlon,nlat,18)     :: vfrac_clm
      real*8, dimension(nlon,nlat,16,12)  :: vlai_clm
      real*8, dimension(nlon,nlat,16)     :: vhit_clm
      real*8, dimension(nlon,nlat,12)     :: vlandf
      real*8, dimension(nlon,nlat,2)      :: crop_cal

      real*8, dimension(nim, njm, 16,12)  :: mlai
      real*8, dimension(nim, njm, 18)     :: mfrac
      real*8, dimension(nim, njm, 16)     :: mheight
      real*8, dimension(nim, njm, 16)     :: mlaix
      real*8, dimension(nim, njm)         :: mcrop
      real*8, dimension(nim, njm)         :: mc4crop
      real*8, dimension(nim, njm, 16)     :: mfrac_isl
      real*8, dimension(nim, njm, 18)     :: mfrac_clm
      real*8, dimension(nim, njm, 16,12)  :: mlai_clm
      real*8, dimension(nim, njm, 16)     :: mhit_clm
      real*8, dimension(nim, njm, 12)     :: mlandf
      real*8, dimension(nim, njm, 16)     :: mhit_yibs
      Integer, dimension(nim, njm)        :: mvegid
      Integer, dimension(nim, njm,2)      :: mcrop_cal
      end module vegfor_com


      subroutine getlaix_m16(laimax,hdata,I0,I1,J0,J1,do_clmpft,lm,lms)

      use vegfor_com
      implicit none

      integer, intent(in) :: I0, I1, J0, J1
      logical, intent(in) :: do_clmpft
      type(lm_exchange) , intent(in) :: lm
      type(lm_state) , intent(in) :: lms
      real*8, dimension(18,I0:I1,J0:J1),intent(out) :: laimax
      real*8, dimension(18,I0:I1,J0:J1),intent(out) :: hdata
      real*8 :: long(360), latg(181), dlon, dlat
      integer i,j,k,ii,jj

      dlon = 360.d0/dble(nIM)
      dlat = 180.d0/dble(nJM-1)
      do i = 1,nIM
         long(i) = -180.0d0 + dble(i-1)*dlon + 1.d0/2.d0
      enddo
      do j = 1, nJM
         latg(j) = -90.0d0 + dble(j-1)*dlat
      enddo
      
      laimax(:,:,:) = 0.d0
      hdata(:,:,:)  = 0.d0
      do j = J0, J1
      do i = I0, I1 
      do k = 1, 16
#ifdef ACTIVE_GROWTH

#ifdef RESTART_HEIGHT
         hdata(k,i,j)=mhit_yibs(ii,jj,k)

#else
         !!hdata(k,i,j)  = canopy_ht(k)
         hdata(k,i,j)=lms%vheight(i,j,k)
#endif

#else
         if (do_clmpft) then
           hdata(k,i,j)=lms%clmhit(i,j,k)
         else
            laimax(k,i,j)=lms%vlaix(i,j,k)
            hdata(k,i,j)=lms%vheight(i,j,k)
         endif
#endif
      enddo
      enddo
      enddo

      return
      end subroutine getlaix_m16

  
      subroutine getlai_m16(jday, lai,I0,I1,J0,J1,do_clmpft,lm,lms)

      use vegfor_com
      implicit none

      integer, intent(in)  :: I0, I1, J0, J1
      integer, intent(in)  :: jday
      logical, intent(in)  :: do_clmpft
      type(lm_exchange) , intent(in) :: lm
      type(lm_state) , intent(in) :: lms
      real*8, dimension(16,I0:I1,J0:J1),intent(out) :: lai
      integer, dimension(0:JMperY+1) :: startday
      integer :: jmon, offset, itd, totdays, k, i, j, nm1, nm2, ii, jj
      real*8  :: alpha, beta, tempnm1, tempnm2

      startday(0:JMperY+1)=jdMIDofM(0:JMperY+1)-1

      jmon  = 1
      totdays = 0
      do while (totdays < jday)
         totdays = totdays + mdays(jmon)
         jmon = jmon + 1
      end do
      jmon = jmon - 1

      if(jday < startday(jmon))then
        offset=-1
        nm1 = mod(jmon+10, 12) + 1 
        nm2 = jmon
      else
        offset=0
        nm1 = jmon
        nm2 = mod(jmon,12) + 1
      endif
      itd = startday(jmon+1+offset) - startday(jmon+offset)

      beta  = dble(jday-startday(jmon+offset)) / dble(itd)
      alpha = 1.d0 - beta

      lai(:,:,:) = 0.d0
      do j = J0, J1
      do i = I0, I1
      
      do k = 1, 16
      If (do_clmpft) then
         tempnm1 = lms%clmlai(i, j, k, nm1)
         tempnm2 = lms%clmlai(i, j, k, nm2)
         if(tempnm1.gt.10) tempnm1 = 0
         if(tempnm2.gt.10) tempnm2 = 0
         lai(k,i,j) = alpha*tempnm1
     &              + beta*tempnm2
      else
         tempnm1 = lms%vlai(i, j, k, nm1)
         tempnm2 = lms%vlai(i, j, k, nm2)
         if(tempnm1.gt.10) tempnm1 = 0
         if(tempnm2.gt.10) tempnm2 = 0
         lai(k,i,j) = alpha*tempnm1
     &              + beta*tempnm2
      Endif
      enddo
      enddo
      enddo

      end subroutine getlai_m16


      subroutine get_landf(jday, frac,I0,I1,J0,J1,lm,lms)

      use vegfor_com
      implicit none
      integer, intent(in)  :: I0, I1, J0, J1
      integer, intent(in)  :: jday
      type(lm_exchange) , intent(in) :: lm
      type(lm_state) , intent(in) :: lms
      real*8, dimension(I0:I1,J0:J1),intent(out) :: frac
      integer, dimension(0:JMperY+1) :: startday
      integer :: jmon, offset, itd, totdays, k, i, j, nm1, nm2,ii,jj
      real*8  :: alpha, beta, tempnm1, tempnm2
      real*8 :: long(360), latg(181), dlon, dlat

      dlon = 360.d0/dble(nIM)
      dlat = 180.d0/dble(nJM-1)
      do i = 1,nIM
         long(i) = -180.0d0 + dble(i-1)*dlon + 1.d0/2.d0
      enddo
      do j = 1, nJM
         latg(j) = -90.0d0 + dble(j-1)*dlat
      enddo
      
      startday(0:JMperY+1)=jdMIDofM(0:JMperY+1)-1

      jmon  = 1
      totdays = 0
      do while (totdays < jday)
         totdays = totdays + mdays(jmon)
         jmon = jmon + 1
      end do
      jmon = jmon - 1

      if(jday < startday(jmon))then
        offset=-1
        nm1 = mod(jmon+10, 12) + 1 
        nm2 = jmon
      else
        offset=0
        nm1 = jmon
        nm2 = mod(jmon,12) + 1
      endif
      itd = startday(jmon+1+offset) - startday(jmon+offset)

      beta  = dble(jday-startday(jmon+offset)) / dble(itd)
      alpha = 1.d0 - beta

      frac(:,:) = 0.d0
      do j = J0, J1
      do i = I0, I1
         tempnm1 = lms%vlandf(i, j, nm1)
         tempnm2 = lms%vlandf(i, j, nm2)
         if(tempnm1.gt.1) print*,'landf greater than 1: ',i,j,tempnm1
         if(tempnm2.gt.1) print*,'landf greater than 1: ',tempnm2
         frac(i,j) = alpha*tempnm1 + beta*tempnm2
      enddo
      enddo

      end subroutine get_landf


      subroutine readlai_m16(year,long,latg,IM,JM)
!@sum READDLAI read in leaf area indicies from selected file
      use vegfor_com
      use filemanager, only: openunit,closeunit
      use yibs_const, only: COVER_SAND, COVER_DIRT, N_COVERTYPES
      use mod_mppparam
      implicit none
      include 'netcdf.inc'
      
      integer, intent(in) :: year
      integer, intent(in) :: IM, JM
      real*8, intent(in) :: long(IM), latg(JM)
      integer :: k,iunit, n, i, j, ii, jj
      integer ncid, nstat, varid, lenc
      real*4  :: vlai1(nlon,nlat), vhit1(nim,njm), varc4(nim,njm)
      real*4, dimension(nim, njm, 18)     :: mtmp
      
      real*8 lons(nlon), lats(nlat), lon0, dlon, dlat, s, a
      character*80 :: fbin, head
      character*120:: fname
      character*2  :: cmonth

! read first month's LAI's:
      If (nim .ne. im .or. njm .ne. jm) 
     &   call stop_model("read_veg_forcing: Incorrect dim",255)

      vlai(:,:,:,:) = 0.0d0
      do n = 1,12 
      write(cmonth, '(I2.2)') n
      fbin=trim(dirter)//pthsep//'LAIM16_'//cmonth
      call openunit(trim(fbin),iunit,.true.,.true.)
      do k=1,16
        read(iunit) head, vlai1
        vlai(:,:,k,n) = vlai1
      enddo
      call closeunit(iunit)
      enddo

      vfrac(:,:,:) = 0.0d0
      fbin=trim(dirter)//pthsep//'VEG16'
      call openunit(trim(fbin),iunit,.true.,.true.)
      do k=1,16
        read(iunit) head, vlai1
        vfrac(:,:,k) = vlai1
      enddo
      call closeunit(iunit)
       
      vfrac0(:,:,:) = 0.0d0
      fbin=trim(dirter)//pthsep//'VEG8'
      call openunit(trim(fbin),iunit,.true.,.true.)
      do k=1,10
        read(iunit) head, vlai1
        vfrac0(:,:,k) = vlai1
      enddo
      call closeunit(iunit)
      vfrac(:,:,COVER_DIRT) = vfrac0(:,:,10)
      vfrac(:,:,COVER_SAND) = vfrac0(:,:,1)

! Read soil fraction
      vlandf(:,:,:) = 0.0d0
      fbin=trim(dirter)//pthsep//'LANDF'
      call openunit(trim(fbin),iunit,.true.,.true.)
      do k=1,12
        read(iunit) head, vlai1
        vlandf(:,:,k) = vlai1/100.d0
      enddo
      call closeunit(iunit)

! Read CLM cover and LAI
      fbin=trim(dirter)//pthsep//'CLMVEG'
      call openunit(trim(fbin),iunit,.true.,.true.)
      do k=1,11
        read(iunit) head, vlai1
        If (clmid(k) .gt. 0) vfrac_clm(:,:,clmid(k)) = vlai1
      enddo
      call closeunit(iunit)   

      fbin=trim(dirter)//pthsep//'CLMLAI'
      call openunit(trim(fbin),iunit,.true.,.true.)
      do n = 1,12 
      do k=1,9
        read(iunit) head, vlai1
        If (clmid(k) .gt. 0) vlai_clm(:,:,clmid(k),n) = vlai1
      enddo
      enddo
      call closeunit(iunit)   

      fbin=trim(dirter)//pthsep//'CLMHIT'
      call openunit(trim(fbin),iunit,.true.,.true.)
      do k=1,9
        read(iunit) head, vlai1
        If (clmid(k) .gt. 0) vhit_clm(:,:,clmid(k)) = vlai1
      enddo
      call closeunit(iunit)   

#ifdef RESTART_HEIGHT
      
      fbin=trim(dirter)//pthsep//'YIBSHT'
      call openunit(trim(fbin),iunit,.true.,.true.)
      do k=1,16
        read(iunit) head, vhit1
        mhit_yibs(:,:,k) = vhit1
      enddo
      call closeunit(iunit)   

#endif

      ! make sure that veg fractions are reasonable
      do j=1,nlat
          if ( myid == italk ) then
            write(*,*) 'xiaoxie lat at : ', j
          end if
        do i=1,nlon
          do k=1,N_COVERTYPES
            ! get rid of unreasonably small fractions
            if ( vfrac(i,j,k) < 1.d-4 ) vfrac(i,j,k) = 0.d0
          enddo
          s = sum( vfrac(i,j,:) )
          if ( s > .9d0 ) then
            vfrac(i,j,:) = vfrac(i,j,:)/s
          else if ( s < .1d0 ) then 
            print *, "missing veg data at ",i,j,"assume bare soil"
            vfrac(i,j,: ) = 0.d0
            vfrac(i,j,COVER_SAND) = 1.d0
          else
            print *,'modify fraction at ',i,j,s
            print *, vfrac(i,j,:)
            vfrac(i,j,COVER_SAND) = vfrac(i,j,COVER_SAND)+1.d0-s
          endif
        enddo
      enddo
 
      vheight(:,:,:) = 0.0d0
      fbin=trim(dirter)//pthsep//'VHT16'
      print *, "reading VHT16 file from : ", trim(fbin)
      call openunit(trim(fbin),iunit,.true.,.true.)
      do k=1,16
        read(iunit) head, vlai1
        vheight(:,:,k) = vlai1
      enddo
      call closeunit(iunit)
       
      vlaix(:,:,:) = 0.0d0
      fbin=trim(dirter)//pthsep//'VLAIX16'
      print *, "reading VLAIX16 file from : ", trim(fbin)
      call openunit(trim(fbin),iunit,.true.,.true.)
      do k=1,16
        read(iunit) head, vlai1
        vlaix(:,:,k) = vlai1
      enddo
      call closeunit(iunit)

      fbin=trim(dirter)//pthsep//'CROPS_TYPE'
      nstat    = NF_OPEN(trim(fbin),NF_NOWRITE,ncid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat    = NF_INQ_VARID(ncid,'C4_annuals',varid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat    = NF_GET_VAR_REAL(ncid,varid,varc4)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat=NF_CLOSE(ncid)
      IF(nstat .NE. NF_NOERR) Call handle_err(nstat)
      mc4crop = varc4*1.0d0

      fbin=trim(dirter)//pthsep//'ISLSCP'
      nstat    = NF_OPEN(trim(fbin),NF_NOWRITE,ncid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat    = NF_INQ_VARID(ncid,'type',varid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat    = NF_GET_VAR_INT(ncid,varid,mvegid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat    = NF_INQ_VARID(ncid,'cover',varid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat    = NF_GET_VAR_REAL(ncid,varid,mtmp)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat=NF_CLOSE(ncid)
      IF(nstat .NE. NF_NOERR) Call handle_err(nstat)
      mfrac_isl(:,:,:) = 0.d0
      do j = 1, njm
      do i = 1, nim
         a = 0.d0
         do k = 1, 18
            a = a + dble(mtmp(i,j,k))/100.d0
         enddo
         mtmp(i,j,:) = mtmp(i,j,:)/a
         do k = 1, 18
            if (id2id(k) .gt. 0) 
     &          mfrac_isl(i,j,id2id(k)) = mfrac_isl(i,j,id2id(k)) 
     &                                  + dble(mtmp(i,j,k))/100.d0
         enddo
      enddo
      enddo

      vcrop(:,:) = 0.0d0
      call readcrop_m16(year,vcrop)
 
      do i = 1, nlon
         lons(i) = -178.75d0 + dble(i-1)*2.5d0
      enddo
      do j = 1, nlat
         lats(j) = -89.0d0 + dble(j-1)*2.0d0
      enddo

      dlon = lons(10) - lons(9)
      dlat = lats(10) - lats(9)
      do j = 1, njm
      do i = 1, nim
         lon0 = long(i)
         If (lon0 .gt. 180.0d0) lon0 = lon0 - 360.0d0
         ii  = nint((lon0-lons(1))/dlon)+1
         jj  = nint((latg(j)-lats(1))/dlat)+1
         ii  = max( 1, min( ii, nlon) )
         jj  = max( 1, min( jj, nlat) )
         mcrop(i,j) = vcrop(ii,jj)
         mcrop_cal(i,j,:) = crop_cal(ii,jj,:)
         do k = 1, 18
            mfrac(i,j,k)     = vfrac(ii,jj,k)
            mfrac_clm(i,j,k) = vfrac_clm(ii,jj,k)
         enddo
         do k = 1, 12
            mlandf(i,j,k)    = vlandf(ii,jj,k)
         enddo
         do k = 1, 16
            mheight(i,j,k)  = vheight(ii,jj,k)
            mhit_clm(i,j,k) = vhit_clm(ii,jj,k)
            mlaix(i,j,k)    = vlaix(ii,jj,k)
            do n = 1, 12
               mlai(i,j,k,n)     = vlai(ii,jj,k,n)
               mlai_clm(i,j,k,n) = vlai_clm(ii,jj,k,n)
            enddo
         enddo
         mvegid(i,j) = id2id(mvegid(i,j)+1)
      enddo
      enddo

      return
      end subroutine readlai_m16


      subroutine readcrop_m16(year,cropdata)
!@sum Read in crop cover fractions from file.
!@+   Calculates crop fraction for given year.
      use vegfor_com
      use FILEMANAGER, only : openunit,closeunit,nameunit
      implicit none
      include 'netcdf.inc'

      integer, intent(in) :: year
      real*8, intent(out) :: cropdata(nlon,nlat)
      integer i, nstat, ncid, id1, id2
      !----------
      integer :: iu_CROPS
      integer :: year1, year2
      real*4 crop4(nlon,nlat)
      real*8 wt, crop1(nlon,nlat), crop2(nlon,nlat)
      real*4 date1(nlon,nlat), date2(nlon,nlat)
      character*80 title,fbin
      
      !* Calculate fraction for given gcmtime:  interpolate between years*/
        
      year1 = -32768 ; crop1(:,:) = 0.d0
      year2 = -32767 ; crop2(:,:) = 0.d0
      wt = 1.d0
          
      fbin=trim(dirter)//pthsep//'VCROPS'
      call openunit(trim(fbin),iu_CROPS,.true.,.true.)
      do while( year2 < year )
        year1 = year2
        crop1(:,:) = crop2(:,:)
        read (iu_CROPS,end=10) title , crop4
        read(title,*) year2 !Read year integer out of character array title
        crop2 = crop4
      enddo
      wt = (year-year1)/(real(year2-year1,kind=8))
 10   continue
      call closeunit(iu_CROPS)

      cropdata(:,:) = max(0.d0, crop1(:,:)
     &     + wt * (crop2(:,:) - crop1(:,:))) 

      fbin=trim(dirter)//pthsep//'CROPS_CAL'
      nstat=NF_OPEN(trim(fbin),NCNOWRIT,ncid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat=NF_INQ_VARID(ncid,'crop_plant_day',id1)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat=NF_GET_VAR_REAL(ncid,id1,date1)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat=NF_INQ_VARID(ncid,'crop_harvest_day',id2)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat=NF_GET_VAR_REAL(ncid,id2,date2)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      nstat=NF_CLOSE(ncid)
      If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      crop_cal(:,:,1)=date1
      crop_cal(:,:,2)=date2

      return
      end subroutine readcrop_m16

