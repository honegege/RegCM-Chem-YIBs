
      module geos_com
      use mod_dynparam

      implicit none
      save

      integer, parameter :: nlon=360
      integer, parameter :: nlat=180
      integer, parameter :: nymax = 1000
      integer :: myear(nymax)
      integer :: nco2y
      real*8  :: co2c(nymax)
      real*8, allocatable, dimension(:,:) :: elev
      
      end module geos_com


      subroutine read_co2
       
      use geos_com
      implicit none
      
      integer ighg, i, j, nyy
      character*80  title
      character*30 filename
      real*8 temp_ghg(6)
      
      filename=trim(dirter)//pthsep//'GHG'
      open(ighg, file=filename)
      do i=1,5; read(ighg,'(a80)') title; enddo
      i = 1
      do while (i .lt. nymax)
         read(ighg,*,end=101) nyy,(temp_ghg(j),j=1,6)
         myear(i) = nyy
         co2c(i)  = temp_ghg(1)
         i = i + 1
      enddo
101   continue
      close(ighg)
      nco2y = i - 1

      end subroutine read_co2


      subroutine get_co2(year, co2v)
      
      use geos_com
      use mod_mppparam
      implicit none
      
      integer year, year2
      real*8  co2v
      
      year2 = Min(Max(year, myear(1)), myear(nco2y))
      co2v  = co2c(year2 - myear(1) + 1)
      if ( myid == italk ) then
      print*, '[CO2] = ', co2v, ' ppm in ', year2
      end if       

      end subroutine get_co2


      subroutine get_elev(IM, JM, I0, I1, J0, J1, long, latg, lms)

      use geos_com
      use mod_regcm_types, only: lm_state
      implicit none
      include 'netcdf.inc'

      integer, intent(in) :: IM, JM, I0, I1, J0, J1
      type(lm_state) , intent(in) :: lms
      real*8, intent(in)  :: long(IM), latg(JM)
      real*4, dimension(nlon,nlat) :: elev4
      real*8, dimension(nlon,nlat) :: elev0
      real*4 lon0(nlon), lat0(nlat)
      real*8 lon(nlon), lat(nlat), dlon, dlat
      integer ncid, nstat, varid, i, j, ii, jj
      integer lenc
      character*120 dire, fele, filename

      !fele = trim(dirter)//pthsep//'etopo.1b1.nc'
     
      !print*, 'Read elevation from '//fele(1:lenc(fele))
      !filename = fele(1:lenc(fele))

      !nstat    = NF_OPEN(filename(1:lenc(filename)),NF_NOWRITE,ncid)
      !If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      !nstat    = NF_INQ_VARID(ncid,'lat',varid)
      !If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      !nstat    = NF_GET_VAR_REAL(ncid,varid,lat0)
      !If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      !nstat    = NF_INQ_VARID(ncid,'lon',varid)
      !If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      !nstat    = NF_GET_VAR_REAl(ncid,varid,lon0)
      !If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      !nstat    = NF_INQ_VARID(ncid,'elev',varid)
      !If(nstat .NE. NF_NOERR) Call handle_err(nstat)
      !nstat    = NF_GET_VAR_REAL(ncid,varid,elev4)
      !If(nstat .NE. NF_NOERR) Call handle_err(nstat)

      !do j = 1, nlat
      !   lat(j) = lat0(nlat-j+1)
      !enddo
      !do i = 1, nlon
      !   lon(i) = lon0(i)
      !enddo
      !do j = 1, nlat
      !do i = 1, nlon
       !  elev0(i,j) = elev4(i,nlat-j+1)
      !enddo
      !enddo
      !dlon = lon(10) - lon(9)
      !dlat = lat(10) - lat(9)

      allocate(elev(I0:I1,J0:J1))

      Do i = I0, I1
      Do j = J0, J1
         elev(i,j)=lms%mvegid(i,j)
      Enddo
      Enddo

      end subroutine get_elev

      SUBROUTINE HANDLE_ERR(NSTAT)
C     *****************
C     *****************

      IMPLICIT NONE
      include 'netcdf.inc'
      INTEGER NSTAT
      IF (NSTAT .NE. NF_NOERR) THEN
        PRINT*, NF_STRERROR(NSTAT)
      STOP
      ENDIF

      RETURN
      END subroutine HANDLE_ERR

      INTEGER FUNCTION LENC(C)
      IMPLICIT NONE
      CHARACTER*(*) C
      INTEGER I
      I=1
  12  IF(C(I:I) .NE. ' ') THEN
      I=I+1
      GOTO 12 
      ENDIF
      LENC=I-1
      RETURN 
      END
