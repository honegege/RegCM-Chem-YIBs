
      module gcm_com

      implicit none
      save

      integer, parameter :: nlon=144
      integer, parameter :: nlat=90
      integer, parameter :: nim =144
      integer, parameter :: njm =90
      integer, parameter :: mdays(12) = (/
     &    31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/)
      character*3, parameter :: sdays(12) = (/
     &    'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', 
     &    'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC'/)
      character*8, parameter :: svar1(12) = (/
     &    'VSAT', 'VTCAN', 'VQCAN', 'VPRES', 'VCA', 'VO3', 
     &    'VCH', 'VWIND', 'VPART', 'VPARD', 'VCOSZ', 'VFW'/)
      character*8, parameter :: svar2(4) = (/
     &    'VSOILT', 'VSOILM', 'VSOILP', 'VSOILI'/)

      real*4, allocatable, dimension(:,:,:)     :: var1
      real*8, allocatable, dimension(:,:,:,:)   :: var2d   ! hourly variables
      real*8, allocatable, dimension(:,:,:,:,:) :: var3d   ! hourly variables

      Integer IG0, IG1, JG0, JG1
      real*8 long(nim), latg(njm), lons(nlon), lats(nlat)
      real*8 dlon, dlat
      real*8 Pfrac(nim,njm,12)  
      real*4 par1(nim,njm)
      Integer start(3), ncount(3)
      logical :: ifirst = .true.
      Integer :: nd_old = -1

      end module gcm_com


      subroutine read_gcm(IM,JM,I0,I1,J0,J1,nm,nd,nh, 
     &      tag, tcang, qfg, pmg, cag, chg, wdg, ptg, prg, coszg,
     &      fwg, o3s, soiltg, soilmg, soilpg, soilig,mddom)

      use gcm_com
      use filemanager, only: openunit,closeunit
      use mod_regcm_types, only: domain

      implicit none
      include 'netcdf.inc'
      type(domain) , intent(in):: mddom
      integer, intent(in) :: IM, JM, I0, I1, J0, J1, nm, nd, nh
      real*8, dimension(I0:I1,J0:J1), intent(out) :: tag
      real*8, dimension(I0:I1,J0:J1), intent(out) :: tcang
      real*8, dimension(I0:I1,J0:J1), intent(out) :: qfg
      real*8, dimension(I0:I1,J0:J1), intent(out) :: pmg
      real*8, dimension(I0:I1,J0:J1), intent(out) :: cag
      real*8, dimension(I0:I1,J0:J1), intent(out) :: chg
      real*8, dimension(I0:I1,J0:J1), intent(out) :: wdg
      real*8, dimension(I0:I1,J0:J1), intent(out) :: ptg
      real*8, dimension(I0:I1,J0:J1), intent(out) :: prg
      real*8, dimension(I0:I1,J0:J1), intent(out) :: coszg
      real*8, dimension(I0:I1,J0:J1), intent(out) :: fwg
      real*8, dimension(I0:I1,J0:J1), intent(out) :: o3s
      real*8, dimension(6,I0:I1,J0:J1), intent(out) :: soiltg
      real*8, dimension(6,I0:I1,J0:J1), intent(out) :: soilmg
      real*8, dimension(6,I0:I1,J0:J1), intent(out) :: soilpg
      real*8, dimension(6,I0:I1,J0:J1), intent(out) :: soilig
   
      integer yyncid, nstat, varid, i, j, k, ii, jj
      integer lenc, nv, ni, nn, iunit
      character*120 dirg, fgcm, filename
      character*80  head
      character*8 cvar, cvar2
      character*2 cidx
      real*8 lon0
      real tempv(144,90),tempx
	  
      If (nd_old .ne. nd) then              ! new day
         nd_old = nd
         dirg  = '/store/test/RegCM-4.6.0/'//
     &           'run/input/'
         fgcm  = sdays(nm)//'2000.EVEGSUBDD.nc'

         do nv = 1, 12
         cvar  = svar1(nv)
         print*, 'Read GCMF : ',dirg(1:lenc(dirg))//cvar(1:lenc(cvar))//
     &        fgcm(1:lenc(fgcm))
         filename = dirg(1:lenc(dirg))//cvar(1:lenc(cvar))
     &            //fgcm(1:lenc(fgcm))
         nstat  = NF_OPEN(filename(1:lenc(filename)),NF_NOWRITE,yyncid)
          If(nstat .NE. NF_NOERR) Call handle_err(nstat)
        
         If (ifirst) then
         nstat    = NF_INQ_VARID(yyncid,'lat',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VAR_DOUBLE(yyncid,varid,lats)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_INQ_VARID(yyncid,'lon',varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat    = NF_GET_VAR_DOUBLE(yyncid,varid,lons)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         dlon = lons(10) - lons(9)
         dlat = lats(10) - lats(9)
         lats(1) = -89.0d0
         lats(nlat) = 89.0d0

         do i = 1, 144
            long(i) = -180.0d0 + dble(i-1)*2.5d0 + 1.25
         enddo 
         do j = 1, 90
            latg(j) = -90.0d0 + dble(j-1)*2.5d0
         enddo

         lon0 = long(I0)
         If (lon0 .gt. 180.0d0) lon0 = lon0 - 360.0d0
         IG0  = nint((lon0-lons(1))/dlon)+1
         JG0  = nint((latg(J0)-lats(1))/dlat)+1
         IG0  = IG0 -1 
         JG0  = JG0 -1
         IG0  = max( 1, min( IG0, nlon) )
         JG0  = max( 1, min( JG0, nlat) )
         lon0 = long(I1)
         If (lon0 .gt. 180.0d0) lon0 = lon0 - 360.0d0
         IG1  = nint((lon0-lons(1))/dlon)+1
         JG1  = nint((latg(J1)-lats(1))/dlat)+1
         IG1  = IG1 +1 
         JG1  = JG1 +1
         IG1  = max( 1, min( IG1, nlon) )
         JG1  = max( 1, min( JG1, nlat) )

         allocate( var1(144,90,24) )
         allocate( var2d(144,90,24,12) )
         allocate( var3d(144,90,24,6,4) )

         start(1)   = 1
         start(2)   = 1
         ncount(1)  = 144
         ncount(2)  = 90
         ncount(3)  = 24

!        PAR dir-to-total ratio
         Pfrac(:,:,:) = 0.0d0
         call openunit('/store/test/RegCM-4.6.0/run/sim2015/input/
     &PFRAC',iunit,.true.,.true.)
         do k=1,12
           read(iunit) head, par1
           Pfrac(:,:,k) = par1
         enddo
         call closeunit(iunit)

         ifirst     = .false.

         Endif

         start(3)   = (nd-1)*24+1

         nstat = NF_INQ_VARID(yyncid,cvar(1:lenc(cvar)),varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat = NF_GET_VARA_REAL(yyncid,varid,start,ncount,var1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat=NF_CLOSE(yyncid)
         IF(nstat .NE. NF_NOERR) Call handle_err(nstat)
         var2d(:,:,:,nv) = var1

         enddo

         do nv = 1, 4
         cvar  = svar2(nv)
         do ni = 1, 6
         write(cidx, '(i1)') ni
         cvar2 = cvar(1:lenc(cvar))//cidx(1:lenc(cidx))
         print*,'Read GCMF : '//cvar2(1:lenc(cvar2))//fgcm(1:lenc(fgcm))
         filename = dirg(1:lenc(dirg))//cvar2(1:lenc(cvar2))
     &            //fgcm(1:lenc(fgcm))
         nstat  = NF_OPEN(filename(1:lenc(filename)),NF_NOWRITE,yyncid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat = NF_INQ_VARID(yyncid,cvar2(1:lenc(cvar2)),varid)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat = NF_GET_VARA_REAL(yyncid,varid,start,ncount,var1)
         If(nstat .NE. NF_NOERR) Call handle_err(nstat)
         nstat=NF_CLOSE(yyncid)
         IF(nstat .NE. NF_NOERR) Call handle_err(nstat)
         var3d(:,:,:,ni,nv) = var1
         enddo   
         enddo
      
      Endif


      Do i = I0, I1
      Do j = J0, J1
      ii  = nint((mddom%xlon(i,j)-long(1))/dlon)+1
      jj  = nint((mddom%xlat(i,j)-latg(1))/dlat)+1
      ii  = max( 1, min( ii, nim) )
      jj  = max( 1, min( jj, njm) )
         nn  = 1
         tag(i,j)=var2d(ii,jj,nh+1,nn)
         nn  = nn + 1
         tcang(i,j)=var2d(ii,jj,nh+1,nn)
         nn  = nn + 1
         qfg(i,j)=var2d(ii,jj,nh+1,nn)
         nn  = nn + 1
	     pmg(i,j)=var2d(ii,jj,nh+1,nn)
         nn  = nn + 1
         cag(i,j)=var2d(ii,jj,nh+1,nn)
         nn  = nn + 1
         o3s(i,j)=var2d(ii,jj,nh+1,nn)
         nn  = nn + 1
         chg(i,j)=var2d(ii,jj,nh+1,nn)
         nn  = nn + 1
         wdg(i,j)=var2d(ii,jj,nh+1,nn)
         nn  = nn + 1
         ptg(i,j)=var2d(ii,jj,nh+1,nn)
         nn  = nn + 1
         prg(i,j)=var2d(ii,jj,nh+1,nn)
         nn  = nn + 1
         coszg(i,j)=var2d(ii,jj,nh+1,nn)
         nn  = nn + 1
         fwg(i,j)=var2d(ii,jj,nh+1,nn)

         do ni = 1, 6
         soiltg(ni,i,j)=var3d(ii,jj,nh+1,ni,1)
         soilmg(ni,i,j)=var3d(ii,jj,nh+1,ni,2)
         soilpg(ni,i,j)=var3d(ii,jj,nh+1,ni,3)
         soilig(ni,i,j)=var3d(ii,jj,nh+1,ni,4)
         enddo
      enddo
      enddo
      print*,'there',long,latg,nIM,nJM,
     & mddom%xlon(15,15),mddom%xlat(15,15)
      print*,'here',IG0,IG1,JG0,JG1,pmg(15,15),var2d(120,48,nh+1,4)
      end subroutine read_gcm
	  

