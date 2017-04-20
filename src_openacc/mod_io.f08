!====================================
! PKUGCM MODULE "io"
! -----------------------------------
! Open, collect, and write outputs
! Xinyu Wen, Peking Univ, Apr/13/2017
!====================================

module io

   implicit none

   integer, parameter :: Nvar    = 8

   integer, dimension(Nvar)           :: lev = (/ 1   , 10  , 10  , 10  , 10  , 10  ,  10 ,  11 /)
   integer, dimension(Nvar)           :: fid = (/ 901 , 902 , 903 , 904 , 905 , 906 , 907 , 908 /)
   character(len=3),  dimension(Nvar) :: vid = (/"ps ","u  ","v  ","div","vor","t  ","z  ","zh "/)
   character(len=30), dimension(Nvar) :: desc= (/"Surface Pressure (mb)         ", &
                                                 "U Wind (m/s)                  ", &
                                                 "V Wind (m/s)                  ", &
                                                 "Div (m/s^2)                   ", &
                                                 "Vor (m/s^2)                   ", &
                                                 "Temperature (K)               ", &
                                                 "Geopotential (m2/s2)          ", &
                                                 "Geopotential at half (m2/s2)  " /)
   integer(kind=8) :: writecount

end module io


!=====================
! io_write_output
!=====================

subroutine io_write_output
   use io
   use pumamod
   implicit none

   ! local
   real, dimension(NLON,NLAT,NLEV)  :: x3d, y3d
   real, dimension(NLON,NLAT,NLEV+1):: x3dh,y3dh
   real, dimension(NLON,NLAT)       :: x2d, y2d
   real, dimension(NLAT)            :: scal
   real                             :: offs
   real, dimension(NLEV+1)          :: sigmahalf
   real, dimension(NLEV)            :: alphar
   integer :: jlon,jlat,jlev



   !--- Surface Pressure ---
   call sp2fc(sp,x2d)
   call fc2gp(x2d,NLON,NLAT)
   call alt2reg(x2d,1)
   x2d = exp(x2d)*psurf    ! psurf=1011mb
   
   write(fid(1)) x2d

   !---  U and V ---
   x3d   = reshape(gu(:,:),(/NLON,NLAT,NLEV/))     ! gu(NLONxNLAT,NLEV)
   y3d   = reshape(gv(:,:),(/NLON,NLAT,NLEV/))     ! gv(NLONxNLAT,NLEV)
   scal  = cv/sqrt(csq(:))
   offs  = 0.0

   do jlev = 1, NLEV
      do jlon = 1, NLON
         x3d(jlon,:,jlev) = x3d(jlon,:,jlev) * scal(:) + offs
         y3d(jlon,:,jlev) = y3d(jlon,:,jlev) * scal(:) + offs
      end do
   end do

   call alt2reg(x3d(:,:,:),NLEV)
   call alt2reg(y3d(:,:,:),NLEV)

   write(fid(2)) x3d
   write(fid(3)) y3d

   !--- Div and Vor ---
   do jlev = 1, NLEV
      ! div
      call sp2fc(sd(:,jlev),x3d(:,:,jlev))
      call fc2gp(x3d(:,:,jlev),NLON,NLAT)
      ! vor
      call sp2fc(sz(:,jlev),y3d(:,:,jlev))
      call fc2gp(y3d(:,:,jlev),NLON,NLAT)
   end do

   offs = 0.0
   x3d = x3d * ww_scale + offs
   y3d = y3d * ww_scale + offs
         
   call alt2reg(x3d,NLEV)
   call alt2reg(y3d,NLEV)

   write(fid(4)) x3d
   write(fid(5)) y3d

   !--- Temperature ---
   do jlev = 1, NLEV
      call sp2fc(st(:,jlev),x3d(:,:,jlev))
      call fc2gp(x3d(:,:,jlev),NLON,NLAT)
      x3d(:,:,jlev) = x3d(:,:,jlev) * ct + t0(jlev)*ct
   end do

   call alt2reg(x3d,NLEV)

   write(fid(6)) x3d

   !--- Geopotential at half ---
   sigmahalf(1)         = 0.0
   sigmahalf(2:NLEV+1)  = sigmh
   alphar(1:NLEV)       = 1.0-sigmahalf(1:NLEV)/dsigma(1:NLEV)*    &
                          log(sigmahalf(2:NLEV+1)/sigmahalf(1:NLEV)) ! Note: 最顶层sigmahalf==0 在log中做分母会被除零

   x3dh(:,:,NLEV+1) = 0.0       ! 这是地表位势 以后需要求精
   do jlev = NLEV, 1, -1
      x3dh(:,:,jlev) = x3dh(:,:,jlev+1)+dsigma(jlev)*x3d(:,:,jlev)   ! zh
      y3d (:,:,jlev) = x3dh(:,:,jlev+1)+alphar(jlev)*x3d(:,:,jlev)   ! z
   end do

   write(fid(7)) y3d
   write(fid(8)) x3dh

   !--- count ---
   writecount = writecount + 1

end subroutine


!=====================
! io_open_output
!=====================

subroutine io_open_output
   use io
   implicit none

   ! local
   integer :: i,f
   character(len=50) :: fn

   do i = 1, Nvar
      f  = fid(i)
      fn = "output_"//trim(vid(i))//".bin"
      open(unit=f, file=trim(fn), access="stream", form="unformatted", status="replace")
   end do

   writecount = 0
end subroutine io_open_output


!=====================
! io_close_output
!=====================

subroutine io_close_output
   use io
   use pumamod, only: NLON,NLAT,NLEV
   implicit none

   ! local
   integer :: i,f
   character(len=50) :: fn
   real :: dd, halfdd

   dd     = 360.0/NLON
   halfdd = dd/2.0

   do i = 1, Nvar
      f = fid(i)
      ! close binary data file
      close(f)
      ! creat GrADS header file
      fn = "output_"//vid(i)
      open(unit=f, file=trim(fn)//".ctl", access="sequential", form="formatted", status="replace")
         write(f,*) "DSET  ^"//trim(fn)//".bin"
         write(f,*) "TITLE "//trim(vid(i))
         write(f,*) "OPTIONS yrev zrev"
         write(f,*) "UNDEF -999"
         write(f,*) "XDEF  ",NLON," LINEAR ",halfdd," ",dd
         write(f,*) "YDEF  ",NLAT," LINEAR ",(-90+halfdd)," ",dd
         if (vid(i)=="ps  ") then
            write(f,*) "ZDEF  1 LEVELS 1000"
         else
            !write(f,*) "ZDEF ",NLEV," LEVELS 1000 900 800 700 600 500 400 300 200 100"
            write(f,*) "ZDEF ",NLEV," LEVELS 950 850 750 650 550 450 350 250 150 50"
         end if
         write(f,*) "TDEF ",writecount," LINEAR 1jan2000 1dy"
         write(f,*) "VARS   1" 
         write(f,*) trim(vid(i))//"   ", lev(i), "   99   "//trim(desc(i))
         write(f,*) "ENDVARS"
      close(fid(i))
   end do
end subroutine io_close_output

