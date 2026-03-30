
      real    :: SizeX,SizeY
      real    :: f0,beta
      real    :: x0,y0,dx,dy
      integer :: i,j
      real :: h_slope = hslope0

      ! Additional vars that are needed
      !real :: max_depth
      !real :: min_depth
      !real :: slope
      !real :: max_h_shelf

      !bathymetry params
      !max_depth = 50
      !min_depth= 2     
      !slope=0.015 

      !tuned so you get 1 m less than  min depth nearshore (which will be masked)
      !max_h_shelf = (min_h_shelf-1) + LLm*h_slope*dx_grd
      SizeX = LLm*dx_grd   !! Domain size in x-direction [m]
      SizeY =  MMm*dy_grd  !! Domain size in y-direction [m]

      f0 = f0_grd;
      beta = 0

      dx = SizeX/gnx   !! grid size in x-direction
      dy = SizeY/gny

      
# ifdef MPI
      x0=dx*dble(iSW_corn)             ! Coordinates of south-west
      y0=dy*dble(jSW_corn)             ! corner of MPI subdomain
# else
      x0=0. ; y0=0.
# endif

      do j=-1,ny+2          ! Extended ranges for x,y arrays
        do i=-1,nx+2
          xr(i,j)=x0+dx*(dble(i)-0.5D0)
          yr(i,j)=y0+dy*(dble(j)-0.5D0)

          pm(i,j)=1./dx
          pn(i,j)=1./dy
        enddo
      enddo


      x0=SizeX/2.   ! Define center of the domain
      y0=SizeY/2.
      do j=-1,ny+2          ! Extended ranges for x,y arrays
        do i=-1,nx+2
          f(i,j)=f0+beta*( yr(i,j)-y0 )
# if defined NONTRAD_COR
!         feta(i,j) = f0*cos(pi/4)
!         fxi(i,j)  = f0*sin(pi/4)
# endif
        enddo
      enddo


      do j=-1,ny+2
        do i=-1,nx+2
          !Linear slope 
	  h(i,j) = max_h_shelf - h_slope*xr(i,j) 

	  !Enforce a min and max depth
          if (h(i,j)<min_h_shelf) then
	      h(i,j)=min_h_shelf
	  endif 
	  if (h(i,j)>max_h_shelf) then
	     h(i,j) = max_h_shelf
	  endif
        enddo
      enddo

# ifdef MASKING
      do j=-1,ny+2
        do i=-1,nx+2
           if (h(i,j)<=min_h_shelf) then
	     !land
	     rmask(i,j) = 0
	   else
	     !water
             rmask(i,j) = 1
	   endif
        enddo
      enddo
# endif
