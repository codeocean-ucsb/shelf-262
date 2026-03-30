
      real    :: SizeX,SizeY
      real    :: f0,beta
      real    :: x0,y0,dx,dy
      integer :: i,j
      !real    :: y0_headland

      integer :: n_smooth, iter
      real :: h_smooth(GLOBAL_2D_ARRAY)
      real :: h_slope
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
          !headland system made by making a guassian in the slope in y (along-shore)
	  !y0 is center of domain
	  h_slope = hslope0 + (hslope0 / 2.0) * exp(-((yr(i,j) - y0)**2) / (2.0 * headland_width**2))
	  h(i,j) = shelf_extent - h_slope*xr(i,j) 
	  !if (h(i,j)<0) then
	  !  h(i,j) =0
	  !endif
          if (h(i,j)>max_h_shelf) then
            h(i,j) = max_h_shelf
          endif
        enddo
      enddo
      !Smoothing
      !n_smooth = 3  ! number of smoothing iterations (tweak as needed)
      ! Make a copy of h into h_smooth initially
      !h_smooth = h
      !do iter = 1, n_smooth
      ! Apply smoothing pass
      !	  do j = 1, ny
      !	    do i = 1, nx
!	      h_smooth(i,j) =  (h(i,j) + h(i+1,j) + h(i-1,j) + h(i,j+1) + h(i,j-1) ) / 5.0
!	    enddo
!	  enddo
	  ! Update h with smoothed values for next iteration
!	  h = h_smooth
!	enddo
# ifdef MASKING
      do j=-1,ny+2
        do i=-1,nx+2
           if (h(i,j)<min_h_shelf) then
	     !land
	     h(i,j) = 0.1 !h can't be zero (which it can be with the headland function)
	     rmask(i,j) = 0
	   else
	     !water
             rmask(i,j) = 1
	   endif
        enddo
      enddo
# endif

