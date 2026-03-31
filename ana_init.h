 
      ! Everything after the implicit none
      ! Replace with something less trivial when needed

      integer :: i,j,k


      real :: g_alpha, alpha, b0
      real :: b, Nb, temp0
      real :: rand1, rand2, noise, normal
 
      ! necessary for FIRST_TIME_STEP flag to work, as per set_global_definitions.h
      ! this had no value for analytical examples previously.
#ifdef EXACT_RESTART
      forw_start=ntstart
#endif

      ! Calculate temperature from buoyancy, lin EOS
      g_alpha = g*Tcoef/rho0
      alpha = Tcoef/rho0

    
      temp0 = 1 !bottom temp
      b0 = temp0 * g * alpha

      b = 0 
      do k=1,nz
        do j=-1,ny+2
          do i=-1,nx+2
	     !Constant strat
	     b = b0 + N2back * (z_r(i,j,k) + max_h_shelf)   

             t(i,j,k,1,itemp) = b / (g*alpha)
             t(i,j,k,2,itemp) = t(i,j,k,1,itemp)
#ifdef MASKING
	     t(i,j,k,1,itemp) = t(i,j,k,1,itemp)*rmask(i,j)
             t(i,j,k,2,itemp) = t(i,j,k,1,itemp)*rmask(i,j)
#endif

# ifdef SALINITY
            t(i,j,k,1,isalt)=36.
            t(i,j,k,2,isalt)=t(i,j,k,1,isalt) 
# endif
	 
#ifdef DYE	
	  t(i,j,k,1,iDye)=0.
	  t(i,j,k,2,iDye)=t(i,j,k,1,iDye)
#endif

!#ifdef TRACER_RESTORE
!	    temp_restore(i,j,k) = t(i,j,k,1,itemp)
!#endif


          enddo
        enddo
      enddo

      !Initialize zero velocity and zero free-surface
      do j=-1,ny+2
        do i=-1,nx+2
          ubar(i,j,1)=0.*umask(i,j)
          ubar(i,j,2)=ubar(i,j,1) *umask(i,j)
          vbar(i,j,1)=0 *vmask(i,j)
          vbar(i,j,2)=vbar(i,j,1) * vmask(i,j)
          zeta(i,j,1)=0 * rmask(i,j)
          zeta(i,j,2)=zeta(i,j,1) * rmask(i,j)
          do k=1,nz
            u(i,j,k,1)=0. *umask(i,j)
            u(i,j,k,2)=u(i,j,k,1) * umask(i,j)
            v(i,j,k,1)=0.*vmask(i,j)
            v(i,j,k,2)=v(i,j,k,1)*vmask(i,j)
          enddo
        enddo
      enddo

