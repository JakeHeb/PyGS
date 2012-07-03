
module GaussDensity

 implicit none
 
 !=============================================================================
 !0. Module-Level Variables
 !=============================================================================
 	 !---Fundamental Parameters -----------------------------------------------
 	 	 real(8), parameter		:: pi = 3.141592653589
 	 	 
 	 !---Saved Variables ------------------------------------------------------
 	 	 real(8), allocatable	:: density_contrast(:,:)
 	 	 
	real(8)			::bin						!Bin size
	real(8), allocatable	::coord_x(:),coord_y(:),coord_z(:) 		!Co-ordinates of the density calcs
	integer			::Q,n_bins					!Number of co-ords
	real(8),allocatable	::delta(:)					!Density Contrast
	real(8),allocatable	::xr(:),yr(:),zr(:)				!Spherical Co-ords of density matrix
	real(8)			::c_inc						!Spacing of the density matrix
	real(8)			::mass,r					!Weighted mass at a point
	real(8),allocatable	::DenHist(:)					!Array of summed Densities in distance bins
	integer, allocatable	::NumHist(:)					!Array of Numbers of lattice points in radial bins
	real(8),allocatable	::xc(:),yc(:),zc(:)
	complex(8),allocatable	::delta_kr(:)					!Fourier Array
 
 !=============================================================================
 !1. GaussianLattice - defines a density lattice based on Gaussian kernel
 	 
 	 !This subroutine sets up a rectilinear lattice and calculates the density
 	 !at each node based on the number of galaxies close by, convolved with a 
 	 !Gaussian kernel. The output is the density contrast, which is the ratio
 	 !of the difference of the density at each point and the mean density to the
 	 !mean density.
 	 	 
 	 !This method is slow and probably gives non-robust results. It either needs
 	 !to be extended and refined, or a better method needs to be used.
 !=============================================================================
	subroutine GaussianLattice(pos, grid_res)
		real(8), intent(in)		:: pos(*,*)	!Galaxy positions in spherical coords
		integer, intent(in)		:: grid_res	!The number of grid points on each axis
		
		real(8)					:: dx	!The spatial separation between grid nodes
		
		integer					::i,j,k
		real(8)					::x,y,z	!Co-ordinates of lattice points.
 	 
  !------------------------Calculate Density Contrast--------------------------	 
  	 write(*,*)'  Calculating Density Contrast Matrix...'
  	 write(999,*)'  Calculating Density Contrast Matrix...'
  	 write(*,*)'  Gaussian Width (Sigma):',sigma, 'Mpc/h'
  	 write(999,*)'  Gaussian Width (Sigma):',sigma, 'Mpc/h'
  	  
   !Evaluate the co-ordinates that receive a density value
    !Set up a Cartesian Lattice
    	xmax = maxval(pos(1,:))
  	  	dx = 2.d0*xmax/grid_res
  	  	
  	  	
  	  	do i=1,grid_res
  	  		do j=1,grid_res
  	  			do k = 1,grid_res
  	  				x = -xmax + (i-1)*dx
  	  				y = -xmax + (j-1)*dx
  	  				z = -xmax + (k-1)*dx

  	  	  
   !Evaluate delta at each point
  	  	  	  allocate(delta(Q))

    !Repeat the following for every lattice point	
  	  	  	  do i=1,Q
  	  	  	  
     !Write how much has been done
  	  	  	  if(floor(i/100.)==(i/100.))then
  	  	  	  write(*,*)'  ',100.0*real(i)/real(Q), '% done'
  	  	  	  end if
  	  	  	  
     !Change the centre to the lattice point
  	  	  	  allocate(coord_x(N),coord_y(N),coord_z(N))
  	  	  	  call Cart(z,ra,dec,coord_x,coord_y,coord_z)
  	  	  	  
  	  	  	  coord_x = coord_x - xc(i)
  	  	  	  coord_y = coord_y - yc(i)
  	  	  	  coord_z = coord_z - zc(i)

  	  	  	
  	  	  	  call Spherical(coord_x,coord_y,coord_z,z,ra,dec)
  	  	  	  
     !Evaluate Weighted mass at the lattice point (using Gaussian)
  	  	  	  mass = 0.0d0
  	  	  	  do j=1,N
  	  	  	  mass = mass + wt(i)*Gauss(z(j),sigma)
  	  	  	  end do
  	  	  	  
     !Divide by the volume and then get density contrast.
  	  	  	  delta(i) = (mass/(2*pi)**(1.50d0))/(sigma**3)
  	  	  	  
  	  	  	  open(unit=155, file=trim(PATH)//prefix//'Density.dat')
  	  	  	  write(155,'(es20.8)')delta(i)
  	  	  	  
     !Change Centre back to original frame
  	  	  	  call Cart(z,ra,dec,coord_x,coord_y,coord_z)
  	  	  	  
  	  	  	  coord_x = coord_x + xc(i)
  	  	  	  coord_y = coord_y + yc(i)
  	  	  	  coord_z = coord_z + zc(i)
  	  	  	  
  	  	  	  call Spherical(coord_x,coord_y,coord_z,z,ra,dec)
  	  	  	  deallocate(coord_x,coord_y,coord_z)
  	  	  	  end do
  	  	  	  
  	  	  	  !do i=1,Q
  	  	  	  !delta(i) = (delta(i)-nb)/nb
  	  	  	  !end do
  	 open(unit=142,file=trim(PATH)//prefix//'DensityContrast.dat')
 	 
 	 do i=1,Q
 	 write(142,9)xc(i),yc(i),zc(i),delta(i)
 	 end do
 	 close(142)
 	 close(155)
 
  !------------------------Non-Centred 3D Transform---------------------------
 	 
   if(opt_all/=0)then
     write(*,*)'  Performing Pair-Separation Analysis...'	
     write(999,*)'  Performing Pair-Separation Analysis...'
     
     n_bins = NINT(Dc*(max_z-min_z)/dr)
     
     allocate(DenHist(n_bins))
     allocate(NumHist(n_bins))
     DenHist = 0.0d0
     NumHist = 0
     
     do j=1,Q-1
       write(*,*)100.0d0*real(j)/real(Q-1), ' % done'
       do i_b=1,Q-j
 	 r=Sqrt((xc(j)-xc(mod(j+i_b,Q)))**2+(yc(j)-yc(mod(j+i_b,Q)))**2&
 	 	&+(zc(j)-zc(mod(j+i_b,Q)))**2)
 	 	 
 	 	i = Floor(r/dr)+1
	 	 
 	 	if(i.LT.n_bins)then
 	 		DenHist(i) = DenHist(i)+delta(j)+delta(i_b)
 	 		NumHist(i) = NumHist(i) + 2
 	 	end if
       end do
     end do  
     
     open(unit=81, file=trim(PATH)//'NonCentre\\'//prefix//'DenHist.dat')
     open(unit=82, file=trim(PATH)//'NonCentre\\'//prefix//'NumHist.dat')
     
     do i=1,n_bins
       write(81,14)dr/2.0d0+(i-1)*dr, DenHist(i)
       write(82,17)dr/2.0d0+(i-1)*dr, NumHist(i)
     end do
     
     do i=1,n_bins
       if(NumHist(i)/=0)then
         DenHist(i) = DenHist(i)/NumHist(i)
       end if
     end do
     
     deallocate(NumHist)
     
     open(unit=83, file=trim(PATH)//'NonCentre\\'//prefix//'AvgHist.dat')
     
     do i=1,n_bins
       write(83,14)dr/2.0d0+(i-1)*dr, DenHist(i)
     end do
     
     
     close(81)
     close(82)
     close(83)
  !-----------------------DFT of the Density Histogram------------------------  
    write(*,*)'  Performing DFT on Density Histogram...'
    write(999,*)'  Performing DFT on Density Histogram...'
    allocate(xr(n_bins))
    do i=1,n_bins
      xr(i) = dr/2.0d0+(i-1)*dr
    end do
    
    allocate(delta_kr(n_r))
    call denDFT(DenHist,xr,r_0,end_r,n_r,delta_kr,'3DFT',PATH,prefix)
    deallocate(delta_kr)
    
    write(*,*)'  Performing DFT on Density Contrast Histogram...'
    write(999,*)'  Performing DFT on Density Constrast Histogram...'
    DenHist = (DenHist-nb)/nb
    
    allocate(delta_kr(n_r))
    call denDFT(DenHist,xr,r_0,end_r,n_r,delta_kr,'3DFT\\Contrast',PATH,prefix)
    deallocate(delta_kr)
    
    deallocate(DenHist) 
    end if
  !------------------------1DFT of Densities---------------------------------- 
  	 write(*,*)'  Performing 1DFT on Density Array...'
  	 write(999,*)'  Performing 1DFT on Density Array...'
  	  
  	  allocate(xr(Q),yr(Q),zr(Q))
  	  call Spherical(xc,yc,zc,xr,yr,zr)
  	  deallocate(xc,yc,zc)
  	  allocate(delta_kr(n_r))
  	  call denDFT(delta,xr,r_0,end_r,n_r,delta_kr,'R',PATH,prefix)
  	  deallocate(delta_kr)
  	  
  	  call cpu_time(tf_sect)
  	  write(*,*)'  Section Time =', tf_sect-ti_sect
  	  write(999,*)'  Section Time =', tf_sect-ti_sect
  	  
  !-----------------------Centre Changing ------------------------------------
  	  !allocate(xc(Q),yc(Q),zc(Q))
  	  !n_bins = NINT(Dc*(max_z-min_z)/dr)
!  	  
  	  !do i = -n_change/2,n_change/2
  	    !allocate(DenHist(n_bins),NumHist(n_bins))
  	    !call Cart(xr,yr,zr,xc,yc,zc)
  	    !xc = xc - r_change*i
!  	  
  	    !call Spherical(xc,yc,zc,xr,yr,zr)
!  	  
  	    !do i_b=1,Q
  	    !j = Floor(xr(i_b)/dr)+1
!	 	 
 	 	!if(j.LT.n_bins)then
 	 		!DenHist(i) = DenHist(i)+delta(i_b)
 	 		!NumHist(i) = NumHist(i) + 2
 	 	!end if
 	  !end do
 	  !end do
 	  
  	  
 
  !------------------------2DFT of Densities R,RA----------------------------- 
  	 !write(*,*)'  Performing 2DFT on Density Array...'
  	 !write(999,*)'  Performing 2DFT on Density Array...'
  	 
  	 !allocate(delta_kr2(n_r,n_th))
  	 !call den2DFT(
  	

   
 
  	  	  
  	  	  
  	 deallocate(z,ra,dec,wt)
 	 77 continue 	 
 	 end do
 	 
 	 call cpu_time(tf_all)
 	 write(*,*)'Total CPU Time:', tf_all-ti_all
 	 write(999,*)'Total CPU Time:', tf_all-ti_all
 	 
 	 8 format(es20.8)
 	 9 format(4es20.8)
 	 14 format(2es20.8)
 	 16 format(3es20.8)
 	 17 format (es20.8, i10)
 	 18 format (f18.14)
 	 20 format (5es20.8)
 	 21 format (i10,es20.8,es20.8)
 	 
end program GaussDensity	  	  

