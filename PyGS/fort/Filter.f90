module Filters
 implicit none


 !============================================================================
 !0. Module Level Variables
 !============================================================================ 
 	 complex(8), parameter	:: ci = (0.0d0,1.d0)	! Complex i
 	 integer, allocatable	:: remains(:)	!Remaining galaxies after filtering

 !============================================================================
 !1. 1D Filter - filters galaxies based on resonance with a chosen 1D 
 	 			 !Fourier peak
 !============================================================================
 	 subroutine filter_1d(quantity,kmin,kmax,n_k,min_phase)
 	 	real(8), intent(in)	:: quantity(*)
 	 	real(8), intent(in)	:: kmin, kmax 	!Minimum and maximum wavenumbers in peak.
 	 	integer, intent(in)	:: n_k			!Number of wavenumbers to integrate.
 	 	real(8), intent(in)	:: min_phase	!The phase cut-off for acceptance/rejection
 	 	
 	 	
 	 	complex(8), allocatable	:: fil(:)	!Filter phases.
 	 	real(8)					:: inc		!Increments between wavenumbers.
 	 	integer					:: nn 		!Number of input galaxies.
 	 	real(8), allocatable	:: filr(:)	!Real part of filter phases.
 	 	integer					:: i,ii, TN
 	 	nn = size(quantity)
 	 	inc = (kmax-kmin)/(n_k-1)
 	 
 	 	allocate(fil(nn))
 	 
 	 	fil(:) = (0.0d0,0.0d0)
 	 
 	 	!Calculate phase values for all galaxies.
 	 	do ii = 1, nn
 	 	 do i=1,n_k
 	 	 	k= 2.0d0*pi/(kmin+(i-1)*inc)
 	 		fil(ii)=fil(ii)+Exp(ci*k*dat(ii))
 	 	 end do
 	 	end do
 	 
 	 	allocate(filr(nn))
 	 	filr = real(fil)/n_k
 	 	deallocate(fil)

 	 	!Determine which galaxies are 'in phase' with the peak.
   		TN=0
   		do i=1,nn
   		 if(filr(i).GT.min_phase)then
   		 	 Tn = TN + 1
   		 end if
   		end do
   		allocate(remains(TN))
   		
   		TN=0
   		do i=1,nn
   		 if(filr(i).GT.min_phase)then
   		 	 TN = TN + 1
   		 	 remains(TN) = i-1
   		 end if
   		end do
   	end subroutine
 
    
 !============================================================================
 !2. 3D Filter
 !============================================================================
 	 subroutine filter_3d(quantity,kmin,kmax,n_k,min_phase)
 	 	real(8), intent(in)	:: quantity(3,*)
 	 	real(8), intent(in)	:: kmin(3), kmax(3) 	!Minimum and maximum wavenumbers in peak.
 	 	integer, intent(in)	:: n_k(3)				!Number of wavenumbers to integrate.
 	 	real(8), intent(in)	:: min_phase			!The phase cut-off for acceptance/rejection
 	 	
 	 	
 	 	complex(8), allocatable	:: fil(:)	!Filter phases.
 	 	real(8)					:: inc(3)	!Increments between wavenumbers.
 	 	integer					:: nn 		!Number of input galaxies.
 	 	real(8), allocatable	:: filr(:)	!Real part of filter phases.
 	 	integer					:: i,ii, TN
 	 	real(8)					:: k1,k2,k3
 	 	nn = size(quantity)
 	 	inc = (kmax-kmin)/(n_k-1)
 	 
 	 	allocate(fil(nn))
 	 
 	 	fil(:) = (0.0d0,0.0d0)
 	 
 		nn = size(quantity(1,:))
 		inc(1) = (kmax(1)-kmin(1))/(n_k(1)-1)
 		inc(2) = (kmax(2)-kmin(2))/(n_k(2)-1)
 		inc(3) = (kmax(3)-kmin(3))/(n_k(3)-1)
 		 
 		allocate(fil(nn))
 		
 		fil(:) = (0.0d0,0.0d0)
 		
 		do ii=1,nn
 		 do jj=1,n_k(3)
 		write(*,*)'  ',((real(ii)+(real(jj)/real(n_k3)))/real(nn))*100.,'% done'
 		  k3= 2.0d0*pi/(kmin(3)+(jj-1)*inc(3))
 		  do j =1,n_k(2)
 		   k2= 2.0d0*pi/(kmin(2)+(j-1)*inc(2))
 		   do i=1,n_k(1)
 		 	k1= 2.0d0*pi/(kmin(1)+(i-1)*inc(1))
 			fil(ii)=fil(ii)+Exp(ci*(k1*quantity(1,ii)+k2*quantity(2,ii)+&
 			&k3*quantity(3,ii)))
 		   end do
 		  end do
 		 end do
 		end do
 	 
 		allocate(filr(nn))
 		filr = real(fil)/(n_k(1)*n_k(2)*n_k(3))
 		deallocate(fil)
 		 
 	 	!Determine which galaxies are 'in phase' with the peak.
   		TN=0
   		do i=1,nn
   		 if(filr(i).GT.min_phase)then
   		 	 Tn = TN + 1
   		 end if
   		end do
   		allocate(remains(TN))
   		
   		TN=0
   		do i=1,nn
   		 if(filr(i).GT.min_phase)then
   		 	 TN = TN + 1
   		 	 remains(TN) = i-1
   		 end if
   		end do
   	end subroutine
 
 	 
end module
  	  
