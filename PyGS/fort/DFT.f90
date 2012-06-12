module DFT
 implicit none

 !============================================================================
 !MODULE DESCRIPTION
 !============================================================================
 !  This module contains various routines for performing Discrete Fourier 
 !  Transforms (DFT's) for various numbers of dimensions and for various types
 !  of input data. The necessity of the module lies in the fact that generally
 !  faster FFT algorithms (such as FFTW) perform the DFT in frequency space,
 !  which is not so nice for plotting in wavelength space. This module has
 !  SLOWER routines, but they are much more customizable and produce 'nicer'
 !  results.
 !============================================================================
   
 !============================================================================
 !Module Data
 !============================================================================ 
   !Note that the reason for the inclusion of the module data is for ease-of-
   !use from python, while it makes the module slightly more difficult to 
   !understand.
   !
   !The variables are merely all arrays that would be allocatable in the
   !associated subroutines (which cannot be done with f2py). They must be 
   !allocated and filled from python (or however) before the corresponding
   !subroutine is called.
   !  
   !The 'dat' arrays are the values used for the DFT. The 'dist' arrays are 
   !the co-ordinate vectors of the dat arrays. Thus, the dat arrays must have
   !shape of (size(dist1),size(dist2),...). The wavelength arrays contain
   !the values of wavelengths to be calculated in the power spectrum.
     
   real(8), parameter	  ::pi= 3.141592653589  !Value of Pi
   complex(8), parameter::ci= (0.0d0,1.0d0)   !The Complex Unit (i)
 	 
 contains
 
 !============================================================================
 !1.Phase DFT in 1D
 !============================================================================ 
  subroutine PhaseDFT_One(n,nn,dist1,wvl1,ps)

   !Requires dist1 and wvl1 to be defined.
	 integer,intent(in)  ::n        !The size of the wavelength array
	 integer             ::nn       !The size of the data
	 real(8)             ::dist1(nn)!Co-ordinates for 1D transform
	 real(8)             ::wvl1(n)  !Co-ordinates for 1D transform
	 
 	 real(8), intent(out)::ps(n)    !The power-spectrum
 	 
 	 complex(8)   	     ::camp(n)  !The complex amplitudes of the DFT
 	 integer				     ::i, j     !Iterators            
 	 
 	 
     ! Initialize variables.
     camp = (0.0,0.0)
     
     !Perform the grunt-work of the DFT
     do j = 1, n
      do i=1,nn
      camp(j)=camp(j)+Exp(ci*(2.d0*pi/wvl1(j))*dist1(i))
      end do
     end do
     
     ! Normalize the amplitudes.
     camp = (nn**(-.50))*camp
     
     ! Calculate the Power
     do i=1,n
      ps(i) = (abs(camp(i)))**2
     end do

 	 
  end subroutine
 	 
 !============================================================================
 !3. Real data (eg. from un-normalized histogram) DFT 1D
 !============================================================================ 	
  subroutine DFT_One(dat1,nn,dist1,n,wvl1,ps)

   !Requires dist1 and wvl1 to be defined.
	 integer,intent(in)  ::n        !The size of the wavelength array
	 integer,intent(in)  ::nn       !The size of the data
	 real(8),intent(in)  ::dist1(nn)!Co-ordinates for 1D transform
	 real(8),intent(in)  ::wvl1(n)  !Co-ordinates for 1D transform
	 real(8),intent(in)  ::dat1(nn) !Values for the FT 
	 
 	 real(8), intent(out)::ps(n)    !The power-spectrum
 	 
 	 complex(8)   	     ::camp(n)  !The complex amplitudes of the DFT
 	 integer				     ::i, j     !Iterators            
 	 
 	 
     ! Initialize variables.
     camp = (0.0,0.0)
     
     !Perform the grunt-work of the DFT
     do j = 1, n
      do i=1,nn
      camp(j)=camp(j)+dat1(i)*Exp(ci*(2.d0*pi/wvl1(j))*dist1(i))
      end do
     end do
     
     ! Normalize the amplitudes.
     camp = (nn**(-.50))*camp
     
     ! Calculate the Power
     do i=1,n
      ps(i) = (abs(camp(i)))**2
     end do

 	 
  end subroutine
 !============================================================================
 !4. 2DFT Phase
 !============================================================================ 
 	subroutine PhaseDFT_Two(nn,dist1,dist2,no_1,wvl1,no_2,wvl2,ps)
 	 integer, intent(in)  ::nn                          !Number of data-points
 	 integer, intent(in)  ::no_1,no_2                   !Sizes of wavelength arrays 
 	 real(8)              ::dist1(nn),dist2(nn)         !Co-ordinates for 2D transform
 	 real(8)              ::wvl1(no_1),wvl2(no_2)       !Wavelengths for 2D transform
 	 
 	 real(8),  intent(out)::ps(no_1,no_2)   !Power Spectrum
 	 
 	 complex(8)           ::camp(no_1,no_2) !Complex FT amplitudes
 	 integer				      ::i, ii, j,jj                 !Iterators
 	 
 	 real(8)              ::k_1,k_2                     !Wavenumbers
 	 
 	 !Initialize Variables
 	 camp(:,:) = 0.0d0


 	 ! Perform 2DFT
 	 do jj =1,no_2
 	  write(*,*)'  ',(real(j)/real(no_2))*100.,'% done'
 	  k_2 = 2.d0*pi/wvl2(jj)
 	  do j =1,no_1
 	    k_1 = 2.d0*pi/wvl1(jj)
 	   do ii=1,nn
 	    do i = 1,nn
 	      camp(j,jj)=camp(j,jj) + exp(ci*(dist1(i)*k_1+dist2(ii)*k_2))
 	    end do
 	   end do
 	  end do
 	 end do


 	 camp(:,:) = camp(:,:)/nn

 	 ! Calculate the Power
 	 do i=1,no_2
 	  do j=1,no_1
      ps(j,i) = (abs(camp(j,i)))**2
    end do
 	 end do
 	 
 	end subroutine
 	
 !============================================================================
 !5. 2DFT Lattice
 !============================================================================ 
 	subroutine DFT_Two(dat2,nn_1,dist1,nn_2,dist2,no_1,wvl1,no_2,wvl2,ps)
 	 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 	   !This subroutine takes as the input values a 2D array. The co-ordinates
 	   !of each array point are specified by the 1D vectors dist1 and dist2. These
 	   !may be irregular, but must be rectilinear. To input completely irregular
 	   !data, use pointDFT_Two.
 	 integer, intent(in)  ::nn_1,nn_2                   !Number of data-points
 	 integer, intent(in)  ::no_1,no_2                   !Sizes of wavelength arrays 
 	 real(8), intent(in)  ::dist1(nn_1),dist2(nn_2)         !Co-ordinates for 2D transform
 	 real(8), intent(in)  ::wvl1(no_1),wvl2(no_2)       !Wavelengths for 2D transform
 	 real(8), intent(in)  ::dat2(nn_1,nn_2)             !Values for the FT
 	 
 	 real(8),  intent(out)::ps(no_1,no_2)   !Power Spectrum
 	 
 	 complex(8)           ::camp(no_1,no_2) !Complex FT amplitudes
 	 integer				      ::i, ii, j,jj                 !Iterators
 	 
 	 real(8)              ::k_1,k_2                     !Wavenumbers
 	 
 	 !Initialize Variables

 	 camp(:,:) = 0.0d0

 	 ! Perform 2DFT
 	 do jj =1,no_2
 	  write(*,*)'  ',(real(j)/real(no_2))*100.,'% done'
 	  k_2 = 2.d0*pi/wvl2(jj)
 	  do j =1,no_1
 	    k_1 = 2.d0*pi/wvl1(j)
 	   do ii=1,nn_1
 	    do i = 1,nn_2
 	      camp(j,jj)=camp(j,jj) + &
 	             &dat2(ii,i)*exp(ci*(dist1(ii)*k_1+dist2(i)*k_2))
 	    end do
 	   end do
 	  end do
 	 end do


 	 camp(:,:) = ((nn_1*nn_2)**(-.50))*camp(:,:)

 	 ! Calculate the Power
 	 do i=1,no_2
 	  do j=1,no_1
      ps(j,i) = (abs(camp(j,i)))**2
    end do
 	 end do
 	 
 	end subroutine
 	
 !============================================================================
 !6. 3DFT Phase
 !============================================================================ 
 	subroutine PhaseDFT_Three(nn,dist1,dist2,dist3,no_1,wvl1,no_2,wvl2,&
 	                          &no_3,wvl3,ps)
 	 integer, intent(in)  ::nn                          !Number of data-points
 	 integer, intent(in)  ::no_1,no_2,no_3              !Sizes of wavelength arrays 
 	 real(8), intent(in)  ::dist1(nn),dist2(nn),dist3(nn)!Co-ordinates for 2D transform
 	 real(8), intent(in)  ::wvl1(no_1),wvl2(no_2),wvl3(no_3)           !Wavelengths for 2D transform
 	 
 	 real(8),  intent(out)::ps(no_1,no_2,no_3)   !Power Spectrum
 	 
 	 complex(8)           ::camp(no_1,no_2,no_3) !Complex FT amplitudes
 	 integer				      ::i, ii, j,jj,jjj,iii         !Iterators
 	 
 	 real(8)              ::k_1,k_2,k_3                 !Wavenumbers
 	 
 	 !Initialize Variables
 	 camp(:,:,:) = 0.0d0
 	 
 	 ! Perform DFT
 	 do jjj=1,no_3
 	  k_3 = wvl3(jjj)
 	  do jj =1,no_2
      write(*,*)'  ',((real(jjj)+real(jj)/real(no_2))/(real(no_3)))*100.,'%done'
 	   k_2 = wvl2(jj)
 	   do j =1,no_1
 	    k_1 = wvl1(j)
 	    do iii=1,nn
 	      do ii=1,nn
 	        do i=1,nn
             camp(j,jj,jjj)=camp(j,jj,jjj) + &
             & exp(ci*(dist1(i)*k_1+dist2(ii)*k_2+dist3(iii)*k_3))
          end do
        end do
      end do
 	   end do
 	  end do
 	 end do

 	 !Normalize
 	 camp(:,:,:) = (nn**(-1.50))*(camp(:,:,:))
 	 
 	 ! Calculate the Power
 	 do i=1,no_3
 	  do j=1,no_2
 	    do ii = 1,no_1
 	      ps(ii,j,i) = (abs(camp(ii,j,i)))**2
 	    end do
 	  end do
 	 end do
 	end subroutine

 !============================================================================
 !7. 3DFT
 !============================================================================ 
 	subroutine DFT_Three(dat3,nn_1,dist1,nn_2,dist2,nn_3,dist3,&
 	                    &no_1,wvl1,no_2,wvl2,no_3,wvl3,ps)

 	 integer,intent(in)   ::no_1,no_2,no_3 !Size of wavelength arrays
 	 integer,intent(in)   ::nn_1,nn_2,nn_3 !Column length of data    
 	 real(8), intent(in)  ::dist1(nn_1),dist2(nn_2),dist3(nn_3)!Co-ordinates for 2D transform
 	 real(8), intent(in)  ::wvl1(no_1),wvl2(no_2),wvl3(no_3)           !Wavelengths for 2D transform
 	 real(8), intent(in)  ::dat3(nn_1,nn_2,nn_3)  !Values for the FT
 	 
 	 real(8), intent(out)	::ps(no_1,no_2,no_3)    !The Power Spectrum
 	 complex(8)           ::camp(no_1,no_2,no_3)  !The complex amplitudes of the FT
 	 
 	 
 	 integer				::i, ii, j, jj,iii,jjj   !Iterators
 	 
 	 real(8)        ::k_1,k_2,k_3    !Wavenumbers
 	 
 	 !Initialize Variables
 	 camp(:,:,:) = 0.0d0
 	 
 	 ! Perform DFT
 	 do jjj=1,no_3
 	  k_3 = 2.d0*pi/wvl3(jjj)
 	  do jj =1,no_2
      write(*,*)'  ',((real(jjj)+real(jj)/real(no_2))/(real(no_3)))*100.,'%done'
 	   k_2 = 2.d0*pi/wvl2(jj)
 	   do j =1,no_1
 	    k_1 = 2.d0*pi/wvl1(j)
 	    do iii=1,nn_3
 	      do ii=1,nn_2
 	        do i=1,nn_1
             camp(j,jj,jjj)=camp(j,jj,jjj) + dat3(i,ii,iii)*&
             &exp(ci*(dist1(i)*k_1+dist2(ii)*k_2+dist3(iii)*k_3))
          end do
        end do
      end do
 	   end do
 	  end do
 	 end do

 	 !Normalize
 	 camp(:,:,:) = ((nn_1*nn_2*nn_3)**(-.50))*(camp(:,:,:))
 	 
 	 ! Calculate the Power
 	 do i=1,no_3
 	  do j=1,no_2
 	    do ii = 1,no_1
 	      ps(ii,j,i) = (abs(camp(ii,j,i)))**2
 	    end do
 	  end do
 	 end do
 	end subroutine 	
 
end module DFT
