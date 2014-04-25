      program smoothmodel

c
c     convolve model with quick and dirty psf 
c

      real*8 vlight
      parameter(maxi1=91901,maxx=2048,nxtra=100,
     .          mmax=8192,ma=7,nca=7,nmax=199,nover=2,
     .          atmpres=0.7895,vlight=2.99792458d5)
      parameter(ncols=4098)              
    
      real*8 modw(maxi1),mods(maxi1),y2mod(maxi1),           !stellar model arrays
     .       wsp(maxx),sp(maxx),yfit(maxx),spdiff(maxx),     !target spectrum and model arrays
     .       spectrum(ncols),slamscl(ncols)    !HRS target       
      real*8 work0(maxx), !work and rv post-proc arrays
     .     terwwork(mmax),terswork(mmax),work(mmax),tsmoo(mmax),
     .     tlam(mmax),xout(mmax)
      real*8 rms,mjd,x,dwl,wl0,ttpsf,trsh,xx,hh,y,          !scalar variables
     .       bccor,spsf,chisq,ave,var,rvold,wave,wvar,cmax
      double complex cpsf(mmax),cdata(mmax),fsp(mmax)   !arrays for FFT convolution
      character dateobs*10,ut*11,ra*11,dec*11,equinox*7,gascell*3,
     .          targname*18,HRSdata*80,junk*80,wat*8
      
      real*8 params(ma),psf0(nmax)                        !Marquardt param and error arrays
      

      open(42,file='coelho_models/4000_45_p00p00_red.dat')  !model stellar spectrum
      open(43,file='coelho_models/smoothM0.dat')

      write(6,*) 'Loading model stellar spectrum...'
      do i=1,maxi1
         read(42,*) modw(i),mods(i)
      enddo

      call dspline(modw,mods,maxi1,1.d30,1.d30,y2mod)!derivative array for stellar tpl


c    want 6700-6750

      npx1=35039 
      npx2=37535

      npixels=npx2-npx1+1

      dwl=(modw(npx2)-modw(npx1))/(npixels-1) !linear dispersion
      wl0=modw(npx1)
     
      npixels1=npixels+2*nxtra  !extra space for spline 
      write(*,*) npixels1
c      bccor=9.5100159d0*5.5d0   !shift to match observations, for plots
c      z=bccor/vlight
      z=0.d0

c resample and multiply
      x=wl0-nxtra*dwl
      do i=1,npixels1                                    
        xx=x/(1.d0+z)                           !Doppler shift stellar template
        call dsplint(modw,mods,y2mod,maxi1,xx,y) !resample at *nxfac resolution
        tlam(i)=x                                !put result in work arrays
        work(i)=y

        x=x+dwl
        cdata(i)=dcmplx(work(i),0.d0)!fill array for convolution
        xout(i)=x
      enddo
      
     
c make psf
             params(1)=0.d0     !intial rv offset, psf and wl-scale
             params(2)=wl0      !initial zero
             params(3)=dwl      !initial step
             params(4)=0.3d0
             params(5)=0.3d0
             params(6)=-0.02d0
             params(7)=-0.02d0

            write(6,*) params

            call psfgen(params,npixels1,hh,cpsf)         !generate psf for template smoothing
        
            write(*,*) 'done calling psfgen'

            call convolve(npixels1,cdata,cpsf,fsp) !do smoothing
            write(*,*) 'done calling convolve'

            do i=1,npixels1       !copy result to work array
               tsmoo(i)=dreal(fsp(i)/(npixels*hh))
            enddo

            do i=1,npixels
               write(43,*) xout(i+nxtra-1),tsmoo(i+nxtra-1) 
            enddo
 
 



      stop
      end

C *************************************************************************

      subroutine psfgen(params,npts,hh,cpsf)

c     Routine to generate an 5-component Gaussian psf following
c     the Butler et al. (1996) prescription. The central peak is fixed
c     at unit height and sigma = 1.5 pix. The other peaks have widths
c     fixed at sigma = 0.75 pix and adjustable hights. They are spaced
c     at the 1.5 pixel intervals.
c 
c     params(4,5): height of the inner pair of sidelobes 
c     ..................................................
c     params(6,7): height of the outer pair of sidelobes
c     w1-w5: offsets of sidelobes from the main peak (fixed)
c
c     ONLY USING 2 pairs of side lobes!

      parameter(nmax=199,lpsf=30,npar=7,nover=2,mmax=8192)
      real*8 params(npar),psf(nmax)
      real*8 j,jm,jp,x,y1,z1,y2,z2,y3,z3,y4,z4,y5,z5,
     .       alpha,beta,w1,w2,w3,w4,w5,norm,hh
      parameter(alpha=1.d0/(nover*2.02d0)**2,
     .          beta=1.d0/(nover*1.06d0)**2)
      parameter(w1=1.50d0*nover,w2=3.00d0*nover,w3=4.50d0*nover,
     .          w4=6.00d0*nover,w5=7.70d0*nover)
      double complex cpsf(mmax)
      logical badpsf
      common /psfflag/ badpsf

      badpsf=.false.
c      hh=2.*0.8326*12.0d0              !approx FWHM (norm factor after convolution)***
      hh=2.*0.8326*2.70d0
      n2=nmax/2-lpsf
      write(*,*) 'n2 = ',n2
      n=2*n2+1                                  !# pts in psf
      n0=(n+1)/2                                !peak of psf

      if(n.gt.nmax) then        !set flag, if psf too wide
         write(*,*) 'psf too wide'
      endif

      y1=-beta*w1*w1                            !calculate norm factor
      y2=-beta*w2*w2
c      y3=-beta*w3*w3
c      y4=-beta*w4*w4
c      y5=-beta*w5*w5
      norm=1.d0+params(4)*dexp(y1)+params(5)*dexp(y1)+
     .          params(6)*dexp(y2)+params(7)*dexp(y2)!+

c     .          params(5)*dexp(y3)+params(6)*dexp(y3)+
c     .          params(7)*dexp(y4)+params(8)*dexp(y4)+
c     .          params(9)*dexp(y5)+params(10)*dexp(y5)

      do i=1,n                                  !generate psf
        j=dfloat(i-n0)
        x=-alpha*j*j
        jm=j-w1
        jp=j+w1
        y1=-beta*jm*jm
        z1=-beta*jp*jp
        jm=j-w2
        jp=j+w2
        y2=-beta*jm*jm
        z2=-beta*jp*jp
c        jm=j-w3
c        jp=j+w3
c        y3=-beta*jm*jm
c        z3=-beta*jp*jp
c        jm=j-w4
c        jp=j+w4
c        y4=-beta*jm*jm
c        z4=-beta*jp*jp
c        jm=j-w5
c        jp=j+w5
c        y5=-beta*jm*jm
c        z5=-beta*jp*jp
        psf(i)=(dexp(x)+params(4)*dexp(y1)+params(5)*dexp(z1)+
     .                  params(6)*dexp(y2)+params(7)*dexp(z2))/norm!+
c
c     .                  params(5)*dexp(y3)+params(6)*dexp(z3)+
c     .                  params(7)*dexp(y4)+params(8)*dexp(z4)+
c     .                  params(9)*dexp(y5)+params(10)*dexp(z5))
c     .                  /norm
      enddo


      do i=1,npts                              !reset work array
         cpsf(i)=(0.d0,0.d0)         
      enddo      
      nwrk=1-n0
      do i=n0,n                                     !fill work array for convolution
        cpsf(i+nwrk)=dcmplx(psf(i),0.d0)
      enddo
      nwrk=npts-n2
      do i=1,n2
        cpsf(nwrk+i)=dcmplx(psf(i),0.d0)
      enddo

      
      return
      end



c***********************************************************************************

      subroutine convolve(n,sp,psf,fsp)

c     Routine to convolve an input spectrum with a PSF

      include '/sw/include/fftw3.f' 
      integer*8 plan
      double complex sp(n),fsp(n),psf(n),fpsf(n)

      call dfftw_plan_dft_1d(plan,n,sp,fsp,                      !forward ft of spectrum
     .                       FFTW_FORWARD,FFTW_ESTIMATE)
      call dfftw_execute(plan)
      call dfftw_destroy_plan(plan)
      call dfftw_plan_dft_1d(plan,n,psf,fpsf,                    !forward ft of PSF
     .                       FFTW_FORWARD,FFTW_ESTIMATE)
      call dfftw_execute(plan)
      call dfftw_destroy_plan(plan)

      do i=1,n                                                   !multiply transforms
        fpsf(i)=fsp(i)*fpsf(i)
      enddo
      call dfftw_plan_dft_1d(plan,n,fpsf,fsp,                    !backward ft of product
     .                     FFTW_BACKWARD,FFTW_ESTIMATE)
      call dfftw_execute(plan)
      call dfftw_destroy_plan(plan)

      return
      end




C *************************************************************************


      include '../../fortran.lib/dNrF77/dspline.f'
      include '../../fortran.lib/dNrF77/dsplint.f'

