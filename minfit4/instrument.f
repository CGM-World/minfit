c.......................................................................
c
c  routine instrument  stuffs the instrumental spread function array; this
c                      called once
c
c  routine phi         computes the relative instrumental response 
c                      (assumes a Gaussian)
c
c routine initconv      sets up the padding arrays for the FFTS; this 
c                       routine must be called everytime M changes
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      SUBROUTINE          instrument(flag)
c
c
c     must be called after routine GETLIST
c     given the instrumental profile sigma in velocity units
c
c     these profiles are loaded in wrap-around order for the convolution
c     zero spatial information is the first index
c
c     flag = 0 ; set up; do not communicate
c
c     flag = 1 ; set up; communicate
c
c.......................................................................
c
      include             'minfit.h'
      include             'const.dek'
      integer             i,j,iplus,iminus,flag,pixi
      double precision    dv,xdv,norm,phi
      double precision    wave,vel1,vel
      character*80        ion_str

c     scan the data and grab the mean pixel velocity width, dv we do not
c     account for bad or masked pixels since the pixel information is
c     not corrupt, only the flux values


      WRITE(SCREEN,*) '************************************************'
      WRITE(SCREEN,*) ' '
      write(SCREEN,*) 'DATA CHARACTERISTICS'
      WRITE(SCREEN,*) ' '
      WRITE(STDOUT,*) '************************************************'
      WRITE(STDOUT,*) ' '
      write(STDOUT,*) 'DATA CHARACTERISTICS'
      WRITE(STDOUT,*) ' '


      pixi = 0
      dv   = 0.0d0

      DO 07 j=1,ions
       ion_str = ion_name(j)
       OPEN(unit=1,file=ion_str,err=999,status='old')
       read(1,*,end=06) wave,vel1
       pixi = pixi + 1
       do 05 i=2,maxpix
         read(1,*,end=06) wave,vel
         pixi = pixi + 1
         dv = dv + abs(vel1 - vel)
         vel1 = vel
 05    continue
 06    CLOSE(unit=1)
 07   CONTINUE

      dv = dv/float(pixi)

c     compute the instrumental resolution in velocity units [km/s] sigma
c     in km/s given by FWHM/2.35

      profile = ckms/(2.35*R_fac)  ! km/s

c     dv is the velocity sampling of the pixels km/s/pixel...  the
c     number of pixels per resolution element = profile/dv

      pixpres = 2.35 * profile/(dv*slit)
      hdv = dv/resfac
      nresponse  = 2*int(conwindo*profile/hdv) + 1

c     now stuff the response function in wrap around order

      response(1) = phi(0.0d0,profile)
      norm        = response(1)
      do 11 i=1,int(nresponse/2)
       xdv              = real(i)*hdv 
       iplus            = i + 1
       iminus           = nresponse - (i-1)
       response(iplus)  = phi(xdv,profile)
       response(iminus) = response(iplus)
       norm = norm + response(iplus) + response(iminus)
 11   continue
  
c     for the convolution integral, the response function needs to be
c     normalized or flux is not conserved...  unit width is used because
c     data are discretized by pixel indices in the convolution

      do 13 i=1,nresponse
       response(i) = response(i)/norm
 13   continue
  

c     coomunicate?

      IF (flag.eq.1) then 
       DO 21 i=1,2
        IF (i.eq.1) j = SCREEN
        IF (i.eq.2) j = STDOUT
C        write(j,'(a,f8.0)') 
C     @  ' -Spectrograph R     [lam/dlam]     =',R_fac
C        write(j,'(a,f8.3)') 
C     @  ' -Effective Slit Width [arcsec]     =',slit
        write(j,'(a,f8.3)') 
     @  ' -Instrumental Sigma         [km/s] =',profile
        write(j,'(a,f8.3)') 
     @  ' -Pixel Resolution           [km/s] =',dv
        write(j,'(a,f8.3)') 
     @  ' -Pix/ResEl                  [FWHM] =',pixpres
        write(j,'(a,i5)')   
     @  ' -Convolution Response Length [pix] =',nresponse
        WRITE(j,*) ' '
 21    CONTINUE
      END IF

c     OK return

      RETURN

c     error on read

 999  WRITE(SCREEN,*) ' ERROR(instrument) - cannot read'
      WRITE(SCREEN,*) ' data for ',ion_str
      WRITE(STDOUT,*) ' ERROR(instrument) - cannot read'
      WRITE(STDOUT,*) ' data for ',ion_str
      STOP

      END

c  
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  
  
      double precision function phi(dv,width)
  
c  
c     given the sigma in velocity units this routine computes the
c     relative Gaussian value of the instrumental profile.  called by
c     routine instrument iteratively
c.......................................................................
c  

      include             'minfit.h'
      include             'const.dek'
      double precision     dv,z,width


c     dv is value at which the instrumental response is to be evaluated

      z   = dv/width
      phi = dexp(-(z**2)/2.0d0)

      return
      end

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      SUBROUTINE          initconv(m)
c
c
c     the convolution requires the number of elements in the array being
c     convolved is a power of 2; thus, we compute the number of array
c     elements including the interpolation (RESFAC) and then use a look
c     up table to determine the next lowest power of 2 and grab that
c
c.......................................................................
c
      include             'minfit.h'
      include             'const.dek'
      integer             i,m,np2
      parameter           (np2 = 11)
      integer             pwrsof2(np2)

      data pwrsof2 /16,32,64,128,256,512,1024,2048,4096,8192,16384/


c     multiply by resfac in order to increase the sampling rate (for a
c     smoother convolution); compute NFFT

      ncondat = int(resfac)*(m-1) - 1
      nfft   = ncondat + (nresponse-1)/2 + 1 

c     replace NFFT with a power of 2 that is just larger

      do 15 i=1,np2
       if (nfft.le.pwrsof2(i)) then ! grab first occurance
        nfft = pwrsof2(i)
        GOTO 17
       end if
 15   continue

c     this is seen only if we exceed the powers of 2 table

      WRITE(SCREEN,*) ' ERROR:(initconv): failed'
      write(SCREEN,*) ' M         = ',m
      write(SCREEN,*) ' NCONDAT   = ',ncondat
      write(SCREEN,*) ' NFFT      = ',nfft
      write(SCREEN,*) ' NFFT MAX  = ',maxcon
      write(STDOUT,*) ' NRESPONSE = ',nresponse
      write(STDOUT,*) ' M         = ',m
      write(STDOUT,*) ' NCONDAT   = ',ncondat
      write(STDOUT,*) ' NFFT      = ',nfft
      write(STDOUT,*) ' NFFT MAX  = ',maxcon
      if (nfft.gt.maxcon) then
       STOP ' NFFT > MAXCON'
      end if

c     double sanity check (what the hell?)

 17   if (nfft.gt.maxcon) then
       STOP ' ERROR:(initconv): NFFT > MAXCON'
      end if

C      write(SCREEN,*) ' CONVOLUTION INFORMATION -'
C      write(SCREEN,*) ' M         = ',m
C      write(SCREEN,*) ' NCONDAT   = ',ncondat
C      write(SCREEN,*) ' NFFT      = ',nfft
C      write(SCREEN,*) ' NFFT MAX  = ',maxcon



c     return

      RETURN

      END

c  
c     eof
