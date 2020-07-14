c.......................................................................
c  

      DOUBLE PRECISION FUNCTION fluxpix(a,ion,cenpix)

c  
c     return the model value at a single pixel after convolving (by
c     hand) over a small range of the model, centered on
c     WRKFLX(ion,cenpix); this routine called repeatedly by routine
c     errors, it is used for estimating the uncertainties in the fitting
c     parameters- it is used to compute the curvature matrix; to do so,
c     a model must be convolved for each test, but we only need to treat
c     the small chunk of spectrum surrounding the pixel being examined
c     corresponding to the small adjustment in the target parameter
c
c     we don't get here unless FEATURES(ION) = .true.
c
c     *** DO NOT CHANGE the value of CENPIX ****
c
c.......................................................................
c  

      include            'minfit.h'
      include            'const.dek'
      integer            i,k,ion,cenpix,pix1,pix2,pixc,pixi
      integer            window,ndat
      double precision   a(maxcoeffs),f,df,norm,phi,dv
      double precision   xa(maxpix),ya(maxpix),x,y,y2a(maxpix)
      double precision   isf


c     some notes...
c
c     NRESPONSE, the response length is already known from the call to
c     routine INSTRUMENT
c
c     the convolution WINDOW can be odd or even; here we define the
c     pixel range, do a sanity check, see further comments below

      window =  nresponse/int(resfac) 
      ndat   = 2*window + 1
      pix1   = cenpix - window
      pix2   = cenpix + window
      if (pix1.lt.1) pix1 = 1
      if (pix2.gt.ndata(ion)) pix2 = ndata(ion)

c     stuff the chunk of the WRKFLX array

      CALL modelchunk(a,ion,pix1,pix2)

c     set the pixel counter and prepare for the interpolation...  pad
c     the data if we are reach the pixel limits

      i = 0

c     pad lower edge?

      if (pix1.ne.(cenpix-window)) then
       do 25 pixi=1,pix1-(cenpix-window)
        i = i + 1
        xa(i)  = real(i)
        ya(i)  = 1.0d0
        y2a(i) = 0.0d0
 25    continue
      end if

c     stuff WRKFLX chunk

      do 23 pixi=pix1,pix2
       i = i + 1
       xa(i)  = real(i)
       ya(i)  = wrkflx(ion,pixi)
       y2a(i) = 0.0d0
 23   continue

c     pad upper edge?

      if (pix2.ne.(cenpix+window)) then
       do 27 pixi=1,(cenpix+window)-pix2
        i = i + 1
        xa(i)  = real(i)
        ya(i)  = 1.0d0
        y2a(i) = 0.0d0
 27   continue
      end if

c     determine the interpolation coefficients for the change in pixel
c     sampling (done to make the convolution smoother)

      call spline(xa,ya,ndat,y2a)

c     we only need the result for the target pixel, at CENPIX, thus, we
c     need only convolve about this pixel and we can loop over both +/-
c     windows simultaneously; this loop computes the unnormalized
c     incremental change to be applied to the target pixel, DF; it also
c     computes the normalization, NORM

c     redefine the central pixel in terms of the above indexing and
c     stuff the unconvolved flux at CENPIX; for the convolution, the
c     "pixel" sampling is increased by factor RESFAC (entered in the
c     minfit.par file by the user)

      pixc = window + 1
      f    = ya(pixc)
      norm = phi(0.0d0,profile)
      df   = 0.0d0
      k    = window*int(resfac)

c     manually convolve the model in this chunk of spectrum' loop over
c     the pixels in the convolution window and pluck off the
c     interpolated flux values; manually convolve around the pixel; note
c     the call to phi, which returns the ISF

      do 29 i=1,k

       dv   = real(k-(i-1))*hdv
       isf  = phi(dv,profile)
       norm = norm + 2.0d0*isf

c     lower window convolution contribution

       x = 1.0d0 + real(i-1)/resfac
       call splint(xa,ya,y2a,ndat,x,y)
       df = df + y*isf

c     upper window convolution contribution

       x = real(ndat) - real(i-1)/resfac
       call splint(xa,ya,y2a,ndat,x,y)
       df = df + y*isf

 29   continue

c     add the "correction", the normalized flux increment

      fluxpix = (f + df)/norm

c     return

      return
      end

c
c..............................................................................
c  

      SUBROUTINE        modelchunk(a,ion,pix1,pix2)

c  
c     returns a small chunk of workflux, for ion ION from PIX1 to PIX2
c     called by routine convpix
c  
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      include           'const.dek'
      integer           linei,speciesi,ion,pixi,in,iz,ib
      integer           pix1,pix2,voigtwin,Naction,fidx
      double precision  a(maxcoeffs),Tp,satNscale
      double precision  w0,b1,b2,b3,w,x,y,u,v,tcon


c     set Naction for the call satNscale

      Naction = 0 ! code units

c     set the a array vector pointer

      voigtwin = 2.0d0*conwindo

c     initialize the wrkflx array (this is a critical step)

       do 01 pixi=pix1,pix2
        wrkflx(ion,pixi) = 1.0d0
 01    continue

c     loop over lines, loop over species of this ion

      do 05 linei=1,lines
       do 07 speciesi=1,species

         fidx = fitindex(speciesi)
         if (spec_idx(ion).eq.fidx) then      

c     compute the indices and yank parameters

          w0   = lambda0(ion)
          in   = idx_an(speciesi,linei)
          ib   = idx_ab(linei)
          iz   = idx_az(linei)

c     if we have defined saturated regions and we have a line that
c     classifies as being saturated, then use the column density
c     scaling; if not saturated then grab the coumn density directly;
c     there is a bit of a sign issue here because the columns can go
c     negative during the iterations; in order to preserve the negative
c     sign we make the scaling in two steps

          if (satflag(speciesi,linei)) then
           b1 = satNscale(a(in),acol(fidx),bcol(fidx),Naction)
           b1 = b1*sign(1.0d0,a(in))
          else
           b1   = a(in)
          end if

          b2   = w0 * (1.0d0+a(iz))/(1.0d0+zabs)
          Tp   = a(ib)
          b3   = w0/c * sqrt(2.0d0*kerg*Tp/spec_mass(fidx))
          y    = con2(ion) / b3
          tcon = con1(ion) * (b1/b3)

c     loop over the narrow pixel range and stuff the work flux with the
c     model

          do 11 pixi=pix1,pix2
           w = lambda(ion,pixi)/(1.0d0+zabs) 
           x = (w-b2)/b3 
           if (abs(x).le.voigtwin) then
            CALL voigt(x,y,u,v)
            wrkflx(ion,pixi) = wrkflx(ion,pixi) * exp(-tcon*u)
           end if
 11       continue
         end if
 07    continue
 05   continue

c     return

      return

      end

c     eof
