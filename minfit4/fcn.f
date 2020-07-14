c     Module: fcn.f
c
c     ROUTINE fcn(iflag,m,n,a,fvec,dum,idum)
c     ROUTINE model(a,m)
c
c     fcn       - the chi-square for the LSF engine
c     model     - computes the model
c
c.......................................................................
c  

      SUBROUTINE         fcn(iflag,m,n,a,fvec,dum,idum)

c  
c     Another very important routine.  This routine is called by the
c     dnls1e fitting machine of SLATEC; you will not find a call to it
c     amongst the subroutines of MIINFIT.
c
c     We compute the model and communicate the progresss of the fit
c
c.......................................................................
c  

      include            'minfit.h'
      include            'const.dek'
      integer            k,m,n,iflag,idum
      integer            ioni,pixi,linei
      double precision   a(maxcoeffs),fvec(maxvec),dum(idum,1)

c
c     iflag = 0 is called every NPRINT iterations during the fit; it is
c     used for communication purposes, if so desired

      if (iflag.eq.0) then

c     update the zlines array so that the component ticks are current,
c     compute the fit statistics, and update the plotting window and
c     return

       k = 0
       do 100 linei=1,lines
        k = idx_az(linei)
        zline(linei) = a(k)
 100   continue
       CALL fitstats(a,n,m,iflag)
       if (plotting) CALL plt
       if (dodebug) CALL modelcomm(a,n,2) ! we are debugging so print the As
       return

      end if

c     OK, so we are computing the vector for which the quadrature sum is
c     being minimized by SLATEC's dnls1 routine; as written this is the
c     Chi square that all the fuss is about

c     stuff the wrkflux array with the call to model 

      CALL model(a,m)

c     stuff fvec; which is the chi-square; it is fvec that is minimized
c     note: dnls1e provides for an analytical Jacobian if one can be
c     writen out; however, as hard as I have tried to think about how to
c     do this (even approximations via expansion series), I have found
c     that the numerical calculation of the Jacobian from dnls1e does a
c     superior job and is not highly time consuming.  So, we do not
c     compute the Jacobian here, it is done by the fitting machine

c     use only unmasked pixels (quality=1); final k will equal m

      k = 0
      do 31 ioni=1,ions
       if (features(ioni)) then
        do 33 pixi=1,ndata(ioni)
         if (quality(ioni,pixi).eq.1) then
          k = k + 1
          fvec(k) = (data(ioni,pixi)-wrkflx(ioni,pixi)) 
     @            / abs(sigma(ioni,pixi))
         end if
 33     continue
       end if
 31   continue

c     return

      return
      end

c
c..............................................................................
c  
  
      SUBROUTINE         model(a,m)
  
c     compute the model (Merit Function); we call this routine everytime
c     we want to stuff the wrkflx array with the model
c  
c     Some notes:
c
c     b1  = column density
c     b2  = rest frame central wavelength
c     b3  = rest frame Doppler width
c     y   = natural broadening term in doppler units
c     x   = wavelength difference in doppler units
c     w0  = rest frame wavelength for transition
c     w   = rest frame wavelength for observed wavelength
c     tau = optical depth
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      include           'const.dek'
      integer           in,iz,ib
      integer           m,linei,speciesi,sidx,fidx,ioni,pixi,Naction
      double precision  a(maxcoeffs),Tp,satNscale,sngb
      double precision  w0,b1,b2,b3,w,x,y,u,v,tcon


c     set Naction for the call to satNscale

      Naction = 0  ! code units

c     clear the wrkflx array (this is a critical step)

      do 03 ioni=1,ions
       do 01 pixi=1,ndata(ioni)
        wrkflx(ioni,pixi) = 1.0d0
 01    continue
 03   continue
  
c     loop over lines, loop over species, do all ions of the species

      do 05 linei=1,lines

       do 07 speciesi=1,species

        do 09 ioni=1,ions

         sidx = spec_idx(ioni)
         fidx = fitindex(speciesi)

         if (features(ioni).AND.(sidx.eq.fidx)) then      

c     compute the indices and yank parameters; we have to take care on a
c     few issues; thus the conditionals
            
          w0   = lambda0(ioni)
          in   = idx_an(speciesi,linei)
          ib   = idx_ab(linei)
          iz   = idx_az(linei)

c     if we have defined saturated regions and we have a line that
c     classifies as being saturated, then use the column density
c     scaling; if not saturated then grab the coumn density directly;
c     there is a bit of a sign issue here because the columns can go
c     negative during the iterations; in order to preserve the negative
c     sign we make the scaling in two steps

c     the column density

          if (satflag(speciesi,linei)) then
           b1 = satNscale(a(in),acol(fidx),bcol(fidx),Naction)
           b1 = b1*sign(1.0d0,a(in))
          else
           b1   = a(in)
          end if

c     the line center

          b2   = w0 * (1.0d0+a(iz))/(1.0d0+zabs)

c     the doppler parameter

          sngb = sign(1.0d0,a(ib))
          Tp   = abs(a(ib))
          b3   = w0/c * sqrt(2.0d0*kerg*Tp/spec_mass(fidx))
          b3   = sngb*b3

c     loop over the data and stuff the workflux; experience shows that
c     we do not need to run through the full set of pixels; we go +/- 6
c     Doppler widths; this works for undamped lines; remove the
c     conditional if you are doing damped lines

c     if the data point is masked, do not include it in the model

          y    = con2(ioni) / b3
          tcon = con1(ioni) * (b1/b3)

          do 11 pixi=1,ndata(ioni)
           if (quality(ioni,pixi).eq.1) then
            w = lambda(ioni,pixi)/(1.0d0+zabs) 
            x = (w-b2)/b3 
            if (abs(x).le.50.d0) then
             CALL voigt(x,y,u,v)
             wrkflx(ioni,pixi) = wrkflx(ioni,pixi) * exp(-tcon*u)
            end if
           end if
 11       continue

         end if

 09     continue

 07    continue

 05   continue

c     we are GO for the convolve; if the logical is set high



      if (convolving) CALL convolve(m,0)

c     return

      return
      end
c  
c..............................................................................
c     eof
