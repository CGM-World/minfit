c  
c.......................................................................
c  

      SUBROUTINE        fitstats(a,n,m,iflag)
  
c     compute the chisq, data variance, and fit variance; exclude masked
c     pixels with QUALITY=0
c
c     IFLAG=-1; this is a call from routine FTEST and we do not
c     communicate anything to the user
c
c     IFLAG=0; this is a call from routine FCN and we communicate the
c     intermediate result
c
c     IFLAG=1; this is a call from routine FITREGION or routine
c     TESTLINES and it is the current the LSF result; communicate this
c     to the use
c  
c     IFLAG=2; this is a final call from routine FITREGION and it is the
c     final LSF result for the current region; communicate this and
c     store the values
c  
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      integer           n,m,ioni,pixi,iflag
      double precision  a(maxcoeffs),nu,
     @                  vardata,chisqnu,varfit
 

c     compute the model, stuffs the WRKFLX array with the current model

      CALL model(a,m)

c     zero the variables

      chisq   = 0.0d0
      vardata = 0.0d0

c     compute the fit chisq and variances; include only unmasked pixels
c     (QUALITY=1); we compute the degrees of freedom in this loop by
c     accounting for the number of unmasked pixels (NU = number of
c     quality functions less number of fitting parameters); the number
c     of functions is equal to the number of quality=1 pixels from all
c     transitions; the value of M containes the number of functions, the
c     valur N containes the number of fitting parameters

      nu = real(m-n)

      do 21 ioni=1,ions
       if (features(ioni)) then 
        do 23 pixi=1,ndata(ioni)
         if (quality(ioni,pixi).eq.1) then
           chisq = chisq + (data(ioni,pixi)-wrkflx(ioni,pixi))**2 
     @                    / sigma(ioni,pixi)**2
           vardata = vardata + 1.0d0/(sigma(ioni,pixi)**2)
         end if
 23     continue
       end if
 21   continue

      vardata = real(m)/vardata   
      chisqnu = chisq/nu          
      varfit  = chisqnu*vardata

c     output the fitted parameters and then the fit statsistics
c     iflag = -1 (actually any value !=0,1,2) no communication

c     iflag = 0 called from FCN, communicate intermediate result

      if (iflag.eq.0) then
      if (oldchi2.gt.0.0d0) delchi2 = oldchi2-chisqnu
      write(SCREEN,102) chisqnu,delchi2
      write(STDOUT,102) chisqnu,delchi2
      end if

c     iflag = 1 called from FITREGION or TESTLINES

      if (iflag.eq.1) then
       WRITE(SCREEN,600)
       write(SCREEN,100)
       write(SCREEN,101) nu,chisq,chisqnu,vardata,varfit
       WRITE(STDOUT,600)
       write(STDOUT,100)
       write(STDOUT,101) nu,chisq,chisqnu,vardata,varfit
      end if

c     iflag = 2 call from driver MINFIT (final model for region)

      if (iflag.eq.2) then
       f_nu(regioni)      = nu
       f_npar(regioni)    = float(n)
       f_chisq(regioni)   = chisq
       f_chisqnu(regioni) = chisqnu
       f_vardata(regioni) = vardata
       f_varfit(regioni)  = varfit
       write(SCREEN,103)
       write(SCREEN,101) nu,chisq,chisqnu,vardata,varfit
       write(STDOUT,103)
       write(STDOUT,101) nu,chisq,chisqnu,vardata,varfit
      end if


c     store the chi-square so that we can compare for the next call

      oldchi2 = chisqnu

c     we are done, return

      return

c     formats

 600  FORMAT(1x,'+++++++++++++++++++++++++++++++++++++++++++++++++++')
 100  format(1x,/,1x,t20,'LSF STATISTICS',/,/,
     @       1x,t5,'nu',t16,'ChiSq',t27,'ChiSqnu',
     @       t38,'Var Data',t49,'Var Fit')
 101  format(1x,1p5e11.3)
 102  format(1x,'Reduced Chi-Sq = ',1pe11.3,
     @       4x,'Del Chi-Sq =',1pe11.3)
 103  format(1x,/,1x,'LSF FIT MODEL STATISTICS',/,/,
     @       1x,t5,'nu',t16,'ChiSq',t27,'ChiSqnu',
     @       t38,'Var Data',t49,'Var Fit')

      end

c  
c.......................................................................
c  
  
      SUBROUTINE        inputmod(m)

c  
c     compute the full model for all ions on input, called once before
c     the LSF
c  
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      include           'const.dek'
      integer           m,linei,speciesi,ioni,pixi
      double precision  w0,b1,b2,b3,w,x,y,u,v,tcon

  
c     some notation for reading ease
c
c     b1  = column density
c     b2  = rest frame central wavelength
c     b3  = rest frame Doppler width
c     y   = natural broadening term in doppler units
c     x   = wavelength difference in doppler units
c     w0  = rest frame wavelength for transition
c     w   = rest frame wavelength for observed wavelength
c     tau = optical depth

c     clear the wrkflx array (this is a critical step)

      do 03 ioni=1,ions
       do 01 pixi=1,ndata(ioni)
        wrkflx(ioni,pixi) = 1.0d0
 01    continue
 03   continue
  
c     if this is the call before the LSF, then we use the input
c     parameters before folding them into the fitting array A and SIGA;
c     if this is the final call, we have already unpacked the A and SIGA
c     arrays and computed the limits

c     loop over lines, loop over species, do all ions of the species,
c     but avoid those components that have only limits

      do 05 linei=1,lines
       do 07 speciesi=1,total

c     skip if the component for this species is a limit

        if (dnline(speciesi,linei).ne.-1.0) then

c     loop over the ions

         do 09 ioni=1,ions

c     skip unless this ion is of the current species

          if (spec_idx(ioni).eq.speciesi) then      

c     OK, got one, now compute the rad tran inputs

           w0   = lambda0(ioni)
           b1   = 10.0d0**(nline(speciesi,linei)-13.0d0)
           b2   = w0 * (1.0d0+zline(linei))/(1.0d0+zabs)
           b3   = w0 * bline(speciesi,linei)/ckms
           y    = con2(ioni) / b3
           tcon = con1(ioni) * (b1/b3)

c     loop over the data and stuff the workflux

           do 11 pixi=1,ndata(ioni)
            w = lambda(ioni,pixi)/(1.0d0+zabs) 
            x = (w-b2)/b3 
            CALL voigt(x,y,u,v)
            wrkflx(ioni,pixi) = wrkflx(ioni,pixi) * exp(-tcon*u)
 11        continue

          end if

 09      continue

        end if

 07    continue
 05   continue

c     and now convolve

      if (convolving) CALL convolve(m,0)

c     done, return

      return
      end

c.......................................................................
c  
  
      SUBROUTINE        outputmod(m)

c  
c     compute the full model for all ions; called once after all regions
c     have completed LSF
c
c     this routine called once from routine PLTOUTPUT, which is also
c     called once from driver MINFIT after all regions are fit.  
c
c     this routine is is called prior to the call to routine FINALCOMM,
c     which is important because this means the the final model of all
c     regions is stuffed in the WRKFLUX array; thus we can compute the
c     final overall chi-square from WRKFLUX here
c  
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      include           'const.dek'
      integer           m,linei,speciesi,ioni,pixi
      double precision  w0,b1,b2,b3,w,x,y,u,v,tcon

  
c     some notation for reading ease
c
c     b1 = column density
c     b2 = rest frame central wavelength
c     b3 = rest frame Doppler width
c     y = natural broadening term in doppler units
c     x = wavelength difference in doppler units
c     w0 = rest frame wavelength for transition
c     w = rest frame wavelength for observed wavelength
c     tau = optical depth

c     clear the wrkflx array (this is a critical step)
c     clear the modflx array (this is also a critical step)

      do 03 ioni=1,ions
       do 01 pixi=1,ndata(ioni)
        wrkflx(ioni,pixi) = 1.0d0
C        do 02 linei=1,vpcomps
C         modflx(ioni,linei,pixi) = 1.0d0
C 02     continue
 01    continue
 03   continue
  

c     loop over all lines of all regions, loop over all species, do all
c     ions of the species, but avoid those components that have only
c     limits

      do 05 linei=1,vpcomps
       do 07 speciesi=1,total

c     skip if the component for this species is a limit

        if (dNfinal(speciesi,linei).ne.-1.0d0) then

c     loop over the ions

         do 09 ioni=1,ions

c     skip unless this ion is of the current species

          if (spec_idx(ioni).eq.speciesi) then      

c     OK, got one, now compute the rad tran inputs

           w0   = lambda0(ioni)
           b1   = 10.0d0**(Nfinal(speciesi,linei)-13.0d0)
           b2   = w0 * (1.0d0+zfinal(linei))/(1.0d0+zabs)
           b3   = w0 * bfinal(speciesi,linei)/ckms
           y    = con2(ioni) / b3
           tcon = con1(ioni) * (b1/b3)

c     loop over the data and stuff the workflux

           do 11 pixi=1,ndata(ioni)
            w = lambda(ioni,pixi)/(1.0d0+zabs) 
            x = (w-b2)/b3 
            CALL voigt(x,y,u,v)
            wrkflx(ioni,pixi)       = wrkflx(ioni,pixi)*exp(-tcon*u)
C            modflx(ioni,linei,pixi) = modflx(ioni,linei,pixi)
C     &                                *exp(-tcon*u)
 11        continue

          end if

 09      continue

        end if

 07    continue
 05   continue

c     and now convolve

      IF (convolving) CALL convolve(m,1)

c     return

      return
      end

c.......................................................................
c  

      SUBROUTINE        storefits
  
c
c     store the final results for this region, increment the total
c     number of VP components, these values used in routine OUTPUTMOD
c  
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      integer           i,k,vpknt
 

c     the counter VPKNT is local and starts at the last value of VPCOMPS
c     which is global

      DO 13 i=1,lines
       DO 11 k=1,total

         vpknt = vpcomps + i

         Rfinal(vpknt)    = regioni
         Nfinal(k,vpknt)  = nline(k,i)
         dNfinal(k,vpknt) = dnline(k,i)
         bfinal(k,vpknt)  = bline(k,i)
         dbfinal(k,vpknt) = dbline(k,i)
         zfinal(vpknt)    = zline(i)
         dzfinal(vpknt)   = dzline(i)

 11    CONTINUE
 13   CONTINUE

      vpcomps = vpknt

      RETURN
      END


c.......................................................................
c  

      SUBROUTINE        finalcomm
  
c
c     communicate the final model in all its glory
c  
c     this routine is called once following the single call to routine
c     PLTOUTPUT which calls routine OUTPUTMOD once.  Thus, we have the
c     final WRKFLUX array in hand.
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      include           'const.dek'
      integer           i,j,ireg
      double precision  vel
      character*80      label(4)
      data label        /' logN','dlogN','  b  ','  db '/


c     write the header 

      WRITE(SCREEN,*) ' '
      WRITE(SCREEN,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(SCREEN,*) '+                                              +'
      WRITE(SCREEN,*) '+      FINAL LSF RESULTS (ALL REGIONS)         +'
      WRITE(SCREEN,*) '+                                              +'
      WRITE(SCREEN,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(SCREEN,*) ' '

      WRITE(STDOUT,*) ' '
      WRITE(STDOUT,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(STDOUT,*) '+                                              +'
      WRITE(STDOUT,*) '+      FINAL LSF RESULTS (ALL REGIONS)         +'
      WRITE(STDOUT,*) '+                                              +'
      WRITE(STDOUT,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(STDOUT,*) ' '

      write(SCREEN,600) (spec_name(j),spec_name(j),j=1,total)
      write(SCREEN,601) (label(1),label(2),label(3),
     @                   label(4),j=1,total)

      write(STDOUT,600) (spec_name(j),spec_name(j),j=1,total)
      write(STDOUT,601) (label(1),label(2),label(3),
     @                   label(4),j=1,total)

c     write params

      do 53 i=1,vpcomps
       ireg = int(Rfinal(i))
       vel = ckms * (zfinal(i)-zabs)/(1.0d0+zabs)
       write(SCREEN,602) ireg,zfinal(i),(Nfinal(j,i),dNfinal(j,i),
     @               bfinal(j,i),dbfinal(j,i),j=1,total)
       write(STDOUT,602) ireg,zfinal(i),(Nfinal(j,i),dNfinal(j,i),
     @               bfinal(j,i),dbfinal(j,i),j=1,total)
 53   continue

c
c     communicate the LSF stats for each region and the total fit
c

      WRITE(SCREEN,700)
      DO 33 i=1,regions
       WRITE(SCREEN,701) i,f_nu(i),f_chisq(i),f_chisqnu(i),
     &                   f_vardata(i),f_varfit(i)
 33   CONTINUE

      CALL finalstats      

      i = regions + 1
      WRITE(SCREEN,702) 
      WRITE(SCREEN,703)  f_nu(i),f_chisq(i),f_chisqnu(i),
     &                   f_vardata(i),f_varfit(i)

      write(SCREEN,*) ' '
      write(SCREEN,*) '************************************************'
      write(STDOUT,*) ' '
      write(STDOUT,*) '************************************************'


      RETURN

 600  format(1x,t22,10a15)
 601  format(1x,'Reg',t7,'redshift',t18,10(a7,a8))
 602  format(1x,i3,2x,f8.6,10(2x,f6.2,1x,f6.2))

 700  format(1x,/,1x,'LSF STATISTICS - BY REGION',/,/,
     @       1x,'Reg',t8,'nu',t19,'ChiSq',t30,'ChiSqnu',
     @       t41,'Var Data',t52,'Var Fit')
 701  format(1x,i3,1p5e11.3)
 702  format(1x,/,1x,'FULL LSF STATISTICS - ALL REGIONS',/,/,
     @       1x,t8,'nu',t19,'ChiSq',t30,'ChiSqnu',
     @       t41,'Var Data',t52,'Var Fit')
 703  format(1x,3x,1p5e11.3)

      END

c  
c  
c.......................................................................
c  

      SUBROUTINE        finalstats
  
c     compute the chisq, data variance, and fit variance; exclude masked
c     pixels with QUALITY=0
c
c  
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      integer           j,ioni,pixi,npix,npar,idx
      double precision  vardata,chisqnu,varfit,nu
      double precision  w0,z1,z2,z 

c     zero the variables

      chisq   = 0.0d0
      vardata = 0.0d0
      npix    = 0
      npar    = 0

c     compute the fit chisq and variances; include only unmasked pixels
c     (QUALITY=1); we compute the degrees of freedom in this loop by
c     accounting for the number of unmasked pixels (NU = number of
c     quality functions less number of fitting parameters); the number
c     of functions is equal to the number of quality=1 pixels from all
c     transitions; the value of M containes the number of functions, the
c     valur N containes the number of fitting parameters

      do 22 j=1,regions
       z1 = zlim(j,1)
       z2 = zlim(j,2)
       npar = npar + int(f_npar(j))
       do 21 ioni=1,ions
        if (featreg(ioni,j)) then 
         w0 = lambda0(ioni)
         do 23 pixi=1,ndata(ioni)
          z = lambda(ioni,pixi)/w0 - 1.0d0
          if ((z.ge.z1).AND.(z.le.z2)) then ! we are in the region J
           if (quality(ioni,pixi).eq.1) then
             chisq = chisq + (data(ioni,pixi)-wrkflx(ioni,pixi))**2 
     @                     / sigma(ioni,pixi)**2
             vardata = vardata + 1.0d0/(sigma(ioni,pixi)**2)
             npix = npix + 1
           end if
          end if
 23      continue
        end if
 21    continue
 22   continue

      nu      = real(npix-npar)
      vardata = real(npix)/vardata   
      chisqnu = chisq/nu          
      varfit  = chisqnu*vardata

c     we use the final stats arrays, setting the index to regions+1

      idx = regions + 1

      f_nu(idx)      = nu
      f_chisq(idx)   = chisq
      f_chisqnu(idx) = chisqnu
      f_vardata(idx) = vardata
      f_varfit(idx)  = varfit

c     we are done, return

      return

      end

c  
c     eof
