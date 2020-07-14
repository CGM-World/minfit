c..............................................................................
c

      INTEGER FUNCTION findfeat(ioni)

c
c     find the regions in the spectrum that are objective features to
c     some predetermined sigma level using the aperture method described
c     in LTW87
c
c     ioni = ion identification 
c     N_sigma = number of sigma for a feature detection
c
c     this routine called twice....
c
c     (1) upon input of the data to assertain if any NSIGMA features are
c     in the spectrum for the current region, for which case the WRKFLX
c     array is set to unity; called from routine readdata
c
c     (2) following the final fit to ascertain if any residuals are
c     found at the input significance level, NSIGMA.  For this case
c     WRKFLX is the final model result; called from routine fitregion
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      include           'minfit.h'
      integer            i,j,k,ioni,nwidth,buff,npix
      parameter          (nwidth = 3)
      double precision   lile(maxpix),lils(maxpix),bige,bigs,dw

c     assumed guilty until proven innocent; translation- assume no
c     features to start

      findfeat = 0

c     set the flag and grab the ion index

      j    = ioni
      buff = nwidth + 1
      npix = ndata(j)

c     compute the detection spectrum, be careful on the ends because the
c     lambda interval is an average over adjacent pixels

c     1st pixel

      i  = 1
      dw = 0.5 * abs(lambda(j,i)-lambda(j,i+1))
      IF (quality(j,i).eq.1) then
       lile(i) = (wrkflx(j,i)-data(j,i)) * dw
      ELSE
       lile(i) = 0.0d0
      END IF
      lils(i) = sigma(j,i) * dw

c     last pixel

      i       = npix
      dw      = 0.5 * abs(lambda(j,i-1)-lambda(j,i))
      IF (quality(j,i).eq.1) then
       lile(i) = (wrkflx(j,i)-data(j,i)) * dw
      ELSE
       lile(i) = 0.0d0
      END IF
      lils(i) = sigma(j,i) * dw

c     embedded pixels

      do 13 i=2,ndata(ioni)-1
       dw      = 0.5 * abs(lambda(j,i-1)-lambda(j,i+1))
       IF (quality(j,i).eq.1) then
        lile(i) = (wrkflx(j,i)-data(j,i)) * dw
       ELSE
        lile(i) = 0.0d0
       END IF
       lils(i) = sigma(j,i) * dw
 13   continue

c     if the region is so small that an aperture makes no sense, then
c     simply check that the region is not consistent with noise to the
c     input sigma, N_sigma, note the embedded RETURN statement

      if (npix.le.2*buff) then
       bige = 0.0d0
       bigs = 0.0d0
       do 21 k=1,npix
        bige = bige + lile(k)
        bigs = bigs + lils(k)**2
 21    continue
       if (abs(bige/sqrt(bigs)).ge.N_sigma) findfeat = 1
       return
      end if

c     since we apparently have the requisite number of pixels do the
c     full blown reduced flux detection, which basically checks for
c     N_sigma equivalent widths in each aperture

c     it only takes a single N_sigma aperture for an affirm but lets
c     count the number of "sigificant apertures"

      do 14 i=buff,npix-buff
       bige = 0.0d0
       bigs = 0.0d0
       do 17 k=i-nwidth,i+nwidth
        bige = bige + lile(k)
        bigs = bigs + lils(k)**2
 17    continue
       if (abs(bige/sqrt(bigs)).ge.N_sigma) then
        findfeat = findfeat + 1
       end if
 14   continue

c     return

      return
      end

c     eof
