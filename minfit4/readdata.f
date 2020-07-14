c.........................................................................
c

      SUBROUTINE         readdata(flag,error)

c
c     reads in data in the following format (compatible with PROFIT):
c     spectral info (file called ion_name):
c
c     wavelength velocity flux noise
c
c     flag=0   we have been called before or after the fit, so that
c              we do not want to remove an ion because of no features
c
c     flag=1   we have been called for fitting the region and we need
c              to discard regions where no features are found
c
c     error    logical flag that says we found no data in this region
c              so that the calling routine can bounce to the next region
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include            'minfit.h'
      logical            error
      integer            i,j,k,ioni,pixi,pixq,findfeat,flag,qdata
      double precision   wave,vel,flux,noise1,noise2,cont,z,dz
      character*80       ion_str,ion_mask,ch_stat


c     white-wash previous region, set the logical features flag high,
c     and set the logical error flag low

      if (flag.eq.1) then
      WRITE(SCREEN,*) '************************************************'
      WRITE(SCREEN,*) ' '
      write(SCREEN,*) 'DATA INPUT'
      WRITE(SCREEN,*) ' '
      WRITE(STDOUT,*) '************************************************'
      WRITE(STDOUT,*) ' '
      write(STDOUT,*) 'DATA INPUT'
      WRITE(STDOUT,*) ' '
      end if

      error = .false.

      do 02 ioni=1,maxions
       features(ioni)        = .true.
       featreg(ioni,regioni) = .true.
       do 01 pixi=1,maxpix
        lambda(ioni,pixi)  = 0.0d0
        data(ioni,pixi)    = 0.0d0
        sigma(ioni,pixi)   = 0.0d0
        wrkflx(ioni,pixi)  = 1.0d0
        quality(ioni,pixi) = 1
 01    continue
 02   continue

c     PHASE 1: read in the data within the current redshift region do
c     array accounting on the fly

      if (flag.eq.1) then
       write(SCREEN,602) regioni,zlim(regioni,1),zlim(regioni,2)
       write(STDOUT,602) regioni,zlim(regioni,1),zlim(regioni,2)
      end if

c     loop over the ions (ionic transitions), open their data files and
c     read in the arrays in the current redshift region

      do 03 ioni=1,ions

       ion_str     = ion_name(ioni)
       pixi        = 0
       ndata(ioni) = 0

c     open the file containing the spectra, assume 6 columns; we the
c     file does not exist, then we terminate (999)

       OPEN(unit=1,file=ion_str,err=999,status='old')

       do 05 i=1,maxpix
        read(1,*,end=06) wave,vel,flux,noise1,noise2,cont
        z = wave/lambda0(ioni) - 1.0d0
        if ((z.ge.zlim(regioni,1)).AND.(z.le.zlim(regioni,2))) then
         pixi = pixi + 1
         lambda(ioni,pixi) = wave  
         data(ioni,pixi)   = flux/cont
         sigma(ioni,pixi)  = noise1 * sqrt(fudge) / cont
        end if
  05   continue

  06   CLOSE(unit=1)
       ndata(ioni) = pixi
       nqual(ioni) = pixi
       nmask(ioni) = 0

c     read in the quality spectrum for statistics accounting of the LSF
c     if the mask exists, then read it in; if not skip it and all pixels
c     are set with quality high (=1) ; no sanity check is performed; it
c     is expected that if the user has masks that they are in good order
c     masked pixels are created in program MASKORDERS

       CALL fappend(ion_str,'mask',ion_mask)
       OPEN(unit=3,file=ion_mask,err=03,status='old')  
       pixq = 0
       pixi = 0
       DO 07 i=1,maxpix
        READ(3,*,end=08) wave,qdata
        z = wave/lambda0(ioni) - 1.0d0
        if ((z.ge.zlim(regioni,1)).AND.(z.le.zlim(regioni,2))) then
         pixi = pixi + 1
         quality(ioni,pixi) = qdata
         if (qdata.eq.1) pixq = pixq + 1
        end if
 07    CONTINUE

 08    CLOSE(unit=3)
       nqual(ioni) = pixq
       nmask(ioni) = ndata(ioni) - pixq

c     sanity check on the pixel book keeping with masking involved

       IF (pixi.ne.ndata(ioni)) THEN
        WRITE(6,*) ' ERROR(readdata)  ',ion_name(ioni)
        WRITE(6,*) ' Number of pixels in mask not equal to data'
        STOP ' MINFIT terminated'
       END IF

c     next ion transition

  03  continue

c     PHASE 2: two checks...  if either case below is true, then the
c     COMMON logical flag called "features" is set low
c
c     (1) check that each ion has data for the entire region some ions
c     may have data in most regions, but not in all regions.  we want to
c     use the data where we can, but not bomb out on a region if the
c     data DNE

      do 11 ioni=1,ions

c     did we obtain any data in the region?

       if (ndata(ioni).eq.0) then
        features(ioni)        = .false.
        featreg(ioni,regioni) = .false.
        write(SCREEN,603) ion_name(ioni)
        write(STDOUT,603) ion_name(ioni)
        goto 11 ! next ion
       end if

c     did we obtain data for the entire region, within 2 pixels of the
c     lower/upper zlimit (arbitrarily chosen buffer window, the constant
c     2.0 in dz?

c     check region lower limit

       wave = lambda(ioni,1)
       z    = wave/lambda0(ioni) - 1.0d0
       dz   = abs(2.0*(wave - lambda(ioni,2)))
       if (z-dz.gt.zlim(regioni,1)) then
        features(ioni)        = .false.
        featreg(ioni,regioni) = .false.
        write(SCREEN,604) ion_name(ioni),
     @                    nint((z-zlim(regioni,1))/dz)
        write(STDOUT,604) ion_name(ioni),
     @                    nint((z-zlim(regioni,1))/dz)
        goto 11 ! next ion
       end if
       
c     check region upper limit

       wave = lambda(ioni,ndata(ioni))
       z    = wave/lambda0(ioni) - 1.0d0
       dz   = abs(2.0*(wave - lambda(ioni,ndata(ioni)-1)))
       if (z+dz.lt.zlim(regioni,2)) then
        features(ioni)        = .false.
        featreg(ioni,regioni) = .false.
        write(SCREEN,605) ion_name(ioni),
     @                     nint((zlim(regioni,2)-z)/dz)
        write(STDOUT,605) ion_name(ioni),
     @                     nint((zlim(regioni,2)-z)/dz)
        goto 11 ! next ion
       end if

 11   continue

c     (2) now check that a given ion has significant features, otherwise
c     we are fitting noise and this reduces leverage in the chi square
c     landscape; only do this if we are called for purposed of fitting

      if (flag.eq.1) then

       write(SCREEN,600)
       write(STDOUT,600)
       j = 0

       do 15 ioni=1,ions
         ch_stat = '1 fitting: will provide VP comps and limits'
         i = findfeat(ioni)
         j = j + i
         if (i.eq.0) then 
           features(ioni)        = .false.
           featreg(ioni,regioni) = .false.
           ch_stat = '0 omiting: will provide limits only'
         end if 
         write(SCREEN,601) ioni,ion_name(ioni),spec_idx(ioni),
     @                     ndata(ioni),i,nmask(ioni),ch_stat
         write(STDOUT,601) ioni,ion_name(ioni),spec_idx(ioni),
     @                     ndata(ioni),i,nmask(ioni),ch_stat
 15    continue

c     do we have any detections at all?

       if (j.eq.0) then
        write(SCREEN,607) 
        write(STDOUT,607) 
        error = .true.
        return
       end if

      end if ! close flag.eq.1 loop

c     OK... if we removed an ion from the list for which there are other
c     species, then we do not have much of a chore with the book
c     keeping; we simply utilize the features logical flag in the
c     appropriate places in the code. If we have removed an entire
c     species, then things are more difficult...for this reason we keep
c     a parallel index called fitindex for each ion.  the fitting then
c     is done based upon the fitindex, but all communication is in terms
c     of spec_idx, as input by the user.  we also make the distinction
c     between "total" species and "species", which is the number of
c     species being fitted for that region

c     did we drop a species?  book keeping loop

      species = total
      k = 0
      do 21 j=1,total
       do 19 ioni=1,ions
        if (features(ioni).AND.(j.eq.spec_idx(ioni))) then
         k           = k + 1         
         fitindex(k) = j      
         GOTO 21
        end if
 19    continue
 21   continue
      if (k.gt.0) species = k

c     return
      return

c     error trapping

 999  write(SCREEN,*) ' NOT FOUND- data : ',ion_str(1:30)
      write(STDOUT,*) ' NOT FOUND- data : ',ion_str(1:30)
      STOP ' *** MINFIT(readdata): terminated ***'

c     formats

 600  format(1x,t3,'Idx',t10,'Ion/Tran',t20,'Ion',t25,'Npix',
     @          t31,'Ndet',t37,'Nmask',t45,'Status')
 601  format(1x,i4,4x,a10,i3,3i6,4x,a45)
 602  format(1x,'REGION:',i3,
     @          '   z=(',f8.6,',',f8.6,')',/)
 603  format(1x,'  --> no data in region for ',a10,/,
     @       1x,'    > dropping transition in this region')
 604  format(1x,'  --> partial data in region for ',a10,/,
     @       1x,'    > transition begins',i4,' pixels redward',/,
     @       1x,'    > dropping transition in this region')
 605  format(1x,'  --> partial data in region for ',a10,/,
     @       1x,'    > transition ends',i4,' pixels blueward',/,
     @       1x,'    > dropping transition in this region')
 607  format(1x,'   *** NO DETECTIONS IN THIS REGION ***',/,
     @       1x,'   no LSF, advancing to the next region')

      end

c     eof

