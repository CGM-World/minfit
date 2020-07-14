c.........................................................................
c

      SUBROUTINE          packa(a,siga,n,m)

c      
c     pack the coefficent vector
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      include           'const.dek'
      integer           i,k,n,m,speciesi,fidx,linei,in,ib,iz
      double precision  a(maxcoeffs),siga(maxcoeffs)


c     check array bounds; if we are out of bounds terminate
     
      if (lines.gt.maxlines) stop ' MAXLINES exceeded '
      if (species.gt.maxspecies) stop ' MAXSPECIES exceeded '
      if (regions.gt.maxregions) stop ' MAXREGIONS exceeded '

c     zero the coefficient arrays

      do 01 n=1,maxcoeffs
       a(n)   = 0.0d0
       siga(n) = 0.0d0
 01   continue

c     initialize this region

      delchi2 = 0.0d0
      oldchi2 = 0.0d0
      n       = 0
      m       = 0

c     set up the index arrays for the fitting corfficients; this routine
c     stuffs the common integer variables IDX_AN, IDX_AB, and IDX_AZ,
c     which are indexed by species number and line number

      CALL setidx(n)

c     fold in the columns densities

      do 11 k=1,species
        fidx = fitindex(k)
        do 15 i=1,lines
         in    = idx_an(k,i)
         a(in) = 10.0d0**(nline(fidx,i)-13.0d0)
 15     continue
 11   continue

c     fold in the b parameters and the line redshifts
c     we assume that the first ion in the list has representative b pars

      k = 1 ! species number
      do 17 i=1,lines
       ib = idx_ab(i)
       iz = idx_az(i)
       a(ib) = spec_mass(k)*1.0d10*(bline(k,i))**2 
     &       / (2.0d0*kerg)
       a(iz) = zline(i)
 17   continue

c     account for tied saturated lines

      k = 0
      CALL jukea(a,n,k)

c     set the number of functions; check array bounds, if we are out of
c     bounds terminate; only use unmasked pixels for the number of
c     functions

      do 19 k = 1,ions
       if (ndata(k).gt.maxpix) stop ' MAXPIX exceeded'
       if (features(k)) m = m + nqual(k)
 19   continue
      if (m.gt.maxvec) stop ' MAXVEC exceeded '

c     return

      return
      end

c

c.........................................................................
c

      SUBROUTINE          setidx(n)
c
c     the folding <--> unfolding formulae are
c
c     define:
c 
c        k       = species number   ( 1 <= k <= species )
c        i       = cloud number     ( 1 <= i <= lines )
c        lines   = number of current clouds being fit
c        species = number of elemental species currently being fit
c     
c
c
c     OPTION 1. 
c     - column densities independent between species 
c     - temperatures independent between species for clouds with tied redshifts
c     - redshifts tied for clouds 
c
c     idx_an(k,i) = lines*(k-1) + i
c     idx_ab(k,i) = lines*(species+k-1) + i
c     idx_az(i)   = 2*species*lines + i
c
c
c     OPTION 2. 
c     - column densities independent between species 
c     - temperatures tied between clouds with tied redshifts 
c     - redshifts tied for clouds 
c
c     idx_an(k,i) =  lines*(k-1) + i
c     idx_ab(i)   =  lines*species + i
c     idx_az(i)   =  lines*(species+1) + i
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      include           'const.dek'

      integer            n,i,j,k,pass

      

c     intitilize the counter

      j = 0

c     allow for dynamic indexing, formalism is that the index is cloud
c     ordered first and species ordered second

      do 07 k=1,species
       do 05 i=1,lines
          j = j + 1
          idx_an(k,i) = j
 05    continue
 07   continue

c     stuff the b parameter indices

      do 09 i=1,lines
         j = j + 1
         idx_ab(i) = j
 09   continue

c     stuff the redshift indices

      do 11 i=1,lines
         j = j + 1
         idx_az(i) = j
 11   continue

      n = j

      if (n.gt.maxcoeffs) stop ' MAXCOEFFS exceeded '


c     return

      RETURN

      END

c.........................................................................
c

      SUBROUTINE          unpacka(a,siga)
c
c     unpack the coefficent vector
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      include           'const.dek'

      integer           i,jj,k,fidx,ioni,in,ib,iz,Naction
      double precision  a(maxcoeffs),siga(maxcoeffs),Tp,satNscale


c     set Naction for the call to satNscale

      Naction = 1 ! log10 units

c     unfold the columns, b parameters, and redshifts; we have clean up
c     all the messy values (insignificant column densities and b
c     parameters) when we find the N limits on the ions, so here we can
c     perform a clean unpacking (big improvement over older versions)

c     when unpacking column densities; avoid logs of negative numbers
C     (old bpar scaling) Tp = ionpot(j)*a(ib)/ionpot(1)

      do 13 i=1,lines
       ib = idx_ab(i)
       do 11 k=1,species
         fidx = fitindex(k)

c     unpack column densities, flag negative columns

         in = idx_an(k,i)
         IF (a(in).le.0.0d0) then
          nline(fidx,i)  = sign(1.0d0,a(in))*(log10(abs(a(in)))+13.0d0)
          dnline(fidx,i) = -1.0d0
         ELSE
          if (satflag(k,i)) then
          nline(fidx,i) = satNscale(a(in),acol(fidx),bcol(fidx),Naction)
C          nline(fidx,i)  = acol(fidx)*(13.+log10(a(in)))+bcol(fidx)
          else
           nline(fidx,i)  = 13.0d0+log10(a(in))
          end if
          dnline(fidx,i) = 0.4343*siga(in)/a(in)
          if (dnline(fidx,i).gt.999.99d0) dnline(fidx,i) = 999.99 ! avoid bad format
         END IF

c     unpack b parameters

         Tp             = a(ib)
         bline(fidx,i)  = sqrt(2.0d0*kerg*Tp/spec_mass(fidx))
         dbline(fidx,i) = kerg*siga(ib)/(spec_mass(fidx)*bline(fidx,i))
         bline(fidx,i)  = 1.0d-5*bline(fidx,i)    ! convert to km/s
         dbline(fidx,i) = 1.0d-5*dbline(fidx,i)   ! convert to km/s
         if (dbline(fidx,i).gt.999.99d0) dbline(fidx,i) = 999.99 ! avoid bad format

 11    continue

c     unpack redshifts

       iz        = idx_az(i)
       zline(i)  = a(iz)
       dzline(i) = siga(iz)

 13   continue

c     for transitions that have no FITINDEX (meaning no detectable
c     features) we want to have the best bpar for computing limits, not
c     the one that was input by the user, so we repack the BLINE array
c     of the transitions with no features before we compute limits; the
c     NLINE array for these transitions will be replaced in getNlimits,
c     so we do not modify it; 

c     note, whereas SPEC_MASS, NLINE, etc are usually reference by
c     FITINDEX, here we use the counter K, which loops through the total
c     species; this is because we are accessing directly (using K)
c     rather than indirectly (which uses FITINDEX)

      do 25 k=1,total  
       do 21 ioni=1,ions
         fidx  = spec_idx(ioni)
         IF (fidx.eq.k) then 
          IF (features(ioni)) then ! skip to next species 
          GOTO 25
          END IF
         END IF
 21    continue
c     if we make it here replace
       do 23 i=1,lines
         ib          = idx_ab(i)
         Tp          = a(ib)
         bline(k,i)  = sqrt(2.0d0*kerg*Tp/spec_mass(k))
         bline(k,i)  = 1.0d-5*bline(k,i)    ! convert to km/s
         dbline(k,i) = -1.0d0
 23    continue
 25   continue

c     sometimes lines get swapped in redshift space; sort them for good
c     measure

      CALL sortlines

C     TEST UPACKING
C
C      do 53 i=1,lines
C       write(SCREEN,602) i,zline(i),(nline(j,i),dnline(j,i),
C     @               bline(j,i),dbline(j,i),j=1,total)
C 53    continue
C 602  format(1x,i3,2x,f8.6,10(2x,f6.2,1x,f6.2))

c     return

      return
      end

c
c.........................................................................
c

      SUBROUTINE          zeroall(fvec,wa,iw)

c
c     zero the working arrays before the call to LSF the monster
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include          'minfit.h'
      integer           i,iw(maxcoeffs)       
      double precision  fvec(maxvec),wa(lwa)


c  zero the integer working array
      
      do 01 i=1,maxcoeffs
       iw(i) = 0
 01   continue

c     zero the function vector

      do 02 i=1,maxvec
       fvec(i) = 0.0d0
 02   continue

c     zero the massive working array (we use this array for a multitude
c     of uses throughout the program); it is also required for the
c     SLATEC libary routine dnls1e, which is the monster fitter (the
c     engine under the hood)

      do 03 i=1,lwa
       wa(i) = 0.0d0
 03   continue

c     return

      return
      end

c  EOF
