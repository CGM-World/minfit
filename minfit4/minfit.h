c................................................... minfit.h .................

c     this file contains all global variables in COMMON blocks
c
c     the COMMON blocks are broken by function and data type; the types are
c     (1) logical variables, (2) integer variables, (3) double precision
c     variables, and (4) character variables; by keeping the COMMON block
c     by data type, we ensure clean use of the stack and avoid bizarre 
c     Segmentation Faults
c
c     the "implicit none" insures that all declarations in a routine are 
c     explicitly declared


      implicit none


c     the following parameters set physical array sizes; they should be larger
c     than any anticipated logical sizes, but not too large that we use up 
c     machine memory and corrupt the data segments in the machine, which will
c     result in a Segmentation Fault or our of memory statement crash
c
c     maxspecies --- max number of species to be fit
c     maxions    --- max number of ion transitions to be fit
c     maxpix     --- max number of pixels for a given ion transition
c     maxlines   --- max number of Voigt profiles to be fit
c     maxcoeffs  --- max number of coefficients to be fit

c     integer MAXVEC is the vector length for the LSF engine dnls1e; it is 
c     the maximum number of functions being fit

c     integer LWA is the size of a working array, called WA, required by the
c     LSF engine dnls1e; WA is declared in the main routine and is passed via 
c     subroutine calls;  we use exploit it throughout when it is dormant


      integer            maxspecies,maxions,maxpix,maxlines,maxcoeffs,
     &                   maxregions,maxvec, lwa
      parameter          (maxspecies =                           7   ,
     @                    maxions    =                          20   ,  
     @                    maxpix     =                        2000   ,
     @                    maxlines   =                          50   ,
     @                    maxregions =                           7   ,
     @                    maxcoeffs  = 2*maxspecies*maxlines+maxlines, 
     @                    maxvec     =              maxions*maxpix   ,
     @                    lwa        = maxcoeffs*(maxvec+5)+maxvec   ) 
      integer            SCREEN,STDIN,STDOUT
      parameter          (SCREEN = 6, STDIN = 5, STDOUT = 91)


c     block for variables pertaining to the statistics of the LSF

      double precision   chisq,delchi2,oldchi2,tol
      COMMON/statblck/   chisq,delchi2,oldchi2,tol


c     blocks for variables pertaining to the atomic constants and general 
c     book keeping, and conditional decision making during the LSF


      logical            features,adjusting,doerrors,satflag
      logical            featreg,dodebug
      integer            species,ions,lines,ndata,nqual,regions,
     @                   spec_idx,fitindex,regioni,nprint,
     @                   nsatreg,nmask
      integer            total,lngth,flagion,adjflag,quality
      integer            idx_an,idx_ab,idx_az,sat_spec
      double precision   con1,con2,lambda0,zlim,zabs,spec_mass
      double precision   zsat,acol,bcol,ionpot,fosc
      double precision   f_nu,f_chisq,f_chisqnu,f_vardata,f_varfit,
     @                   f_npar
      character*30       profitfile
      character*80       ion_name,spec_name
      character*80       ionfile

      COMMON/featblck/features(maxions),satflag(maxspecies,maxlines),
     @                featreg(maxions,maxregions)
      COMMON/logiblck/adjusting,doerrors,dodebug
      COMMON/atmrblck/con1(maxions),con2(maxions),lambda0(maxions),
     @                zlim(maxregions,2),spec_mass(maxspecies),
     @                acol(maxspecies),bcol(maxspecies),
     @                ionpot(maxspecies),fosc(maxions)
      COMMON/redzblck/zsat(maxregions,2),zabs
      COMMON/atomchar/ion_name(maxions),spec_name(maxspecies),
     @                profitfile,ionfile
      COMMON/accntint/species,ions,lines,regions,regioni,nsatreg,
     @                nprint,ndata(maxions),nqual(maxions),
     @                lngth(maxions),flagion(maxions),
     @                adjflag,quality(maxions,maxpix),
     @                nmask(maxions)
      COMMON/idxblck/spec_idx(maxions),fitindex(maxions),
     @               idx_an(maxspecies,maxlines),
     @               idx_ab(maxlines),idx_az(maxlines),
     @               sat_spec(maxregions,2),total
      COMMON/statblk/f_nu(maxregions),f_npar(maxregions),
     @ 	             f_chisqnu(maxregions),f_chisq(maxregions),
     @		     f_vardata(maxregions),f_varfit(maxregions)


c     blocks for variables pertaining to the data arrays and parameter arrays

      logical            plotting
      integer            iturb,vpcomps
      double precision   nline,zline,bline,dnline,dzline,dbline
      double precision   Rfinal,Nfinal,bfinal,zfinal
      double precision   dNfinal,dbfinal,dzfinal
      double precision   lambda,data,sigma,wrkflx,fudge,N_sigma,
     @                   badness,confidence,fberr,bmin,bmax,Nmax
C      double precision   modflx
      character*80       ionlabel
c
      COMMON/pltblck/plotting
      COMMON/fitblck/lambda(maxions,maxpix),data(maxions,maxpix),
     @                sigma(maxions,maxpix),wrkflx(maxions,maxpix),  
     @                fudge,N_sigma,
     @                badness,confidence,fberr,bmin,bmax,Nmax
C     @                modflx(maxions,maxlines,maxpix)
      COMMON/turbblck/iturb,vpcomps
      COMMON/testblck/ionlabel(maxions)
      COMMON/parblck/nline(maxspecies,maxlines),zline(maxlines),
     @                bline(maxspecies,maxlines),
     @                dnline(maxspecies,maxlines),dzline(maxlines),
     @                dbline(maxspecies,maxlines),Rfinal(maxlines),
     @                Nfinal(maxspecies,maxlines),zfinal(maxlines),
     @                bfinal(maxspecies,maxlines),
     @                dNfinal(maxspecies,maxlines),dzfinal(maxlines),
     @                dbfinal(maxspecies,maxlines)



c     blocks pertaining to the convolution of the model with the instrumental
c     spread function (ISF)

c     ***********************************************************************
c                  ***** IMPORTANT NOTE FOR CONVOLUTION *****
c     ***********************************************************************

c     Parameter MAXCON defines the working space for the FFTs, it must 
c     be a power of 2.

c     NOTE: there is a DATA statement in routine initconv in the instruments.f
c     module that contains a list of powers of 2.  It goes up to 16384.  
c     There is parameter NMAX in the routine convlv in fft.f module that 
c     defines the size of the complex variables for the FFT.  It is also 
c     set to 16384.  

c     If you increase MAXCON to > 16384 (by a power of 2) then you must edit 
c     initconv in instruments.f so that the data statements to go up to the 
c     corresponding power of 2 you are increasing MAXCON to in order to be 
c     consistent.  You must also edit the paramtter NMAX in the convlv 
c     routine of fft.f to equal MAXCON. This is all to ensure that the value 
c     of MAXPIX and MAXCON should meet the condition 
c                       MAXCON > MAXPIX*RESFAC 
c     or there will not be enough physical space to perform the convolution;
c     parameter RESFAC is set in the minfit.par file; MINFIT checks
c     this condition each time it sets up the convolution arrays

      logical           convolving
      integer           maxcon
      parameter         (maxcon = 16384)
      integer           nresponse,ncondat,nfft
      double precision  response,convdata,pixpres
      double precision  R_fac,slit,conwindo,resfac,profile,hdv

      COMMON/convlblck/convolving
      COMMON/convrblck/response(maxpix),convdata(maxcon),pixpres
      COMMON/conviblck/nresponse,ncondat,nfft
      COMMON/instblck/R_fac,slit,conwindo,resfac,profile,hdv

c............................................... end minfit.h .................
