c..............................................................................
c

      SUBROUTINE          fitregion(fcn,iopt,m,n,a,siga,fvec,toler,
     @                              io_iter,info,iw,wa,ilwa)
 
c
c     this routine performs the first pass LSF of the input model.  it
c     works one region at a time, called from the main loop in driver
c     MINFIT that steps through REGIONI one at a time.  The current
c     region is defined by global variable REGIONI
c     
c     this routine is divided into four main parts: 
c
c     PART 1: SET UP FOR THE INITIAL FIT ITERATIONS
c
c     sets up the LSF, reads in the data (routine READATA) and initial
c     model (routine GETLINES), and book keeps the VP parameter arrays
c     (routine PACKA) and sets up the instrumental convolution arrays
c     for the FFT convolution of the instrumental spread function
c     (routine INITCONV) for the current REGIONI
c
c     PART 2: THE INITIAL FIT ITERATION LOOP 
c
c     this is the first pass LSF stage of the current REGIONI; after the
c     working arrays are nulled (routine ZEROALL), the LSF is performed
c     on the initial user input model (routine DNLS1E); if integer flag
c     ADJFLAG is 1 or 2, then the resulting LSF model is checked for
c     insignificant and unconstrained lines (routine CHECKA); an
c     insignificant line is one for which the FLAGION or primary ion has
c     a negative column density; the primary ion is the one that has the
c     been designated as IDX=1 in the "ions.table" file.  If the VP
c     component of the FLAGION or primary ion is negative of highly
c     uncertain then it is unceremoniously deleted from the model (see
c     routine CHECKA).  If a line is deleted from the model then the LSF
c     is repeated.  In this routine, this process is repeated until all
c     lines are significant and/or constrained
c
c     PART 3: ERRORS AND STATISTICAL TESTING OF LINES
c
c     this part is executed only if the user has defined logicals
c     DOERRORS=.TRUE. and ADJUSTING=.TRUE. (set by the user in the
c     "minfit4.par" file).  if DOERRORS=.TRUE. then the uncertainties,
c     SIGA, in the VP parameters are computed (routine ERRORS).  if
c     ADJUSTING=.TRUE. then the routine TESTLINES is called and the LSF
c     is further refined.  See the module "adjustfit.f" for the
c     description of routine TESTLINES. The LSF model returned by
c     routine TESTLINES is the final LSF for the current REGIONI
c
c     PART 4: CLOSE UP SHOP FOR CURRENT REGIONI
c
c     communicates the final model for the current REGIONI, performs
c     some checks to see if the LSF has any significant residuals and
c     communicates; unfolds the A and SIGA arrays into the intermediate
c     results arrays (routine UNPACKA), and computes the column density
c     limits for unconstrained lines (routine GETNLIMITS), and then
c     stores the these final results in final storage arrays (routine
c     STOREFITS).  returns to the main loop in driver MINFIT
c
c     further explanations are commented throughout this routine
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include          'minfit.h'
      logical           error,nukeline
      external          fcn
      integer           i,iopt,io_iter,ioni,findfeat,loop
      integer           m,mm,n,info,iw(maxcoeffs),ilwa,nknt
      integer           blwa(maxlines)
      double precision  a(maxcoeffs),siga(maxcoeffs),
     @                  fvec(maxvec),wa(lwa),toler


c     intiailization of local logical

      error    = .false.

c     ---------------------------------------------
c     PART 1: SET UP FOR THE INITIAL FIT ITERATIONS
c     ---------------------------------------------

c     read in the data, set fitting flags and adjust the current number
c     of species to be fit based upon a detection screening for each
c     transition; also reads in the mask spectrum for each transition;
c     if there is an error, return with error set high

      CALL readdata(1,error)
      IF (error) RETURN

c     get the FLAGION for REGIONI; the flag ion, if defined, is used for
c     determining if lines in the model become insignificant; if a flag
c     ion is set, then we use this to determine if a given VP component
c     is significant simply by checking if its column desnity goes
c     negative... if so, we remove the VP component in routine CHECKA
c     and restart the LSF. also, if FLAGION is specified then in routine
c     TESTLINES only the errors in the flag ion parameters will be used
c     for constraining the VP model.  if a FLAGION is not specified for
c     the current REGIONI, then the first species in the ions.table
c     input list is used for CHECKA and the errors for the components
c     for all species are invoked for constraining the LSF in routine
c     TESTLINES

      CALL getflagions

c     get the user intial model VP parameters for the current REGIONI;
c     this reads in the initial VP model parameters; it retains only
c     those that are in the curent REGIONI

      CALL getlines

c     pack the VP parameter array, A, and the uncertainty array, SIGA;
c     initilaize the number of LSF parameters, N, and the number of LSF
c     functions, M, which is basically the total number of pixels for
c     all transitions (minus any masked pixels)

      CALL packa(a,siga,n,m)

c     if CONVOLVING=.TRUE. then the VP model is convolved with the
c     instrumental spread function, thus we set up the convolution
c     arrays for the FFTs; further details provided in comments of
c     routine INITCONV; the use of MM in place of M here is to account
c     for the full length of all the pixels included in the LSF, M book
c     keeps the number of unmasked pixels (fitting functions) in the
c     LSF, whereas MM book keeps the total number of pixels in the data.

      mm = 0
      if (convolving) then
       DO 99 i=1,ions
        IF (features(i)) mm = mm + ndata(i) 
 99    CONTINUE
       CALL initconv(mm)
      end if
c
c     END PART 1
c

      WRITE(SCREEN,*) ' '
      WRITE(STDOUT,*) ' '

c     --------------------------------------
c     PART 2: THE INITIAL FIT ITERATION LOOP 
c     --------------------------------------

c     obtain the LSF; upon completion we now have an LSF model, but it
c     may be a bit messy, on output, INFO contains the convergence
c     status, which we communicate, if INFO=0 then an error has
c     occured... attempt to exit gracefully

c     first we communicate the loop iteration, we then null the M
c     elements of FVEC (the function being minimuzed, see routine FCN),
c     and the working arrays used by the LSF engine; after the LSF, we
c     then compute the fit statistics and communicate (in routine
c     FITSTATS), the LSF convergence is communicated (routine INFORM),
c     and the model parameters are communicated (routine MODELCOMM).
c     We then perform the checks on the LSF parameters (see comments
c     below)

      loop     = 1       ! initialize the loop counter

c     code line 01 is the loop starting point; initalize logicals flags
c     and NKNT, where NKNT provides the number of insignificant lines
c     removed from the LSF model before iterating the loop on the next
c     model

 01   nukeline = .false.
      error    = .false.
      nknt     = 0

c     communicate the loop iteration

      WRITE(SCREEN,600) loop
      WRITE(STDOUT,600) loop
 600  FORMAT(1x,/,
     &       1x,'+++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     &       1x,'*** FITREGION **** iteration loop =',i4,/,
     &       1x,'+++++++++++++++++++++++++++++++++++++++++++++++++++')

c     initilize and perform LSF 

      CALL zeroall(fvec,wa,iw)
      CALL dnls1e(fcn,iopt,m,n,a,fvec,toler,
     &            io_iter,info,iw,wa,ilwa)

c     communicate the LSF result

      IF (plotting) call plt
      CALL fitstats(a,n,m,1)
      CALL inform(info,error)
      CALL modelcomm(a,n,1) ! 1 = an LSF result
      IF (error) RETURN

c     this is the CHECKA segment of the loop; if we remove a line or
c     several lines in routine CHECKA then NUKELINE=.TRUE. is returned
c     (the array book keeping is done by routine JUKEA called in routine
c     CHECKA)

      IF (adjflag.ne.0) THEN
        CALL checka(a,n,nukeline,blwa,nknt)
        WRITE(SCREEN,*) ' '
        WRITE(STDOUT,*) ' '
      END IF

c     depending upon the values of NUKELINE we communicate,
c     book keep the iteration number in LOOP and jump to the loop
c     starting point (code line 01)

      IF (nukeline) then
        WRITE(SCREEN,*) '--- RESTARTING FIT W/O NUKED LINE(S) ---'
        WRITE(STDOUT,*) '--- RESTARTING FIT W/O NUKED LINE(S) ---'
        WRITE(SCREEN,*) ' '
        WRITE(STDOUT,*) ' '
        CALL modelcomm(a,n,0) 
        loop = loop + 1
        GOTO 01               ! go back to the loop start
      ELSE
        WRITE(SCREEN,*) '--- FITREGION COMPLETE ---'
        WRITE(STDOUT,*) '--- FITREGION COMPLETE ---'
        WRITE(SCREEN,*) ' '
        WRITE(STDOUT,*) ' '
      END IF

c
c     END PART 2
c

c     -----------------------------------------------
c     PART 3: ERRORS AND STATISTICAL TESTING OF LINES
c     -----------------------------------------------

c     do we compute the errors?; this routine can be somewhat time
c     intensive and a bit of a bottleneck.  logical flag ERROR is set
c     .TRUE. if the curvature matrix was singular, in which case there
c     will no erros and no rejection criterion

 02   IF (doerrors) CALL errors(a,siga,n,m,error)

c     aremed with uncertainties in the VP parameters we now determine
c     the statistical significance of these lines; if the user has set
c     logical ADJUSTING=.TRUE. then we call routine TESTLINES, which
c     runs iterative F-tests on the targeted lines based upon their
c     fractional errors; the routine is quite intensive and is the true
c     heart and soul of MINFIT; upon returning we are done with the
c     current REGIONI
  
      IF (adjusting) THEN
       IF (.not.error) CALL testlines(fcn,iopt,m,n,a,siga,fvec,
     &                      toler,io_iter,info,iw,wa,ilwa)
      END IF 

c
c     END PART 3
c

c     -------------------------------------
c     PART 4: CLOSE UP SHOP FOR THIS REGION
c     -------------------------------------


c     communicate this is it

      WRITE(SCREEN,*) ' '
      WRITE(SCREEN,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(SCREEN,*) '+                                              +'
      WRITE(SCREEN,*) '+      FITTING OF THIS REGION IS COMPLETED     +'
      WRITE(SCREEN,*) '+                                              +'
      WRITE(SCREEN,*) '++++++++++++++++++++++++++++++++++++++++++++++++'

      WRITE(STDOUT,*) ' '
      WRITE(STDOUT,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE(STDOUT,*) '+                                              +'
      WRITE(STDOUT,*) '+      FITTING OF THIS REGION IS COMPLETED     +'
      WRITE(STDOUT,*) '+                                              +'
      WRITE(STDOUT,*) '++++++++++++++++++++++++++++++++++++++++++++++++'


c     for good measure, inform use if there are any left-over
c     significant features in the residuals; use the FINDFEAT function;
c     the 3 is a hardwire that should be set to the number of pixels per
c     resolution element

      DO 11 ioni=1,ions
        i = 0
        IF (features(ioni)) i = findfeat(ioni)
        IF (i.ge.3) THEN
          WRITE(SCREEN,601) ion_name(ioni),i,N_sigma
          WRITE(STDOUT,601) ion_name(ioni),i,N_sigma
        END IF
 11   CONTINUE

 601  FORMAT(1x,'*WARNING* - a check of the residuals for ',a30,/,
     @  1x,' reveals ',i3,' pixels with ',f4.1,' sigma equivalent',/,
     @  1x,' widths.  Perhaps an additional component is needed?')

      WRITE(SCREEN,*) ' '
      WRITE(STDOUT,*) ' '

c     make sure the last view is the last model!

       IF (plotting) THEN
         CALL model(a,m)
         DO 100 i=1,lines  ! update tick locations for plotting
          zline(i) = a(idx_az(i))
 100     CONTINUE
         CALL plt
       END IF

c     unfold the coefficent vector back into the NLINE, BLINE, and ZLINE
c     arrays

      CALL unpacka(a,siga)

c     get the column density limits in the bad spots; we have not set up
c     communication if an error occurs (NEED TO ADD)

      if (doerrors) CALL getNlimits(error)

c     store the VP parameters before moving to the next region

      CALL storefits

c
c     END PART 4
c

c     congratulations- I hope

      RETURN
      END

c
c..............................................................................
c

      LOGICAL FUNCTION ifnan(x)

c  
c     cheap way to check if the value is an IEEE floating point
c     exception
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      logical             error
      double precision    x,y,value2
      character*20        ch_x


      error = .false.
      ifnan = .false.

      WRITE(ch_x,'(1pe12.5)') x

      y = value2(ch_x,error)      

      IF (error) ifnan = .true.

      RETURN
      END


c..............................................................................
c..
      DOUBLE PRECISION FUNCTION value2(string,err)
      implicit none
c..
c..this routine takes the character string number and converts 
c..it to a real. if trouble is encountered during the conversion, this
c..routine returns the logical err as .true.
c..
c..required routines: none
c..
c..date code: 28may90    Frank Timmes
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c..
c..declare
      logical          pflag,err
      character*(*)    string
      character*1      plus,minus,decmal,blank,se,sd,se1,sd1
      integer          noblnk,long,ipoint,power,psign,iten,j,z,i
      double precision x,thesign,factor,rten,temp
      parameter        (plus = '+'  , minus = '-' , decmal = '.'   ,
     1                  blank = ' ' , se = 'e'    , sd = 'd'       ,
     2                  se1 = 'E'   , sd1 = 'D'   , rten =  10.0,
     3                  iten = 10                                   )
c..
c..initialize
      err    =  .false.
      x      =  0.0d0
      thesign   =  1.0d0
      factor =  rten
      pflag  =  .false.
      noblnk =  0
      power  =  0
      psign  =  1
      long   =  len(string)
c..
c..remove leading blanks and the sign
      do 10 z = 1,7
       noblnk = noblnk + 1
       if ( string(noblnk:noblnk) .eq. blank) then
        if (noblnk .gt. 6 ) goto 1000
       else
        if (string(noblnk:noblnk) .eq. plus) then
         noblnk = noblnk + 1
        else if (string(noblnk:noblnk) .eq. minus) then
         noblnk = noblnk + 1
         thesign =  -1.0d0
        end if
        go to 100
       end if
10    continue
c..
c..main number conversion loop
100   do 200 i = noblnk,long
       ipoint = i + 1
c..
c..if blank character then we are done
       if ( string(i:i) .eq. blank ) then
        x = x * thesign
        value2 = x 
        return
c..
c..
c..if it is an exponent process it
       else if (string(i:i).eq.se  .or. string(i:i).eq.sd .or.
     1          string(i:i).eq.se1 .or. string(i:i).eq.sd1   ) then
        if (x .eq. 0.0d0 .and. ipoint.eq.2)     x = 1.0d0
        if (thesign .eq. -1.0d0 .and. ipoint.eq.3) x = 1.0d0
        if (string(ipoint:ipoint) .eq. plus) ipoint = ipoint + 1
        if (string(ipoint:ipoint) .eq. minus) then
         ipoint = ipoint + 1
         psign = -1
        end if
        do 150 z = ipoint,long
         if (string(z:z) .eq. blank)  then
          x = thesign * x * rten**(power*psign)
          value2 = x
          return
         else
          j = ichar(string(z:z)) - 48
          if ( (j.lt.0) .or. (j.gt.9) ) go to 1000
          power= (power * iten)  + j
         end if
150     continue
c..
c..if it is a number process it
       else if (string(i:i) .ne. decmal) then
        j = ichar(string(i:i)) - 48
        if ( (j.lt.0) .or. (j.gt.9) ) go to 1000
        if (.not.(pflag) ) then
         x = (x*rten) + j
        else
         temp   = j
         x      = x + (temp/factor)
         factor = factor * rten
         go to 200
        end if
c..
c..must be a decimal point
       else
        if (pflag) go to 1000
        pflag = .true.
       end if
200   continue
c..
c..errors
1000  err = .true.
      return
      end


c..............................................................................
c eof
