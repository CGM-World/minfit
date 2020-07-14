c..............................................................................
c

      SUBROUTINE        testlines(fcn,iopt,m,n,a,siga,fvec,toler,
     @                           io_iter,info,iw,wa,ilwa)

c
c     the point of MINFIT is to obtain a statistically meaningful
c     chi-square minimized model of the spectra... statistically
c     meaningful means that we use as few parameters as possible for the
c     given umber of degrees of freedom.  Thus, we loop over the LSF fit
c     parameters and flag those with total fractional errors > BADNESS,
c     fracional errors in the b parameter > FBERR, and/or b parameters
c     less than BMIN, where BADNESS, FBERR, and BMIN are user defined in
c     the minfit.par file.
c 
c     If we have some "bad lines", we temporarily remove the "bad" line,
c     being sure to store all the intial parameters for safe keeping; we
c     then re-run the fitting of the full model without the bad line.
c
c     Once we have the new fit (minus the bad line), we compute the
c     fitting statistics and perform an F-test to see if the original
c     line was statistically significant.  If it was, we put it back and
c     move to the next bad line in the list.  If it was not, then we
c     retain the new model and run the tests on this new model.  The
c     process is iterated until there are no remaining insignficant/bad
c     lines.
c     
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include          'minfit.h'
      include          'const.dek'

      logical           error,nukeline,replace,fixline
      logical           satflagsave(maxspecies,maxlines)
      external          fcn
      integer           i,j,k,kk,m,n,linei,badline(maxlines)
      integer           iopt,io_iter,n_save,lines_save,loop,nunique
      integer           info,iw(maxcoeffs),ilwa,maxloops,nknt
      integer           idxnsave(maxspecies,maxlines)
      integer           idxbsave(maxlines),idxzsave(maxlines)
      integer           blwa(maxlines)
      double precision  asave(maxcoeffs,2)
      double precision  a(maxcoeffs),siga(maxcoeffs),wa(lwa)
      double precision  fvec(maxvec),toler,fracerr
      double precision  ftest,fprob,var1,nu1,var2,nu2

c     this routine is structured into a outer and inner loop.  

c     The outer loop is the main driver that takes a given model and
c     tests it.  The first step is the identification of badlines.
c     These badlines are stored and then ordered by line index number.
c     They are then scouted for which lines are the most offensive
c     (dominant transitions with worse errors) and ranked in order of
c     "most critical" to "least critical". Then in this order, one by
c     one, the badlines are removed from the fit and the new model is
c     tested.  If FLAGION is not set by the user, then all transitions
c     are checked for badness, but the routine is species and line index
c     driven.

c     the inner loop iterates on the test model until all lines in the
c     test model are fit with reasonable values

C
C     BEGIN OUTER LOOP (code line 01)
C

c     the loop to test suspect bad lines starts here at 01; there is
c     some intialization required for the start of each looping


 01   error       = .false.
      nukeline    = .false.
      fixline     = .false.
      nknt        = 0

C
C     OUTER LOOP STEP 1: IDENTIFY AND PRIORITY SORT THE BADLINES FOR THE
C                        CURRENT MODEL
C

c     announce that we are testing the current model

      DO 1000 kk=1,2
       IF (kk.eq.1) i = SCREEN
       IF (kk.eq.2) i = STDOUT
       WRITE(i,401) badness,confidence,bmin,bmax,Nmax
 1000 CONTINUE

 401  FORMAT(1x,/,
     &       1x,'++++++++++++++++++++++++++++++++++++++++++++++++',/,
     &       1x,'+                                              +',/, 
     &       1x,'+        SIGNIFICANCE TESTING OF MODEL         +',/, 
     &       1x,'+                                              +',/, 
     &       1x,'+        fractional error limit = ',f7.3,'      +',/,
     &       1x,'+        confidence level       = ',f7.3,'      +',/,
     &       1x,'+        minimum b parameter    = ',f7.3,'      +',/,
     &       1x,'+        maximum b parameter    = ',f7.3,'      +',/,
     &       1x,'+        maximum N [log10]      = ',f7.3,'      +',/,
     &       1x,'+                                              +',/, 
     &       1x,'++++++++++++++++++++++++++++++++++++++++++++++++',/,/,
     &       1x,'--- identifying and prioritizing lines with',/,
     &       1x,'    poorly constrained VP components ...',/)


c     determine the number the bad components, then determine the unique
c     line numbers (some bad components will be the same line number but
c     for different species), then prioritize the unique line numbers
c     for F-testing: on return from routine TAGLINES, NUNIQUE contains
c     the number of unique lines and array BADLINE contains the line
c     numbers in priority order

      CALL taglines(m,n,a,siga,iw,wa,ilwa,nunique,badline)

c     if all the lines are good, then we are done and return to the
c     calling routine, FITREGION

      IF (nunique.eq.0) RETURN

c     otherwise, continue on with the pain and suffering, first
c     communicate the priority order that the lines for the current
c     model will be tests; of course, at any point if a line is deemed
c     insignifican, then the current BADLINE list is no longer relevant
c     - we will have to test the new model and re-establish the BADLINE
c     list for it

      WRITE(SCREEN,*) '--- LINE PRIORITY TESTING ORDER ...'
      WRITE(STDOUT,*) '--- LINE PRIORITY TESTING ORDER ...'

      WRITE(SCREEN,'(60i4)') (badline(i),i=1,nunique)
      WRITE(STDOUT,'(60i4)') (badline(i),i=1,nunique)
      WRITE(SCREEN,*) ' '
      WRITE(STDOUT,*) ' '

C
C     OUTER LOOP STEP 2: PERFORM SIGNIFICANCE TESTS
C     (this step includes the inner loop)
C

c     loop through the unique prioritized badlines one at a time

      DO 999 j=1,nunique

c     communicate the current line for which we are significance testing

       WRITE(SCREEN,501) badline(j)
       WRITE(STDOUT,501) badline(j)

 501  format(1x,'*** CHECKING SIGNIFICANCE OF LINE =',i4,/)

c     save the current model (prior to removing the target badline; this
c     is a book keeping step; also save the previous fit statistics, the
c     A and ASIG vectors, array index pointers, and the length of the
c     logicel size of these vectors (N and LINES)

c     save the current model fit statistics for the F-Test

       CALL fitstats(a,n,m,-1) ! do not communicate
       nu1   = real(m-n)
       var1  = chisq/nu1

c     save the current model fitting parameters

       DO 21 i=1,n
        asave(i,1) = a(i)
        asave(i,2) = siga(i)
 21    CONTINUE

c     save the current model indices and flags

       DO 25 i=1,lines
        idxbsave(i) = idx_ab(i)
        idxzsave(i) = idx_az(i)
        DO 23 k=1,species
         idxnsave(k,i)    = idx_an(k,i)
         satflagsave(k,i) = satflag(k,i)
 23     CONTINUE
 25    CONTINUE

       n_save      = n
       lines_save  = lines

c     routine JUKEA sets up the new test model by creating the new A,
c     SIGA, IDX_AN, IDX_AB, IDX_AZ, SATFLAG arrays; the bad line being
c     removed, N and LINES is decremented (within routine JUKEA)

       CALL jukea(a,n,badline(j))

c
c     BEGIN INNER LOOP (code line 02)
c

c     We now enter the inner loop (starts at code line 02).  initialize
c     logical flags, zero everything and obtain the new LFS; code line
c     02 will be restarted if the LSF of the test model after lines that
c     are adjusted in routines CHECKA and/or FIXA (see comments below)
  
c     using MAXLOOPS, we limit the number of iterations from CHECKA and
c     FIXA to 3

       loop     = 1
       maxloops = 3       ! should move this to minfit.par read in

 02    fixline  = .false.
       nukeline = .false.
       error    = .false.
       nknt     = 0

      WRITE(SCREEN,600) loop
      WRITE(STDOUT,600) loop
 600  FORMAT(1x,/,
     &       1x'+++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     &       1x,'*** ADJUSTFIT **** iteration loop =',i4,/,
     &       1x,'+++++++++++++++++++++++++++++++++++++++++++++++++++')

c     initialize and fit the test model

       CALL zeroall(fvec,wa,iw)     
       CALL dnls1e(fcn,iopt,m,n,a,fvec,toler,
     @             io_iter,info,iw,wa,ilwa)  

c     communicate test model fit results

       IF (plotting) CALL plt
       CALL fitstats(a,n,m,1)
       CALL inform(info,error)               
       CALL modelcomm(a,n,1) ! 1 = an LSF result
                                             
C     ERROR? (NB. INCOMPLETE... FIX CODE HERE) did we screw up or go
C     into index space where no minimizer has gone before? I WANT TO
C     EVENTUALY REMOVE THE STOP COMMAND AND HAVE MINFIT RESET 

       IF ((info.eq.0).OR.(error)) THEN
        WRITE(SCREEN,*) 
     @   ' PROBLEM - parameters incorrectly posed for LSF...',
     @   '           restoring previous fit and terminating',
     @   '           all further LSF adjustments of this region.'
        WRITE(STDOUT,*) 
     @   ' PROBLEM - parameters incorrectly posed for LSF...',
     @   '           restoring previous fit and terminating',
     @   '           all further LSF adjustments of this region.'
        error = .true.
        STOP
       END IF

c     be sure to attend to poorly constrained lines for the test fit.

c     logic for restarting inner loop (code line 02): if a line is
c     removed in routine CHECKA, then logical flag NUKELINE is set high.
c     CHECKA calls routine JUKEA, so the model vector is once again
c     modified (the first iteration of the test model is lost).  if an
c     adjustment has been made in FIXA, then logical flag FIXLINE is set
c     high.

c     If NUKELINE is high, then we have a totally new model and we reset
c     the looping counter at the start again and continue on.  If
c     FIXLINE is high, then we increment the loop counter and minimize
c     the model again in the hopes that the adjustment improves things;
c     if both NUKELINE and FIXLINE are low, then we are ready for the
c     F-TEST (happy day!). if the iteration of the model reaches
c     MAXLOOPS then we go ahead and do the F-test because the
c     adjustments in FIXA are not helping the convergence

        IF (loop.eq.maxloops) THEN
          WRITE(SCREEN,599)
          WRITE(STDOUT,599)
          GOTO 03  ! do the F-Test, we are floundering
        END IF
 599    FORMAT(1x,/,' *** MAXIMUM ITERATIONS REACHED - MOVING ON ***',/)

        CALL checka(a,n,nukeline,blwa,nknt)  ! remove <0 column in flagion
        CALL fixa(a,n,fixline,blwa,nknt)     ! check all a values

        WRITE(SCREEN,*) ' '
        WRITE(STDOUT,*) ' '

        IF (fixline.OR.nukeline) THEN

          IF (fixline) then
            WRITE(SCREEN,*) '--- RESTARTING FIT W/ FIXED PARAMETERS ---'
            WRITE(STDOUT,*) '--- RESTARTING FIT W/ FIXED PARAMETERS ---'
            IF (.not.nukeline) then
              WRITE(SCREEN,*) ' -> stepping loop iterations '
              WRITE(STDOUT,*) ' -> stepping loop iterations '
              WRITE(SCREEN,*) ' '
              WRITE(STDOUT,*) ' '
              loop = loop + 1
            END IF
          END IF

          IF (nukeline) then
            WRITE(SCREEN,*) '--- RESTARTING FIT W/O NUKED LINE ---'
            WRITE(STDOUT,*) '--- RESTARTING FIT W/O NUKED LINE ---'
            WRITE(SCREEN,*) ' -> resetting loop iterations '
            WRITE(STDOUT,*) ' -> resetting loop iterations '
            WRITE(SCREEN,*) ' '
            WRITE(STDOUT,*) ' '
            loop = 1
          END IF

          CALL modelcomm(a,n,0) ! to be input into the LSF
          GOTO 02  ! re minimize the adjusted model

        END IF


C
C     END OF INNER LOOP

c
c     F-TEST THE ADJUSTED MODEL, DOFTEST returns the logical flag
c     REPLACE, if it is high, then we keep the adjusted model, if it is
c     low, we restore the ORIGINAL saved model and move to the next bad
c     line in the priority lies
c

 03    CALL doFtest(n,m,a,var1,nu1,var2,nu2,ftest,fprob,replace)

c     if replacing (keeping the test model), then communicate, compute
c     the errors in the new fitting parameters and then check if the
c     errors worked out.  if they did not then communicate that we are
c     restoring the original model and reset logical REPLACE low for
c     downstream actions

       IF (replace) THEN

        WRITE(SCREEN,*) ' ' 
        WRITE(SCREEN,*) '  ***                                  ***'
        WRITE(SCREEN,*) '  ***                                  ***'
        WRITE(SCREEN,602) ftest,(1.0d0-fprob)
        WRITE(SCREEN,*) '  ***      REJECTING TEST LINE(S)      ***'
        WRITE(SCREEN,*) '  ***        ADOPTING NEW MODEL        ***'
        WRITE(SCREEN,*) '  ***                                  ***'
        WRITE(SCREEN,*) '  ****************************************'
        WRITE(SCREEN,*) '  ***                                  ***'
        WRITE(SCREEN,*) '  ***                                  ***'
        WRITE(SCREEN,*) ' ' 

        WRITE(STDOUT,*) ' ' 
        WRITE(STDOUT,*) '  ***                                  ***'
        WRITE(STDOUT,*) '  ***                                  ***'
        WRITE(STDOUT,602) ftest,(1.0d0-fprob)
        WRITE(STDOUT,*) '  ***      REJECTING TEST LINE(S)      ***'
        WRITE(STDOUT,*) '  ***        ADOPTING NEW MODEL        ***'
        WRITE(STDOUT,*) '  ***                                  ***'
        WRITE(STDOUT,*) '  ****************************************'
        WRITE(STDOUT,*) '  ***                                  ***'
        WRITE(STDOUT,*) '  ***                                  ***'
        WRITE(STDOUT,*) ' ' 

        CALL fitstats(a,n,m,1) 
        CALL modelcomm(a,n,1) ! current result of LSF

        WRITE(SCREEN,*) ' ' 
        WRITE(SCREEN,*) '  ***                                  ***'
        WRITE(SCREEN,*) '  ***   PROCEEDING TO TEST NEW MODEL   ***'
        WRITE(SCREEN,*) '  ***                                  ***'
        WRITE(STDOUT,*) ' ' 
        WRITE(STDOUT,*) '  ***                                  ***'
        WRITE(STDOUT,*) '  ***   PROCEEDING TO TEST NEW MODEL   ***'
        WRITE(STDOUT,*) '  ***                                  ***'


        CALL errors(a,siga,n,m,error)
        IF (error) THEN
         WRITE(SCREEN,*) ' Trouble with the errors in the F-Test model'
         WRITE(SCREEN,*) '     restoring Pre F-Test Model'
         WRITE(STDOUT,*) ' Trouble with the errors in the F-Test model'
         WRITE(STDOUT,*) '     restoring Pre F-Test Model'
         replace = .false.
        END IF
       END IF

c     now that we have checked that the errors worked out, we can go
c     ahead and finalize the replacement, assuming logical REPLACE is
c     still set high; communicate the new fit statistics; update the
c     ZLINE array (for plotting purposes), and plot

c     OK... so we have a new model with reasonable parameters and
c     computed errors.  this new model now needs to be checked if it is
c     to survive as the final model. so back to the beginning of the
c     OUTER LOOP (code line 01) to start the whole game over again; this
c     means that we have to examine each line in the new model to see if
c     it is not a bad line

       IF (replace) THEN
        DO 103 linei=1,lines
         i = idx_az(linei)
         zline(linei) = a(i)
 103    CONTINUE
        IF (plotting) CALL plt
        GOTO 01 ! start all over again with this new model!        
       END IF

C
C     END OF THE OUTER LOOP (which is embedded in DO 999 CONTINUE LOOP)
C

C
C     IF WE ARE HERE THEN WE RETAINED THE CURRENT MODEL; WE BOOK KEEP
C     AND THEN MOVE ON TO THE NEXT BAD LINE
C

c     the tested line was significant or the erros crapped out; put the
c     line back by reloading the A and SIGA vectors, and restoring LINES
c     and the number of coefficients N; replot the original model

      IF (.not.replace) THEN

c     and finally, if we are not replacing, we continue on with the
c     current bad line list until it is exhausted

        lines = lines_save
        n     = n_save
        DO 31 i=1,n
         a(i)    = asave(i,1)
         siga(i) = asave(i,2)
 31     CONTINUE

        DO 35 i=1,lines
         idx_ab(i) = idxbsave(i) 
         idx_az(i) = idxzsave(i) 
         DO 33 k=1,species
          idx_an(k,i)  = idxnsave(k,i) 
          satflag(k,i) = satflagsave(k,i)
 33      CONTINUE
 35     CONTINUE


 602    FORMAT(1x,'  ****************************************',/,
     &         1x,'  ***                                  ***',/,
     &         1x,'  ***         F-TEST RESULTS           ***',/,
     &         1x,'  ***                                  ***',/,
     &         1x,'  ***        F = ',f9.5,'             ***',/,
     &         1x,'  ***       CL = ',f9.5,'             ***',/,
     &         1x,'  ***                                  ***')     

        WRITE(SCREEN,*) '  ***                                  ***'
        WRITE(SCREEN,*) '  ***                                  ***'
        WRITE(SCREEN,602) ftest,(1.0d0-fprob)
        WRITE(SCREEN,*) '  *** RESTORING PRE F-TEST MODEL       ***'
        WRITE(SCREEN,*) '  ***                                  ***'
        IF (j.lt.nunique) then
        WRITE(SCREEN,*) '  *** ADVANCING TO THE NEXT BAD LINE   ***'
        WRITE(SCREEN,*) '  ***                                  ***'
        END IF
        WRITE(SCREEN,*) '  ****************************************'
        WRITE(SCREEN,*) '  ***                                  ***'
        WRITE(SCREEN,*) '  ***                                  ***'
        WRITE(SCREEN,*) ' ' 
        WRITE(STDOUT,*) '  ***                                  ***'
        WRITE(STDOUT,*) '  ***                                  ***'
        WRITE(STDOUT,602) ftest,(1.0d0-fprob)
        WRITE(STDOUT,*) '  *** RESTORING PRE F-TEST MODEL       ***'
        WRITE(STDOUT,*) '  ***                                  ***'
        IF (j.lt.nunique) then
        WRITE(STDOUT,*) '  *** ADVANCING TO THE NEXT BAD LINE   ***'
        WRITE(STDOUT,*) '  ***                                  ***'
        END IF
        WRITE(STDOUT,*) '  ****************************************'
        WRITE(STDOUT,*) '  ***                                  ***'
        WRITE(STDOUT,*) '  ***                                  ***'
        WRITE(STDOUT,*) '  '
        CALL fitstats(a,n,m,1)  
        CALL modelcomm(a,n,1) ! this is an LSF model
        DO 100 linei=1,lines  ! update tick locations for plotting
         i = idx_az(linei)
         zline(linei) = a(i)
 100    CONTINUE
        IF (plotting) call plt

       END IF

c     go to the next bad line in the current model

 999  CONTINUE


c     return

      RETURN

      END

c
c..............................................................................
c

      SUBROUTINE        checka(a,n,nukeline,blwa,nknt)

c
c     this routine will unceremoniously nukes multiple lines from the
c     fit if the column density goes negative (in the flag ion only), or
c     the line goes outside the redshift window of the region.  
c
c     Since some weak transitions can have this happen but the line is
c     critical in the stronger lines it is important to examine only the
c     strongest line, or we will adversely effect the F-testing process
c     be skewing the statistics
c
c     if this routine nukes a line during an F-test run, then if and only
c     if the pre-F-test model is restored will the line be restored.
c     this is why it is important to NOT nuke lines based upon weak
c     transition, because that will artificually elevate the
c     significance of tested lines 
c
c     after nuking a line, the number of lines has decreased so we must
c     adjust by resetting the book keeping
c
C     some old ideas
C     IF (a(ib).lt.0.0d0) a(ib) = abs(a(ib))   ! make T positive
C     Tbmin = spec_mass(speciesi)*1.0d10*(bmin**2) / (2.0d0*kerg)
C     IF (a(ib).le.Tbmin) a(ib) = 3.0d0*Tbmin   ! make T 3*Minimum allowed
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include          'minfit.h'
      include          'const.dek'
      logical           nukeline,ifnan,flagN,flagz
      integer           speciesi,linei,in,ib,n,iz,badline,flgline,iflag
      integer           nknt,blwa(maxlines)
      double precision  a(maxcoeffs),siga(maxcoeffs)
      double precision  fna,Tp,btest
      character*6       mess


c     if we nuke lines then the line numbers change, but we want to keep
c     track of the old line numbers for communicated the correct
c     original line number to the user, the number of nuked lines is
c     counted in NKNT and the array BLWA stores the nuked line numbers;
c     we later unfold the information in routine FIXA using the array
c     ORIGLINE

c     null NKNT and zero the BLWA array

      nknt = 0   
      DO 09 linei=1,maxlines
       blwa(linei) = 0
 09   CONTINUE

c     we only want to check the most significant ion in the user
c     list. if there is no flagion set by the user, we then assume that
c     the first ion in the input list, with spec_idx(ioni)=1, is the
c     most significant ion, in which case we set FLGLINE=1

      IF (flagion(regioni).eq.0) then
         flgline = 1
      ELSE
         flgline = flagion(regioni)
      END IF

c     communicate which transition we are checking

      WRITE(SCREEN,401) spec_name(flgline)
      WRITE(STDOUT,401) spec_name(flgline)

 401  FORMAT(1x,/,1x,'--- EXAMINING ',a7,' PARAMETERS FOR NUKE ... ')


c
c     MAIN PART OF THE SUBROUTINE
c

c     loop over the species and then each line for the flagion species
c     we see this loop as many times as we need to in order to nuke them
c     all

 01   DO 11 speciesi=1,species

        IF (fitindex(speciesi).eq.flgline) then ! flagion/first species only 

        flagN = .false.  ! will communicate bad N when set .true.
        flagz = .false.  ! will communicate bad z when set .true.

        do 15 linei=1,lines
         in = idx_an(speciesi,linei)
         ib = idx_ab(linei)
         iz = idx_az(linei)

c     here is the check for negative column densities in the flagion,
c     nuke only if it is not in a saturated region. -- if the column is
c     negative and in a saturated region we hold off on any action here;
c     the column will be adjusted in FIXA

         IF ((a(in).le.0.0d0).AND.(.not.satflag(speciesi,linei))) THEN
           nukeline = .true.
           badline  = linei
           flagN    = .true.
           GOTO 21 ! jump out of loops to book keep and start over
         END IF

c     here is the check for bad redshifts; if the LSF has sent a line
c     outside the redshift range of the current working redshift region,
c     then we nuke the line

         IF ((a(iz).lt.zlim(regioni,1)).OR.
     &       (a(iz).gt.zlim(regioni,2))) THEN
           nukeline = .true.
           badline  = linei
           flagz    = .true.
           GOTO 21 ! jump out of loops to book keep and start over
         END IF

 15     CONTINUE  ! next linei

       END IF

 11   CONTINUE  ! next speciesi

c     if we get here we are done nuking, return; if we did not nuke any
c     lines then communicate

      IF (.not.nukeline) WRITE(SCREEN,*) ' -> all lines significant'
      IF (.not.nukeline) WRITE(STDOUT,*) ' -> all lines significant'

      RETURN

C
C     JUMP POINT
C

c     yeah, if we found that a single line model went bad there is
c     nothing much we can do here in MINFIT, so if this is the only
c     line, then the user is SOL and we STOP

 21   IF (lines.eq.1) THEN
        WRITE(SCREEN,*) ' DEATH(checka): single line model went bad'
        WRITE(STDOUT,*) ' DEATH(checka): single line model went bad'
        STOP 
      END IF

c     communicate line nuke

      fna = sign(1.0d0,a(in))

      IF (flagN) mess = 'bad N'
      IF (flagz) mess = 'bad z'

      WRITE(SCREEN,602) mess,(linei+nknt),a(iz),
     @                  fna*(13.0+log10(abs(a(in)))),a(ib)
      WRITE(STDOUT,602) mess,(linei+nknt),a(iz),
     @                  fna*(13.0+log10(abs(a(in)))),a(ib)

 602  FORMAT(1x,' -> Nuking Line: (',a6,') line=',i3,' z =',f9.6,
     &          '  logN =',f7.2,'  T =',1pe12.5)


c     nuke and book keep the line nuke; update the A array

c     first increment the nuke counter NKNT and store the line number so
c     that we communicate the correct line number from the original
c     model (at the time of the call to this routine)

      nknt = nknt + 1
      blwa(nknt) = linei + nknt - 1

      CALL jukea(a,n,badline)

c
c     REPEAT WITH NEW MODEL - get all the bad lines
c     

C      RETURN
      GOTO 01

c     we do not need a RETURN here because it is the trap after the
c     SPECIESI and LINEI loops above when no bad lines are found

      END

c
c..............................................................................
c

      SUBROUTINE        fixa(a,n,fixline,blwa,nknt)

c
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include          'minfit.h'
      include          'const.dek'
      logical           fixline,ifnan,limflag
      integer           i,j,k,speciesi,fidx,linei,in,ib,n,iz,badline
      integer           blwa(maxlines),nknt,hit,origline(maxlines)
      double precision  a(maxcoeffs),Tbmin,Tbmax


c     if FIXLINE is set high, the we will be iterating the fit upon
c     return to routine FITREGION or CHECKLINES


c     the paramter nknt tracks how many lines were nuked in the
c     proceeding call to CHECKA, we use it for communication only

      WRITE(SCREEN,401) 
      WRITE(STDOUT,401) 

 401  FORMAT(1x,/,1x,'--- EXAMINING ALL PARAMETERS FOR ADJUSTMENT ... ')

      limflag = .false.
      fixline = .false.

c     first, for communication purposes of the original line numbers,
c     store the original line number prior to line nuking from CHECKA
c     NKNT contains the number of lines nuked in CHECKA and BLWA
c     contains the original line numbers that were nuked.  this is a bit
c     convoluted logic, but it works

       linei = 0
       DO 71 j=1,lines+nknt   ! over the original number of lines
        hit = 0
        DO 72 k=1,j
         DO 73 i=1,nknt
          IF (blwa(i).eq.j) then
           hit = hit + 1
           IF (k.eq.j) GOTO 71 ! if the current line then do not increment
          END IF
 73      CONTINUE
 72     CONTINUE
        linei           = linei + 1
        origline(linei) = j
C        WRITE(6,*) origline(linei),linei
 71    CONTINUE

c     get to work, loop over the species and lines

 01   CONTINUE  ! start... if we nuke a line start over

      DO 11 speciesi=1,species

        fidx = fitindex(speciesi)

        DO 15 linei=1,lines

C         write(6,*) ' DEBUG: checking species ',spec_name(fidx)
C         write(6,*) ' DEBUG: line             ',linei

         in = idx_an(speciesi,linei)
         ib = idx_ab(linei)
         iz = idx_az(linei)

C         write(6,*) 'IN,IB,IZ = ',in,ib,iz
C         write(6,*) 'N,b,z    = ',sign(1.0d0,a(in))*
C     &                            (13.+log10(abs(a(in)))),
C     &                            a(ib),a(iz)

c     check that the column density is not an IEEE floating point
c     execption; if so, remove it from the fit in routine JUKEA

         IF (ifnan(a(in))) THEN
           WRITE(SCREEN,603) spec_name(fidx),linei,origline(linei)
           WRITE(STDOUT,603) spec_name(fidx),linei,origline(linei)
           badline = linei
           CALL jukea(a,n,badline)
           fixline    = .true.
           nknt       = nknt + 1
           blwa(nknt) = badline
           GOTO 01  ! lines will have been changed and the linei loop 
         END IF     ! will be corrupt - start over

c     check that the b parameter is not an IEEE floating point
c     execption; if so, remove it from the fit in routine JUKEA

         IF (ifnan(a(ib))) THEN
           WRITE(SCREEN,604) spec_name(fidx),linei,origline(linei)
           WRITE(STDOUT,604) spec_name(fidx),linei,origline(linei)
           badline = linei
           CALL jukea(a,n,badline)
           fixline    = .true.
           nknt       = nknt + 1
           blwa(nknt) = badline
           GOTO 01  ! lines will have been changed and the linei loop 
         END IF     ! will be corrupt - start over

c     check that the b parameter is not negative; if so, take its
c     absolute value and iterate the fit

         IF (a(ib).le.0.0d0) THEN
           a(ib)   = abs(a(ib))
           fixline = .true.
           WRITE(SCREEN,602) spec_name(fidx),linei,origline(linei)
           WRITE(STDOUT,602) spec_name(fidx),linei,origline(linei)
         END IF

c     check that the b parameter is not less than bmin; if so, increase
c     to bmin

         Tbmin = spec_mass(fidx)*1.0d10*(bmin**2) / (2.0d0*kerg)
         IF (a(ib).lt.Tbmin) THEN
           a(ib)   = Tbmin   
           fixline = .true.
           WRITE(SCREEN,605) spec_name(fidx),linei,origline(linei)
           WRITE(STDOUT,605) spec_name(fidx),linei,origline(linei)
         END IF

c     check that the b parameter is not less than bmin; if so, increase
c     to bmin

         Tbmax = spec_mass(fidx)*1.0d10*(bmax**2) / (2.0d0*kerg)
         IF (a(ib).gt.Tbmax) THEN
           a(ib)   = Tbmax   
           fixline = .true.
           WRITE(SCREEN,606) spec_name(fidx),linei,origline(linei)
           WRITE(STDOUT,606) spec_name(fidx),linei,origline(linei)
         END IF

c     we call FIXN last to make sure we have a fixed bpar
c
c     check that columnn density is not negligible; if so then call
c     routine FIXN, which returns an estimate of the column density
c     based upon the equivalent width limit in the spectrum or the AOD
c     based upon the flux values in the spectrum

         IF (a(in).le.1.0d-5) THEN  !(A<1.e-5 is logN<8)
           WRITE(SCREEN,601) spec_name(fidx),linei,origline(linei)
           WRITE(STDOUT,601) spec_name(fidx),linei,origline(linei)
           CALL fixN(speciesi,fidx,linei,a(in),a(ib),a(iz),limflag)
           fixline = .true.
         END IF

 15     CONTINUE
 11   CONTINUE

c     return

      IF (.not.fixline) WRITE(SCREEN,*) ' -> all parameters OK ' 
      IF (.not.fixline) WRITE(STDOUT,*) ' -> all parameters OK ' 

      RETURN


 601  FORMAT(1x,' -> Fix Line: ',a6,' new line=',i3,' old line=',i3,
     &        ' N insignificant')
 602  FORMAT(1x,' -> Fix Line: ',a6,' new line=',i3,' old line=',i3,
     &          ' b negative - sign changed')
 603  FORMAT(1x,' -> Fix Line: ',a6,' new line=',i3,' old line=',i3,
     &          ' N nan (line removed)')
 604  FORMAT(1x,' -> Fix Line: ',a6,' new line=',i3,' old line=',i3,
     &          ' b nan (line removed)')
 605  FORMAT(1x,' -> Fix Line: ',a6,' new line=',i3,' old line=',i3,
     &          ' b < bmin - reset to bmin')
 606  FORMAT(1x,' -> Fix Line: ',a6,' new line=',i3,' old line=',i3,
     &          ' b > bmax - reset to bmax')

      END
c
c..............................................................................
c

      SUBROUTINE        jukea(a,n,badline)

c
c     adjust the coefficient array; this is my masterpiece!
c
c     this is a very carefully crafted routine that is called under two
c     conditions
c
c     condition 1: badline=0; this call is made from routine defcoeffs,
c     which sets up the initial a array and index arrays prior to
c     consideration of tied saturated lines; under this condition,
c     badline is ignored and only the conditionals for tieing saturated
c     lines are examined; when a saturated tie is to occur, the a array
c     and the index array are dynamically adjusted
c
c     condition 2: badline is non-zero; if routine checka flags a
c     badline then this routine is used during the iterative F-test to
c     check the significance of the line to the LSF. 
c
c     on output the original a and index arrays are destroyed (but, they
c     are saved in the calling routine in the case of the F-test in case
c     we restore the pre-F-test model).
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include          'minfit.h'
      include          'const.dek'

      logical           satzone
      integer           i,j,k,k2,n,badline,badlinesave
      integer           idxnsave(maxspecies,maxlines)
      double precision  a(maxcoeffs),zcl
      double precision  nsave(maxspecies,maxlines),
     &                  bsave(maxlines),zsave(maxlines)


      IF (dodebug) then 
c     PRINT INPUT ARRAY FOR TESTING PURPOSES ONLY
      write(6,*) ' ------ '
      write(6,*) ' BADLINE = ',badline
      do 91 k=1,species
       write(6,600) (idx_an(k,i),sign(1.0d0,a(idx_an(k,i)))*
     &               (13.+log10(abs(a(idx_an(k,i))))),
     &               i=1,lines)
 91   continue
      write(6,600) (idx_ab(i),a(idx_ab(i)),i=1,lines)
      write(6,601) (idx_az(i),a(idx_az(i)),i=1,lines)
      write(6,*) ' ------ '
      END IF
C     END PRINT TEST


c     save the N indices and N,b,z groups; logical satflag is used in
c     routine model to employ the scaling between tied columns and in
c     routine unpacka to also employ the scaling

      do 11 i=1,lines
       bsave(i)    = a(idx_ab(i))
       zsave(i)    = a(idx_az(i))
       do 13 k=1,species
        idxnsave(k,i) = idx_an(k,i)
        nsave(k,i)    = a(idx_an(k,i))
        satflag(k,i)  = .false.
 13    continue
 11   continue

c     decrement the number of lines if one of the lines is being tossed;
c     if not set badline to larger than the number of lines so that we
c     avoid the line tossing logic; save badline so that it is not
c     changed in the calling routine on output

       badlinesave = badline
       If (badline.ne.0) then
         lines = lines - 1
       else
         badline = lines + 1
       end if

c     set the dynamic a index counter (J) to zero

      j      = 0 

c     this is where the killer logic is applied; we telescope the column
c     densities part of the a and index array accounting for saturated
c     ties and we also keep track of badline removal; K is the species
c     being examined, K2 is the species that it is scaled to if there is
c     saturated

      do 21 k=1,species
       do 23 i=1,lines
        if (i.lt.badline) then ! this is pre badline 
         zcl = zsave(i)
         if (satzone(k,k2,zcl)) then ! sat zone
          write(SCREEN,700) i,spec_name(k),spec_name(k2)
          write(STDOUT,700) i,spec_name(k),spec_name(k2)
          idx_an(k,i)   = idxnsave(k2,i)
          satflag(k,i)  = .true.
          satflag(k2,i) = .true.
         else                                                  ! NOT sat zone
          j = j + 1
          a(j) = nsave(k,i)
          idx_an(k,i) = j 
         end if
        else                   ! this is badline or post badline
         zcl = zsave(i+1)
         if (satzone(k,k2,zcl)) then ! sat zone
          write(SCREEN,700) i,spec_name(k),spec_name(k2)
          write(STDOUT,700) i,spec_name(k),spec_name(k2)
          idx_an(k,i)  = idxnsave(k2,i)
          satflag(k,i)  = .true.
          satflag(k2,i) = .true.
         else                                                  ! NOT sat zone
          j = j + 1
          a(j) = nsave(k,i+1)
          idx_an(k,i) = j 
         end if
        end if
 23    continue
 21   continue

c     now telescope the b parameters; remove badline

      do 31 i=1,lines
       j = j + 1
       if (i.lt.badline) then ! this is pre badline 
        a(j) = bsave(i)
       else                   ! this is badline or post badline
        a(j) = bsave(i+1)
       end if
       idx_ab(i) = j 
 31   continue

c     now telescope the redshifts; remove badline

      do 41 i=1,lines
       j = j + 1
       if (i.lt.badline) then  ! this is pre badline 
        a(j) = zsave(i)
       else                    ! this is badline or post badline
        a(j) = zsave(i+1)
       end if
       idx_az(i) = j 
 41   continue


c     reset the number of fitted coefficients (the length of a)

      n = j

c     we are done; restore badline
      
      badline = badlinesave

      IF (dodebug) then
C     PRINT ADJUSTED ARRAY FOR TESTING PURPOSES ONLY
       write(6,*) ' N = ',n
       do 99 k=1,species
        write(6,600) (idx_an(k,i),sign(1.0d0,a(idx_an(k,i)))*
     &                (13.+log10(abs(a(idx_an(k,i))))),
     &                i=1,lines)
 99    continue
       write(6,600) (idx_ab(i),a(idx_ab(i)),i=1,lines)
       write(6,601) (idx_az(i),a(idx_az(i)),i=1,lines)
       write(6,*) ' ------ '
      END IF
C     END TEST PRINT

c     return

      return

c     formats

 600  FORMAT(1x,40(2x,i3,':',1pe10.2))
 601  FORMAT(1x,40(2x,i3,':',2x,f8.6))
 700  FORMAT(2x,'saturated line=',i2,': ',a6,
     @          ' column density scaled to ',a6)

      end

c
c
c..............................................................................
c

      LOGICAL FUNCTION satzone(speciesi,j2,zcl)

c
c     logical function to determine if cloud redshift is in a saturation
c     zone defined in the sat_regions.dat file; the input data for
c     saturation zones is put in the file sat_regions.dat; if the file
c     DNE then nsatreg=0; we determine if this is a species and cloud
c     redshift for which saturation applies, if so we set logical
c     function satzone high and return the anchor transition for the
c     scaled column density (j2)
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      include           'const.dek'

      integer           i,speciesi,j1,j2
      double precision  zcl,z1,z2


c     initialize to false

      satzone = .false.

      do 09 i=1,nsatreg
        
       j1 = sat_spec(i,1) 
       j2 = sat_spec(i,2)
       z1 = zsat(i,1)
       z2 = zsat(i,2)

       if ((speciesi.eq.j1).AND.
     &     (zcl.ge.z1).AND.(zcl.le.z2)) then
        satzone = .true.
        return
       end if

 09   continue


c     return

      return

      end
c
c..............................................................................
c

      DOUBLE PRECISION FUNCTION satNscale(Nvp,a,b,Naction)
c  
c     the saturation scaling is log10(x2) = A*log10(x1) + B
c
c     identically, we have x2 = x1^a * 10^b
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
C
      implicit none 
      integer           Naction
      double precision  Nvp,a,b

c     Naction is defined in the calling routine

c     if Naction=0, is called from routines MODEL and MODELCHUNK, and
c     indicates that we are working in code units for the column density
c     (1.e-13 atoms cm^-2)

      If (Naction.eq.0) satNscale = (abs(Nvp)**a)*(10.0d0**b)

c     if Naction=1, is called from routines UNPACKA and INFORM and
c     indicates we are working in log10 units for the column density, so
c     we add 13.0 to convert from code units and scale using the linear
c     relation of logarithms

      If (Naction.eq.1) satNscale =  b + a*(13.0d0+log10(abs(Nvp)))

      RETURN
      END 
c eof




