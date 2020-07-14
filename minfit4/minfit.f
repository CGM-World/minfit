c..............................................................................
c
      PROGRAM           minfit

c     author: Chris Churchill (cwc@nmsu.edu)
c     started at Penn State, updated and expanded at New Mexico State
c     - if you use this code, or any of its support codes, please reference
c       "Churchill, C. W. 1997, Ph.D. Thesis, University of California, 
c        Santa Cruz"
c
c     This program performs least square fitting and chi2 minimization
c     of Voigt profiles to absorption lines in spectra
c
c
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     "to the kitchen...."
c

      include           'minfit.h'
      external          fcn
      integer           j,jj,iopt,io_iter,imongo
      integer           m,n,info,iw(maxcoeffs),ilwa
      real*4            xpos,ypos
      double precision  a(maxcoeffs),siga(maxcoeffs),
     @                  fvec(maxvec),wa(lwa),toler
      character*80      ionlist
      character*12      adate,atime


c     COMMUNICATION

c     communication to the user: we write all communication to the
c     screen (unit=SCREEN) and to the "minfit.runlog" file
c     (unit=STDOUT), the unit values SCREEN and STDOUT are defined 
c     in the "minfit.h" module

c     grab the machine time and date      

      CALL today(adate,atime)

c     open the "minfit.runlog" file

      OPEN(unit=STDOUT,file='minfit.runlog',status='unknown')

c     communicate header

      DO 901 j=1,2
       IF (j.eq.1) jj = SCREEN
       IF (j.eq.2) jj = STDOUT
       WRITE(jj,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
       WRITE(jj,*) '+                                              +'
       WRITE(jj,*) '+                 MINFIT V4.3                  +'
       WRITE(jj,*) '+                 March 2011                   +'
       WRITE(jj,*) '+                                              +'
       WRITE(jj,*) '+         VOIGT PROFILE MINIMIZATION           +'
       WRITE(jj,*) '+                                              +'
       WRITE(jj,*) '+                  RUN STAMP                   +'
       WRITE(jj,900) adate,atime
       WRITE(jj,*) '+                                              +'
       WRITE(jj,*) '++++++++++++++++++++++++++++++++++++++++++++++++'
 901  CONTINUE

 900  FORMAT(1x,'+',17x,a12,17x,'+',/,1x,'+',18x,a12,16x,'+')


c     initialize the array sizes

      total    = 0          ! the total number of species read in
      species  = 0          ! the number of species in the LSF in REGIONI
      ions     = 0          ! the number of transitions read in
      lines    = 0          ! the number of VP components in REGIONI
      regions  = 0          ! the number of fitting regions
      vpcomps  = 0          ! the total number of VP components all regions
      dodebug  = .false.    ! logical flag for extra print statements

c     COMMAND LINE INPUTS (see comments in routine GETCOMMLINE)

      CALL getcommline(ionlist)

c     set up the LSF

      CALL getlist(ionlist) ! loads the ion/transition list
      CALL getparams        ! loads the run parameters 'minfit.par'
      CALL getatomic        ! loads the atomic data
      CALL instrument(1)    ! sets up the instrumental spread function
      CALL getregions       ! loads the wavelength/redshift regions
      CALL getsatreg        ! lads info on saturation (in velocity space)

c     are we plotting?  if so, initialize and plot input model

      if (plotting) then
       CALL intro(11)       ! open the plotting window and communicate
       CALL sleep(3)        ! pause so communication can be read
       CALL pltinput        ! plot the input VP model (loaded within)
      end if

c     intialize the SLATEC parameters (must follow routine GETPARAMS)

      iopt     = 1          ! DNL1SE computes the Jacobian
      toler    = tol        ! tolerance for convergence
      ilwa     = lwa        ! maximum size of working arrays
      io_iter  = nprint     ! frequency at which updates are communicate



c
c     MAIN ROUTINE: THE MINIFIT DRIVER LOOP
c

c     now we begin the fitting loop- region by region; first, open the
c     output files; then loop over regions, calling routine fitregions,
c     which calls the slatec LSF.  If we return unhappily (info=0), then
c     kick out, and do not write to files; if (info.neq.0), [it can have
c     several values], then write the results for the region.  when
c     done, close the files.

c     open the data files that we will write to region by region

      CALL openfiles

c     loop over the wavelength/redshift regions

      do 21 regioni=1,regions

c     a little communication nicety

      write(SCREEN,*) ' '
      write(SCREEN,*) '************************************************'
      write(STDOUT,*) ' '
      write(STDOUT,*) '************************************************'

c     if plotting, then we query the user to click on the window to
c     start each new region

      if (plotting) then
       write(SCREEN,600) regioni         ! query to advance to the next region
       CALL mongohairs(imongo,xpos,ypos) ! wait for the click
      end if

 600  FORMAT(1x,' Click on plotting window ',/,
     &       1x,' to proceed with REGION :',i3,/)

c     do the LSF for the region; if an error close up shop and exit

      CALL fitregion(fcn,iopt,m,n,a,siga,fvec,tol,
     @               io_iter,info,iw,wa,ilwa)
      if (info.eq.0) goto 1000

c     write data to the files and store the LSF model stats for the
c     current region

      CALL writefiles
      CALL fitstats(a,n,m,2) 

 21   continue

c     close the files

 1000 CALL closefiles

c     ERROR TRAP: successful or unsuccessful? do we have a solution, if
c     not abort (bad solution communicated already).

      if (info.eq.0) STOP ' aborted MINFIT'


c
c     END MINIFIT DRIVER LOOP
c

c     CLOSE UP SHOP

c     a little communication nicety

      write(SCREEN,*) ' '
      write(SCREEN,*) '************************************************'
      write(STDOUT,*) ' '
      write(STDOUT,*) '************************************************'

c     if we are plotting, then plot the final solution! we query the use
c     to click on the window to finally present the overall results; we
c     we are not plotting we still communicate the final numbers to the
c     screen and the "minfit.runlog" file

      if (plotting) then
       write(SCREEN,601)                  ! query to view overall model
       CALL mongohairs(imongo,xpos,ypos)  ! wait for the click
       CALL pltoutput                     ! plt the overall final model
       CALL finalcomm                     ! communicate final model params
       write(SCREEN,602)                  ! query for minfit quit
       CALL mongohairs(imongo,xpos,ypos)  ! wait for the click
      else
       CALL finalcomm                     ! communicate final model params
      end if

 601  FORMAT(1x,' Click on plotting window ',/,
     &       1x,' to view final LSF model',/)
 602  FORMAT(1x,' Click on plotting window ',/,
     &       1x,' to QUIT MINFIT',/)

c     since we have all the VP components stored, we write the model
c     spectra to individual files for each transition; inlcudes the
c     residuals as well

      CALL vpmods

c     stamp the end time 

      CALL today(adate,atime)

      write(SCREEN,*) ' ' 
      write(SCREEN,*) ' '
      write(SCREEN,*) ' MINFIT end ',adate,' @ ',atime
      write(STDOUT,*) ' ' 
      write(STDOUT,*) ' '
      write(STDOUT,*) ' MINFIT end ',adate,' @ ',atime

c     close the "minfit.runlog" file

      CLOSE(unit=STDOUT)

c     successful terminate

      STOP ' normal termination (minfit)'

      END

c
c..............................................................................
c

      SUBROUTINE          today(adat,atim)   

c
c     this routine gets the date and time out of a machine.  the output
c     format is
c
c     adat dd mmmyy
c     atim hh:mm:ss
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      character*3        amon(12) 
      character*12       atim,adat
      character*10       adigit   
      parameter          (adigit='0123456789')
      integer*4          idat(3),itim(3)  

      data amon/ 'Jan' , 'Feb' , 'Mar' , 'Apr' , 'May' , 'Jun' ,    
     +           'Jul' , 'Aug' , 'Sep' , 'Oct' , 'Nov' , 'Dec' /

      adat=' '  
      atim=' '  

c     grab the time and "write" into variable atim

      CALL itime(itim)  
      write (atim,113) itim  

c     grab the date and write into variable adat

      CALL idate(idat)  
      write (adat,114) idat(1),amon(idat(2)),idat(3)

 113  format(i2,':',i2,':',i2)  
 114  format(i2,'-',a3,i4)  

      return

      end   

c eof
