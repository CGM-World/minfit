c------------------------------------------------------------------------------
c

        SUBROUTINE        getparams

c
c       Read in the parameter file, which is named "minfit4.par".  The
c       comments below give the order in which the parameters must
c       appear in the file.
c
c
C       input params for minfit - see "helpfile.f" for explanation 
C***************** input parameters for minfit ********************************
c
C "/home/matrix/cwc/Data/Atomic/atoms.dat"        Atomic data file
C 5       nprint - frequency at which iterations are communicated
C 1.e-7   tol - the tolerance of the fit
C 1.0     fudge - scale sigma vector by sqrt of fudge
C 4.25    N_sigma - feature detection threshold using LTW87
C 1.0     bmin - minimum b par allowed for adjustment testing
C 20.0    bmax - maximum b par allowed for adjustment testing
C 16.0    Nmax - maximum log10 column density for adjustment testing
C 10.     fberr - fractional error in b par for N limit computations
C 1       iturb - integer flag: 1=pure turbulance, 0=pure thermal
C 1       pltflag - integer flag: any=plotting OFF  1=plotting ON
C 2       adjflag - integer flag: adjust/toss -> 0=OFF/OFF  1=ON/OFF 2=ON/ON
C 11.5    badness - minimum fractional error in "suspect" line
C 0.97    confidence - confidence level at which to reject a "suspect" line
C 1       errflag - integer flag: any=compute errors OFF  1=compute errors ON
C 1       conflag - integer flag: any=convolution OFF  1=convolution ON
C 45000.  R_fac - the spectrograph resolution R=lambda/(Delta lambda)
C 0.861   slit - the spectrograph slit width in arcseconds
C 3.      conwindo - # of inst. profile sigma to include in convolution 
C 3.      resfac - (odd integer) sample rate factor for convolution arrays 
C******************************************************************************
c

        include           'minfit.h'
        integer           k,i,pltflag,conflag,errflag
        integer           access,istat
        character*80      mf_path,def_file,par_file

c       set some flags low

        plotting   = .false.
        convolving = .false.
        adjusting  = .false.
        doerrors   = .false.

c       check current directory for parameter file, if DNE goto error
c       trap, which grabs the default file which lives where the user
c       has set their MFDIR environment variable

c       stat the default files; istat=0 means it is good to GO

        def_file = 'minfit4.par'
	istat    = access(def_file,'r') 

	IF (istat.ne.0) then
         def_file = '/home/matrix/cwc/Programs/Defaults/minfit4.par'
         istat    = access(def_file,'r') 
         IF (istat.ne.0) GOTO 999
        END IF

	OPEN(unit=3,file=def_file,err=999,status='old')

	WRITE(SCREEN,600) def_file(1:50)
	WRITE(STDOUT,600) def_file(1:50)

 01	READ(3,*)                 ! by pass the header
	READ(3,*) ionfile
        READ(3,*) nprint
        READ(3,*) tol
        READ(3,*) fudge
        READ(3,*) N_sigma
        READ(3,*) bmin
        READ(3,*) bmax
	READ(3,*) Nmax
        READ(3,*) fberr
	READ(3,*) iturb

        READ(3,*) pltflag
        if (pltflag.eq.1) plotting = .true.

c       set if adjusting, then doerrors must be high otherwise, we can
c       compute errors but not adjust the fit

        READ(3,*) adjflag
        READ(3,*) badness
        READ(3,*) confidence
        READ(3,*) errflag
        if (adjflag.eq.2) then
         adjusting = .true.
         doerrors  = .true.
        end if
        if (errflag.eq.1) doerrors = .true.

        READ(3,*) conflag
        if (conflag.eq.1) convolving = .true.

        READ(3,*) R_fac
        READ(3,*) slit
        READ(3,*) conwindo
        READ(3,*) resfac        ! must be odd integer

	CLOSE(3)

c       get the redshift

c       READ in the systemtic redshift in file "zabs.dat" and
c       communicate the file name

        OPEN(unit=1,file='zabs.dat',err=998,status='old')
        READ(1,*) zabs
        CLOSE(unit=1)

        WRITE(SCREEN,601) 
        WRITE(STDOUT,601) 

c       communicate the run parameters

      DO 11i=1,2
      IF (i.eq.1) k = SCREEN
      IF (i.eq.2) k = STDOUT

      WRITE(k,*) '************************************************'
      WRITE(k,*) ' '
      WRITE(k,'(a,f10.7)') ' SYSTEMIC REDSHIFT          = ',zabs

      WRITE(k,*) ' '            
      WRITE(k,*) 'LSF CONSTRAINT PARAMETERS'
      WRITE(k,'(a,1pe10.2)') ' -Fitting tolerance         =',tol
      WRITE(k,'(a,f5.2)') ' -Noise scaling             = ',fudge
      WRITE(k,'(a,f5.2)') ' -Significance level        = ',N_sigma
      WRITE(k,'(a,f5.2)') ' -Maximum N [log10]         = ',Nmax
      WRITE(k,'(a,f5.2)') ' -Minimum b parameter       = ',bmin
      WRITE(k,'(a,f5.2)') ' -Maximum b parameter       = ',bmax
      WRITE(k,'(a,f5.2)') ' -Max frac error in b       = ',fberr
      IF (iturb.eq.1) then
      WRITE(k,*) '-Turb/Therm criterion      = 100% TURBULENCE'
      ELSE 
      WRITE(k,*) '-Turb/Therm criterion      = 100% THERMAL'
      END IF

      WRITE(k,*) ' '
      WRITE(k,*) 'LSF ITERATION PARAMETERS'
      WRITE(k,*) '-Parameter errors (T/F)    = ',doerrors
      WRITE(k,*) '-Examine poor components   = ',adjusting
      WRITE(k,'(a,1pe10.2)') ' -Max frac error (badness)  =',badness
      WRITE(k,'(a,1pe10.2)') ' -Confidence Level (F-Test) =',confidence

      WRITE(k,*) ' '
      WRITE(k,*) 'INSTRUMENT/CONVOLUTION PARAMETERS'
      WRITE(k,'(a,1pe10.2)') ' -Resolution                =',R_fac
      WRITE(k,'(a,1pe10.2)') ' -Slit Width [arcsec]       =',slit
      WRITE(k,*) '-Convolving ISF (T/F)      = ',convolving
      IF (convolving) then
      WRITE(k,'(a,f5.2)') ' -Convolution Window        = ',conwindo
      WRITE(k,'(a,f5.2)') ' -Interoplation rate        = ',resfac
      END IF

      WRITE(k,*) ' '
      WRITE(k,*) 'PLOTTING (T/F)             = ',plotting
      WRITE(k,*) ' '

 11   CONTINUE


c       normal return

      RETURN


c       error traps

c       error opening local parameter file; check for the default
c       settings get the Unix environment variable MFDIR, which contains
c       the path to where the file lives; try to open it; on fail, abort

 998    write(SCREEN,*) ' NOT FOUND- file: zabs.dat'
	write(STDOUT,*) ' NOT FOUND- file: zabs.dat'
        STOP ' *** MINFIT(getparms): terminated ***'

 999    write(SCREEN,*) ' NOT FOUND- file: minfit.par'
	write(STDOUT,*) ' NOT FOUND- file: minfit.par'
        STOP ' *** MINFIT(getparms): terminated ***'

 1999	write(SCREEN,*) ' NOT FOUND- minfit.par: parameter listing'
	write(STDOUT,*) ' NOT FOUND- minfit.par: parameter listing'
        STOP ' *** MINFIT(getparams): terminated ***'

c       formats

 600	FORMAT(1x,'INPUT PARAMETER FILE  : ',a50)
 601	FORMAT(1x,'SYSTEM REDSHIFT FILE : zabs.dat')

	end


c------------------------------------------------------------------------------
c
	SUBROUTINE getcommline(ionlist)

c
c     MINFIT requires two input lists: (1) the input model of the VP
c     components, and (2) the list of ions (actually transition) which
c     is used to read in the spectra and book keep ties between species;
c     since one may wish to experiment with the information in these
c     files, non-default names can be placed on the command line; if no
c     comman line arguments are found, then we use the default file
c     names (1) "ivpfit.ivpars' and (2) 'ions.table'; we allow a third
c     option for debuggin purposes (see below)

c     grab the ionlist file, load the ion_name array, load atomic
c     constants, load feature regions, set input parameters
c
c..............................................................................

      include           'minfit.h'
      integer           istat,access
      character*80      ionlist,input1,input2
      

      CALL getarg(1,input1)
      CALL getarg(2,input2)

c     stat the default files; istat=0 means it is good to GO

c     for debuging we wish to print out a lot of intermediate results
c     for examination, if DEBUG is the first argument on the command
c     line, then set it high (DODEBUG is global); we must allow that the
c     remaining command line arguments may contain useful info

      IF (input1.eq.'DEBUG') then
        dodebug = .true.
        CALL getarg(2,input1)  ! check if there is a 2nd argument
        CALL getarg(3,input2)  ! check if there is a 3rd argument
      end if

c     if the command line is blank, then use the defaults

      IF (input1.eq.'') then
       profitfile = 'ivpfit.ivppars'
       istat = access(profitfile,'r')
       IF (istat.ne.0) profitfile = 'ivpfit.last'
      ELSE
       profitfile = input1
      END IF

      write(SCREEN,*) '************************************************'
      write(SCREEN,*) 'INITIAL VP MODEL FILE : ',profitfile(1:30)
      write(STDOUT,*) '************************************************'
      write(STDOUT,*) 'INITIAL VP MODEL FILE : ',profitfile(1:30)

      IF (input2.eq.'') then
       ionlist = 'ions.table'
      ELSE
       ionlist = input2
      END IF

      write(SCREEN,*) 'LIST OF TRANSITIONS FILE : ',ionlist(1:30)
      write(STDOUT,*) 'LIST OF TRANSITIONS FILE : ',ionlist(1:30)

      RETURN
      END


c       eof

