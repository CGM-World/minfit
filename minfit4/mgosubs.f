c     routines for profile.f
c.......................................................................
c                               mgosubs.f
c.......................................................................
c
c  autovp fitting and general Mongo routines
c  routine pltgraf    initializes the box and labels around the plot
c  routine pltdata    p lots the data in histogram format thin solid line
c  routine pltreject  plots the 'x' on masked data
c  routine pltfitt    plots the overall model with thick solid line
c  routine pltoutput  plots the final model
c  routine pltinput   plots the input model
c
c  code: Christopher W. Churchill  April 1994 - July 2010
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c

      SUBROUTINE         intro(iterm)
c
c
c     put up the intro message, intialize Mongo, and intialize model
c     workspace to unity
c
c.......................................................................

      implicit none
      integer            iterm
      real*4             pxmin,pymin,pxmax,pymax


      pxmin = 0.0
      pymin = 0.0
      pxmax = 10000.
      pymax = 10000.

c     starting Mongo... initialize the plotting device wait a little for
c     the window to be mapped

      CALL device(iterm)
      CALL tsetup
      CALL setphysical(pxmin,pymin,pxmax,pymax)
      CALL sleep(2)
      CALL erase
      CALL getphysical(pxmin,pymin,pxmax,pymax)
      pxmax = pxmax * .45
      pymax = pymax * .95
      CALL setphysical(pxmin,pymin,pxmax,pymax)
      CALL erase
      CALL window(1,1,1)

c     define some colors, the first parameter in the list is the color
c     index in the color table: current definition: 0=background=black,
c     1=foreground=white, 2=red, 3=green 4=blue 5=yellow

      CALL makecolor(0,0.,0.,0.)
      CALL makecolor(1,1.,1.,1.)
      CALL makecolor(2,1.,0.,0.)
      CALL makecolor(3,0.,1.,0.)
      CALL makecolor(4,0.,0.,1.)
      CALL makecolor(5,1.,1.,0.)

c     make a unit sized invisible box 
      CALL setlvis(0)
      CALL setlim(0.,0.,1.,1.)
      CALL rect(-1,-1,-1,-1)

c     make a large announcement
      CALL setlweight(2.0)
      CALL setexpand(5.0)
      CALL setcolor(5)
      CALL relocate(0.5,0.8)
      CALL putlabel(6,'MINFIT',5) 
      CALL setexpand(3.0)
      CALL relocate(0.5,0.65)
      CALL putlabel(23,'Voigt Profile Minimizer',5)
      CALL setexpand(2.0)
      CALL setlweight(1.0)
      CALL relocate(0.5,0.35) 
      CALL putlabel(24,'Christopher W. Churchill',5)
      CALL relocate(0.5,0.25) 
      CALL putlabel(16,'V4.3 : July 2010',5)
      CALL setexpand(1.0)
      CALL setlweight(1.0)
      CALL tidle

      CALL setexpand(1.0)
      CALL setlweight(0.5)
      CALL setcolor(1)
      CALL tidle

      return
      end

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c

      SUBROUTINE         plt

c
c     quick and dirty subroutine during the fitting to plot the
c     intermediate results
c
c.......................................................................

       CALL erase
       CALL pltgraf
       CALL pltdata
       CALL pltfitt
       CALL pltreject

       CALL tidle

       return
       end 

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c

      SUBROUTINE         pltgraf

c
c     plot the box, limits, and labels
c
c.......................................................................

c     declarations
      include            'minfit.h'
      integer            nywin,ioni
      real*4             xxlo,xxhi,yylo,yyhi
      character*80       pltstr


c     loop over the number of elements define current window limits

      xxlo  = zlim(regioni,1)
      xxhi  = zlim(regioni,2)
      yylo  = -0.10
      yyhi  =  1.30
      nywin = ions

c     plot current window box

 12   CALL submargins(0.5,0.0)
      CALL xlabel(80,'Redshift\\e')
      CALL ylabel(80,'Normalized Flux\\e')

      do 11 ioni=1,nywin
       pltstr = ionlabel(ioni)
       CALL window(1,nywin,ioni)              
       CALL setlim(xxlo,yylo,xxhi,yyhi)
       CALL setlweight(0.5)
       if (ioni.eq.1) then
        CALL abox(1,2,-1,0)
       else
        CALL abox(-1,2,-1,0)
       end if
       CALL setltype(1)
       CALL relocate(xxlo,0.0)
       CALL draw(xxhi,0.0)
       CALL setltype(0)

       CALL setcolor(5) ! yellow
       CALL relocate(xxhi,0.0)
       CALL setexpand(1.2)
       CALL setlweight(1.0)
       CALL putlabel(lngth(ioni),pltstr,7)
       CALL setlweight(0.5)
       CALL setexpand(1.0)
       CALL setcolor(1)

       CALL setexpand(1.0)

 11   continue

      CALL tidle

      return
      end
c..
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c
      SUBROUTINE         pltdata
c
c
c     plot the data in histogram format
c
c.......................................................................

c     declarations
      include            'minfit.h'
      integer            i,nplt,ioni,nywin
      real*4             xx(maxpix),yy(maxpix)

c     set the number of sub windows
      nywin = ions
      CALL setlweight(0.5)
      CALL setltype(0)

c     loop over ions
      do 03 ioni=1,ions
       nplt = ndata(ioni)

c     covert wavelength to redshift and stuff data arrays
       do 11 i=1,nplt
        xx(i) = lambda(ioni,i)/lambda0(ioni) - 1.0d0
        yy(i) = data(ioni,i)
 11    continue

c     plot histogram style
       CALL window(1,nywin,ioni)              
       CALL histogram(xx,yy,nplt)
       CALL tidle

 03   continue ! next ion

c..reset line width and bail
      CALL setlweight(0.5)

      return
      end

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c
 
      SUBROUTINE pltreject

c
c     plot the rejected points in redshift space
c
c.......................................................................

      include            'minfit.h'
      integer            i,j,ioni,nywin,nplt
      real*4             xx(maxpix),yy(maxpix),style

  
c     set the line weight and line type (thin solid)

  
       CALL setlweight(0.5)
       CALL setltype(0)
       CALL setcolor(5) ! yellow

c     set number of sub windows

      nywin = ions

c     loop over ions

      DO 03 ioni=1,ions

c     count the number of masked pixels, stored as NPLT
       nplt = 0
       j    = 0

c     stuff the plotting arrays
       DO 11 i=1,ndata(ioni)
        IF (quality(ioni,i).eq.0) then 
         j = j + 1
         xx(j) = lambda(ioni,i)/lambda0(ioni) - 1.0d0
         yy(j) = data(ioni,i)
        END IF
 11    CONTINUE

c     set NPLT, and plot the rejected points
       nplt = j
       IF (nplt.gt.0) then
        style = 42.0
        CALL window(1,nywin,ioni)              
        CALL setangle(45.)
        CALL points(style,1,xx,yy,nplt)
        CALL setangle(0.)
        CALL tidle
       END IF

 03   CONTINUE ! next ion

c     back to defaults, clear and clean

      CALL setltype(0)
      CALL setcolor(1)
      CALL setlweight(0.5)
      CALL tidle
c

      RETURN
      END
 
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c

      SUBROUTINE         pltfitt

c
c     grab the new model and plot; the components are not plotted, only
c     the resulting fit
c
c.......................................................................
c
c     declarations
      include            'minfit.h'
      integer            i,ioni,nplt,nywin,linei
      real*4             xx(maxpix),yy(maxpix),zcen


c     set the number of sub windows
       nywin = ions

c     set line weight very thick
      CALL setltype(0)
      CALL setlweight(0.5)
      CALL setcolor(3)

c     convert wavelength to redshift, stuff the data arrays 
      do 03 ioni=1,ions
       if (features(ioni)) then
        nplt     = ndata(ioni)
        do 11 i=1,nplt
         xx(i) = lambda(ioni,i)/lambda0(ioni) - 1.0d0
         yy(i) = wrkflx(ioni,i)  
 11     continue

c     plot histogram style
        CALL window(1,nywin,ioni)
        CALL connect(xx,yy,nplt)
        CALL tidle

c     mark the components
        do 12 linei=1,lines
         zcen = zline(linei)
         CALL relocate(zcen,1.15)
         CALL draw(zcen,1.30)
         CALL tidle
 12     continue
       end if
 03   continue

c     bail

      return
      end

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c

      SUBROUTINE         pltoutput

c
c     plot the final model
c
c.......................................................................
c
      include           'minfit.h'
      include           'const.dek'
      logical           error
      integer           pixi
      integer           i,j,n,m,ioni,nplt,nywin,linei,imongo
      real*4            xx(maxions*maxpix),yy(maxions*maxpix),vcen
      real*4            vmin,vmax,zmin,zmax,yylo,yyhi,xpos,ypos,z
      real*4            style
      double precision  wave,vel,flux,noise1,noise2,cont
      double precision  a(maxcoeffs),siga(maxcoeffs)
      character*80      pltstr,ion_str



c     an important step is the use of SPECIES, which must be reset to
c     the TOTAL number of species

      species = total


c     obtain the overall plotting limits, and number of windows

      vmin = ckms * (zlim(regions+1,1)-zabs)/(1.0d0+zabs) - 20.0
      vmax = ckms * (zlim(regions+1,2)-zabs)/(1.0d0+zabs) + 20.0
      yylo = -0.10
      yyhi =  1.30
      nywin = ions


c     plot the whole kabootle, first the windows and labels loop over
c     the number of ions (including those not being fit)

      CALL erase
      CALL submargins(0.5,0.0)
      CALL xlabel(80,'Rest-Frame Velocity\\e')
      CALL ylabel(80,'Normalized Flux\\e')

      do 03 ioni=1,nywin
       CALL setcolor(1)
       CALL window(1,nywin,ioni)              
       CALL setlim(vmin,yylo,vmax,yyhi)
       CALL setlweight(0.5)
       if (ioni.eq.1) then
        CALL abox(1,2,-1,0)
       else
        CALL abox(-1,2,-1,0)
       end if
       CALL setltype(1)
       CALL relocate(vmin,0.0)
       CALL draw(vmax,0.0)
       CALL setltype(0)
       pltstr = ionlabel(ioni)
       CALL setcolor(5)
       CALL relocate(vmax,0.1)
       CALL setexpand(1.2)
       CALL setlweight(1.0)
       CALL putlabel(lngth(ioni),pltstr,7)
       CALL setlweight(0.5)
       CALL setexpand(1.0)
       CALL setcolor(1)
       CALL setexpand(1.0)
       CALL tidle
 03   continue

c     read in the data, set fitting flags and adjust the current number
c     of species within, if need be species if no detected features then
c     return (should not happen!!!!)

       CALL readdata(0,error)
       if (error) GOTO 99

c     reset the number of functions (we are done fitting)

       m = 0
       DO 19 ioni=1,ions
       m = m + ndata(ioni)
 19    CONTINUE

c     set up the FFT

       IF (convolving) CALL initconv(m)

c     stuff the WRKFLX array, using the results routine

       CALL outputmod(m)


c     NOW- DO THE PLOTTING OF THE DATA, SIGMA, MODEL, AND MASKS

c     loop over ions and plot

       do 07 ioni=1,ions

c     plot the windows
        CALL setcolor(1) ! white on black
        CALL window(1,nywin,ioni)              
        nplt = ndata(ioni)

c     plot the data
         do 25 i=1,nplt
          z     = lambda(ioni,i)/lambda0(ioni) - 1.0d0
          xx(i) = ckms * (z-zabs)/(1.0d0+zabs)
          yy(i) = data(ioni,i)
 25      continue
         CALL histogram(xx,yy,nplt)
         CALL tidle

c     plot the model
         CALL setcolor(3) ! green
         do 27 i=1,nplt
          yy(i) = wrkflx(ioni,i)
 27      continue
         CALL connect(xx,yy,nplt)
         CALL tidle

c     plot the cloud positions as ticks
         do 29 linei=1,vpcomps
          vcen = ckms * (zfinal(linei)-zabs)/(1.0d0+zabs)
          CALL relocate(vcen,1.15)
          CALL draw(vcen,1.30)
          CALL tidle
 29      continue

c     plot the sigma spectrum
         CALL setcolor(2) ! red
         do 31 i=1,nplt
          yy(i) = sigma(ioni,i)
 31      continue
         CALL histogram(xx,yy,nplt)
         CALL tidle

c     plot the masked data if they exists

       j    = 0

       DO 11 i=1,ndata(ioni)
        IF (quality(ioni,i).eq.0) then 
         j = j + 1
         z     = lambda(ioni,i)/lambda0(ioni) - 1.0d0
         xx(j) = ckms * (z-zabs)/(1.0d0+zabs)
         yy(j) = data(ioni,i)
        END IF
 11    CONTINUE

       IF (j.gt.0) then
        CALL setcolor(5) ! yellow
        style = 42.0
        CALL window(1,nywin,ioni)              
        CALL setangle(45.)
        CALL points(style,1,xx,yy,j)
        CALL setangle(0.)
        CALL tidle
       END IF

 07    continue  ! next ion

 99   continue   ! next region


c     return

      return
      end
c

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c

      SUBROUTINE         pltinput

c
c   plot the input model
c
c.......................................................................
c
c
      include           'minfit.h'
      include           'const.dek'
      logical           error
      integer           i,j,n,m,ioni,nplt,nywin,linei,imongo
      real*4            xx(maxions*maxpix),yy(maxions*maxpix),vcen
      real*4            vmin,vmax,zmin,zmax,yylo,yyhi,xpos,ypos,z
      real*4            style
      double precision  a(maxcoeffs),siga(maxcoeffs)
      character*80      pltstr




c     the TOTAL number of species

      species = total

c     obtain the overall plotting limits, and number of windows

      vmin = ckms * (zlim(regions+1,1)-zabs)/(1.0d0+zabs) !- 20.0
      vmax = ckms * (zlim(regions+1,2)-zabs)/(1.0d0+zabs) !+ 20.0
      yylo = -0.10
      yyhi =  1.30
      nywin = ions
c
c
c      
c  plot the whole kabootle, first the windows and labels
c  loop over the number of ions (including those not being fit)
c
      CALL erase
      CALL submargins(0.5,0.0)
      CALL xlabel(80,'Rest-Frame Velocity\\e')
      CALL ylabel(80,'Normalized Flux\\e')
c
      do 03 ioni=1,nywin
       CALL setcolor(1)
       CALL window(1,nywin,ioni)              
       CALL setlim(vmin,yylo,vmax,yyhi)
       CALL setlweight(0.5)
       if (ioni.eq.1) then
        CALL abox(1,2,-1,0)
       else
        CALL abox(-1,2,-1,0)
       end if
       CALL setltype(1)
       CALL relocate(vmin,0.0)
       CALL draw(vmax,0.0)
       CALL setltype(0)
       pltstr = ionlabel(ioni)
       CALL setcolor(5)
       CALL relocate(vmax,0.1)
       CALL setexpand(1.2)
       CALL setlweight(1.0)
       CALL putlabel(lngth(ioni),pltstr,7)
       CALL setlweight(0.5)
       CALL setexpand(1.0)
       CALL setcolor(1)
       CALL setexpand(1.0)
       CALL tidle
 03   continue


c     now the windows are in place, loop over the regions and plot them
c     one by one

      do 99 regioni=1,regions

c     read in the data, set fitting flags and adjust the current number
c     of species within, if need be species if no detected features then
c     return (should not happen!!!!)

       CALL readdata(0,error)
       if (error) GOTO 99

c     get the lines, this is the call before the LSF starts
 
       CALL getlines

c     compute the number of funtions

       m = 0
       DO 19 ioni=1,ions
       m = m + ndata(ioni)
 19    CONTINUE

c     set up the FFT?

       IF (convolving) CALL initconv(m)

c     stuff the WRKFLX array, using the inputmod routine

       CALL inputmod(m)

c     loop over ions and plot

       do 07 ioni=1,ions
        if (features(ioni)) then
        CALL setcolor(1)
        CALL window(1,nywin,ioni)              
         nplt = ndata(ioni)
         do 25 i=1,nplt
          z     = lambda(ioni,i)/lambda0(ioni) - 1.0d0
          xx(i) = ckms * (z-zabs)/(1.0d0+zabs)
          yy(i) = data(ioni,i)
 25      continue
         CALL histogram(xx,yy,nplt)
         CALL tidle

c     plot the user's initial model

         CALL setcolor(3) ! green
         do 27 i=1,nplt
          yy(i) = wrkflx(ioni,i)
 27     continue
         CALL connect(xx,yy,nplt)
         CALL tidle

c     plot the line center as ticks

         do 29 linei=1,lines
           vcen = ckms * (zline(linei)-zabs)/(1.0d0+zabs)
           CALL relocate(vcen,1.15)
           CALL draw(vcen,1.30)
           CALL tidle
 29       continue

c     plot the sigma spectra

         CALL setcolor(2) ! red
         do 31 i=1,nplt
          yy(i) = sigma(ioni,i)
 31      continue
         CALL histogram(xx,yy,nplt)
         CALL tidle
        end if

c     plot the masked data if they exists

       j    = 0

       DO 11 i=1,ndata(ioni)
        IF (quality(ioni,i).eq.0) then 
         j = j + 1
         z     = lambda(ioni,i)/lambda0(ioni) - 1.0d0
         xx(j) = ckms * (z-zabs)/(1.0d0+zabs)
         yy(j) = data(ioni,i)
        END IF
 11    CONTINUE

       IF (j.gt.0) then
        CALL setcolor(5) ! yellow
        style = 42.0
        CALL window(1,nywin,ioni)              
        CALL setangle(45.)
        CALL points(style,1,xx,yy,j)
        CALL setangle(0.)
        CALL tidle
       END IF

 07    continue ! next ion

 99   continue  ! next region


c     return

      return
      end
c
c     eof



