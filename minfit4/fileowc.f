c.........................................................................
c

      SUBROUTINE        openfiles

c
c     this routine opens the files that will be written to save the
c     results of the LSF
c
c     MINFIT.ZDAT file; contains the VP model parameters as a function
c     of redshift, it has a header
c
c     MINFIT.VDAT file; contains the VP model parameters as a
c     function of velocity, it has no header but the contents are the
c     same format as the minfit.zdat file
c
c     TICKS.DAT file contains the VP components as a function of
c     velocity, wavelength, and redshift, the first column is the
c     constant 1.25; this file is used for plotting the ticks at the VP
c     component positions over the spectra and VPMOD spectra
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      integer           j
      character*80      label(4)
      data label        /' logN','dlogN','  b  ','  db '/



c     open the ZDAT, VDAT, TICKS, and TDAT files

      OPEN(unit=70,file='minfit4.zdat',status='unknown')      
      OPEN(unit=71,file='minfit4.vdat',status='unknown')      
      OPEN(unit=72,file='ticks.dat',status='unknown')      
      OPEN(unit=73,file='minfit4.Tdat',status='unknown')      

c     write the header to the ZDAT file 

      write(70,600) (spec_name(j),spec_name(j),j=1,total)
      write(70,601) (label(1),label(2),label(3),
     @               label(4),j=1,total)

c     write the header to the TICKS file 

      write(72,602) (ion_name(j),j=1,ions)

c     formats

 600  format(1x,t22,10a15)
 601  format(1x,t7,'redshift',t18,10(a7,a8))
 602  format(1x,t3,'y',t12,'v',t22,'z',t33,20a12)

c     return

      return
      end

c
c.........................................................................
c

      SUBROUTINE        writefiles

c
c     write to designated files and to STDOUT; the files are explained
c     in the header comments to routine openfiles
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      include           'const.dek'

      integer           i,j
      double precision  vel,temp,dtemp
      character*80      label(4)
      data label        /' logN','dlogN','  b  ','  db '/



c     we are going to communicate to STDOUT as well; write the headers

      write(SCREEN,603) regioni,zlim(regioni,1),zlim(regioni,2)
      write(SCREEN,600) (spec_name(j),spec_name(j),j=1,total)
      write(SCREEN,601) (label(1),label(2),label(3),
     @                  label(4),j=1,total)
      write(STDOUT,603) regioni,zlim(regioni,1),zlim(regioni,2)
      write(STDOUT,600) (spec_name(j),spec_name(j),j=1,total)
      write(STDOUT,601) (label(1),label(2),label(3),
     @                  label(4),j=1,total)

c     write data to the .ZDAT, .VDAT, and TICKS.DAT file; the observed
c     wavelengths are output in the same order of ions in the fitlist
c     file


      do 53 i=1,lines
       vel = ckms * (zline(i)-zabs)/(1.0d0+zabs)
       write(SCREEN,602) i,zline(i),(nline(j,i),dnline(j,i),
     @               bline(j,i),dbline(j,i),j=1,total)
       write(STDOUT,602) i,zline(i),(nline(j,i),dnline(j,i),
     @               bline(j,i),dbline(j,i),j=1,total)
       write(70,602) i,zline(i),(nline(j,i),dnline(j,i),
     @               bline(j,i),dbline(j,i),j=1,total)
       write(71,604) i,vel,(nline(j,i),dnline(j,i),
     @               bline(j,i),dbline(j,i),j=1,total)
       write(72,605) 1.25,vel,zline(i),
     @               (lambda0(j)*(1.0+zline(i)),j=1,ions)
       temp = 1.0d10*spec_mass(spec_idx(1))*bline(1,i)**2 / (2.0d0*kerg)
       dtemp =  1.0d10*spec_mass(spec_idx(1))*bline(1,i)*dbline(1,i) 
     @          / (2.0d0*kerg)
       write(73,606) i,zline(i),vel,temp,dtemp,dtemp/temp
 53   continue

      WRITE(SCREEN,*) ' '
      WRITE(STDOUT,*) ' '

c     formats

 600  format(1x,t22,10a15)
 601  format(1x,t7,'redshift',t18,10(a7,a8))
 602  format(1x,i3,2x,f8.6,10(2x,f6.2,1x,f6.2))
 603  format(1x,/,1x,'LSF FIT MODEL --- REGION:',i2,
     @       ' z=(',f8.6,',',f8.6,')',/)
 604  format(1x,i3,2x,f9.2,10(2x,f6.2,1x,f7.3))
 605  format(1x,f4.2,1x,f10.4,1x,f10.6,1x,20f12.4)
 606  format(1x,i3,2x,f8.6,2x,f9.2,3(1pe12.3))

c     return

      return
      end

c
c.........................................................................
c

      SUBROUTINE        closefiles

c
c     this routine called after all regions are completed; close the
c     files that were opened in routine openfiles, and written to in
c     routine writefiles
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c


c     close the ZDAT, VDAT, TICKS, and TDAT files

      CLOSE(unit=70)
      CLOSE(unit=71)
      CLOSE(unit=72)
      CLOSE(unit=73)

      RETURN

      END

c
c.........................................................................
c

      SUBROUTINE        vpmods

c
c     VPMOD files; we create a file with the VP model spectrum and its
c     residuals about the dat for each ion/transition for which we have
c     obtained a fit; i.e., those that have detected 
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      include           'const.dek'

      integer           ioni,pixi,linei
      double precision  wdum,vdum,one,zip
      double precision  w,w0,vel,res,ressig
      character*80      ion_str,outfile



      one = 1.0d0
      zip = 0.0d0

c     compute pre-padding

      wdum = 0.0
      vdum = -9999.9999d0

c     open the VPMOD files; one for each ion/transition; write
c     pre-padding

      do 07 ioni=1,ions
        ion_str = ion_name(ioni)
        call fappend(ion_str,'vpmod',outfile)
        OPEN(unit=10+ioni,file=outfile,status='unknown')
        WRITE(10+ioni,101) wdum,vdum,one,zero,zip
C     @                  (1.0d0,linei=1,vpcomps)         
 07   continue

c     compute the residuals and the ratio of the residual to the
c     uncertainty; write the results

      do 11 ioni=1,ions
        w0  = lambda0(ioni)*(1.0d0+zabs)
        do 09 pixi=1,ndata(ioni)
         w      = lambda(ioni,pixi) 
         vel    = ckms * (w-w0)/w0 
         res    = data(ioni,pixi)-wrkflx(ioni,pixi)
         ressig = abs(res/sigma(ioni,pixi))
         write(10+ioni,101) w,vel,wrkflx(ioni,pixi),res,ressig
C     @                  (modflx(ioni,linei,pixi),linei=1,vpcomps)         
 09     continue
 11   continue

c     compute and write the post-padding; close the VPMOD files one by
c     one

      wdum = 99999.9999d0
      vdum =  9999.9999d0

      do 15 ioni=1,ions
        WRITE(10+ioni,101) wdum,vdum,one,zero,zero
C     @                  (1.0d0,linei=1,vpcomps)         
        CLOSE(unit=10+ioni)
 15   continue

c     format

 101  format(1x,f10.4,1x,f10.4,1x,1p3e12.4,1p40e12.4)


      RETURN
      END
