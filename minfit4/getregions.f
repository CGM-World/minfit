c------------------------------------------------------------------------------
c

      SUBROUTINE        getregions

c
c     read in the ew_regions.dat file and store the redshifts of the
c     regions where features have been identified
c
c------------------------------------------------------------------------------
c

      include           'minfit.h'
      integer           i,j
      double precision  y,w0,w1,w2,reflam,zmin,zmax,vmin,vmax
      include           'const.dek'


      write(SCREEN,*) '************************************************'
      write(SCREEN,*) ' ' 
      write(SCREEN,*) 'FITTING REGIONS '
      write(SCREEN,*) ' ' 
      write(SCREEN,600)

      write(STDOUT,*) '************************************************'
      write(STDOUT,*) ' ' 
      write(STDOUT,*) 'FITTING REGIONS '
      write(STDOUT,*) ' ' 
      write(STDOUT,600)


      zmax = -1.0d0
      zmin = 1000.

      OPEN(unit=1,file='ew_regions.dat',err=999,status='old')
      do 11 i=1,maxlines
       read(1,*,end=6) j,y,w1,w2,w0,reflam
       zlim(i,1) = w1/reflam - 1.0d0
       zlim(i,2) = w2/reflam - 1.0d0
       vmin      = ckms*(zlim(i,1)-zabs)/(1.0d0+zabs)
       vmax      = ckms*(zlim(i,2)-zabs)/(1.0d0+zabs)
       regions   = regions + 1
       write(SCREEN,601) regions,zlim(regions,1),zlim(regions,2),
     &                   vmin,vmax
       write(STDOUT,601) regions,zlim(regions,1),zlim(regions,2),
     &                   vmin,vmax
       zmin = min(zmin,zlim(i,1))
       zmax = max(zmax,zlim(i,2))
 11   continue
 6    CLOSE(unit=1)

      write(SCREEN,*) ' ' 
      write(SCREEN,'(a,i3)') ' NUMBER OF REGIONS = ',regions
      write(SCREEN,*) ' ' 

      write(STDOUT,*) ' ' 
      write(STDOUT,'(a,i3)') ' NUMBER OF REGIONS = ',regions
      write(STDOUT,*) ' ' 

c     determine the full region spread; for initial and final plotting
c     purposes, we buffer approximately 90 km/s for the final plotting

      zlim(regions+1,1) = zmin - 0.0003d0/(1.0d0+zabs)
      zlim(regions+1,2) = zmax + 0.0003d0/(1.0d0+zabs)

c     return

      return

c     error trap

 999  write(SCREEN,*) ' cannot locate file : ew_regions.dat'
      write(STDOUT,*) ' cannot locate file : ew_regions.dat'

      STOP ' *** MINFIT(getregions): terminated ***'

 600  format(1x,t3,'Reg',t12,'z_lo',t22,'z_hi',
     &       t35,'v_lo',t45,'v_hi [km/s]')
 601  format(1x,i4,3x,2f10.6,2f10.2)

      end
c
c------------------------------------------------------------------------------
c

      SUBROUTINE        getsatreg

c
c     read in the sat_region.dat file and store the redshifts of the
c     regions and species where the data are saturated; if the file DNE
c     null all counters; the file format is
c
c
c     Tiespec  Anchor   vmin   vmax
c
c     where  tiespec = species number (in ions.table file) 
c            anchor  = species number to anchor column density to (usually 1)
c            vmin    = lower velocity bound for this tiespec/anchor combo
c            vmax    = upper velocity bound for this tiespec/anchor combo
c
c     there can be multiple tiespec/anchor combos, one for each line in
c     the file even if they are in the same or similar velocity windows
c
c------------------------------------------------------------------------------
c

      include           'minfit.h'
      include           'const.dek'
      integer           i,j1,j2
      double precision  v1,v2

c     communicate

      write(SCREEN,*) '************************************************'
      write(SCREEN,*) ' '
      write(SCREEN,*) 'SATURATION REGIONS'
      write(SCREEN,*) ' '

      write(STDOUT,*) '************************************************'
      write(STDOUT,*) ' '
      write(STDOUT,*) 'SATURATION REGIONS'
      write(STDOUT,*) ' '

c     initialize

      nsatreg = 0

      do 03 i=1,maxregions
       sat_spec(i,1) = 0
       sat_spec(i,2) = 0
       zsat(i,1) = 0.0d0
       zsat(i,2) = 0.0d0
 03   continue

c     try to open the file, if DNE trap, communicate, and return

      OPEN(unit=1,file='sat_regions.dat',err=999,status='old')

      write(SCREEN,600)
      write(STDOUT,600)

c     OK, read in the file

      do 11 i=1,maxregions
       read(1,*,end=6) j1,j2,v1,v2
       sat_spec(i,1) = j1
       sat_spec(i,2) = j2
       zsat(i,1) = v1*(1.0d0+zabs)/ckms + zabs
       zsat(i,2) = v2*(1.0d0+zabs)/ckms + zabs
       nsatreg   = nsatreg + 1
       write(SCREEN,601) nsatreg,spec_name(j1),spec_name(j2),
     &                   zsat(i,1),zsat(i,2),acol(j1),bcol(j1)
       write(STDOUT,601) nsatreg,spec_name(j1),spec_name(j2),
     &                   zsat(i,1),zsat(i,2),acol(j1),bcol(j1)
 11   continue

 6    CLOSE(unit=1)

c     communicate

      write(SCREEN,*) ' ' 
      write(SCREEN,'(a,i3)') ' NUMBER OF SAT REGIONS = ',nsatreg
      write(SCREEN,*) ' ' 
      write(SCREEN,*) '************************************************'

      write(STDOUT,*) ' ' 
      write(STDOUT,'(a,i3)') ' NUMBER OF SAT REGIONS = ',nsatreg
      write(SCREEN,*) ' ' 
      write(STDOUT,*) '************************************************'

      return

c     error trap

 999  nsatreg  = 0

      write(SCREEN,*) ' ' 
      write(SCREEN,*) 'NO SAT REGIONS'
      write(SCREEN,*) ' ' 
      write(SCREEN,*) '************************************************'

      write(STDOUT,*) ' ' 
      write(STDOUT,*) 'NO SAT REGIONS'
      write(STDOUT,*) ' ' 
      write(STDOUT,*) '************************************************'


      return

 600  format(1x,t3,'Idx',t10,'Tie',t16,'Anchor',
     &          t25,'z_low',t35,'z_high',t45,'Acol',t55,'Bcol')
 601  format(1x,i4,2x,a6,2x,a6,4f10.6)

      end
c
c------------------------------------------------------------------------------
c

      SUBROUTINE        getflagions

c
c
c------------------------------------------------------------------------------
c

      include           'minfit.h'
      integer           i,j,idx,regi
      double precision  dum


      write(SCREEN,*) ' '
      write(SCREEN,*) '************************************************'
      write(SCREEN,*) ' '
      write(SCREEN,*) 'CONSTRAINT IONS FOR REGION',regioni
      write(SCREEN,*) ' '

      write(STDOUT,*) ' '
      write(STDOUT,*) '************************************************'
      write(STDOUT,*) ' '
      write(STDOUT,*) 'CONTRAINT IONS FOR REGION',regioni
      write(STDOUT,*) ' '


c     null the flagion 

      flagion(regioni) = 0

c     open the FITREGIONS file, if it DNE skip to line 30

      OPEN(unit=1,file='fitregions',err=30,status='old')

c     we have a successful read, now stuff the flagion array

c     the file format of file FITREGIONS is "regioni indexi", if not all
c     regions are present the file ends early and we jump out of the
c     loop and set IDX=0

      do 11 i=1,regions
       read(1,*,end=20,err=991) regi,idx
       IF (regi.eq.regioni) then
         flagion(regioni) = idx
         GOTO 20
       END IF
 11   continue

 20   CLOSE(unit=1)

c     perform sanity check that this ion has representative transitions
c     (i.e., that are being included in the fit!)

      DO 13 i=1,species
       IF (fitindex(i).eq.flagion(regioni)) then ! all is well
        write(SCREEN,600) regioni,spec_name(idx)
        write(STDOUT,600) regioni,spec_name(idx)
        RETURN
       END IF
 13   CONTINUE

c     user choked, this ion is not being used in the fitting, null
c     FLAGION

      flagion(regioni) = 0

c     file does not exists, FLAGION already nulled above

 30   write(SCREEN,601) regioni
      write(STDOUT,601) regioni

      RETURN

c     error trap

 991  write(SCREEN,*) ' ERROR(getflagions): cannot read',
     &                ' file fitregions'
      write(STDOUT,*) ' ERROR(getflagions): cannot read',
     &                ' file fitregions'
      STOP 'MINFIT error- program terminated early'

c     formats

 600  format(1x,'REGION:', i3,' -- fit constrained by ion ',a20)
 601  format(1x,'REGION:', i3,' -- fit constrained by all ions')

      end

c     eof
