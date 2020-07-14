c------------------------------------------------------------------------------
c

      SUBROUTINE          getatomic
c      
c     check current directory for ions.dat, if DNE goto error trap,
c     which grabs the default file; loops over ions and stuffs the
c     atomic data arrays
c
c     IMPORTANT: the constants con1 and con2 are used to compute the
c     absorption coefficients (via the Voigt profile).  con1 is the
c     constant multiplying the Voigt profile integral and con2 is the
c     lorenztian width of the ion/transition.  Throughout, we carry the
c     column density in units of 1.e13 and wavelength in units of
c     angstroms.  Thus, there is a factor of 1.e5 scaling con1, which is
c     the product of 1.e13 for the column density and 1.e-8 for the
c     conversion of angstroms into cemtimeters.  Similarly, there is the
c     factor 1.e-8 in front of con2.
c
c..............................................................................
c

      include           'minfit.h'
      include           'const.dek'
      integer           i,ioni
      double precision  w,f,gamma,amu,amu1,ip,gam0(maxions)
      character*80      atom_list,mf_path,def_file,ion_str



c     we assume a the ionfile is given in the minfit.par file; it
c     contains the atomic data; if that does not exist, then we bail
      
      def_file  = ionfile

      write(SCREEN,*) '************************************************'
      write(SCREEN,*) ' '
      write(SCREEN,*) 'ATOMIC DATA'
      write(SCREEN,*) ' '
      write(SCREEN,600)

      write(STDOUT,*) '************************************************'
      write(STDOUT,*) ' '
      write(STDOUT,*) 'ATOMIC DATA'
      write(STDOUT,*) ' '
      write(STDOUT,600)

      OPEN(unit=1,file=def_file,err=999,status='old')

c     scan the atomic data file for each transition; we asssume that
c     file atomic data file is no longer then 1000 lines; if we exceed
c     that length, send a warning; we loop through the atomic data one
c     at a time and then for each entry, and check if it matches one of our
c     ions

 01   do 11 i=1,1000

       read(1,*,end=18) ion_str,w,f,gamma,amu,ip

       do 13 ioni=1,ions

c     found it, grab the wavelength, and create con1 and con2; close up
c     the whole shop

        if (ion_name(ioni).eq.ion_str) then

         lambda0(ioni) = w
         fosc(ioni)    = f
         gam0(ioni)    = 1.0d-8*gamma
         con1(ioni)    = 1.0d+5 * f * w**2 * sqrt(pi)*e*e / (me*c*c)
         con2(ioni)    = 1.0e-8 * gamma * w**2 / (4.0d0*pi*c)

         spec_mass(spec_idx(ioni)) = mp*amu 
         ionpot(spec_idx(ioni))    = ip

         GOTO 11 ! next entry in the atomic data file

        end if

 13    continue

 11   continue


c     if you are here, you are in trouble, the atomic data file is
c     longer than 1000 lines
      
      WRITE(SCREEN,*) ' WARNING(getatomic): atomic data file too long'
      WRITE(SCREEN,*) ' --- suggest you shorten it and begin again...'
      WRITE(STDOUT,*) ' WARNING(getatomic): atomic data file too long'
      WRITE(STDOUT,*) ' --- suggest you shorten it and begin again...'

 18   CLOSE(unit=1)


c     communicate, adopt the tubulence/thermmal condition by assigning
c     either the same mass to all ions or assigning them the real
c     masses; if turbulence is to be used, then we use the mass of the
c     first ion in the list

      DO 09 i=1,ions

        If (iturb.eq.1) spec_mass(spec_idx(i)) = spec_mass(spec_idx(1))

        WRITE(SCREEN,601) ion_name(i),lambda0(i),fosc(i),gam0(i),
     &                    (spec_mass(spec_idx(i))/mp),
     &                    ionpot(spec_idx(i))  
        WRITE(STDOUT,601) ion_name(i),lambda0(i),fosc(i),gam0(i),
     &                    (spec_mass(spec_idx(i))/mp),
     &                    ionpot(spec_idx(i))  

C       spec_mass(spec_idx(ioni)) = mp*amu1  ! keep all the same  

 09   CONTINUE

      write(SCREEN,*) ' '
      write(STDOUT,*) ' '

c     return

      RETURN

c  error traps

 999  WRITE(6,*) ' ERROR(getatomic): cannot find atomic data file'
      WRITE(6,*) ionfile
      STOP

c     formats

 600  format(1x,t5,'Ion/Tran',t19,'lambda',t30,'f-val',t38,'Gam8',
     &       t47,'mass',t56,'IPot')
 601  format(1x,3x,a10,f12.4,f8.4,f8.4,1x,f8.4,1x,f10.6)

      end

c     eof
