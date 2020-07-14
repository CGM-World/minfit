


c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  
  
      SUBROUTINE         inform(info,error)

c  
c.......................................................................
c

      include            'minfit.h'
      logical            error
      integer            info


c     the following messages are paraphrased from the dnl1se instruction
c     header, commented in the dnl1se.f source code file
  
       write(SCREEN,*) ' '
       write(STDOUT,*) ' '


      if (info.eq.0) then
       write(SCREEN,*) ' improper input parameters.'
       write(STDOUT,*) ' improper input parameters.'
       error = .true.
      end if
  
      if (info.eq.1) then
       write(SCREEN,*) ' **** 1st convergence condition satisfied ****'
       write(STDOUT,*) ' **** 1st convergence condition satisfied ****'
      end if
  
      if (info.eq.2) then
       write(SCREEN,*) ' **** 2nd convergence condition satisfied ****'
       write(STDOUT,*) ' **** 2nd convergence condition satisfied ****'
      end if
  
      if (info.eq.3) then
       write(SCREEN,*) ' * 1st & 2nd convergence condition satisfied *'
       write(STDOUT,*) ' * 1st & 2nd convergence condition satisfied *'
      end if
  
      if (info.eq.4) then
       write(SCREEN,*) ' **** 3rd convergence condition satisfied ****'
       write(STDOUT,*) ' **** 3rd convergence condition satisfied ****'
      end if
  
      if (info.eq.5) then
       write(SCREEN,*) ' *** convergence conditions not satisfied ***'
       write(SCREEN,*) '     INFO = 5; stopped at maximum iterations '
       write(STDOUT,*) ' *** convergence conditions not satisfied ***'
       write(STDOUT,*) '     INFO = 5; stopped at maximum iterations '
      end if
  
      if (info.eq.6) then
       write(SCREEN,*) ' *** convergence conditions not satisfied ***' 
       write(SCREEN,*) '     INFO = 6; ftol too tight '
       write(STDOUT,*) ' *** convergence conditions not satisfied ***' 
       write(STDOUT,*) '     INFO = 6; ftol too tight '
      end if
  
      if (info.eq.7) then
       write(SCREEN,*) ' *** convergence conditions not satisfied ***' 
       write(SCREEN,*) '     INFO = 7; xtol too tight '
       write(STDOUT,*) ' *** convergence conditions not satisfied ***' 
       write(STDOUT,*) '     INFO = 7; xtol too tight '
      end if

      if (info.eq.8) then
       write(SCREEN,*) ' *** convergence conditions not satisfied ***' 
       write(SCREEN,*) '     INFO = 8; gtol too tight '
       write(STDOUT,*) ' *** convergence conditions not satisfied ***' 
       write(STDOUT,*) '     INFO = 8; gtol too tight '
      end if
  
       write(SCREEN,*) ' '
       write(STDOUT,*) ' '

c     return
  
      return
      end

c
c..............................................................................
c

      SUBROUTINE         modelcomm(a,n,icall)

c  
c     communicate the fit results prior to getting the limits
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      include           'const.dek'

      integer           n,i,k,jj(maxspecies),fidx
      integer           in,ib,iz,icall
      character*80      label(2)
      double precision  id_isnan,Naction
      double precision  a(maxcoeffs),sngb,satNscale
      double precision  col(maxspecies,maxlines),b(maxlines,maxlines)
      double precision  z(maxlines),Tp

      data label        /'logN','Dopp'/


c     set Naction for call to satNscale

      Naction = 0 ! log10 units

c     since COL, B, and Z are local array and we are using them only to
c     communicate, we can index them directly by 1->SPECIES, thus, we
c     just use streight indexing; but we do need to keep track of the
c     FITINDEX for certain unfolding information

      do 13 i=1,lines
       do 11 k=1,species

         fidx  = fitindex(k)
         jj(k) = fidx
         in    = idx_an(k,i)
         ib    = idx_ab(i)
         iz    = idx_az(i)

c     unpack column densities; avoid logs of negative numbers

         if (satflag(k,i)) then
           col(k,i) = satNscale(a(in),acol(fidx),bcol(fidx),Naction)
C          col(k,i) = acol(fidx)*(13.+log10(abs(a(in))))+bcol(fidx)
           col(k,i) = col(k,i)*sign(1.0d0,a(in))
         else
          col(k,i) = 13.0d0+log10(abs(a(in)))
          col(k,i) = col(k,i)*sign(1.0d0,a(in))
         end if

c     unpack b parameters, convert to km/s and carry the sign

         sngb   = sign(1.0d0,a(ib))
         Tp     = abs(a(ib))
         b(k,i) = sqrt(2.0d0*kerg*Tp/spec_mass(fidx))
         b(k,i) = 1.0d-5*sngb*b(k,i) ! km/s 

 11    continue

c     unpack redshifts

       z(i) = a(iz)

 13   continue


c     icall=0 means we are reporting the model being sent to the LSF
c     engine

      IF (icall.eq.0) then 
      write(SCREEN,603) regioni,zlim(regioni,1),zlim(regioni,2)
      write(STDOUT,603) regioni,zlim(regioni,1),zlim(regioni,2)
      END IF

c     icall=1 means we are reporting the resulting LSF model after we
c     have returned from the LSF engine

      IF (icall.eq.1) then 
      write(SCREEN,604) regioni,zlim(regioni,1),zlim(regioni,2)
      write(STDOUT,604) regioni,zlim(regioni,1),zlim(regioni,2)
      END IF

c     icall=2 means we are reporting the model every NPRINT iterations
c     from the LSF engine

      IF (icall.eq.2) then 
      write(SCREEN,605) regioni,zlim(regioni,1),zlim(regioni,2)
      write(STDOUT,605) regioni,zlim(regioni,1),zlim(regioni,2)
      END IF

c     now report the model parameters, first the header, and then the
c     parameters themselves

      write(SCREEN,600) (spec_name(jj(k)),k=1,species)
      write(SCREEN,601) (label(1),label(2),k=1,species)
      write(STDOUT,600) (spec_name(jj(k)),k=1,species)
      write(STDOUT,601) (label(1),label(2),k=1,species)
      do 17 i=1,lines
        write(SCREEN,602) i,z(i),(col(k,i),b(k,i),k=1,species)
        write(STDOUT,602) i,z(i),(col(k,i),b(k,i),k=1,species)
 17   continue

c     we are done, return

      return

c     formats

 600  format(1x,t22,10a15)
 601  format(1x,t7,'redshift',t19,10(a7,a8))
 602  format(1x,i3,2x,f8.6,10(2x,f6.2,1x,f6.2))
 603  format(1x,/,1x,'--- CURRENT INPUT MODEL --- REGION:',i3,
     @          '   (',f8.6,',',f8.6,')')
 604  format(1x,/,1x,'--- CURRENT LSF MODEL --- REGION:',i3,
     @          '   (',f8.6,',',f8.6,')')
 605  format(1x,/,1x,'--- REPORT OF CURRENT LSF MODEL --- REGION:',i3,
     @          '   (',f8.6,',',f8.6,')')

      end

c
c..............................................................................
c     eof
