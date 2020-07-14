c------------------------------------------------------------------------------
c

      SUBROUTINE        getlines

c
c     read in fit lines from the AUTOVP / IVPFIT program output
c
c==============================================================================
c
      include           'minfit.h'
      include           'const.dek'
      logical           flag
      integer           ioni,i,j,k,nfits,idx,totlines,linei
      double precision  vline,dumm,col,bpar,z
      character*80      in_file,ion_str,label(2)
      data label        /'logN','Dopp'/


c     zero the arrays for good measure, this way the user can
c     immediately see which lines are held at 0.0, since they are
c     statistically meaningless

       do 03 i=1,maxspecies
        do 05 linei=1,maxlines
         zline(linei)    = 0.0d0
         nline(i,linei)  = 0.0d0
         bline(i,linei)  = 0.0d0
         dzline(linei)   = 0.0d0
         dnline(i,linei) = 0.0d0
         dbline(i,linei) = 0.0d0
 05     continue
 03    continue

      in_file = profitfile

      write(SCREEN,603) regioni,zlim(regioni,1),zlim(regioni,2)
      write(STDOUT,603) regioni,zlim(regioni,1),zlim(regioni,2)

c     the big assumption here is that the .pro files are in tidy order -
c     CHECK this for each system: MOD July 2006 - .pro files now store
c     column densities in log10 but program still works in e13

      OPEN(unit=1,file=in_file,err=999,status='old')

      read(1,*) nfits,dumm,dumm

      do 09 i=1,nfits

       flag = .false.
       read(1,*) ion_str,totlines

       do 11 ioni=1,ions

        if (ion_str(1:4).eq.ion_name(ioni)(1:4)) then
         flag = .true.
         idx   = spec_idx(ioni)
         linei = 0

         do 30 j=1,totlines
          read(1,*) k,col,vline,bpar
          z = zabs + (1.0d0+zabs)*vline/ckms
          if ((z.ge.zlim(regioni,1)).AND.(z.le.zlim(regioni,2))) then
           linei            = linei + 1
           zline(linei)     = z
           nline(idx,linei) = col
           bline(idx,linei) = bpar

c     flag if the transition has no detections for later output

          if (.not.features(ioni)) dnline(idx,linei) = -1.0d0

          end if 

 30      continue

         GOTO 09

        end if

 11    continue

c     read the next series and discard

       if (.not.flag) then      
        do 31 j=1,totlines
         read(1,*) k,col,vline,bpar
 31     continue
       end if

 09   continue

      CLOSE(unit=1)

      lines = linei

c     now put the lines in redshift order for ease they are not
c     necessarily in redshift order from PROFIT

      CALL sortlines

c     communicate the sorted input model output the column density in
c     log10

      write(SCREEN,600) (spec_name(j),j=1,total)
      write(SCREEN,601) (label(1),label(2),j=1,total)
      write(STDOUT,600) (spec_name(j),j=1,total)
      write(STDOUT,601) (label(1),label(2),j=1,total)
      do 23 linei=1,lines
       write(SCREEN,602) linei,zline(linei), 
     &       (nline(j,linei),bline(j,linei),j=1,total)
       write(STDOUT,602) linei,zline(linei), 
     &       (nline(j,linei),bline(j,linei),j=1,total)
 23   continue
      write(SCREEN,*) ' '
      write(SCREEN,'(a,i3)') ' Number of LINES = ',lines
      write(SCREEN,*) ' '
      write(STDOUT,*) ' '
      write(STDOUT,'(a,i3)') ' Number of LINES = ',lines
      write(STDOUT,*) ' '

c     return

      return

c     error traps

 999  write(SCREEN,*) ' NOT FOUND- file: ',in_file(1:30)
      write(STDOUT,*) ' NOT FOUND- file: ',in_file(1:30)
      STOP ' *** MINFIT(getlines): terminated ***'
c     formats

 600  format(1x,t22,10a15)
 601  format(1x,t7,'redshift',t19,10(a7,a8))
 602  format(1x,i3,2x,f8.6,10(3x,f5.2,2x,f5.2))
 603  format(1x,/,1x,'--- STARTING MODEL --- REGION:',i3,
     @          '   (',f8.6,',',f8.6,')')

      end

c..............................................................................
c

      SUBROUTINE        sortlines

c
c     we call this routine twice 
c
c     before LSF: from routine GETLINES(getlines.f)
c
c     after  LSF: from routine UNPACKA(defcoeffs.f)
c
c     after fitting we need to also keep track of the uncertainties.  so
c     we do this too
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h' 
      integer           i,j,k,id,nl,indx(maxlines)
      double precision  col(maxspecies,maxlines),zcl(maxlines),
     &                  dcol(maxspecies,maxlines),dzcl(maxlines),
     &                  b(maxspecies,maxlines),
     &                  db(maxspecies,maxlines)

c     if we have a single line, the no need to sort!  (also avoids seg
c     fault!)

      IF (lines.eq.1) RETURN

c     store all to avoid overwriting, also it is not allowed to pass
c     common variables, so we need to handle that for the call to
c     routine INDEXX

      do 12 j=1,lines
       zcl(j)  = zline(j)
       dzcl(j) = dzline(j)
       do 05 k=1,total
        col(k,j)  = nline(k,j)
        dcol(k,j) = dnline(k,j)
        b(k,j)    = bline(k,j)
        db(k,j)   = dbline(k,j)
 05    continue
 12   continue

c     get the index ordering of the zline array

      nl = lines

      CALL indexx(nl,zcl,indx)

c     now sort by redshift

      do 22 j=1,lines
       id = indx(j)
       zline(j) = zcl(id)
       dzline(j) = dzcl(id)
       do 15 k=1,total
        nline(k,j)  = col(k,id)  
        dnline(k,j) = dcol(k,id) 
        bline(k,j)  = b(k,id)    
        dbline(k,j) = db(k,id)   
 15    continue
 22   continue

c     return

      return
      end

c     eof
