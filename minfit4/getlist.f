c.......................................................................
c 

       SUBROUTINE         getlist(ionlist)

c
c     reads in ionlist from the command line, sets the number of species
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

       include            'minfit.h'
       integer            i,j,k,l,idxnum,jstart
       double precision   atmp,btmp
       character*80       ionlist,ion_str,level(3),pltstr

       data level         /'I','V','X'/



       total = 0

c     open and read in the master fitting list; this file has the format
c     $1=transition name 
c     $2=species index number 
c     $3=log column scaling slope            log10scaled = $3*log(N)+$4
c     $4=log column scaling intercept
c     the log scaling is applied in regions defined to be saturated

       OPEN(unit=1,file=ionlist,err=999,status='old')
       do 05 i=1,maxions
        read(1,*,end=6) ion_str,idxnum,atmp,btmp
        ion_name(i)       = ion_str
        spec_idx(i)       = idxnum
        spec_name(idxnum) = ion_name(i)(1:4)
        acol(idxnum)      = atmp
        bcol(idxnum)      = btmp
        total             = max(total,idxnum)
        ions              = ions + 1
 05    continue
 06    CLOSE(unit=1)


c     construct the species names for pleasant communication to the user
c     most species names are two characters but not all, so search for
c     those (i.e., C=carbon, O=oxygen, N=nitrogen)

       do 07 i=1,total
         if ((spec_name(i)(2:2).eq.level(1)).OR.
     @       (spec_name(i)(2:2).eq.level(2)).OR.
     @       (spec_name(i)(2:2).eq.level(3))) then
           ion_str = spec_name(i)(1:1)//' '   ! if true, then this is 1 character species  
           jstart = 2
         else
           ion_str = spec_name(i)(1:2)//' '   ! else this is 2 character species  
           jstart = 3
         end if
        do 09 j=jstart,20
         if ((spec_name(i)(j:j).eq.level(1)).OR.
     @       (spec_name(i)(j:j).eq.level(2)).OR.
     @       (spec_name(i)(j:j).eq.level(3))) then
          ion_str = ion_str(1:j)//spec_name(i)(j:j)
         end if
 09     continue
        spec_name(i) = ion_str
 07    continue  

c     construct the ion/transition names for labels

       do 17 i=1,ions
        pltstr = ion_name(i)(1:2)//' '
        l = 3
        do 19 j=3,8
         if ((ion_name(i)(j:j).eq.level(1)).OR.
     @       (ion_name(i)(j:j).eq.level(2)).OR.
     @       (ion_name(i)(j:j).eq.level(3))) then
          pltstr = pltstr(1:l)//ion_name(i)(j:j)
          l = l + 1
         else
          pltstr = pltstr(1:l)//' ' 
          l = l + 1
          GOTO 18
         end if
 19     continue
 18     do 11 k=j,80
         if (ion_name(i)(k:k).ne.' ') then
          pltstr = pltstr(1:l)//ion_name(i)(k:k)
          l = l + 1
         else
          pltstr = pltstr(1:l)//' '
          l = l + 1
          GOTO 16
         end if
 11     continue
 16     ionlabel(i) = pltstr(1:l)
        lngth(i)  = l
 17    continue

c     return

        return
       
c     error trapping

 999  write(SCREEN,*) 'NOT FOUND-  : ',ionlist(1:30)
      write(STDOUT,*) 'NOT FOUND-  : ',ionlist(1:30)
      STOP ' *** MINFIT(getlist): terminated ***'

      end

c     eof
