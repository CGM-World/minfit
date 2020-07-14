c..............................................................................
c

      SUBROUTINE taglines(m,n,a,siga,iw,wa,ilwa,nunique,badline)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      include          'minfit.h'
      include          'const.dek'

      integer           maxbl
      parameter         (maxbl=maxlines*maxions+maxlines)
      logical           check1,check2,check3,check4,check5
      integer           i,j,k,kk,m,n,speciesi,linei,fidx,in,ib,iz
      integer           blknt,lastbadline,indx(maxbl),sibad(maxbl)
      integer           in1,in2,ib1,ib2,iz1,iz2,nunique,priority
      integer           bad(maxbl),badline(maxbl),specidx(maxbl)
      integer           rank(maxbl),pal(maxbl),drop(maxbl)
      integer           iw(maxcoeffs),ilwa

      double precision  N13max,fracerr,bpar,dz,ddz,febadmax
      double precision  wa(lwa)
      double precision  fNerr(maxbl),fTerr(maxbl),
     @                  fzerr(maxbl),ferr(maxbl),febad(maxbl)
      double precision  a(maxcoeffs),siga(maxcoeffs)

      character*80      mess


c     initialize BLKNT, NUNIQUE

      blknt   = 0
      nunique = 0

c     loop over the lines for each species; if the line has a total
c     fractional error larger than input parameter BADNESS, then
c     increment the bad line counter (BLKNT) and store the line
c     information

c     If a flagion is not defined, i.e. =0, we will check all species.
c     If flagion is defined then check only that tranistions species
c     lines.

c     we use working array WA here, which we stuff with the line number;
c     this will be so that we can sort the lines in order of increasing
c     line number following this loop

      DO 07 speciesi=1,species

       fidx = fitindex(speciesi)

       DO 05 linei=1,lines

        IF ((flagion(regioni).eq.0).OR.
     &      (fidx.eq.flagion(regioni))) THEN

          in      = idx_an(speciesi,linei)
          ib      = idx_ab(linei)
          iz      = idx_az(linei)

c     compute the total fractional error, FRACERR, and the b parameter,
c     BPAR, in km/s

          N13max  = sign(1.0d0,a(in))*(log10(abs(a(in)))+13.0d0)
          fracerr = sqrt((siga(in)/a(in))**2+(siga(ib)/a(ib))**2
     &                +(siga(iz)/a(iz))**2)
          bpar    = sqrt(2.0d0*kerg*a(ib)/spec_mass(fidx))
          bpar    = 1.0d-5*bpar   ! convert to km/s

C          WRITE(6,*) spec_name(fidx)(1:5),N13max,a(in),bpar 

c     initialize the logical check flags and check parameter constraints:
c     - CHECK1 flags if total frational error exceeds user defined BADNESS
c     - CHECK2 flags if the b parameter fractional error exceeds the
c       maximum allowed specified as FBERR by the user
c     - CHECK3 flags if the b parameter is smaller than the minimum
c       allowed, BMIN, sepcidied by the user

          check1  = .false.
          check2  = .false.
          check3  = .false.
          check4  = .false.
          check5  = .false.

          IF ((fidx.eq.1).OR.((fidx.gt.1).AND.(a(in).gt.0.0d0))) then
           IF (fracerr.gt.badness)      check1 = .true.
           IF (siga(ib)/a(ib).gt.fberr) check2 = .true.
           IF (bpar.lt.bmin)            check3 = .true.
           IF (bpar.gt.bmax)            check4 = .true.
           IF (N13max.gt.Nmax)          check5 = .true.         
          END IF

c     update the book keeping arrays for "bad" lines

          IF (check1.OR.check2.OR.check3.OR.check4.OR.check5) THEN

            blknt          = blknt + 1
            ferr(blknt)    = fracerr
            fNerr(blknt)   = siga(in)/a(in)
            fTerr(blknt)   = siga(ib)/a(ib)
            fzerr(blknt)   = siga(iz)/a(iz)
            badline(blknt) = linei
            specidx(blknt) = fidx
            wa(blknt)      = float(linei)

c     there are 7 possible permutations to why a line is considered
c     "bad"; set the BAD flag for communication purposes

         IF (check1.AND.(.not.check2).AND.(.not.check3)) bad(blknt) = 1 
         IF ((.not.check1).AND.check2.AND.(.not.check3)) bad(blknt) = 2 
         IF ((.not.check1).AND.(.not.check2).AND.check3) bad(blknt) = 3
         IF (check1.AND.check2.AND.(.not.check3))        bad(blknt) = 4
         IF ((.not.check1).AND.check2.AND.check3)        bad(blknt) = 5
         IF (check1.AND.(.not.check2).AND.check3)        bad(blknt) = 6
         IF (check1.AND.check2.AND.check3)               bad(blknt) = 7
         IF (check4)                                     bad(blknt) = 9
         IF (check5)                                     bad(blknt) = 10

          END IF

        END IF

 05    CONTINUE  ! next linei
   
 07   CONTINUE   ! next speciesi

c     another criterion is that two lines are at very similar redshifts
c     as determined by their uncertainties; this can happen in saturated
c     regions easily where the redshift is not well constrained, so
c     don't check saturated line; if we do the check, flag both lines as
c     being bad

      DO 09 linei=1,lines-1 

      IF (.not.satflag(speciesi,linei)) then ! do only if not saturated

       in1 = idx_an(1,linei)    ! use primary species
       in2 = idx_an(1,linei+1)  ! use primary species
       ib1 = idx_ab(linei)
       ib2 = idx_ab(linei+1)
       iz1 = idx_az(linei)
       iz2 = idx_az(linei+1)

       dz  = abs(a(iz1)-a(iz2))
       ddz = sqrt(siga(iz1)**2 + siga(iz2)**2)

        IF (dz.lt.ddz) THEN  ! 1 sigma  (3 sigma too often true)

         fracerr = sqrt((siga(in1)/a(in1))**2+(siga(ib1)/a(ib1))**2
     &                +(siga(iz1)/a(iz1))**2)
         blknt          = blknt + 1
         ferr(blknt)    = fracerr
         fNerr(blknt)   = siga(in1)/a(in1)
         fTerr(blknt)   = siga(ib1)/a(ib1)
         fzerr(blknt)   = siga(iz1)/a(iz1)
         badline(blknt) = linei
         specidx(blknt) = 1             ! primary species for highest priority
         wa(blknt)      = float(linei)
         bad(blknt)     = 8
         pal(blknt)     = badline(blknt) + 1

         fracerr = sqrt((siga(in2)/a(in2))**2+(siga(ib2)/a(ib2))**2
     &                +(siga(iz2)/a(iz2))**2)
         blknt          = blknt + 1
         ferr(blknt)    = fracerr
         fNerr(blknt)   = siga(in2)/a(in2)
         fTerr(blknt)   = siga(ib2)/a(ib2)
         fzerr(blknt)   = siga(iz2)/a(iz2)
         badline(blknt) = linei+1
         specidx(blknt) = 1             ! primary species for highest priority
         wa(blknt)      = float(linei+1)
         bad(blknt)     = 8
         pal(blknt)     = badline(blknt) - 1

        END IF
       END IF

 09    CONTINUE  ! next linei


c     we now have identified and stored all the lines that have suspect
c     fractional errors. 

c     if all lines are "cool", we are done, communicate all lines are
c     good to go; return to routine FITREGION

      IF (blknt.eq.0) then
        WRITE(SCREEN,402)
        WRITE(STDOUT,402) 
        nunique = 0
        RETURN
      END IF

 402  FORMAT(1x,'--- ALL LINES PARAMETERS WELL CONSTRAINED ---',/)

c     if we found some bad lines

c     communicate: announce that we have found bad lines, provide header
c     for line list

      WRITE(SCREEN,403)
      WRITE(STDOUT,403) 
 403  FORMAT(1x,'*** FOUND POORLY CONSTRAINED LINES ***',/)

      WRITE(SCREEN,404)
      WRITE(STDOUT,404)
 404  FORMAT(1x,'SPECIES',1x,'LINE',3x,'ferr',8x,'fNerr',7x,'fberr',
     &       7x,'fzerr',/,
     &       1x,'-------',1x,'----',3x,'----',8x,'-----',7x,'-----',
     &       7x,'-----')

c     sort the index of the WA (contains the line numbers) so that we
c     can communicate and test the lines in line number order; we used
c     WA because routine INDEXX requires a double precision real,
c     whereas the array BADLINE is integer

      IF (blknt.gt.1) THEN
        CALL indexx(blknt,wa,indx)
      ELSE
        indx(1) = 1
      END IF

c     more communication; announce which lines are bad lines; first
c     announce the redshift overlaps; then announce the lines in order
c     that they were found

      DO 101 j=1,blknt
       IF (bad(indx(j)).eq.8) then ! 8 is the redshift overlap index
        mess= '-> redshift overlaps with line ='
        WRITE(SCREEN,405) spec_name(specidx(indx(j))),badline(indx(j)),
     &                    ferr(indx(j)),fNerr(indx(j)),fTerr(indx(j)),
     &                    fzerr(indx(j)),mess,pal(indx(j))
        WRITE(STDOUT,405) spec_name(specidx(indx(j))),badline(indx(j)),
     &                    ferr(indx(j)),fNerr(indx(j)),fTerr(indx(j)),
     &                    fzerr(indx(j)),mess,pal(indx(j))
       END IF
 101  CONTINUE

 405  FORMAT(1x,a6,i4,1p4e12.2,3x,a33,i2)


      DO 102 j=1,blknt
        IF (bad(indx(j)).ne.8) then
        IF (bad(indx(j)).eq.1)  mess= '-> ferr exceeded'
        IF (bad(indx(j)).eq.2)  mess= '-> fberr exceeded'
        IF (bad(indx(j)).eq.3)  mess= '-> bmin exceeded'
        IF (bad(indx(j)).eq.4)  mess= '-> ferr and fberr exceeded'
        IF (bad(indx(j)).eq.5)  mess= '-> fberr and bmin exceeded'
        IF (bad(indx(j)).eq.6)  mess= '-> ferr and bmin exceeded'
        IF (bad(indx(j)).eq.7)  mess= '-> ferr, fberr, & bmin exceeded'
        IF (bad(indx(j)).eq.9)  mess= '-> bmax exceeded'
        IF (bad(indx(j)).eq.10) mess= '-> Nmax exceeded'
        WRITE(SCREEN,406) spec_name(specidx(indx(j))),badline(indx(j)),
     &                    ferr(indx(j)),fNerr(indx(j)),fTerr(indx(j)),
     &                    fzerr(indx(j)),mess
        WRITE(STDOUT,406) spec_name(specidx(indx(j))),badline(indx(j)),
     &                    ferr(indx(j)),fNerr(indx(j)),fTerr(indx(j)),
     &                    fzerr(indx(j)),mess
       END IF
 102  CONTINUE

 406  FORMAT(1x,a6,i4,1p4e12.2,3x,a40)

      WRITE(SCREEN,407)
      WRITE(STDOUT,407)
 407  FORMAT(1x,'-------',1x,'----',3x,'----',8x,'-----',7x,'-----',
     &       7x,'-----')


c
c     SANITY CHECK (cannot F-Test a single line model)
c

      IF ((lines.eq.1).AND.(blknt.ne.0)) then
        WRITE(SCREEN,*) ' >>>>>>>>>>>>>>>> WARNING <<<<<<<<<<<<<<<<'
        WRITE(SCREEN,*) ' Single Line model is not well constrained'
        WRITE(SCREEN,*) ' Keeping this model at this stage- but you'
        WRITE(SCREEN,*) ' might want to reconsider this region...'
        WRITE(SCREEN,*) ' >>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<'
        WRITE(STDOUT,*) ' >>>>>>>>>>>>>>>> WARNING <<<<<<<<<<<<<<<<'
        WRITE(STDOUT,*) ' Single Line model is not well constrained'
        WRITE(STDOUT,*) ' Keeping this model at this stage- but you'
        WRITE(STDOUT,*) ' might want to reconsider this region...'
        WRITE(STDOUT,*) ' >>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<'
        RETURN
      END IF


c
c     DETERMINE THE UNIQUE LINE NUMBERS
c

c     we want to check each line only once, even if it is "bad" in
c     multiple species and transitions; here we identify the unique line
c     numbers; we use the working array WA to store the unique line
c     numbers, which we will transfer back into the BADLINE array once
c     they have been sorted in priority order

      j = 0
      lastbadline = 0
      DO 61 i=1,blknt
        linei = badline(indx(i))
        IF (linei.ne.lastbadline) then
         j           = j + 1          
         wa(j)       = linei
         lastbadline = linei
        END IF
 61   CONTINUE

c     store the number of unique lines

      nunique = j

c
c     PRIORITIZE THE BAD LINES
c

c     WA now contains the unique bad lines by number (whereas BADLINE
c     contains the line number for all components inclusive- it has
c     redundancies); we now want to prioritize the bad lines according to:
c
c     speciesi absorption strength order with decreasing FERR order; the
c     species order is assumed set in the ions.table file by the used
c     (column 2 - the tie values)
c
c     for prioritizing we need to know store the minimum value if the
c     species ID for each unique line (the smaller value will translate
c     to a higher priority, and we store the maximum fractional error (a
c     larger fractional error will translate to a higher priority)

      DO 767 i=1,nunique
        rank(i)  = 0
        drop(i)  = 0
        k        = int(wa(i)) ! the unique line number stored in WA
        febad(i) = -1.0d99
        sibad(i) = 100
        DO 768 kk=1,blknt  ! loop and find min SPECIDX for unique line number
         linei = badline(kk)
         IF (linei.eq.k) sibad(i) = min(specidx(kk),sibad(i))
 768    CONTINUE
        DO 766 kk=1,blknt  ! loop to find max FERR for this SIBAD
         linei = badline(kk)
         IF (linei.eq.k) then ! and store 
          IF (specidx(kk).eq.sibad(i)) febad(i) = max(ferr(kk),febad(i))
         END IF
 766    CONTINUE
C        WRITE(6,*) '----> ',i,int(wa(i)),sibad(i),febad(i) ! DEBUG 
 767  CONTINUE

c     now use the RANK array to stuff the indices in priority order

      priority = 0
      DO 821 i=1,species
       fidx = fitindex(i) 
 1002  febadmax = -1.0d99
       k        = 0
       DO 831 j=1,nunique
        IF (j.ne.drop(j)) then
         IF ((sibad(j).eq.fidx).AND.(febad(j).ge.febadmax)) then
          k = j
          febadmax = febad(j)
         END IF
         END IF
 831   CONTINUE
       IF (k.ne.0) then
        priority = priority + 1
        drop(k)  = k
        rank(k)  = priority
C        WRITE(6,*) '--> ',i,drop(k),rank(k) ! DEBUG 
        GOTO 1002  ! loop through the unique line list again for this species
       END IF
 821  CONTINUE

C DEBUG (loop)
C      DO 51 j=1,nunique
C       WRITE(6,*) j,int(wa(j)),sibad(j),febad(j),rank(j)
C 51   CONTINUE

c     from this point on we disjoint the BADLINE array from the other
c     arrays describing the fractional errors (we no longer need any of
c     the other data on the bad lines), BADLINE now contains the unique
c     line numbers in priority order

      DO 62 i=1,nunique
       badline(rank(i)) = int(wa(i))
C       WRITE(6,*) i,rank(i),int(wa(i)) ! DEBUG
 62   CONTINUE

c     we can now return, the important info we send back is NUNIQUE and
c     BADLINE in prioirty order


      RETURN
      END

c
c..............................................................................
c

      SUBROUTINE doFtest(n,m,a,var1,nu1,var2,nu2,ftest,fprob,replace)
c  
c     we use the reduced-chisq, which is a variance indicator to
c     ascertain the significance. the F-test checks that the variances
c     of two distributions are statistically similar (null hypothesis);
c     so we compute the F statistic on the reduced-chisq. the F
c     statistic probability, P(F) (which provides the probability that
c     the value of F would be as large as it is for the degrees of
c     freedom) needs to be very small to rule out that the distributions
c     are statistically similar; as such P(F) is a 1-CL, where CL is the
c     confidence level that the null hypothesis is ruled out.  The user
c     has supplied the confidence level in the MINFIT.PAR file as
c     parameter CONFIDENCE the probability P(F) is computed using the
c     incomplete beta function; we assume a double sided distribution
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'

      logical           replace
      integer           n,m
      double precision  a(maxcoeffs),var0,var1,var2,nu0,nu1,nu2
      double precision  xx,aa,bb,fprob,ftest,betai


c     set replace low

      replace = .false.

c     we want to save VAR1 and NU1, but we need to do a two-sided F-test
c     so we will be juking the VAR1, VAR2, NU2, and NU2 values, so save
c     VAR1, and NU1 in VAR0 and NU0, which we use for the test

      var0 = var1
      nu0  = nu1

c     we already called fitstats before examining the test model, but
c     call it again for good measure and save

      CALL fitstats(a,n,m,-1)
      nu2   = real(m-n)
      var2  = chisq/nu2

c     the F-test requires VAR0>VAR2

      IF (var0.gt.var2) THEN
        ftest = var0/var2
      ELSE ! requires switching VAR and NU
        ftest = var2/var0
        nu0 = nu2
        nu2 = nu1
      END IF

      xx    = nu2/(nu2+nu0*ftest)
      aa    = 0.5d0*nu2
      bb    = 0.5d0*nu0
      fprob = 2.0d0 * betai(aa,bb,xx)
      IF (fprob.gt.1.0d0) fprob = 2.0d0 - fprob

c     check if the line is significant; set logical REPLACE high if we
c     are replacing the original fit with the F-test fit; the default is
c     .false. (set above)

      IF ((1.0d0-fprob).lt.confidence) replace = .true.

      RETURN
      END
c
c
c..............................................................................
c     the following are Numerical Recipes used for the statistical
c     testing; they are unmodified
c..............................................................................
c
c     "RIPPED OFF" NUMERICAL RECIPE 
c

      DOUBLE PRECISION FUNCTION BETAI(A,B,X)

c
c     this is a Numerical Recipe that computes the beta function, which
c     is the cummulative distribution function for the students-t
c     distribution, it provides the probability that the value of F in
c     the F test has a value as great as what is measure as compared to
c     unity
c
c     the beta function is a special case of the gamma function
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit double precision (a-h,o-z)

c
      IF(X.LT.0..OR.X.GT.1.)
     &  STOP ' DEATH: bad X in routine betai (adjustfit.f)'
      IF(X.EQ.0..OR.X.EQ.1.)THEN
        BT=0.
      ELSE
        BT=EXP(GAMMLN(A+B)-GAMMLN(A)-GAMMLN(B)
     *      +A*LOG(X)+B*LOG(1.-X))
      ENDIF
      IF(X.LT.(A+1.)/(A+B+2.))THEN
        BETAI=BT*BETACF(A,B,X)/A
        RETURN
      ELSE
        BETAI=1.-BT*BETACF(B,A,1.-X)/B
        RETURN
      ENDIF
c
      END
c
c..............................................................................
c
c     "RIPPED OFF" NUMERICAL RECIPE 
c

      DOUBLE PRECISION FUNCTION BETACF(A,B,X)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      implicit double precision (a-h,o-z)      

      PARAMETER (ITMAX=100,EPS=3.E-16)
      AM=1.
      BM=1.
      AZ=1.
      QAB=A+B
      QAP=A+1.
      QAM=A-1.
      BZ=1.-QAB*X/QAP
      DO 11 M=1,ITMAX
        EM=M
        TEM=EM+EM
        D=EM*(B-M)*X/((QAM+TEM)*(A+TEM))
        AP=AZ+D*AM
        BP=BZ+D*BM
        D=-(A+EM)*(QAB+EM)*X/((A+TEM)*(QAP+TEM))
        APP=AP+D*AZ
        BPP=BP+D*BZ
        AOLD=AZ
        AM=AP/BPP
        BM=BP/BPP
        AZ=APP/BPP
        BZ=1.
        IF(ABS(AZ-AOLD).LT.EPS*ABS(AZ))GO TO 1
11    CONTINUE
      STOP ' DEATH: A/B 2big, or ITMAX 2small in betafc (adjustfit.f)'
1     BETACF=AZ

      RETURN
      END
c
c..............................................................................
c
c     "RIPPED OFF" NUMERICAL RECIPE 
c
 
      DOUBLE PRECISION FUNCTION GAMMLN(XX)

c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit double precision (a-h,o-z)
      double precision          COF(6)

      SAVE COF,STP
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/

      X=XX-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE
      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE
      GAMMLN=TMP+LOG(STP*SER)

      RETURN
      END
c


c..............................................................................
c
c     "RIPPED OFF" NUMERICAL RECIPE 
c

      SUBROUTINE        indexx(n,arrin,indx)   

      implicit none
c..
c..indexes an array arrin of length n. outputs the array indx such that
c..arrin(index(j)) (j=1,...n) is in ascending order. the input array
c..arrin is not changed.
c..
c..declare

      integer          n,indx(n),indxt,j,l,ir,i
      double precision arrin(n),q
c..
c..iniitalize the index array with consecutive integers
      do 11 j=1,n   
        indx(j)=j   
11    continue  
c..
c..from here on out its heapsort, but with indirect indexing through
c..indx in all references to arrin. compare to the heapsort routine above.
      l=n/2+1   
      ir=n  
10    continue  
        if(l.gt.1)then  
          l=l-1 
          indxt=indx(l) 
          q=arrin(indxt)
        else
          indxt=indx(ir)
          q=arrin(indxt)
          indx(ir)=indx(1)  
          ir=ir-1   
          if(ir.eq.1)then   
            indx(1)=indxt   
            return  
          endif 
        endif   
        i=l 
        j=l+l   
20      if(j.le.ir)then 
          if(j.lt.ir)then   
            if(arrin(indx(j)).lt.arrin(indx(j+1)))j=j+1 
          endif 
          if(q.lt.arrin(indx(j)))then   
            indx(i)=indx(j) 
            i=j 
            j=j+j   
          else  
            j=ir+1  
          endif 
        go to 20
        endif   
        indx(i)=indxt   
      go to 10  
      end   

c eof
