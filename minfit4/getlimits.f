c     Module: getlimits.f
c
c     ROUTINE getNlimits(error)
c     ROUTINE fixN(speciesi,fidx,linei,col,zcl,limflag)
c
c.......................................................................
c  

      SUBROUTINE         getNlimits(error)

c  
c     UPSHOT: here,we check each component to determine if the VP
c     parameters are consistent with the noise (or the limiting column
c     density).  we then apply tests to determine whether to adopt a
c     column density limit or the VP model fit.  Care must be taken to
c     account for masked pixels.
c
c     this routine is called after we have the final VP parameters for
c     the current region.  in some cases components for some transitions
c     are either non-existant (no feature at the same velocity where
c     there are components in other transitions).  in other cases in
c     which a component was found, its LSF value may be consistent with
c     the column density detection limit
c
c     with regard to masked pixels.  if a VP component in a higher order
c     transition is masked and it is the only transition representing
c     its ion species, then we replace with the limit for the species;
c     we use the criterion that the VP component cannot be assessed if
c     the component centroid is within +/-1 pixel of the masked pixels
c     (the mask is at the line center).  also, in cases where multiple
c     transitions are available, we exclude calculations in regions that
c     have masked pixels; we are more liberal here, if any of the pixels
c     are masked in the region and there are other transitions
c     representing the ion species to use for constraints, we omit the
c     transition that has masked pixels in the region
c
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  

      include            'minfit.h'
      include            'const.dek'
      double precision   Csig
      parameter          (Csig = 3.0d0)
      integer            i,j,iknt
      integer            pixi,linei,speciesi,ioni
      logical            error,flag0,flag1,flag2,flag3,flag4
      logical            action,flag,mflag(maxions)
      double precision   collim(maxions),col
      double precision   getNlim,getNaod,bpar,dbpar,fberror
      double precision   Nlim,Nlim1,Nlim2
      double precision   Naod,Naodsig,aodsig,Nlimave,wave,knt
      double precision   zcl,Nvp,dNvp,Nvpmin,Nvpmax,vpsng
      double precision   Naodmin,Naodmax,aodsng
      character*1        sflag0,sflag1,sflag2,sflag3,sflag4


      DO 01 j=1,2
        if (j.eq.1) i = SCREEN
        if (j.eq.2) i = STDOUT
        WRITE(i,*) '+++++++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(i,*) '+                                               +'
        WRITE(i,*) '+        ASSESSING VP COLUMN DENSITIES          +'
        WRITE(i,*) '+                                               +'
        WRITE(i,*) '+ using VP errors, limits (Nlim), & AOD (Naod)  +'
        WRITE(i,*) '+                                               +'
        WRITE(i,*) '+ Applying five boolean criteria - (12345)      +'
        WRITE(i,*) '+ 1 (T/F) Nvp returned by LSF?                  +'
        WRITE(i,*) '+ 2 (T/F) (Nvp-dNvp)>Nlim?                      +'
        WRITE(i,*) '+ 3 (T/F) Nvp>Nlim & (Naod-dNaod)>Nlim?         +'
        WRITE(i,*) '+ 4 (T/F) Nvp>Nlim & Naod saturated?            +'
        WRITE(i,*) '+ 5 (T/F) (Nvp-dNvp)<0 & (Naod-dNaod)>Nlim?     +'
        WRITE(i,*) '+ -> (FXXXX) = Nlim adopted                     +'
        WRITE(i,*) '+ -> (TFFFX) = Nlim adopted                     +'
        WRITE(i,*) '+ -> (TFXXT) = Nlim adopted (but check it)      +' 
        WRITE(i,*) '+ -> all other combinations = Nvp adopt         +'
        WRITE(i,*) '+                                               +'
        WRITE(i,*) '+ NOTES:                                        +'
        WRITE(i,*) '+ 1) lines have been sorted and renumbered in   +'
        WRITE(i,*) '+    order of increasing redshift               +'
        WRITE(i,*) '+ 2) condition 5 applies only in unsaturated    +'
        WRITE(i,*) '+    regions (as specified by the user)         +'
        WRITE(i,*) '+ 3) If the 5th condition is T, then one        +'
        WRITE(i,*) '+    should carefully check this component      +'
        WRITE(i,*) '+++++++++++++++++++++++++++++++++++++++++++++++++'
        WRITE(i,*) ' '
 01   CONTINUE

c     set error low; it is passes around but never utilized in this routine

      error =.false.

c     MASTER LOOPS (outer: species; inner: lines)
c     SUB LOOP     (ions)

c     loop over all the species (total) one by one and examine the
c     fitting parameters component by component. 

      DO 09 speciesi=1,total

c     loop over the lines.  for each line examine if the VP column
c     density for speciesi should be replace by its column density
c     limits

c     FYI: column densities have been unpacked and are in log10 now
c     give the limit to 3 sigma 

        DO 07 linei=1,lines

         flag0  = .true.   ! set false of Nvp unconstrained by fit
         flag1  = .false.  ! set true if Nvpmin>Nlim
         flag2  = .false.  ! set true if Nvp>Nlim AND Naodmin>Nlim
         flag3  = .false.  ! set true if Nvp>Nlim AND Naodsig=-1.0
         flag4  = .false.  ! set true if Nvplim<0 AND not saturated
         action = .false.  ! set true if we replace line with limit

         sflag0 = 'T'
         sflag1 = 'F'
         sflag2 = 'F'
         sflag3 = 'F'
         sflag4 = 'F'

C     STEP 1: obtain estimation of limits 

c     compute limits at the Csig level; we must adopt the unweighted
c     average constraining limit from all ions for the species based
c     upon signal to noise (sigma spectrum).  we count the number of
c     transitions associated with the current species (SPECIESI) and
c     then use FLAG to book keep whether the result for a given
c     transition is based upon the uncertainty spectrum, or in the case
c     of having masked pixels, based upon the residuals in the data
c     (conservative estimate)

         knt     = 0.
         iknt    = 0
         flag = .false.
         DO 03 ioni=1,ions
          IF (spec_idx(ioni).eq.speciesi) THEN
            iknt = iknt + 1
            bpar = bline(speciesi,linei)  
            zcl  = zline(linei)
            collim(iknt)   = getNlim(ioni,linei,bpar,zcl,flag)
            mflag(iknt) = flag
C            WRITE(6,*) ion_name(ioni),zcl,iknt,collim(iknt),mflag(iknt)
           END IF
 03      CONTINUE

c     now that we have the column density limit estimate for all the
c     transitions associated with SPECIESI, we must select the best
c     estimate of the limit; if MFLAG(i)-.TRUE. then masked pixels were
c     encountered, if .FALSE. then no masked pixels were encountered.
c     the logic is as follows, we adopt the most stringent limit from
c     those transitions with no masked pixels.  If this is null, then we
c     adopt the most stringent limit from those with masked pixels.

         Nlim  = 0.0d0
         Nlim1 = 1.0d99
         Nlim2 = 1.0d99
         iknt  = 0
         DO 04 ioni=1,ions
          IF (spec_idx(ioni).eq.speciesi) THEN
            iknt = iknt + 1
            IF (.not.mflag(iknt)) then ! non masked pixels 
              Nlim1 = min(Nlim1,collim(iknt))
            END IF            
            IF (mflag(iknt)) then ! masked pixels 
              Nlim2 = min(Nlim2,collim(iknt))
            END IF            
          END IF
 04      CONTINUE

         IF (Nlim1.ne.1.0d99) Nlim = Nlim1
         IF (Nlim.eq.0.0d0)   Nlim = Nlim2
    
C        WRITE(6,'(1x,a5,1pe10.2)') 'Adopted:  ',Nlim

C     STEP 2: obtain estimation of AOD columns

c     compute the AOD and 1-sig uncertainty; we chose the largest value
c     of all ions for the species over the region of the bpar for the
c     fit - this is to check that VP cols with large errors are found in
c     regions of significant column density, if the ion data are masked
c     omit from the calculation

         Naod = -1.0d99
         DO 05 ioni=1,ions
          IF (spec_idx(ioni).eq.speciesi) THEN
           bpar = bline(speciesi,linei)  
           zcl  = zline(linei)
           col  = getNaod(speciesi,ioni,linei,bpar,zcl,aodsig)
C           WRITE(6,'(1x,a5,a13,2(1pe10.2))') 'AOD: ',ion_name(ioni),
C     &               col,aodsig
           IF (max(Naod,col).eq.col) then
             Naod    = max(Naod,col)
             Naodsig = aodsig              
           END IF
          END IF
 05      CONTINUE
C         WRITE(6,'(1x,a5,1pe10.2)') 'Adopted:  ',Naod

c     now store important values of the VP component (from the LSF) and
c     of the AOD column density, grab the sign of the minimum columns,
c     which may go negative if the fractional uncertainty is > 1

         Nvp     = 10.0d0**nline(speciesi,linei) 
         dNvp    = Nvp*dnline(speciesi,linei)/0.4343
         bpar    = bline(speciesi,linei)  
         dbpar   = dbline(speciesi,linei) 
         fberror = abs(dbpar/bpar)

         vpsng   = 1.0d0
         Nvpmin  = 0.0d0         
         Nvpmax  = 0.0d0

         IF  (dNvp.gt.0.0d0) THEN
           Nvpmin = Nvp - dNvp
           Nvpmax = Nvp + dNvp
           IF (Nvpmin.lt.0.0d0) vpsng = -1.0d0
         END IF

         aodsng  = 1.0d0
         Naodmin = 0.0d0         
         Naodmax = 0.0d0

         IF (Naodsig.ne.-1.0d0) then
           Naodmin = Naod - 3.0d0*Naodsig         
           Naodmax = Naod + 3.0d0*Naodsig
           IF (Naodmin.lt.0.0d0) aodsng = -1.0d0
         END IF

         IF (Naod.eq.1.0d0) then
           Naodmin = 1.0
           Naodmax = 1.0
         END IF

c     perform the checks

c     0th condition; this is a hardfast REPLACE condition with Nlim. if
c     DNLINE=-1 then the ion had no detected feature; if there was an
c     Nvp and all the transitions are consistent with log(Naod)=0 (which
c     can happen when masking is involved), then also REPLACE with Nlim

         IF ((dnline(speciesi,linei).eq.-1.0).OR.(Naod.eq.1.0d0)) then
           flag0  = .false.
           sflag0 = 'F'
         END IF

c     the following are flags for KEEP CONDITIONS, if any one of them is
c     set high then we keep the VP column (equivalently, if all of them
c     fail, then we apply the limit); the tests are in order of
c     robustness

c     1st condition; Nvpmin > Nlim; if this happens then we believe that
c     the VP column is significant.  this is a pretty robust condition,
c     there is no real getting around this condition for it implies that
c     our VP column is both well constrained and above the limit

         IF (Nvpmin.gt.Nlim/Csig) THEN
           flag1  = .true.
           sflag1 = 'T'
         END IF

c     2nd condition; Nvp>Nlim AND Naodmin>Nlim; this checks that the Nvp
c     is in fact larger than the limit, then if it is we check if the
c     AOD column is conistent with being greater than the limit; this
c     method is least robust and can fail if the b parameter is small

         IF ((Nvp.gt.Nlim/Csig).AND.(Naodmin.gt.Nlim/Csig)) THEN
           flag2  = .true.
           sflag2 = 'T'
         END IF

c     3rd condition; Nvp>Nlim AND Naodsig=-1.0; this checks that the Nvp
c     is in fact larger than the limit, then if it is we check if the
c     AOD column is a lower limit, which suggest that the line is
c     saturated.  we could have a large error in Nvp though

         IF ((Nvp.gt.Nlim/Csig).AND.(Naodsig.eq.-1.0d0)) THEN
           flag3 = .true.
           sflag3 = 'T'
         END IF

C     THIS FLAG FOR THIS CONDITION IS TURNED OFF, BUT THE CONDITION IS COMMUNICATED
c     4th condition: if the 1st condition is .FALSE. and the 2nd
c     condition is .TRUE. then check if the error in Nvp is so large
c     that the VP column density is consistent with Nvp=0; we must not
c     be in a specified saturation region (provided by the user).  If
c     all the above is .TRUE. then we will replace Nvp by the AOD column

         IF ((.not.satflag(speciesi,linei)).AND.(Nvpmin.lt.0.0d0).AND.
     &       (Naodmin.gt.Nlim)) then
C           flag4  = .true.
           sflag4 = 'T'
         END IF

c     ---------------------------------------------
c     all checks complete.  results are as follows:
c     ---------------------------------------------

c     0th test (trump condition)
c     flag0=.true.  POSSIBLE KEEP    Nvp and uncertainty was delevired by LSF
c     flag0=.false. ENFORCE REPLACE  Nvp->Nlim no fit constraint 

c     1st test
c     flag1=.true.  ENFORCE KEEP     Nvpmin>Nlim
c     flag1=.false. POSSIBLE REPLACE Nvp->Nlim if flag2=flag3=.FALSE.
c                                    Nvp->Naod if flag2=.TRUE.

c     2nd test
c     flag2=.true.  ENFORCE KEEP     Nvp>Nlim AND Naodsig=-1.0
c     flag2=.false. POSSIBLE REPLACE Nvp->>Nlim if flag1=flag2=flag3=.FALSE.

c     3rd test
c     flag3=.true.  ENFORCE KEEP     Nvp>Nlim AND Naodmin>Nlim
c     flag3=.false. POSSIBLE REPLACE Nvp->>Nlim if flag1=flag2=flag3=.FALSE.

C     TURNED OFF
c     4th test
c     flag4=.true.  ENFORCE KEEP     Nvpmin>0 (redundant with 1st test)
c     flag4=.false. ENFORCE REPLACE  Nvp->Naod if flag1=.FALSE. & flag2=.TRUE.

c     -> if flag0=.FALSE. then REPLACE with Nlim. 
c     -> if any one of flag1,2,3=.TRUE. then KEEP
c     -> flag0=.TRUE. & flag1,2,3=.FLASE. REPLACE with Nlim
c     -> flag4=.TRUE. & flag1=.FALSE. & flag2=.TRUE., REPLACE with Naod (DEFUNCT)

         IF (.not.flag0) then  ! then flags 1,2,3,4 moot, set the low for
           flag1 = .false.     ! next conditional statement to be true
           flag2 = .false.     ! to enforce the REPLACE
           flag3 = .false.
           flag4 = .false.
           sflag1 = 'X'
           sflag2 = 'X'
           sflag3 = 'X'
           sflag4 = 'X'
         END IF

c     replace with Nlim? if on the condition that flag0=.FALSE., we have
c     set all the remaining flags to .FALSE. so that this first
c     condition will be TRUE

         IF ((.not.flag1).AND.(.not.flag2).AND.(.not.flag3)) then
           action = .true.  ! we are replacing Nvp with Nlim
           nline(speciesi,linei)  = log10(Nlim)
           dnline(speciesi,linei) = -1.0
           dbline(speciesi,linei) = -1.0 
         END IF

C     THIS CONDITION WILL NOT BE MET BECAUSE FLAG4 SETTING HAS BEEN TURNED OFF
c     replace with Naod? if on the condition that flag0=.FALSE., we have
c     set all the remaining flags to .FALSE. so that this second
c     condition will not be TRUE (we do not want to execute it because
c     the flag0=.FALSE. condition is trumps all other conditions; if the
c     below condition TRUE then we are replacing Nvp with Noad, we do
c     not modify the Doppler parameter

         IF ((.not.flag1).AND.(flag2).AND.(flag4)) then
           action = .true.  ! we are replacing Nvp with Nlim
           nline(speciesi,linei)  = log10(Naod)
           dnline(speciesi,linei) = 0.4343*0.50d0*(Naodmax-Naodmin)/Naod
         END IF

c     communicate 

         DO 1000 j=1,2

         if (j.eq.1) i = SCREEN
         if (j.eq.2) i = STDOUT

         IF (.not.flag0) then ! replace with Nlim (trump condition)
          IF (Naodsig.ne.-1.0d0) then ! AOD uncertainties
          WRITE(i,701) spec_name(speciesi),linei,log10(Nlim),
     &      aodsng*log10(abs(Naodmin)),
     &      sign(1.0d0,Naod)*log10(abs(Naod)),
     &      sign(1.0d0,Naodmax)* log10(abs(Naodmax)),
     &      sign(1.0d0,Nvp)*log10(abs(Nvp))
          ELSE                        ! AOD limit
          WRITE(i,702) spec_name(speciesi),linei,log10(Nlim),
     &      sign(1.0d0,Naod)*log10(abs(Naod)),
     &      sign(1.0d0,Nvp)*log10(abs(Nvp))
          END IF
         END IF

         IF (flag0.and.action) then ! replace with Nlim or Naod
          IF (flag4) then ! replace with Naod
          WRITE(i,707) spec_name(speciesi),linei,log10(Nlim),
     &      aodsng*log10(abs(Naodmin)),
     &      sign(1.0d0,Naod)*log10(abs(Naod)),
     &      sign(1.0d0,Naodmax)* log10(abs(Naodmax)),
     &      vpsng*log10(abs(Nvpmin)),
     &      sign(1.0d0,Nvp)*log10(abs(Nvp)),
     &      sign(1.0d0,Nvpmax)*log10(abs(Nvpmax)),
     &      sflag0,sflag1,sflag2,sflag3,sflag4
          END IF
          IF (.not.flag4) then 
          IF (Naodsig.ne.-1.0d0) then ! replace with Nlim, AOD uncertainties
          WRITE(i,703) spec_name(speciesi),linei,log10(Nlim),
     &      aodsng*log10(abs(Naodmin)),
     &      sign(1.0d0,Naod)*log10(abs(Naod)),
     &      sign(1.0d0,Naodmax)* log10(abs(Naodmax)),
     &      vpsng*log10(abs(Nvpmin)),
     &      sign(1.0d0,Nvp)*log10(abs(Nvp)),
     &      sign(1.0d0,Nvpmax)*log10(abs(Nvpmax)),
     &      sflag0,sflag1,sflag2,sflag3,sflag4
          ELSE                        ! replace with Nlim, AOD limits
          WRITE(i,704) spec_name(speciesi),linei,log10(Nlim),
     &      sign(1.0d0,Naod)*log10(abs(Naod)),
     &      vpsng*log10(abs(Nvpmin)),
     &      sign(1.0d0,Nvp)*log10(abs(Nvp)),
     &      sign(1.0d0,Nvpmax)*log10(abs(Nvpmax)),
     &      sflag0,sflag1,sflag2,sflag3,sflag4
          END IF
          END IF
         END IF

         IF (flag0.and.(.not.action)) then ! keep Nvp
          IF (Naodsig.ne.-1.0d0) then ! AOD uncertainties
          WRITE(i,705) spec_name(speciesi),linei,log10(Nlim),
     &      aodsng*log10(abs(Naodmin)),
     &      sign(1.0d0,Naod)*log10(abs(Naod)),
     &      sign(1.0d0,Naodmax)* log10(abs(Naodmax)),
     &      vpsng*log10(abs(Nvpmin)),
     &      sign(1.0d0,Nvp)*log10(abs(Nvp)),
     &      sign(1.0d0,Nvpmax)*log10(abs(Nvpmax)),
     &      sflag0,sflag1,sflag2,sflag3,sflag4
          ELSE                        ! AOD limits
          WRITE(i,706) spec_name(speciesi),linei,log10(Nlim),
     &      sign(1.0d0,Naod)*log10(abs(Naod)),
     &      vpsng*log10(abs(Nvpmin)),
     &      sign(1.0d0,Nvp)*log10(abs(Nvp)),
     &      sign(1.0d0,Nvpmax)*log10(abs(Nvpmax)),
     &      sflag0,sflag1,sflag2,sflag3,sflag4
          END IF
         END IF

 1000    CONTINUE

 07     CONTINUE  ! next line

        IF (lines.gt.1) then 
         WRITE(SCREEN,*) ' '
         WRITE(STDOUT,*) ' '
        END IF 

 09   CONTINUE ! next species

c     return
      WRITE(SCREEN,*) ' '
      WRITE(STDOUT,*) ' '

      RETURN

c     formats
                   
 701     FORMAT(2x,a6,3x,'line=',i2,2x,'Nlim=',f6.2,
     &     ' (Naod-,Naod,Naod+)=(',f6.2,',',f6.2,',',f6.2,')',
     &     ' (Nvp-,Nvp,Nvp+)=(',' .... ',',',f6.2,',',' .... ',')',
     &     ' Nlim adopted - Nvp unconstrained')
 702     FORMAT(2x,a6,3x,'line=',i2,2x,'Nlim=',f6.2,
     &     ' (Naod saturated  )=(',' .... ','>',f6.2,' ',' .... ',')',
     &     ' (Nvp-,Nvp,Nvp+)=(',' .... ',',',f6.2,',',' .... ',')',
     &     ' Nlim adopted - (Nvp unconstrained)')
 703     FORMAT(2x,a6,3x,'line=',i2,2x,'Nlim=',f6.2,
     &     ' (Naod-,Naod,Naod+)=(',f6.2,',',f6.2,',',f6.2,')',
     &     ' (Nvp-,Nvp,Nvp+)=(',f6.2,',',f6.2,',',f6.2,')',
     &     ' Nlim adopted - (',a1,a1,a1,a1,a1,')')
 704     FORMAT(2x,a6,3x,'line=',i2,2x,'Nlim=',f6.2,
     &     ' (Naod saturated  )=(',' .... ','>',f6.2,' ',' .... ',')',
     &     ' (Nvp-,Nvp,Nvp+)=(',f6.2,',',f6.2,',',f6.2,')',
     &     ' Nlim adopted - (',a1,a1,a1,a1,a1,')')
 705     FORMAT(2x,a6,3x,'line=',i2,2x,'Nlim=',f6.2,
     &     ' (Naod-,Naod,Naod+)=(',f6.2,',',f6.2,',',f6.2,')',
     &     ' (Nvp-,Nvp,Nvp+)=(',f6.2,',',f6.2,',',f6.2,')',
     &     ' Nvp  adopted - (',a1,a1,a1,a1,a1,')')
 706     FORMAT(2x,a6,3x,'line=',i2,2x,'Nlim=',f6.2,
     &     ' (Naod saturated  )=(',' .... ','>',f6.2,' ',' .... ',')',
     &     ' (Nvp-,Nvp,Nvp+)=(',f6.2,',',f6.2,',',f6.2,')',
     &     ' Nvp  adopted - (',a1,a1,a1,a1,a1,')')
 707     FORMAT(2x,a6,3x,'line=',i2,2x,'Nlim=',f6.2,
     &     ' (Naod-,Naod,Naod+)=(',f6.2,',',f6.2,',',f6.2,')',
     &     ' (Nvp-,Nvp,Nvp+)=(',f6.2,',',f6.2,',',f6.2,')',
     &     ' Naod adopted - (',a1,a1,a1,a1,a1,')')

      END


c  
c.......................................................................
c

      DOUBLE PRECISION FUNCTION getNlim(ioni,linei,bpar,zcl,flag)

c
c     this routine returns the column density limit spectrum for presnt
c     ion near the region of present component (LINEI)
c
c     the N limit is found by first obtaining the EW limit over
c     apertures defined by the b parameter (not just unresolved lines)
c
c  
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  

      include            'minfit.h'
      include            'const.dek'

      logical             flag
      integer             i,j,k,ioni,npix,upix,lpix,pixj,linei
      double precision    dw,dv,ewlim,wave,bpar,tau0,nres,bres
      double precision    doppb,zcl,z1,z2

      flag = .false.

c     define the ion number as j

      j    = ioni
      npix = ndata(j)

c     compute the optical depth at the line core for the linear COG,
c     requires physical constants input from const.dek; the factor of
c     1.E-8 converts one of the wavelengths to centimeters (not both,
c     because EW is in angstroms)

      tau0  = 1.0d-8*(pi*e*e*lambda0(j)*lambda0(j))/(me*c*c)
      doppb = bpar

c     if this line is under resolved, reset the bpar to the resolution
c     of the spectrograph so that we are getting the limit for an
c     nresolved line at a minimum; if the line is very broad such that
c     its bpar is greater than BMAX input by the user, then reset to
c     BMAX

      bres = ckms/(2.35*R_fac)  ! km/s
      IF (doppb.lt.bres) doppb = bres
      IF (doppb.gt.bmax) doppb = bmax

c     grab the pixel number corresponding to this line center

      wave = lambda0(j)*(1.0d0+zcl)
      pixj = 0
      DO 11 i=2,npix-1
       IF ((wave.gt.lambda(j,i)).AND.
     &     (wave.le.lambda(j,i+1))) THEN
         pixj = i
         GOTO 12
       END IF
 11   CONTINUE

c     TRAP: if we get here, then the pixel wasn't found, we do not want
c     to have DEATH because we are done fitting, communicate and move on

      WRITE(SCREEN,600) linei,ion_name(j),wave
      WRITE(STDOUT,600) linei,ion_name(j),wave

      getNlim = 1.0d0 ! account for log10 conversion to return log10(Nlim)=0.00

      RETURN

c     compute the aperture, avoiding edge effects

 12   dw = 0.0d0
      IF (pixj.eq.1)    dw = lambda(j,2)-lambda(j,1)
      IF (pixj.eq.npix) dw = lambda(j,npix)-lambda(j,npix-1)
      IF (dw.eq.0.0d0)  dw = lambda(j,pixj)-lambda(j,pixj-1)
      dv = (ckms/wave)*dw

c     set the initial pixel window for the calulation based upon the
c     DOPPB, but, we want to do this over 2*FWHM aperture, so include the
c     factor 2.35 (FWHM=2.35*DOPPB)

      lpix  = pixj - nint(2.35d0*doppb/dv)
      upix  = pixj + nint(2.35d0*doppb/dv)

      IF (lpix.lt.1)    lpix = 1
      IF (upix.gt.npix) upix = npix

C     DEBUGGING?
      IF (dodebug) then
       z1 = lambda(j,lpix)/lambda0(ioni)-1.0d0
       z2 = lambda(j,upix)/lambda0(ioni)-1.0d0
       WRITE(6,*) 'Nlim',linei,doppb,lpix,upix,zcl,z1,z2
      END IF

c     compute the equivalent width limit from the uncertainty spectrum,
c     if some of the data then we use the residuals in the spectrum, if
c     not, we use the standard sigma spectrum method; be sure to account
c     for edge effects if lpix=1 or upix=npix

      ewlim = 0.0d0

c     check for masked pixels in the line core.  if the line core is
c     constained then this is considered a potentially constrainable
c     fit, so we use the sigma spectrum.  here, we define the line core
c     as +/-1 pixel from the center pixel PIXJ, if any of these three
c     pixels are masked then we use the data residuals (based upon
c     setting flag=.TRUE.)

      DO 13 i=lpix,upix
       IF (abs(pixj-i).le.1) then 
        IF (quality(j,i).eq.0) flag = .true.
       END IF
 13   CONTINUE

c     if not pixels are masked, compute the limit from the uncertainty spectrum

       IF (.not.flag) then 

        DO 14 i=lpix,upix
         dw = 0.5 * abs(lambda(j,i-1)-lambda(j,i+1))
         IF (i.eq.1) dw = abs(lambda(j,i)-lambda(j,i+1))
         IF (i.eq.npix) dw = abs(lambda(j,i-1)-lambda(j,i))
         ewlim = ewlim + (dw*sigma(j,i))**2
         IF (dodebug) WRITE(6,*) 'Nlim',i,lambda(j,i),dw,
     &                           sqrt(ewlim),flag
 14    CONTINUE
       getNlim = sqrt(ewlim)/tau0
       RETURN

      END IF
   
c     if one or more pixels is masked then compute the limit from the
c     residuals in the data.  this appraoch provides a conservative
c     limit in the case of bad data

       IF (flag) then 

        DO 15 i=lpix,upix
         dw = 0.5 * abs(lambda(j,i-1)-lambda(j,i+1))
         IF (i.eq.1) dw = abs(lambda(j,i)-lambda(j,i+1))
         IF (i.eq.npix) dw = abs(lambda(j,i-1)-lambda(j,i))
         ewlim = ewlim + dw*abs(1.0d0-data(j,i))
         IF (dodebug) WRITE(6,*) 'Nlim',i,lambda(j,i),dw,
     &                           sqrt(ewlim),flag
 15    CONTINUE
       getNlim = ewlim/tau0
       RETURN

      END IF
   
 600  FORMAT(1x,' WARNING(getNlim): cannot find pixel for line ',i3,/,
     &          ' for ',a13,' at lambda = ',f8.2)

      END

c
c.......................................................................
c

      DOUBLE PRECISION FUNCTION getNaod(speci,ioni,linei,bpar,zcl,
     &                                  sigNaod)

c
c     this routine returns the AOD column density from the spectrum for
c     the ion within the region of the presnet component (line); the
c     aperture is taken as the b parameter of the VP fit
c
c     if there is saturation return the value as a lower limit, this is
c     flagged by setting the uncertainty to -1.0
c  
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  

      include            'minfit.h'
      include            'const.dek'
      integer             i,j,k,ioni,npix,upix,lpix,pixj
      integer             speci,linei,nsat
      double precision    Caod,dv,dw,wave,bpar,bres,doppb,sigNaod
      double precision    zcl,flux,dflux,tausum,dtausum
      double precision    taucore,dtaucore,Ncore,sigNcore
      double precision    z1,z2

c     define the ion number as j

      getNaod = 0.0d0   ! good measure
      sigNaod = 0.0d0   ! good measure

c     if the transition has no detection, then there is no point to
c     returning the AOD, set the value to null and return immediately

       IF (.not.features(ioni)) then
         getNaod = 1.0d0 ! avoid -Inf in the log10 (returns 0.00)
         sigNaod = 1.0d0
         RETURN
       END IF

c     otherwise proceed, set up local indices for short hand

      j    = ioni
      npix = ndata(j)

C      WRITE(6,*) 'N_AOD...'
C      WRITE(6,*) ion_name(j),lambda0(j),npix

c     compute the conversion constant, it is assumed that we are
c     integrating over velocity in [km/s]

c     the optical depth constant requires physical constants input from
c     const.dek; the factor of 1.E13 converts the speed of light to km/s
c     (1.E5) and the wavelength to centimeters (1.E8)

      Caod  = 1.0d13*(me*c/(pi*e*e))/(lambda0(j)*fosc(j)) 
      doppb = bpar

c     if this line is under resolved, reset the bpar to the resolution
c     of the spectrograph so that we are getting the AOD for a resolved
c     line at a minimum; if the line is very broad such that its bpar is
c     greater than BMAX input by the user, then reset to BMAX

      bres = ckms/(2.35*R_fac)  ! km/s
      IF (doppb.lt.bres) doppb = bres
      IF (doppb.gt.bmax) doppb = bmax

c     if we are saturated, we need to be sure to cover a large enough
c     range to obtain a meaningul AOD column.  If the bpar is narrow,
c     then the AOD will be understimated, so if this line is saturated
c     reset the bpar to the bamx given by the use in the .par file

      IF (satflag(speci,linei)) doppb = bmax

c     grab the pixel number corresponding to this line

      wave = lambda0(j)*(1.0d0+zcl)
      pixj = 0
      do 11 i=2,npix-1
       if ((wave.gt.lambda(j,i)).AND.
     &     (wave.le.lambda(j,i+1))) THEN
         pixj = i
         GOTO 12
       end if
 11   continue

c     TRAP: if we get here, then the pixel wasn't found, we do not want
c     to have DEATH because we are done fitting, null the AOD column and
c     return

      WRITE(SCREEN,600) linei,ion_name(j),wave
      WRITE(STDOUT,600) linei,ion_name(j),wave

      getNaod = 1.0d0 ! avoid -Inf in the log10 (returns 0.00)
      sigNaod = 0.0d0
      RETURN

c     compute the aperture to be the b parameter for the current line,
c     avoid edge effects

 12   dw = 0.0d0
      IF (pixj.eq.1)    dw = lambda(j,2)-lambda(j,1)
      IF (pixj.eq.npix) dw = lambda(j,npix)-lambda(j,npix-1)
      IF (dw.eq.0.0d0)  dw = lambda(j,pixj)-lambda(j,pixj-1)
      dv = (ckms/wave)*dw

c     set the initial pixel window for the calulation based upon the
c     DOPPB, but, we want to do this over an total aperture that is
c     1.5*FWHM, so the aperture on each side of the line center is
c     (1.5/2)*FWHM which is 2.35*(1.5/2)*FWHM 

      lpix  = pixj - nint(1.763d0*doppb/dv)
      upix  = pixj + nint(1.763d0*doppb/dv)

      IF (lpix.lt.1)    lpix = 1
      IF (upix.gt.npix) upix = npix

C     DEBUGGING?
      IF (dodebug) then
       z1 = lambda(j,lpix)/lambda0(ioni)-1.0d0
       z2 = lambda(j,upix)/lambda0(ioni)-1.0d0
       WRITE(6,*) 'NAOD',linei,doppb,lpix,upix,zcl,z1,z2
      END IF

c     check for masked pixels in the line core.  if the line core is
c     constained then this is considered a potentially constrainable
c     fit.  here, we define the line core as +/-1 pixel from the center
c     pixel PIXJ, if any of these three pixels are masked then we must
c     null the AOD coumn and return

      DO 13 i=lpix,upix
       IF (abs(pixj-i).le.1) then 
        IF (quality(j,i).eq.0) then
         getNaod = 1.0d0 ! avoid -Inf in the log10 (returns 0.00)
         sigNaod = 1.0d0
         RETURN
        END IF
       END IF
 13   CONTINUE

c     now integrate the optical depth and obtain its uncertainty, we
c     will assume that velocity width of all pixels are equal (though
c     they may be very slightly different)

      nsat     = 0
      tausum   = 0.0d0
      dtausum  = 0.0d0 
      taucore  = 0.0d0
      dtaucore = 0.0d0 

c     avoid masked pixels

      do 14 i=lpix,upix
       IF (quality(j,i).eq.1) then 
        flux  = data(j,i)
        dflux = abs(sigma(j,i))
        IF (flux.le.dflux) then
          flux = dflux  ! trap saturation
          nsat = nsat + 1
        END IF
        IF (flux.ne.1.0d0) then
         tausum  = tausum  - dv*sign(1.0d0,flux)*log(abs(flux))
         dtausum = dtausum + (dv*dflux/flux)**2 
         IF (abs(i-pixj).le.1) then
          taucore  = taucore  - dv*sign(1.0d0,flux)*log(abs(flux)) 
          dtaucore = dtaucore + (dv*dflux/flux)**2 
         END IF
        END IF
C     DEBUGGING?
        IF (dodebug) WRITE(6,*) 'NAOD',i,lambda(j,i),flux,
     &                          dflux,tausum,dtausum
       END IF
 14   CONTINUE
   
c     compute and return; if tausum<=0 then we force routine fixN and
c     getlimits to adopt the limit, we do this by nulling getNaod

      IF (tausum.le.0.0d0) then
        getNaod = 1.0d0  ! avoid -Inf in the log10 (returns 0.00)
        sigNaod = 0.0d0
        RETURN
      ELSE ! check saturation, if so key the limit
       getNaod = Caod*tausum
       Ncore   = Caod*taucore
       IF ((nsat.gt.0).OR.satflag(speci,linei)) then 
         sigNaod = -1.0d0
        ELSE
         sigNaod  = Caod*sqrt(dtausum)
         sigNcore = Caod*sqrt(dtaucore)
       END IF
      END IF

C     DEBUGGING?
      IF (dodebug) then
       WRITE(6,*) getNaod,sigNaod,getNaod-sigNaod,getNaod+sigNaod
      END IF

c     this really needs to be done correctly.  NEEDS FIXING SOMEDAY- it
c     is an attempt to handle parts of the spectra where the Naod is
c     dominated by the wing of a neighboring line... sanity check that a
c     discernable AOD is not outside the core, it should be dominated by
c     the core, if it is not, then set getNaod to consistent with 0.00

C      IF ((getNaod-Ncore)/Ncore.le.1.0d0) then
C        getNaod = 1.0d0  ! avoid -Inf in the log10 (returns 0.00)
C        sigNaod = 0.0d0
C      END IF



      RETURN

 600  FORMAT(1x,' WARNING(getNaod): cannot find pixel for line ',i3,/,
     &          ' for ',a13,' at lambda = ',f8.2)

      END

c
c.......................................................................
c

      SUBROUTINE         fixN(speciesi,fidx,linei,col,Temp,zcl,limflag)
c  
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include            'minfit.h'
      include            'const.dek'

      logical            limflag,flag
      integer            nwidth
      parameter          (nwidth = 1)
      integer            fidx,ioni,i,ii,j,linei,speciesi,numtran
      integer            tran,nt(2)
      double precision   N_lim(2),N_aod(2),dN_aod(2)
      double precision   zcl,col,Temp,fmin,fmax,bpar,Naodmin1,Naodmin2
      double precision   getNlim,getNaod,aodsig,fsig


      limflag =  .false. 

c     compute the b parameter for this species (1.0d-5 converts to km/s)

      bpar = 1.0d-5*sqrt(2.0d0*kerg*Temp/spec_mass(fidx))

c     find the ion and fosc for this species that has the weakest and
c     the strongest oscillator strength; if there is only one ionic
c     tranistion for the species, then NT(1)=NT(2) 

      fmin    =  1.0d10
      fmax    = -1.0d10
      DO 09 ioni=1,ions
        j = spec_idx(ioni)
        IF ((j.eq.fidx).AND.(features(ioni))) THEN
          IF (fosc(ioni).lt.fmin) THEN ! weak tranistion
            fmin  = fosc(ioni) 
            nt(1) = ioni
          END IF
          IF (fosc(ioni).gt.fmax) THEN ! strong tranistion
            fmax  = fosc(ioni)
            nt(2) = ioni
          END IF
        END IF
 09   CONTINUE


c     now check if there is only one ionic transition for this species,
c     if so the set numtran to unity

      IF (nt(1).eq.nt(2)) then
        numtran = 1
      ELSE
        numtran = 2  
      END IF

c     now obtain the limits for each of the selected transitions

      flag = .false.
      DO 21 j=1,numtran
       tran = nt(j)
       N_lim(j)  = getNlim(tran,linei,bpar,zcl,flag)
       N_aod(j)  = getNaod(speciesi,tran,linei,bpar,zcl,aodsig)
       dN_aod(j) = aodsig
C     DEBUGGING?
       IF (dodebug) then
        WRITE(6,*) j,ion_name(tran)(1:13),
     &           sign(1.0d0,N_lim(j))*log10(abs(N_lim(j))),
     &           sign(1.0d0,N_aod(j))*log10(abs(N_aod(j))),
     &           0.4343*dN_aod(j)/N_aod(j)
        Naodmin1 = N_aod(j) - dN_aod(j)
        WRITE(6,*) '                    ',
     &             sign(1.0d0,Naodmin1)*log10(abs(Naodmin1))
       END IF
 21   CONTINUE  ! next ionic transition

c
c     THE CHECKS AND TESTS FOLLOW
c

c     check for masked data, which has N_aod=0 and dN_aod=+1

      i  = 0
      ii = 0
      DO 22 j=1,numtran
        IF ((N_lim(j).eq.0.0d0).AND.(dN_aod(j).eq.1.0d0)) then
         i  = j
         ii = ii + 1
        END IF
 22   CONTINUE

c     if the second one is bad, discard, if the first one is bad, grab
c     the second one; if both are bad we are fucked

      IF (ii.eq.1) then  ! one is bad
        numtran = 1
        IF (i.eq.1) then ! if the 1st one, use the second one
          nt(1)     = nt(2)
          N_lim(1)  = N_lim(2)
          N_aod(1)  = N_aod(2)
          dN_aod(1) = dN_aod(2)
        END IF
      END IF

      IF (ii.eq.2) then ! both are bad, return without changing 
        WRITE(SCREEN,597) 
        WRITE(STDOUT,597) 
        RETURN
      END IF

 597  FORMAT(1x,'    data for this species is masked',
     &       1x,' -> cannot modify log(N)')


c     we now have the EW limits and limiting column densities (1-sigma)
c     for the single transition or for the weakest and strongest
c     transitions of the species; we now a perform a few tests determine
c     if the both columns are limits, one is a limit, or neither are
c     limits; limits are determined at the N_sigma level as specied by
c     the user in the MINFIT.PAR file


c     IF NUMTRAN=1 we take either the limit or the AOD column of the
c     single transition and hope for the best (because noise could fuck
c     us up) and return

      IF (numtran.eq.1) THEN

c     if saturation, DN_aod=-1; or if in a saturated line, adopt AOD

        IF ((dN_aod(1).eq.-1.0d0).OR.(satflag(speciesi,linei))) then
          col = 1.0d-13*N_aod(1)
          WRITE(SCREEN,598) ion_name(nt(1)),
     &                      log10(col)+13.0d0
          WRITE(STDOUT,598) ion_name(nt(1)),
     &                      log10(col)+13.0d0
          RETURN
        END IF

 598  FORMAT(1x,'    single tran AOD (saturation) for ',a8,
     &       1x,' -> adopted AOD log(N) = ',f6.2)

c     if AOD is consistent with LIMIT, then take the limit, if not then
c     take the AOD

        Naodmin1 = N_aod(1) - dN_aod(1)
        IF (Naodmin1.le.N_lim(1)) then     ! adopt limit, set LIMFLAG high
          col = 1.0d-13*N_lim(1)
          WRITE(SCREEN,600) ion_name(nt(1)),
     &                      log10(col)+13.0d0
          WRITE(STDOUT,600) ion_name(nt(1)),
     &                    log10(col)+13.0d0
          limflag = .true.
          RETURN
        ELSE                               ! adopt AOD
          col = 1.0d-13*N_aod(1)
          WRITE(SCREEN,601) ion_name(nt(1)),
     &                      log10(col)+13.0d0
          WRITE(STDOUT,601) ion_name(nt(1)),
     &                      log10(col)+13.0d0
          RETURN
        END IF

      END IF

 600  FORMAT(1x,'    single tran limit for ',a8,
     &       1x,' -> adopted LIMIT log(N) = ',f6.2)
 601  FORMAT(1x,'    single tran AOD for ',a8,
     &       1x,' -> adopted AOD log(N) = ',f6.2)

c     SO NUMTRAN=2 IS IT?!, we now must perform a few tests determine if
c     both are limits, one is a limit, or neither are limits; if both
c     are limits, take the most stringent limit, if one is a limit, we
c     will need to do some tests, if neither are limits then we continue
c     on and examine the AOD columns

      Naodmin1 = N_aod(1) - dN_aod(1)
      Naodmin2 = N_aod(2) - dN_aod(2)

c     are both limits?  take the minimum, set LIMFLAG high and return

      IF ((Naodmin1.le.N_lim(1)).AND.(Naodmin2.le.N_lim(2))) then
        col = 1.0d-13*min(N_lim(1),N_lim(2))        
        WRITE(SCREEN,602) ion_name(nt(1)),ion_name(nt(2)),
     &                    log10(col)+13.0d0
        WRITE(STDOUT,602) ion_name(nt(1)),ion_name(nt(2)),
     &                    log10(col)+13.0d0
        limflag = .true.
        RETURN
      END IF

 602  FORMAT(1x,'    limit for both ',a8,' & ',a8,
     &       ' -> adopted smaller LIMIT log(N) = ',f6.2)

c     is the strong transition a limit, but the weak one a detection;
c     this can happen if we have some noise in the weak transition; if
c     so, take the limit from the strong transition (because it provides
c     the most stringet constraint in the absences of noice), set
c     LIMFLAG high and return

      IF ((Naodmin1.gt.N_lim(1)).AND.(Naodmin2.le.N_lim(2))) then
        col = 1.0d-13*N_lim(2)        
        WRITE(SCREEN,603) ion_name(nt(1)),ion_name(nt(2)),
     &                    log10(col)+13.0d0
        WRITE(STDOUT,603) ion_name(nt(1)),ion_name(nt(2)),
     &                    log10(col)+13.0d0
        limflag = .true.
        RETURN
      END IF

 603  FORMAT(1x,'    AOD for ',a8,' & LIMIT for ',a8,
     &       ' -> adopted LIMIT log(N) = ',f6.2)


c     is the weak transition a limit, but the strong one a detection;
c     this can happen if we have a very weak detection (in the absence
c     of bad noise); if so, take the the AOD column of the stonger
c     transitions, return

      IF ((Naodmin1.le.N_lim(1)).AND.(Naodmin2.gt.N_lim(2))) then
        col = 1.0d-13*N_aod(2)
        WRITE(SCREEN,604) ion_name(nt(1)),ion_name(nt(2)),
     &                    log10(col)+13.0d0
        WRITE(STDOUT,604) ion_name(nt(1)),ion_name(nt(2)),
     &                    log10(col)+13.0d0
        RETURN
      END IF

 604  FORMAT(1x,'    limit from ',a8,' ; AOD from ',a8,
     &       ' -> adopted AOD log(N) = ',f6.2)


c     limits no longer play a role in the decision making

c     it is now assumed that both transitions are detections; then we
c     must check the AOD columns of both transitions; check for
c     consistency with saturation; if saturation condition is met, take
c     the AOD column from the weaker transition; if the two are
c     consistent with one another, take the average of the two
c     transitions; if not consistent and not saturated, take the smaller
c     column density; for this to work the weaker transitions must have
c     an AOD column that is larger than the stronger transition; check
c     this too.. if not?

c     check if consistent with saturation?  then check if the columns
c     have the correct sense?  if not, adopt the larger and return; if
c     so take the weaker transition column and return

      IF (((dN_aod(1).eq.-1.0d0).AND.(dN_aod(2).eq.-1.0d0)).OR.
     &    (satflag(speciesi,linei))) then ! saturation zone
      
        IF (N_aod(2).gt.N_aod(1)) THEN    ! incorrect sense
          col = 1.0d-13*max(N_aod(1),N_aod(2))
          WRITE(SCREEN,605) ion_name(nt(1)),ion_name(nt(2)),
     &                      log10(col)+13.0d0
          WRITE(6,*) '++>',N_aod(1),N_aod(2),max(N_aod(1),N_aod(2)),
     &                      (1.0d-13)*max(N_aod(1),N_aod(2)),
     &                      log10(max(N_aod(1),N_aod(2)))
          WRITE(STDOUT,605) ion_name(nt(1)),ion_name(nt(2)),
     &                    log10(col)+13.0d0
          RETURN
        ELSE                              ! correct sense
          col = 1.0d-13*N_aod(1)
          WRITE(SCREEN,606) ion_name(nt(1)),ion_name(nt(2)),
     &                      log10(col)+13.0d0
          WRITE(STDOUT,606) ion_name(nt(1)),ion_name(nt(2)),
     &                      log10(col)+13.0d0
          RETURN
        END IF

      END IF

 605  FORMAT(1x,'    AOD from ',a8,' & ',a8,' wrong sense',
     &       ' -> adopted max AOD log(N) = ',f6.2)
 606  FORMAT(1x,'    AOD from ',a8,' & ',a8,' both saturated',
     &       ' -> adopted AOD log(N) = ',f6.2, 
     &       ' from weaker transition')

c     so, both are not saturated, check if one is saturated and the
c     other not.  the strong transition should be the saturated one, but
c     check

c     correct sense
      IF ((dN_aod(1).ne.-1.0d0).AND.(dN_aod(2).eq.-1.0d0)) then  
        col = 1.0d-13*N_aod(1)
        WRITE(SCREEN,607) ion_name(nt(1)),ion_name(nt(2)),
     &                    log10(col)+13.0d0
        WRITE(STDOUT,607) ion_name(nt(1)),ion_name(nt(2)),
     &                    log10(col)+13.0d0
        RETURN
      END IF

 607  FORMAT(1x,'    AOD from ',a8,' saturated & AOD from ',a8,
     &      ' measured -> adopted measured AOD log(N) = ',f6.2)

c     wrong sense
      IF ((dN_aod(1).eq.-1.0d0).AND.(dN_aod(2).ne.-1.0d0)) then  
        col = 1.0d0-13*N_aod(2)
        WRITE(SCREEN,608) ion_name(nt(1)),ion_name(nt(2)),
     &                    log10(col)+13.0d0
        WRITE(STDOUT,608) ion_name(nt(1)),ion_name(nt(2)),
     &                    log10(col)+13.0d0
        RETURN
      END IF

 608  FORMAT(1x,'    AOD from ',a8,' measured & AOD from ',a8,
     &      ' saturated -> adopted measured AOD log(N) = ',f6.2,
     &      ' ** caution **')

c     so, neither are saturated, adopt the optimal weight

      col  = N_aod(1)/(dN_aod(1)**2) + N_aod(2)/(dN_aod(2)**2)
      fsig = 1.0d0/(dN_aod(1)**2) + 1.0d0/(dN_aod(2)**2)
      col  = 1.0d-13*(col/fsig)
      WRITE(SCREEN,609) ion_name(nt(1)),ion_name(nt(2)),
     &                  log10(col)+13.0d0
      WRITE(STDOUT,609) ion_name(nt(1)),ion_name(nt(2)),
     &                    log10(col)+13.0d0
      RETURN

 609  FORMAT(1x,'    AOD from ',a8,' & ',a8, 
     &       ' -> adopted weighted mean log(N) = ',f6.2)


c     we should never get here

      WRITE(SCREEN,*) ' DEATH(fixN): failed to determine N'
      WRITE(STDOUT,*) ' DEATH(fixN): failed to determine N'
      STOP

      END

c..............................................................................
c     eof

