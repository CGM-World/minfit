c..............................................................................
c  

      SUBROUTINE         errors(a,siga,n,m,error)

c  
c     compute the uncertainties in the fitted coefficients; this is a
c     key component to MINFIT and is one of the bottlenecks in terms of
c     time spent examining the fit
c
c     there is no analytical method for estimating how the model changes
c     with small changes in the fitting corefficients, so I had to write
c     this one by hand; our job is to compute the curvature matrix for
c     the chi^2 and then invert it to the covariance matrix.  We than
c     assume that the diagonals of the covariance matrix are the
c     uncertainties in the fitting parameters. 
c
c     the errors are stored in the array SIGA
c  
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      logical           error
      integer           jj,speciesi,linei,in,ib,iz,pixcnt,pixstart
      integer           i,j,k,ioni,pixi,n,m,indx(maxcoeffs),np
      integer           knt_dfda,knt_n,knt_b,knt_z
      double precision  a(maxcoeffs),siga(maxcoeffs),dparam,d,
     @                  d1mach
      double precision  dfda(maxcoeffs,maxvec),sig2(maxvec),
     @                  alpha(maxcoeffs,maxcoeffs),
     @                  covar(maxcoeffs,maxcoeffs)

      error    = .false.
      knt_dfda = 0
      knt_n   = 0
      knt_b   = 0
      knt_z   = 0

c     communicate

      WRITE(SCREEN,600) regioni,zlim(regioni,1),zlim(regioni,2),
     &                  lines,n,m
      WRITE(STDOUT,600) regioni,zlim(regioni,1),zlim(regioni,2),
     &                  lines,n,m

 600  FORMAT(1x,/,
     &       1x,'++++++++++++++++++++++++++++++++++++++++++++++++',/,
     &       1x,'+                                              +',/, 
     &       1x,'+        COMPUTING PARAMETER UNCERTAINTIES     +',/, 
     &       1x,'+                                              +',/, 
     &       1x,'+        Region              = ',i3,'             +',/,
     &       1x,'+        z_lower             = ',f8.6,'        +',/,
     &       1x,'+        z_upper             = ',f8.6,'        +',/,
     &       1x,'+                                              +',/, 
     &       1x,'+                                              +',/, 
     &       1x,'+       Number of Lines      = ',i5,'           +',/,        
     &       1x,'+       Number of Parameters = ',i5,'           +',/,        
     &       1x,'+       Number of Functions  = ',i5,'           +',/,        
     &       1x,'+                                              +',/, 
     &       1x,'++++++++++++++++++++++++++++++++++++++++++++++++',/)


c     zero the curvature (ALPHA) and derivative (DFDA) matrices

      do 02 j=1,n
       do 01 k=1,n
        alpha(j,k) = 0.0d0
        covar(j,k) = 0.0d0
 01    continue
       do 04 i=1,m
        dfda(j,i)  = 0.0d0
 04    continue
 02   continue

c     set the COVAR matrix to the identity matrix and set the SIGA = 0

      do 03 j=1,n
       covar(j,j) = 1.0d0
       siga(j)    = 0.0d0
 03   continue

c     compute uncertanties; this is the bottleneck for the entire
c     program...  computing the dfda is very important in order to get
c     useful errors; it is difficult to *not* obtain a singular
c     curvature matrix; the routine DPARAM does a good job of computing
c     numerical derivatives that are not 0; but it is not perfect; thus,
c     if DFDA is returned as 0; then it is reset to the machine
c     precision by function D1MACH

      WRITE(SCREEN ,*) ' - computing df/da elements ...'
      WRITE(STDOUT ,*) ' - computing df/da elements ...'



      do 27 linei=1,lines

        pixcnt = 0

       do 25 speciesi=1,species
        do 07 ioni=1,ions

         jj = spec_idx(ioni)
         if (features(ioni).AND.(jj.eq.fitindex(speciesi))) then      

c     for testing purposes
C          WRITE(SCREEN,*) ' ----------------------------------------'
C          WRITE(SCREEN,*) '   line       = ',linei
C          WRITE(SCREEN,*) '   species    = ',spec_name(jj)(1:8)
C          WRITE(SCREEN,*) '   transition = ',ion_name(ioni)(1:15)
C          WRITE(SCREEN,*) '   pixels     = ',ndata(ioni)

          in   = idx_an(speciesi,linei) 
          ib   = idx_ab(linei)
          iz   = idx_az(linei)

          do 06 pixi=1,ndata(ioni)

           i        = pixcnt + pixi
           knt_dfda = knt_dfda + 1
           sig2(i)  = sigma(ioni,pixi)*sigma(ioni,pixi)

           dfda(in,i) = dparam(a,in,speciesi,linei,ioni,pixi)
           if (dfda(in,i).eq.0.0d0) then
             dfda(in,i) = d1mach(4)
             knt_n = knt_n + 1
           end if
           dfda(ib,i) = dparam(a,ib,speciesi,linei,ioni,pixi)
           if (dfda(ib,i).eq.0.0d0) then
             dfda(ib,i) = d1mach(4)
             knt_b = knt_b + 1
           end if
           dfda(iz,i) = dparam(a,iz,speciesi,linei,ioni,pixi)
           if (dfda(iz,i).eq.0.0d0) then
             dfda(iz,i) = d1mach(4)
             knt_z = knt_z + 1
           end if

 06       continue

          pixcnt = i

         end if
 07     continue
 25    continue

 27   continue

      
      WRITE(SCREEN,*) '   -> df/da elements      = ',knt_dfda
      WRITE(SCREEN,*) '   -> dmach precision (N) = ',knt_n
      WRITE(SCREEN,*) '   -> dmach precision (b) = ',knt_b
      WRITE(SCREEN,*) '   -> dmach precision (z) = ',knt_z
      WRITE(SCREEN,*) '   -> number of a     (n) = ',n
      WRITE(SCREEN,*) '   -> number of func  (m) = ',m
      WRITE(SCREEN,*) '   -> number of pixels    = ',i

      WRITE(STDOUT,*) '   -> df/da elements      = ',knt_dfda
      WRITE(STDOUT,*) '   -> dmach precision (N) = ',knt_n
      WRITE(STDOUT,*) '   -> dmach precision (b) = ',knt_b
      WRITE(STDOUT,*) '   -> dmach precision (z) = ',knt_z
      WRITE(STDOUT,*) '   -> number of a     (n) = ',n
      WRITE(STDOUT,*) '   -> number of func  (m) = ',m
      WRITE(STDOUT,*) '   -> number of pixels    = ',i

c     for testing purposes: debug? write out the derivatives
C      do 1001 i=1,m
C       write(40,1101) (dfda(j,i),j=1,n)
C 1001  continue

c     compute the symmetric curvature matrix alpha(j,k)
c     we are assuming first derivatives of the basis functions in the
c     Hessian... if including 2nd derivatives, then we would invoke
c     alpha(j,k) = dfda(j,i)*dfda(k,i)/sig2(i)-fvec(i)*d2yd2a(j,k,i)

c     below the diagonal and diagonal inclusive

      WRITE(SCREEN,*) 
     &  ' - computing curvature matrix from df/da elements...'
      WRITE(STDOUT,*) 
     &  ' - computing curvature matrix from df/da elements...'

      i = 0
      do 09 j=1,n
       do 08 k=1,j
        do 20 i=1,m
C         write(42,*) j,k,i,dfda(j,i),dfda(k,i),sig2(i)
         alpha(j,k) = alpha(j,k) +  dfda(j,i)*dfda(k,i)/sig2(i)
 20     continue
 08    continue
 09   continue

c     above the diagonal

      WRITE(SCREEN,*) '   -> reflecting the diagonal...'
      WRITE(STDOUT,*) '   -> reflecting the diagonal...'

      do 11 j=2,n
       do 10 k=1,j-1
        alpha(k,j) = alpha(j,k)
 10    continue
 11   continue

c     for testing purposes; debug? write out ALPHA
C      WRITE(41,*) 'alpha'
C      do 1002 i=1,n
C       write(41,1101) (alpha(i,j),j=1,n)
C 1002 continue  


c     compute the coefficent uncertainties siga for an assumed unity
c     change in chi-squared (68% confidence level); the covariance
c     matrix is the inverse of the curvature matrix; invert alpha in two
c     steps first LU decompose, then back substitute row by row the sqrt
c     of the covariance diagonals give the estimated "un-correlated
c     normal" errors -- caveats apply

c     do LU decomposition

      WRITE(SCREEN,*) ' - computing the co-variance matrix...'
      WRITE(SCREEN,*) '   -> performing LU decomposition...'
      WRITE(STDOUT,*) ' - computing the co-variance matrix...'
      WRITE(STDOUT,*) '   -> performing LU decomposition...'


      np = maxcoeffs

      call ludcmp(alpha,n,np,indx,d,error)

c     if singular matrix, exit gracefully with siga(i) = 0.00

      if (.not.error) then

c     do back substitution

       WRITE(SCREEN,*) '   -> performing back substitution...'
       WRITE(STDOUT,*) '   -> performing back substitution...'

       do 12 i=1,n
        call lubksb(alpha,n,np,indx,covar(1,i))
        siga(i) = sqrt(abs(covar(i,i)))
 12    continue


c     for testing purposes
C       WRITE(41,*) 'errors'
C       do 13 i=1,n
C        WRITE(41,*) i,a(i),siga(i),siga(i)/a(i)
C 13   continue
C       PAUSE

      else

       write(SCREEN,*) ' ERROR: singular curvature matrix      *****'
       write(SCREEN,*) ' --- fit uncertainty estimates not available'
       write(STDOUT,*) ' ERROR: singular curvature matrix      *****'
       write(STDOUT,*) ' --- fit uncertainty estimates not available'

      end if
  
c     restuff the WRKFLX array, which was corrupted for computational
c     purposes

      WRITE(SCREEN,*) ' - restoring model vector...'
      WRITE(STDOUT,*) ' - restoring model vector...'

      CALL model(a,m)

c     return

      WRITE(SCREEN,*) ' - errors completed.'
      WRITE(SCREEN,*) ' '
      WRITE(STDOUT,*) ' - errors completed.'
      WRITE(STDOUT,*) ' '


      RETURN

c     formats

C 1001 continue  
 1101 format(1x,1p50e9.1)

      END

c
c..............................................................................
c

      DOUBLE PRECISION FUNCTION dparam(a,k,speciesi,linei,ion,pix)

c
c     compute the derivitives of the function wrt to the fit parmaters
c     this technique is based upon Neville's tableau; first order in the
c     Taylor's expansion; basically I did a major hack on a numerical
c     recipe subroutine! - don't tell Bill Press
c
c     basically, we adjust a fitting parameter up or down by a small
c     increment and then investigate how this effects the model fit;
c     this requires that we recompute the model... however, we do this
c     only for the small chunk of spectrum for which the parameter
c     applies, calls routine convpix
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      include           'minfit.h'
      integer           speciesi,linei
      integer           ntab,i,j,k,ion,pix,inmax,ibmax
      double precision  err,con,consq,big,safe,afac,fluxpix
      parameter         (con  = 1.4d0  , consq = con*con,
     @                   big  = 1.0d30 , ntab  =      10, 
     @                   safe = 2.0d0  , afac  =     0.01)
      double precision  errt,fac,hh,a(maxcoeffs),ak,da(ntab,ntab),
     @                  fup,fdn,asave


c  initialize

      asave = a(k)
      err   = big

c     the maximum index for N and b (note that these are the largest
c     indices SPECIES and LINE, not the current SPECIESI and LINEI

      inmax = idx_an(species,lines)
      ibmax = idx_ab(lines)

c     if column density, saturation checked in fluxpix call, set hh 

      if (k.le.inmax) then
        ak = a(k)
        if (satflag(speciesi,linei)) then
          hh  = 0.5*abs(ak)         ! finite difference for N by afac
        else
          hh  = afac*abs(ak)        ! finite difference for N by afac
        end if
      end if

c     if Temp (b param) then set hh for b param

      if ((k.gt.inmax).AND.(k.le.ibmax)) then
       ak = a(k)
       hh = afac*abs(ak)         ! finite difference for T,b by afac
      end if

c     if redshift, then set hh for z

      if (k.gt.ibmax) then
       ak = a(k)
       hh = 1.0d-5               ! finite difference redshift by 3 km/s
      end if

c     for testing purposes
C      WRITE(41,*) 'k=',k,'pix=',pix,'  a(k)=',ak,'   hh=',hh

c     first order derivatives; compute the first estimate

      a(k)    = ak + hh
      fup     = fluxpix(a,ion,pix)
      a(k)    = ak - hh
      fdn     = fluxpix(a,ion,pix)
      da(1,1) = 0.5d0 * (fup-fdn)/hh

c     compute the Neville table; based upon Num Recipe dfridr successive
c     columns in the Neville table have smaller hh and higher orders of
c     extrapolation; if you understand this, you are in good shape

      do 12 i=2,ntab
       hh      = hh/con
       a(k)    = ak + hh
       fup     = fluxpix(a,ion,pix)
       a(k)    = ak - hh
       fdn     = fluxpix(a,ion,pix)
       da(1,i) = 0.5d0 * (fup-fdn)/hh
       fac     = consq
       do 11 j=2,i
        da(j,i) = (fac*da(j-1,i)-da(j-1,i-1))/(fac-1.0d0)       
        fac    = fac * consq
        errt   = max(abs(da(j,i)-da(j-1,i)),abs(da(j,i)-da(j-1,i-1))) 
        if (errt.le.err) then
         err = errt
         dparam = da(j,i)
        end if
 11    continue

c     if higher order is worse by a significant factor (safe), then
c     return as is, reset a(i)

       if (abs(da(i,i)-da(i-1,i-1)).ge.safe*err) then
        a(k) = asave
        return
       end if

c     go to higher order

 12   continue   

c     reset a(k); after all, we really don't want to change the
c     parameter!

      a(k) = asave

c     return

      RETURN
      END

c
c..............................................................................
c
c     The following are slightly modified Numerical Recipes.  NR V.1987
c
c..............................................................................
c

      SUBROUTINE         ludcmp(a,n,np,indx,d,error)

c
c     given an n by n matrix a, with physical dimensions np by np, this
c     routine replaces a by the lu decompostion of a row-wise
c     permutation of itself.  input are a,n,np and a is output in the
c     form of equation 2.3.14 of numerical recipes. also output is indx
c     which records the row permutations effected by the partial
c     pivoting and d which is 1 if the number of interchanges is even,
c     -1 if odd. this routine is used in combination with the following
c     routine lubksb to solve linear equations or invert a matrix.
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c  

      include            'minfit.h'
      logical            error
      integer            n,np,indx(1),i,j,k,imax
      double precision   a,d,tiny,vv(maxcoeffs),aamax,sum,dum
      parameter          (tiny=1.0e-20)
      dimension          a(np,np)



      error = .false.
  
c     vv stores the implicit scaling of each row loop over the rows to
c     get the scaling information

      d = 1.0
      do 12 i=1,n
       aamax = 0.0
       do 11 j=1,n
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11     continue
       if (aamax.eq.0.) then

c     singular matrix, set error high and bail

        error = .true.
        return
       end if
       vv(i)=1./aamax
12    continue
  
c     for each column apply crouts method; see equation 2.3.12

      do 19 j=1,n
       if (j.gt.1) then
        do 14 i=1,j-1
         sum=a(i,j)
         if (i.gt.1)then
          do 13 k=1,i-1
           sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
         end if
14      continue
       end if
  
c     initialize the search for the largest pivot element
       
       aamax = 0.0
       do 16 i=j,n
        sum=a(i,j)
        if (j.gt.1)then
         do 15 k=1,j-1
          sum=sum-a(i,k)*a(k,j)
15       continue
         a(i,j)=sum
        end if
        dum = vv(i)*abs(sum)
        if (dum.ge.aamax) then
         imax=i
         aamax=dum
        end if
16     continue
  
c     if we need to interchange rows

       if (j.ne.imax)then
        do 17 k=1,n
         dum=a(imax,k)
         a(imax,k)=a(j,k)
         a(j,k)=dum
17      continue
        d=-d
        vv(imax)=vv(j)
       end if
  
c     now divide by the pivot element

       indx(j)=imax
       if (j.ne.n) then
        if(a(j,j).eq.0.)a(j,j)=tiny
        dum=1./a(j,j)
        do 18 i=j+1,n
         a(i,j)=a(i,j)*dum
18      continue
       end if
  
c     and go back for another column

19    continue

      if(a(n,n).eq.0.) a(n,n)=tiny

c     return
      return
      end

c  
c..............................................................................
c  

      SUBROUTINE         lubksb(a,n,np,indx,b)

c  
c     solves the set of n linear equations ax=b.  a is input in its lu
c     decomposition form, determined by the routine above ludcmp. indx
c     is input as the permutation vector also returned by ludcmp. b is
c     input as the right hand side vector and returns with the solution
c     vector x.  a,n and np are not modified by this routine and thus
c     can be left in place for successive calls (i.e matrix inversion)
c  
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer            n,np,indx(1),i,ii,j,ll
      double precision   a,b(1),sum
      dimension          a(np,np)

      
c     when ii is > 0 it becomes the index of the first nonzero element
c     of b this is forward substitution of equation 2.3.6, and the
c     unscamble takes place as one goes.

      ii = 0
      do 12 i=1,n
       ll = indx(i)
       sum = b(ll)
       b(ll) = b(i)
       if (ii.ne.0)then
        do 11 j=ii,i-1
         sum = sum - a(i,j) * b(j)
11      continue

c     a nonzero element was found, so from now on one does the sums in
c     the loop abovs

       else if (sum.ne.0.) then
        ii = i
       end if
       b(i) = sum
12    continue

c     back substitution equation 2.3.7

      do 14 i=n,1,-1
       sum = b(i)
       if(i.lt.n)then
        do 13 j=i+1,n
         sum = sum - a(i,j) * b(j)
13      continue
       end if

c     store a component of the solution vector x

       b(i) = sum/a(i,i)
14    continue

c     return

      return
      end

c     eof
