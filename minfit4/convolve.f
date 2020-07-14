c..............................................................................
c  
      SUBROUTINE          convolve(m,iflag)
c  
c     this routine sets up the call to the Num Recipes convolver
c     convolution resolution is higher than data resolution so that it
c     data are "smoother"
c
c     IFLAG=1 means include all transitions in the convolution
c     regardless of whehter there are features
c
c     IFLAG!=1 means include only ions with features
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      include             'minfit.h'
      integer             i,k,m,nmax2,ndat,ioni,pixi,mfft,iflag
      parameter           (nmax2 = maxcon*2)
      double precision    xa(maxvec),ya(maxvec),x,y,y2a(maxvec)
      double precision    ans(nmax2+2),phiov(maxcon)


       ndat    = m
       mfft    = nfft


c     stuff the "smoothing" data array, one only need keep the index as
c     the abscissa; then obtain the spline coefficients

      i = 0

      IF (iflag.eq.1) then   ! include all ions
       do 14 ioni=1,ions
         do 17 pixi=1,ndata(ioni)
          i = i + 1
          xa(i)  = real(i)
          ya(i)  = wrkflx(ioni,pixi)
          y2a(i) = 0.0d0
 17      continue
 14    continue

      ELSE

       do 24 ioni=1,ions   ! include only fitted ions
        if (features(ioni)) then
         do 27 pixi=1,ndata(ioni)
          i = i + 1
          xa(i)  = real(i)
          ya(i)  = wrkflx(ioni,pixi)
          y2a(i) = 0.0d0
 27      continue
        end if
 24    continue

      END IF

      ndat = i

      call spline(xa,ya,ndat,y2a)

c     stuff the end points and the interpolation points

      convdata(1)       = ya(1)
      convdata(ncondat) = ya(ndat)
      do 29 i=2,ncondat-1
       x = 1.0d0 + real(i-1)/resfac
       call splint(xa,ya,y2a,ndat,x,y)
       convdata(i) = y
 29   continue

c     pad the data array with unity, the continuum flux level this
c     assumes that no features are on the edges of the data array

      do 25 i=ncondat+1,nfft
       convdata(i) = 1.0d0
 25   continue

c     the FFT's do violence to the response function- so re-load it into
c     phi(v)

      do 26 i=1,nresponse
        phiov(i) = response(i)
 26   continue

c     we are "go" for the FFT's

      CALL convlv(convdata,mfft,phiov,nresponse,+1,ans)

c     now, restuff the wrkflx array and bail; pick off every resfac
c     element and stuff into idx element of wrkflx

      i = 0

      IF (iflag.eq.1) then

      do 31 ioni=1,ions
        do 33 pixi=1,ndata(ioni)
         i = i + 1
         k = 1 + (i-1)*int(resfac)
         wrkflx(ioni,pixi) = ans(k)
 33     continue
 31   continue

      ELSE

      do 41 ioni=1,ions
       if (features(ioni)) then
        do 43 pixi=1,ndata(ioni)
         i = i + 1
         k = 1 + (i-1)*int(resfac)
         wrkflx(ioni,pixi) = ans(k)
 43     continue
       end if
 41   continue
   
      END IF

c     return

      return
      end
c
c..............................................................................
c

      SUBROUTINE          spline(x,y,n,y2)   

c
c     given arrays x and y of length n containing a tabulated function
c     y=f(x), with the x monotonically increasing and given values for
c     yp1 and ypn, the first derivative at the pints 1 and n
c     respectively, this routine returns the array y2 of length n which
c     contains the second derivatives of the at the tabulated points x.
c
c     if yp1 and/or ypn are set larger than 1e30, the routine sets the
c     boundary condtion for a natural spline (one with zero second
c     derivative) at the boundary.
c
c     this routine is only called once for any given x and y.
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c

      implicit none
      integer             n,i,k,nmax 
      parameter           (nmax=10000) 
      double precision    x(n),y(n),y2(n),yp1,ypn,u(nmax),sig,p,qn,un
      parameter           (yp1=1.0e33, ypn=1.0e33)


c     the lower boundary condition is set to be either natural

      if (yp1 .gt. .99E30) then 
        y2(1)=0.0   
        u(1)=0.0

c     or it has a specified first derivative

      else  
        y2(1)=-0.5  
        u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1) 
      end if

c     this is the decomposition loop of the tridiagonal algorithm. y2
c     and u are used for temporary storage of the decomposition factors.

      do 11 i=2,n-1 
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))   
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.)/p
        u(i)=(6.*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1)) 
     +      /(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p   
11    continue  

 
c     the upper boundary condition is set to be either natural

      if (ypn .gt. .99E30) then 
        qn=0.0  
        un=0.0  
   
c     or it has a specified first derivative

      else  
        qn=0.5  
        un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1))) 
      end if

c     this is the backsubstitution loop of the tridiagonal algorithm

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)  
      do 12 k=n-1,1,-1  
        y2(k)=y2(k)*y2(k+1)+u(k)
12    continue  

c     return

      return
      end   

c
c............................................................................. 
c 

      SUBROUTINE          splint(xa,ya,y2a,n,x,y)

c
c     given the arrays xa and ya of length n, which tablulate a
c     monotonic func and given the array y2a, which is the output of
c     spline (above), and given a value of x this routine returns a
c     cubic spline interpolated value of y.
c
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c   

      include            'minfit.h'
      integer             n,klo,khi,k   
      double precision    xa(n),ya(n),y2a(n),x,y,h,a,b  

c     find the right place in the table by bisection. this is optimal if
c     the sequential calls to this routine are at random values of x.
c     if the sequential calls are in order and closely spaced, one might
c     store the values of klo and khi and test if they remain
c     appropriate on next call

      klo=1 
      khi=n 
1     if (khi-klo .gt. 1) then  
       k=(khi+klo)/2
       if (xa(k) .gt. x) then   
        khi=k   
       else 
        klo=k   
       end if   
       goto 1   
      end if    

c     klo and khi now bracket the input value of x the xa's must be
c     distinct

      h=xa(khi)-xa(klo) 
      if (h .eq. 0.0d0) then
        write(SCREEN,99) khi,xa(khi),klo,xa(klo)
        write(STDOUT,99) khi,xa(khi),klo,xa(klo)
   99   format(1x,'  khi=',i4,'  xa(khi)=',1pe13.5,'  klo=',
     1         i4,'  xa(klo)=',1pe13.5)
        stop 'bad xa input in routine splint' 
      end if
   
c     evaluate the cubic spline

      a=(xa(khi)-x)/h   
      b=(x-xa(klo))/h   
      y=a*ya(klo)+b*ya(khi)+
     +      ((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0

c     return

      return    
      end   
c  
c..............................................................................
c eof
