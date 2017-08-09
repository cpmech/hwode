cc    This module calculates the r.h.s for the KS equation
cc    
cc    d/dt u = - (dx^2+ dx^4)u - u u'
cc    
cc    in its form with fourier modes.
cc    I.e., the system
cc    
cc    dt u_n = (n^2 q^2 - n^4 q^4) u_n 
cc    - (i q n/2) sum_{n_1+n_2=n} u_{n_1} u_{n_2}
cc    
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      
      subroutine FKS(ND,T,Y,F,RPAR,IPAR)
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER (MMM=9)
      PARAMETER (NH=2**MMM,N=2*NH)      
      DIMENSION y(ND),f(ND)
      DIMENSION U(N)
        COMMON/TRANS/QQ,UZERO
      
c ---    copy y to u
        U(1)=UZERO
        U(2)=0.D0
        DO I=3,N
         U(I)=Y(I-2)
        END DO
        CALL REALFT(U,NH,-1)
        DO I=1,N
         U(I)=2.D0*U(I)
        END DO
C ---    SQUARE ------
        DO I=1,N
         U(I)=U(I)**2/2.D0
        END DO
C --- TRANSFORM BACK ---
        CALL REALFT(U,NH,+1)
        AN=N
        DO I=1,N
         U(I)=U(I)/AN
        END DO
      do J=1,NH-1
        DIAG=(QQ*J)**2*(1.D0-(QQ*J)**2)
        F(2*J-1)=DIAG*Y(2*J-1) +QQ*J*U(2*J+2)
        F(2*J  )=DIAG*Y(2*J  ) -QQ*J*U(2*J+1)
      enddo
      end
c--------------------------------------------------------------
      subroutine JKS(ND,x,y,dfy,ldfy,RPAR,IPAR)
      PARAMETER (MMM=9)
      PARAMETER (NH=2**MMM,N=2*NH)      
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION y(ND),DFY(LDFY,ND)
        COMMON/TRANS/QQ,UZERO
      do J=1,NH-1
        DIAG=(QQ*J)**2*(1.D0-(QQ*J)**2)
        DFY(1,2*J-1)=DIAG
        DFY(1,2*J  )=DIAG
      enddo
      end
      
      
      subroutine fft(xr,n,isign)
      integer n
      integer isign
      real*8  xr(2*n)
c     ***************************************************************
c     *  the fft computes the discrete fast fourier transform of a  *
c     *  sequence of n terms.                                       *
c     *  the forward fft computes                                   *
c     *       y(j)= sum (from k=0 to n-1) x(k)*exp(2*pi*i*j*k/n)    *
c     *  the backward fft computes                                  *
c     *       y(j)= sum (from k=0 to n-1) x(k)*exp(-2*pi*i*j*k/n)   *
c     *                                                             *
c     *  x     is a complex array of length n.                      *
c     *  n     is a power of 2. n<=4096.                            *
c     *  isign is the direction of the transform.if isign>=0 then   *
c     *        the fft is forward, otherwise backward.              *
c     *                                                             *
c     *  reference                                                  *
c     *  j.w.cooley, p.a.w.lewis, and p.d.welch, 'the fast fourier  *
c     *  transform and its applications', ieee trans. on education, *
c     *  vol.e-12, no.1, (march 1969),p.29.                         *
c     *  the fft is a modified version of the program described in  *
c     *  this reference.the modification is the use of a table to   *
c     *  store the roots of unity.                                  *
c     ***************************************************************
      integer maxn
      parameter(maxn=8192)
      integer i
      integer icnt,ii,iir,ijr,j,jj,k,m,ntbl
      real*8  cstore(2*maxn),dcos,dsin,datan
      real*8  pi,temp,vr1,vr2,vrr1,vrr2,wr1,wr2,xr1,xr2
      save cstore
      save ntbl
      data ntbl/0/
      
      
c     ****  store table of roots of unity ****
      
c     the roots of unity exp(pi*i*k/j) for j=1,2,4,..,n/2 and
c     k=0,1,2,..,j-1 are computed once and stored in a table.
c     this table is used in subsequent calls of fft with parameter
c     n <= ntbl.
      
      if(n.gt.maxn) then
         print *,'n too big in fft',n
         stop
      endif
      if (n.gt.ntbl) then
         ntbl=n
         pi=datan(1.0d0)*4.0d0
         j=1
         icnt=0
 10      continue   
c     s=pi*(0,1)/j
         do 20 k=0,j-1
            icnt=icnt+1
            cstore(icnt+icnt-1)=dcos(k*(pi/j))
            cstore(icnt+icnt  )=dsin(k*(pi/j))
 20      continue
         j=j+j
         if (j.lt.n) go to 10
      end if
      
c     ****  bit reversal ****
      
c     the x(j) are permuted, in such a way that each new place
c     number j is the bit reverse of the original place number.
      
      j=1
      do 30 i=1,n
         if (i.lt.j) then
            jj=2*j
            ii=2*i
            vrr1=xr(jj-1)
            xr(jj-1)=xr(ii-1)
            xr(ii-1)=vrr1
            vrr2=xr(jj)
            xr(jj)=xr(ii)
            xr(ii)=vrr2
         end if
         m=n/2
 25      continue
         if (j.gt.m) then
            j=j-m
            m=m/2
            if (m.ge.1) go to 25
         else
            j=j+m
         end if
 30   continue
      
c     ****  matrix multiplication ****
      
c     the roots of unity and the x(j) are multiplied.
      
      j=1
      icnt=0
 40   jj=j+j
      do 50 k=1,j
         icnt=icnt+1
         wr1=cstore(2*icnt-1)
         wr2=cstore(2*icnt  )
         if (isign.lt.0) wr2=-wr2
         
         do 50 i=k,n,jj
            iir=i+i
            ijr=(i+j)*2
            xr1=xr(ijr-1)
            xr2=xr(ijr)
            vr1=wr1*xr1
            vr1=vr1-wr2*xr2
            temp=xr(iir-1)
            xr(iir-1)=temp+vr1
            xr(ijr-1)=temp-vr1
            vr2=wr1*xr2+wr2*xr1
            temp=xr(iir)
            xr(ijr)=temp-vr2
            xr(iir)=temp+vr2
 50      continue
         j=jj
         if (j.lt.n) go to 40
       return
      end
c-----------------------------
      
      subroutine realft(x,n,isign)
      integer n
      integer isign
      real*8  x(2*n)
      
      integer i,i1,i2,i3,i4,n2p3,ppp
      real*8 c1,c2,h1r,h1i,h2r,h2i
      real*8 wr,wi,wpr,wpi,wtemp,theta,pi,dsin,datan
      real*8 x1,x2,x3,x4,uu,vv,wr0,wr1,wr2
      c1=0.5
      pi=datan(1.0D0)*4
      theta=pi/ n
      if (isign.eq.1) then
         c2 = -0.5D+00
         call fft(x,n,1)
      else 
         c2=0.5D+00
         theta = -theta
      endif
      wtemp=dsin(0.5*theta)
      wpr = -2.0*wtemp*wtemp
      wpi=dsin(theta)
      wr=1.0+wpr
      wi=wpi
      n2p3=2*n+3
      do i=2,n/2
         ppp=1
         if(ppp.eq.0)then
c     ... normal code
            i1=2*i-1
            i2=i1+1
            i3=n2p3-i2
            i4=i3+1
            h1r=c1*(x(i1)+x(i3))
            h1i=c1*(x(i2)-x(i4))
            h2r = -c2*(x(i2)+x(i4))
            h2i=c2*(x(i1)-x(i3))
            x(i1)=   h1r +wr*h2r-wi*h2i
            x(i3)=   h1r -wr*h2r+wi*h2i
            x(i2)=   h1i +wr*h2i+wi*h2r
            x(i4) = -h1i +wr*h2i+wi*h2r
            wtemp=wr
            wr=wtemp*wpr-wi*wpi+wr
            wi=wi*wpr+wtemp*wpi+wi
         else
c     ... optimized for apollo
            i1=2*i-1
            i3=(n2p3-1)-i1
            x1=x(i1)
            x3=x(i3)
            h2i=c2*(x1-x3)
            uu= wi*h2i
            vv= wr*h2i
            h1r=c1*(x1+x3)
            x2=x(i1+1)
            x4=x(i3+1)
            h2r = -c2*(x2+x4)
            uu=wr*h2r-uu
            vv=wi*h2r+vv
            x(i1)=h1r+uu
            x(i3)=h1r-uu
            h1i =  c1*(x2-x4)
            x(i1+1)=   h1i+vv
            x(i3+1) = -h1i+vv
            wr0=wi*wpi
            wr1=(1+wpr)
            wi=wi*wr1
            wtemp=wr
            wr2=wtemp*wr1
            wi=wtemp*wpi+wi
            wr=wr2-wr0
         endif
      enddo
      if (isign.eq.1) then
         h1r=x(1)
         x(1) = h1r+x(2)
         x(2) = h1r-x(2)
      else 
         h1r=x(1)
         x(1)=c1*(h1r+x(2))
         x(2)=c1*(h1r-x(2))
         call fft(x,n,-1)
      endif
      end
