c----------------------------------------------------------------
c     SNG(mui,muf) 
c----------------------------------------------------------------
      function sng(q,b)
      implicit none
      real*8 b,sng,beta0,q,qb,aa,bb,cc,mui,euler,CAF,alphas
      real*8 u,CA,CF,Qini,zzh,kkt,QQ,pi,bessel0,coll_ff
      real*8 bmax,c0,Q0
      integer iq,iqg

      bmax = 0.8d0

*      bmax = 1.5d0
*      euler = 0.577216
*      qb=b/dsqrt(1d0+(b/1.5d0)**2d0)
*      mui = 2d0*dexp(-euler)/qb
      qb = b/dsqrt(1d0+(b/bmax)**2d0)
      c0 = 1.122919d0
      Q0 = c0/qb

      CA = 3d0
      CF = 4d0/3d0
      pi = 4d0*atan(1d0)
      CAF = CA*CF*(pi**2d0)/3d0
      beta0 = 11d0-2d0*3d0/3d0
      u = dlog(alphas(Q0)/alphas(q))/beta0
*      u = -dlog(1-2d0*beta0*alphas(q)*dlog(mui))/(4d0*pi*beta0)
      aa = 0.85*CA
      bb = 0.86*CA
      cc = 1.33d0
      sng = dexp(-CAF*(u**2d0)*(1+(aa*u)**2d0)
     & /(1+((bb*u)**2d0)**(cc/2d0)))

      return
      end


c--------------------------------------------------------------
c     Perbative evolution kernel with bstar: Revo
c     iq = 1 / 0 for quark or gluon
c--------------------------------------------------------------
      subroutine evolve(iq,b,Q,Q0,Revo)
Cf2py intent(in)  iq
Cf2py intent(in)  b
Cf2py intent(in)  Q
Cf2py intent(in)  Q0
Cf2py intent(out) Revo
      implicit none
      real*8 Revo,b,Q0,Q,beta0,beta1,bmax,bstar,bdiv
      real*8 c0,pi,CF,CA,mb,alam,alam4,alam5,funcsp
      real*8 b0,b1b0,rnf,xq,xup,xlow,xuplog,xlowlog,sp1,sp2,sp
      integer nbel,nbuu,nz,nf,iq
      
      bdiv = 10d0
      nbel = 10 
      nbuu = 5
      nz   = 5

      c0 = 1.122919d0
      pi = 4d0*atan(1d0)
      CF = 4d0/3d0
      CA = 3d0
      mb = 4.5d0
      alam4 = 0.326d0
      alam5 = 0.226d0
      bmax = 1.5d0
      
      bstar = b/dsqrt(1d0+(b/bmax)**2d0)
      Q0 = c0/bstar
      if(Q.le.mb.and.Q0.le.mb) then
         nf = 4
         rnf = 4d0
         beta0 = 11d0-2d0*rnf/3d0
         beta1 = 102d0-38d0*rnf/3d0
         b0 = beta0
         b1b0  = beta1/(beta0**2d0)
         alam = alam4
         xq = 2d0*dlog(Q/alam)
         xup = xq
         xlow = 2d0 * dlog(Q0/alam)
         xuplog = dlog(xup)
         xlowlog = dlog(xlow)
         sp = funcsp(iq,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)
      elseif(Q.ge.mb.and.Q0.ge.mb) then
         nf = 5
         rnf = 5d0
         beta0 = 11d0-2d0*rnf/3d0
         beta1 = 102d0-38d0*rnf/3d0
         b0 = beta0
         b1b0  = beta1/(beta0**2d0)
         alam = alam5
         xq = 2d0*dlog(Q/alam)
         xup = xq
         xlow = 2d0 * dlog(Q0/alam)
         xuplog = dlog(xup)
         xlowlog = dlog(xlow)
         sp = funcsp(iq,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)
      elseif(Q.ge.mb.and.Q0.le.mb) then
         nf = 5
         rnf = 5d0
         beta0 = 11d0-2d0*rnf/3d0
         beta1 = 102d0-38d0*rnf/3d0
         b0 = beta0
         b1b0  = beta1/(beta0**2d0)
         alam = alam5
         xq = 2d0*dlog(Q/alam)
         xup = xq
         xlow =2d0*dlog(mb/alam)
         xuplog = dlog(xup)
         xlowlog = dlog(xlow)
         sp1 = funcsp(iq,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)

         nf = 4
         rnf = 4d0
         beta0 = 11d0-2d0*rnf/3d0
         beta1 = 102d0-38d0*rnf/3d0
         b0 = beta0
         b1b0  = beta1/(beta0**2d0)
         alam = alam4
         xq = 2d0*dlog(Q/alam)
         xup = 2d0*dlog(mb/alam)
         xlow = 2d0*dlog(Q0/alam)
         xuplog = dlog(xup)
         xlowlog = dlog(xlow)
         sp2 = funcsp(iq,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)
         sp = sp1 + sp2
      elseif(Q.le.mb.and.Q0.ge.mb) then
         nf = 4
         rnf = 4d0
         beta0 = 11d0-2d0*rnf/3d0
         beta1 = 102d0-38d0*rnf/3d0
         b0 = beta0
         b1b0  = beta1/(beta0**2d0)
         alam = alam4
         xq = 2d0*dlog(Q/alam)
         xup = xq
         xlow =2d0*dlog(mb/alam)
         xuplog = dlog(xup)
         xlowlog = dlog(xlow)
         sp1 = funcsp(iq,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)

         nf = 5
         rnf = 5d0
         beta0 = 11d0-2d0*rnf/3d0
         beta1 = 102d0-38d0*rnf/3d0
         b0 = beta0
         b1b0  = beta1/(beta0**2d0)
         alam = alam5
         xq = 2d0*dlog(Q/alam)
         xup = 2d0*dlog(mb/alam)
         xlow = 2d0*dlog(Q0/alam)
         xuplog = dlog(xup)
         xlowlog = dlog(xlow)
         sp2 = funcsp(iq,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)
         sp = sp1 + sp2
      endif

      Revo = dexp( -sp )
      
      return
      end

c--------------------------------------------------------------
c     Perturbative Sudakov exponent: integral of \gamma_F
c--------------------------------------------------------------

      function funcsp(iqg,xq,xup,xuplog,xlow,xlowlog,b0,b1b0,nf)
      implicit none
      real*8 funcsp,xq,xup,xuplog,xlow,xlowlog,b0,b1b0
      real*8 f00,f01,f11,f02,f12,f03,f13,f23,f04,f14,f24
      real*8 CF,CA,pi,G0,G1,G0V,sp0,sp0V,sp1,alam4,alam5,mb,c0,bmax
      integer nf,iqg

      c0 = 1.122919d0
      pi = 4d0*atan(1d0)
      CF = 4d0/3d0
      CA = 3d0
      mb = 4.5d0
      alam4 = 0.326d0
      alam5 = 0.226d0
      bmax = 1.5d0
      
      if(iqg.eq.1) then         ! for quark TMD
         G0  =  4d0 * CF
         G0V = -6d0 * CF
         G1  =  4d0*CF*( (67d0/9d0 - pi*pi/3d0)*CA - 10d0/9d0*nf )
      elseif(iqg.eq.0) then     ! for gluon TMD
         G0  =  4d0 * CA
         G0V = -2d0*(11d0/3d0*CA - 2d0/3d0*nf)
         G1  =  4d0*CA*( (67d0/9d0 - pi*pi/3d0)*CA - 10d0/9d0*nf )
      else
         print *, 'iqg should be 0 (1) for g (q): ', iqg
      endif
         
      f00 = xup - xlow
      f01 = xuplog - xlowlog
      f11 = ( xuplog**2d0 - xlowlog**2d0 ) / 2d0
      f02 = 1d0/xlow - 1d0/xup
      f12 = xlowlog/xlow - xuplog/xup + f02
      f03 = (1d0/xlow**2d0 - 1d0/xup**2d0) / 2d0
      f13 = (xlowlog/xlow**2d0 - xuplog/xup**2d0)/2d0 + f03/2d0
      f23 = ((xlowlog/xlow)**2d0 - (xuplog/xup)**2d0 )/2d0 + f13
      f04 = (1d0/xlow**3d0 - 1d0/xup**3d0)/3d0
      f14 = (xlowlog/xlow**3d0 - xuplog/xup**3d0)/3d0 + f04/3d0
      f24 = (xlowlog**2d0/xlow**3d0 - xuplog**2d0/xup**3d0)/3d0 
     >     + f14*2d0/3d0

      sp0 = G0 /(2d0*b0) * ( (xq * f01 - f00) - b1b0*(xq * f12 - f11) )
      sp0V = G0V /(2d0*b0) * ( f01 - b1b0 * f12 )
      sp1 = G1 /(2d0*b0*b0) * ( (xq*f02 - f01) - 2d0*b1b0*(xq*f13-f12)
     >     +b1b0**2d0*(xq*f24-f23) )

      funcsp = sp0 + sp0V + sp1

      return
      end

      
c--------------------------------------------------------------
c     alphas(q) -- NLO strong coupling constant (from CTEQ)
c--------------------------------------------------------------      
      function alphas(q)
      implicit none
      real*8 q,q2,lambda,lambda2,alphas,b0,b1,tt,pi,mb
      integer nf

      pi = atan(1d0)*4d0
      mb = 4.5d0
      nf=4
      b0=11d0-2d0/3d0*nf
      b1=102d0-38d0/3d0*nf
      lambda=0.227506d0

      if (q.le.0.23d0) then
        q=0.23d0
      endif
      q2=q*q
      lambda2=lambda*lambda
      tt=dlog(q2/lambda2)
      alphas=4d0*pi/(b0*tt)*(1d0-b1/(b0*b0)*dlog(tt)/tt)
      pi = atan(1d0)*4d0
      mb = 4.5d0
      if (q.le.mb) then
         nf=4
         lambda=0.326d0
      else
         nf=5
         lambda=0.226d0
      endif
      b0=11d0-2d0/3d0*nf
      b1=102d0-38d0/3d0*nf
*
      q2=q*q
      lambda2=lambda*lambda
      tt=dlog(q2/lambda2)
      alphas=4d0*pi/(b0*tt)*(1d0-b1/(b0*b0)*dlog(tt)/tt)
      return
      end
