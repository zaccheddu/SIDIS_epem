c......
c.... subroutine per FF di pione_2 e lambda_1 
 
c----------------------------------------------------------------------

c     AKK ROUTINES 2008

c----------------------------------------------------------------------      
C Input:

C IH hadron species
C***************************************************************
C*************** CHARGE-SIGN UNIDENTIFIED FFS ******************
C***************************************************************
C  1  pi^+ + pi^-
C  2  K^+ + K^-
C  3  p + anti-p
C  4  K0short
C  5  Lambda+anti-Lambda
C  6  pi^0, calculated from (pi^+ + pi^-)/2
C  7  K0short, calculated from (K^+ + K^-)/2 with u<->d
C  8  n + anti-n, calculated from p + anti-p with u<->d
C  9  h^+ + h^- (= pi^+ + pi^- + K^+ + K^- + p + anti-p)
C***************************************************************
C************** CHARGE-SIGN ASYMMETRY FFS **********************
C***************************************************************
C 10  h^+ - h^- (= pi^+ - pi^- + K^+ - K^- + p - anti-p)
C 11  pi^+ - pi^-
C 12  K^+ - K^-
C 13  p - anti-p

C Z: longitudinal-momentum fraction

C Q: fragmentation scale (in GeV)

C Output:

C DH(i): fragmentation function of parton i
C 0  1    2    3    4    5    6    7    8    9   10
C g  u  u-bar  d  d-bar  s  s-bar  c  c-bar  b  b-bar


       SUBROUTINE AKK(IH,Z,Q,DH)
      
       IMPLICIT REAL*8 (A-H,M-Z)
      DIMENSION DH(0:10),FF(0:5),FF_pion(0:5),
     ^     FF_kaon(0:5),FF_proton(0:5)
      
      Q2=Q*Q
cf2py intent(in) ih
cf2py intent(in) z
cf2py intent(in) q
cf2py intent(out) dh
      

      if (ih.eq.1) then
      call akk_Lambda(z,Q2,FF)
      DH(0)=FF(0)
      DH(1)=FF(1)
      DH(2)=DH(1)
      DH(3)=FF(2)
      DH(4)=DH(3)
      DH(5)=FF(3)
      DH(6)=DH(5)
      DH(7)=FF(4)
      DH(8)=DH(7)
      DH(9)=FF(5)
      DH(10)=DH(9)
      else if (ih.eq.2) then
      call akk_pion(z,Q2,FF)
      DH(0)=FF(0)
      DH(1)=FF(1)
      DH(2)=DH(1)
      DH(3)=FF(2)
      DH(4)=DH(3)
      DH(5)=FF(3)
      DH(6)=DH(5)
      DH(7)=FF(4)
      DH(8)=DH(7)
      DH(9)=FF(5)
      DH(10)=DH(9)
      else if (ih.eq.3) then
      call akk_kaon(z,Q2,FF)
      DH(0)=FF(0)
      DH(1)=FF(1)
      DH(2)=DH(1)
      DH(3)=FF(2)
      DH(4)=DH(3)
      DH(5)=FF(3)
      DH(6)=DH(5)
      DH(7)=FF(4)
      DH(8)=DH(7)
      DH(9)=FF(5)
      DH(10)=DH(9)
      else if (ih.eq.4) then
      call akk_pionCSA(z,Q2,FF)
      DH(0)=0.d0
      DH(1)=FF(1)
      DH(2)=-DH(1)
      DH(3)=FF(2)
      DH(4)=-DH(3)
      DH(5)=0.d0
      DH(6)=0.d0
      DH(7)=0.d0
      DH(8)=0.d0
      DH(9)=0.d0
      DH(10)=0.d0
      else if (ih.eq.5) then
      call akk_kaonCSA(z,Q2,FF)
      DH(0)=0.d0
      DH(1)=FF(1)
      DH(2)=-DH(1)
      DH(3)=0.d0
      DH(4)=0.d0
c..      write(6,*)'Warning for s->K+/- CSA: SU(3) used !!'
      DH(5)=-DH(1)
      DH(6)=-DH(2)
      DH(7)=0.d0
      DH(8)=0.d0
      DH(9)=0.d0
      DH(10)=0.d0
      end if



      end subroutine AKK
c.........
c......
ccccccc

       subroutine akk_Lambda(x,Q2,FFLambda)

      implicit none

      real*8 x,Q2,FFLambda(0:5)
      integer i,j

      integer nx,nQ2
      parameter(nx=500,nQ2=250)

      real*8 FFLambda_a(nQ2,nx,0:5)
      common/FF_interp_Lambda/FFLambda_a
      real*8 bLambda_a(nQ2,nx,0:5),
     .     cLambda_a(nQ2,nx,0:5),
     .     dLambda_a(nQ2,nx,0:5)
      common/splines_Lambda/bLambda_a,
     .     cLambda_a,dLambda_a

      integer precalc_done
      save precalc_done

      integer grid_read
      common/grid_read_com/grid_read

      if(precalc_done.ne.2375)then
         call akk_unit
         open(unit=grid_read,file='set_FF/grids_akk/grid_Lambda.txt',
     .        status='old')
         do i=1,nQ2
            do j=1,nx
               read(grid_read,'(6e14.6)')
     .              FFLambda_a(i,j,0),
     .              FFLambda_a(i,j,1),
     .              FFLambda_a(i,j,2),
     .              FFLambda_a(i,j,3),
     .              FFLambda_a(i,j,4),
     .              FFLambda_a(i,j,5)
            enddo
         enddo
         close(grid_read)
         call spline_fit(FFLambda_a,bLambda_a,
     .        cLambda_a,dLambda_a)
         precalc_done=2375
      endif
      call getFFs(x,Q2,FFLambda,FFLambda_a,
     .     bLambda_a,cLambda_a,dLambda_a)
      
      end subroutine
c------------------------------

      subroutine akk_pion(x,Q2,FFpion)

      implicit none

      real*8 x,Q2,FFpion(0:5)
      integer i,j

      integer nx,nQ2
      parameter(nx=500,nQ2=250)

      real*8 FFpion_a(nQ2,nx,0:5)
      common/FF_interp_pion/FFpion_a
      real*8 bpion_a(nQ2,nx,0:5),cpion_a(nQ2,nx,0:5),
     .     dpion_a(nQ2,nx,0:5)
      common/splines_pion/bpion_a,cpion_a,dpion_a

      integer precalc_done
      save precalc_done

      integer grid_read
      common/grid_read_com/grid_read

      if(precalc_done.ne.2375)then         
         call akk_unit
         open(unit=grid_read,file='set_FF/grids_akk/grid_pion',
     .        status='old')
         do i=1,nQ2
            do j=1,nx
               read(grid_read,'(6e14.6)')
     .              FFpion_a(i,j,0),
     .              FFpion_a(i,j,1),
     .              FFpion_a(i,j,2),FFpion_a(i,j,3),
     .              FFpion_a(i,j,4),FFpion_a(i,j,5)
            enddo
         enddo
         close(grid_read)
         call spline_fit(FFpion_a,bpion_a,cpion_a,
     .        dpion_a)
         precalc_done=2375
      endif
      call getFFs(x,Q2,FFpion,FFpion_a,bpion_a,
     .     cpion_a,dpion_a)
      
      end 
cc.....
c----------------------------------------------------------------
     
       subroutine akk_kaon(x,Q2,FFkaon)

      implicit none

      real*8 x,Q2,FFkaon(0:5)
      integer i,j

      integer nx,nQ2
      parameter(nx=500,nQ2=250)

      real*8 FFkaon_a(nQ2,nx,0:5)
      common/FF_interp_kaon/FFkaon_a
      real*8 bkaon_a(nQ2,nx,0:5),ckaon_a(nQ2,nx,0:5),
     .     dkaon_a(nQ2,nx,0:5)
      common/splines_kaon/bkaon_a,ckaon_a,dkaon_a

      integer precalc_done
      save precalc_done

      integer grid_read
      common/grid_read_com/grid_read

      if(precalc_done.ne.2375)then
         call akk_unit
         open(unit=grid_read,file='set_FF/grids_akk/grid_kaon',
     .        status='old')
         do i=1,nQ2
            do j=1,nx
               read(grid_read,'(6e14.6)')
     .              FFkaon_a(i,j,0),
     .              FFkaon_a(i,j,1),
     .              FFkaon_a(i,j,2),FFkaon_a(i,j,3),
     .              FFkaon_a(i,j,4),FFkaon_a(i,j,5)
            enddo
         enddo
         close(grid_read)
         call spline_fit(FFkaon_a,bkaon_a,ckaon_a,
     .        dkaon_a)
         precalc_done=2375
      endif
      call getFFs(x,Q2,FFkaon,FFkaon_a,bkaon_a,
     .     ckaon_a,dkaon_a)
      
      end 
c------------------------------------------------------------
  
       subroutine akk_pionCSA(x,Q2,FFpionCSA)

      implicit none

      real*8 x,Q2,FFpionCSA(0:5)
      integer i,j

      integer nx,nQ2
      parameter(nx=500,nQ2=250)

      real*8 FFpionCSA_a(nQ2,nx,0:5)
      common/FF_interp_pionCSA/FFpionCSA_a
      real*8 bpionCSA_a(nQ2,nx,0:5),cpionCSA_a(nQ2,nx,0:5),
     .     dpionCSA_a(nQ2,nx,0:5)
      common/splines_pionCSA/bpionCSA_a,cpionCSA_a,dpionCSA_a

      integer precalc_done
      save precalc_done

      integer grid_read
      common/grid_read_com/grid_read

      if(precalc_done.ne.2375)then         
         call akk_unit
         open(unit=grid_read,file='set_FF/grids_akk/grid_pionCSA',
     .        status='old')
         do i=1,nQ2
            do j=1,nx
               read(grid_read,'(6e14.6)')
     .              FFpionCSA_a(i,j,0),
     .              FFpionCSA_a(i,j,1),
     .              FFpionCSA_a(i,j,2),FFpionCSA_a(i,j,3),
     .              FFpionCSA_a(i,j,4),FFpionCSA_a(i,j,5)
            enddo
         enddo
         close(grid_read)
         call spline_fit(FFpionCSA_a,bpionCSA_a,cpionCSA_a,
     .        dpionCSA_a)
         precalc_done=2375
      endif
      call getFFs(x,Q2,FFpionCSA,FFpionCSA_a,bpionCSA_a,
     .     cpionCSA_a,dpionCSA_a)
      
      end
c----------------------------------------

       subroutine akk_kaonCSA(x,Q2,FFkaonCSA)

      implicit none

      real*8 x,Q2,FFkaonCSA(0:5)
      integer i,j

      integer nx,nQ2
      parameter(nx=500,nQ2=250)

      real*8 FFkaonCSA_a(nQ2,nx,0:5)
      common/FF_interp_kaonCSA/FFkaonCSA_a
      real*8 bkaonCSA_a(nQ2,nx,0:5),ckaonCSA_a(nQ2,nx,0:5),
     .     dkaonCSA_a(nQ2,nx,0:5)
      common/splines_kaonCSA/bkaonCSA_a,ckaonCSA_a,dkaonCSA_a

      integer precalc_done
      save precalc_done

      integer grid_read
      common/grid_read_com/grid_read

      if(precalc_done.ne.2375)then
         call akk_unit
         open(unit=grid_read,file='set_FF/grids_akk/grid_kaonCSA',
     .        status='old')
         do i=1,nQ2
            do j=1,nx
               read(grid_read,'(6e14.6)')
     .              FFkaonCSA_a(i,j,0),
     .              FFkaonCSA_a(i,j,1),
     .              FFkaonCSA_a(i,j,2),FFkaonCSA_a(i,j,3),
     .              FFkaonCSA_a(i,j,4),FFkaonCSA_a(i,j,5)
            enddo
         enddo
         close(grid_read)
         call spline_fit(FFkaonCSA_a,bkaonCSA_a,ckaonCSA_a,
     .        dkaonCSA_a)
         precalc_done=2375
      endif
      call getFFs(x,Q2,FFkaonCSA,FFkaonCSA_a,bkaonCSA_a,
     .     ckaonCSA_a,dkaonCSA_a)
      
      end 
c--------------------------------------------------

c......
 
       subroutine akk_unit

      implicit none

      integer unit_

      integer grid_read
      common/grid_read_com/grid_read

      logical opened_

      do unit_=10,300
         inquire (unit=unit_,opened=opened_)
         if(.not.opened_)then
            grid_read=unit_
            return
         endif
      enddo
      write(6,*)'There is no available I/O unit.'
      stop 

      end subroutine
c......
c......
c......

       subroutine getFFs(x,Q2,FF,FF_a,b_a,c_a,d_a)

      implicit none

      integer xpos,Q2pos,i,iabs,ifail

      integer nx,nQ2
      parameter(nx=500,nQ2=250)

      real*8 x_a(nx),Q2_a(nQ2)
      common/x_Q2_grid/x_a,Q2_a

      real*8 x,FF(0:5),tq,dx,FF1,FF2,Q2,
     .     FF_a(nQ2,nx,0:5),b_a(nQ2,nx,0:5),
     .     c_a(nQ2,nx,0:5),d_a(nQ2,nx,0:5)

      call neighbours(x,x_a,nx,xpos,ifail)
      if(ifail.eq.1)then
         write(6,'(A,f4.2,A,f5.3,A,f3.1)')
     .        'Error: z = ',x,
     .        ', must be in the range ',x_a(1),
     .        ' --- ',x_a(nx)
         stop
      endif

      call neighbours(Q2,Q2_a,nQ2,Q2pos,ifail)
      if(ifail.eq.1)then
         write(6,'(A,f5.1,A,f3.1,A,f5.1,A)')
     .        'Error: Q = ',dsqrt(Q2),
     .        ' GeV, must be in range ',dsqrt(Q2_a(1)),
     .        ' --- ',dsqrt(Q2_a(nQ2)),' GeV'
         stop
      endif

      do i=0,5
         iabs=i
         if(iabs.lt.0)iabs=-iabs
         dx=x-x_a(xpos)
         FF1=FF_a(Q2pos,xpos,iabs)
     .        +dx*(b_a(Q2pos,xpos,iabs)
     .        +dx*(c_a(Q2pos,xpos,iabs)
     .        +dx*d_a(Q2pos,xpos,iabs)))
         FF2=FF_a(Q2pos+1,xpos,iabs)
     .        +dx*(b_a(Q2pos+1,xpos,iabs)
     .        +dx*(c_a(Q2pos+1,xpos,iabs)
     .        +dx*d_a(Q2pos+1,xpos,iabs)))
         tq=(dlog(Q2)-dlog(Q2_a(Q2pos)))
     .        /(dlog(Q2_a(Q2pos+1))
     .        -dlog(Q2_a(Q2pos)))
         FF(i)=(1.d0-tq)*FF1+tq*FF2
      enddo

      end subroutine
c.................................................
c................................
ccccccc
       subroutine neighbours(val,val_a,max,j,ifail)

      implicit none

      integer max,j,ifail
      real*8 val,val_a(max)

      ifail=0
      if(val.le.val_a(1).or.val.gt.val_a(max))then
         ifail=1
      else
         do j=1,max-1
            if(val.gt.val_a(j).and.
     .           val.le.val_a(j+1))goto 1
         enddo
 1       continue
      endif

      return

      end subroutine
c......................
c.................
ccccccc
       subroutine spline_fit(FF_a,b_a,c_a,d_a)

      implicit none

      integer i,j,k

      integer nx,nQ2
      parameter(nx=500,nQ2=250)

      real*8 y(nx),FF_a(nQ2,nx,0:5),b(nx),c(nx),
     .     d(nx),b_a(nQ2,nx,0:5),c_a(nQ2,nx,0:5),
     .     d_a(nQ2,nx,0:5)

      real*8 x_a(nx),Q2_a(nQ2)
      common/x_Q2_grid/x_a,Q2_a

      integer precalc_done
      save precalc_done

      if(precalc_done.ne.2375)then
         call akk_x_Q2_grid
         precalc_done=2375
      endif

      do i=0,5
         do j=1,nQ2
            do k=1,nx
               y(k)=FF_a(j,k,i)
            enddo
            call spline(nx,x_a,y,b,c,d)
            do k=1,nx
               b_a(j,k,i)=b(k)
               c_a(j,k,i)=c(k)
               d_a(j,k,i)=d(k)
            enddo
         enddo
      enddo

      end subroutine
c........................
c...................
ccccccc
       SUBROUTINE SPLINE(N,X,Y,B,C,D)
*     ---------------------------------------------------------------------
*     CALCULATE THE COEFFICIENTS B,C,D IN A CUBIC SPLINE INTERPOLATION.
*     INTERPOLATION SUBROUTINES ARE TAKEN FROM
*     G.E. FORSYTHE, M.A. MALCOLM AND C.B. MOLER,
*     COMPUTER METHODS FOR MATHEMATICAL COMPUTATIONS (PRENTICE-HALL, 1977).
*     
*     SUBROUTINE TAKEN FROM AAC GROUP (KUMANO et al.)
*     
      integer N,K,IB,NM1
      REAL*8 X(N),Y(N),B(N),C(N),D(N),T
*     
      NM1=N-1
      IF(N.LT.2) then
         write(6,*)'N too small, =',N
         stop
      else
         D(1)=X(2)-X(1)
         C(2)=(Y(2)-Y(1))/D(1)
         DO  K=2,NM1
            D(K)=X(K+1)-X(K)
            B(K)=2.0D0*(D(K-1)+D(K))
            C(K+1)=(Y(K+1)-Y(K))/D(K)
            C(K)=C(K+1)-C(K)
         ENDDO
         B(1)=-D(1)
         B(N)=-D(N-1)
         C(1)=0.0D0
         C(N)=0.0D0
         IF(N.ne.3) then
            C(1)=C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
            C(N)=C(N-1)/(X(N)-X(N-2))-C(N-2)
     .           /(X(N-1)-X(N-3))
            C(1)=C(1)*D(1)**2.0D0/(X(4)-X(1))
            C(N)=-C(N)*D(N-1)**2.0D0/(X(N)-X(N-3))
         endif
         DO K=2,N
            T=D(K-1)/B(K-1)
            B(K)=B(K)-T*D(K-1)
            C(K)=C(K)-T*C(K-1)
         enddo
         C(N)=C(N)/B(N)
         DO IB=1,NM1
            K=N-IB
            C(K)=(C(K)-D(K)*C(K+1))/B(K)
         enddo
         B(N)=(Y(N)-Y(NM1))/D(NM1)
     1        +D(NM1)*(C(NM1)+2.0D0*C(N))
         DO  K=1,NM1
            B(K)=(Y(K+1)-Y(K))/D(K)
     1           -D(K)*(C(K+1)+2.0D0*C(K))
            D(K)=(C(K+1)-C(K))/D(K)
            C(K)=3.0D0*C(K)
         enddo
         C(N)=3.0D0*C(N)
         D(N)=D(N-1)
         RETURN
      endif
      B(1)=(Y(2)-Y(1))/(X(2)-X(1))
      C(1)=0.0D0
      D(1)=0.0D0
      B(2)=B(1)
      C(2)=0.0D0
      D(2)=0.0D0
      RETURN
      END SUBROUTINE
c......................
c....................
ccccccc
       subroutine akk_x_Q2_grid

      implicit none

      integer i

      integer nx,nQ2
      parameter(nx=500,nQ2=250)

      real*8 x_a(nx),Q2_a(nQ2)
      common/x_Q2_grid/x_a,Q2_a

      integer grid_read
      common/grid_read_com/grid_read

      open(unit=grid_read,file='set_FF/grids_akk/Q2_grid.txt',
     .     status='old')
      do i=1,nQ2-4,5
         read(grid_read,'(5e14.6)')Q2_a(i),
     .        Q2_a(i+1),Q2_a(i+2),
     .        Q2_a(i+3),Q2_a(i+4)
      enddo
      close(grid_read)

      open(unit=grid_read,file='set_FF/grids_akk/x_grid.txt',
     .     status='old')
      do i=1,nx-4,5
         read(grid_read,'(5e14.6)')x_a(i),
     .        x_a(i+1),x_a(i+2),
     .        x_a(i+3),x_a(i+4)
      enddo
      close(grid_read)

      end subroutine



















