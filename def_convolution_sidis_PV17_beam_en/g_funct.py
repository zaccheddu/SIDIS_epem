# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
import lhapdf






class g_funct():

	CA= 3.
	CF= 4/3
	nf= 3
	beta0 = 11/3*CA -2/3*nf
	beta1 = 34/3*CA**2 - 10/3*CA*nf - 2*CF*nf

	gmk1=16*CF
	gmk2=2*CA*CF*(536/9 - 8*pi**2/3) -160/9*CF*nf
	gmD =6*CF

	pdf=lhapdf.mkAlphaS("NNFF10_PIp_nlo")

	def __init__(self):
		
		self.bmax = 1.
		self.mu= 10.58

		self.a = 0.05 ## GeV**2
		

		self.bmin = 2*e**(-euler_gamma)/10.58	


	def bt_str(self,bt):  #define bt_star
	
		bm= self.bmax	
		
		b_new = sqrt(bt**2 + self.bmin**2)

		bstar = b_new/sqrt(1+b_new**2/bm**2)			
		#print('in bt_star:' + str(bt))		
		return bstar

	def mu_b(self,bt):  #define mu_b

		bstr=self.bt_str(bt)

		mub=2*e**(-euler_gamma)/bstr
		#print('in mu_b:' + str(bt))
		return mub

	def xx(self,bt): #define x
		
		alpha_s = self.pdf.alphasQ(self.mu)
		n = alpha_s/4/pi
		#print(n)
		#print(self.mu)
		ll = log(self.mu/self.mu_b(bt))
		ll = n*ll
		
		return ll #, self.mu/self.mu_b(bt)

	def gK(self,bt):

		gg = self.a
		gg = gg*bt**2

		return gg



	def g1(self,bt):	

		x=self.xx(bt)

		#print('valore x= ' + str(x))
		#x=bt
		bx = 2*self.beta0*x
		#print(2*beta0)
		gg = self.gmk1/4/self.beta0
		#print(gg)	

		gg = gg*(1 + log(1-bx)/(bx) )
		#gg = gg*(log(1-bx)/(bx) )
		return gg

	def g2(self,bt):

		x=self.xx(bt)
		#print('valore x= ' + str(x))

		bx = 2*self.beta0*x

		g1 = self.gmk1/4/self.beta0
		g1 = g1*self.beta1/self.beta0
		g1 = g1*( x/(1-bx) + ( log(1-bx) + ( log(1-bx) )**2/2  )/2/self.beta0  )

		g2 = - self.gmk2/8/self.beta0**2
		g2 = g2*( bx/(1-bx) + log(1-bx) )

		g3 = -self.gmD/2/self.beta0
		g3 = g3*log(1-bx)

		#print(2*beta0)
		#print(gmk1/4/beta0*beta1/beta0)
		#print(gmk2/8/beta0)
		#print(-gmD/2/beta0)

		gg= g1 + g2 + g3

		return gg

	def g2K(self,bt):

		x=self.xx(bt)
		#print('valore x= ' + str(x))

		bx = 2*self.beta0*x

		gg = self.gmk1/2/self.beta0
		gg = gg*log(1-bx)

		return gg

	def g3K(self,bt):		

		x=self.xx(bt)
		#print('valore x= ' + str(x))

		bx = 2*self.beta0*x

		gg = x**2/(1-bx)
		gg = gg*(self.gmk1*self.beta1/self.beta0 - self.gmk2)

		return gg



