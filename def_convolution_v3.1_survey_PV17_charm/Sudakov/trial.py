# -*- coding: utf-8 -*-
from pylab import*
import numpy as np
import lhapdf
import scipy.integrate as integrate

class Soft:

	CA= 3.
	CF= 4/3
	nf= 3	#active flavors
	Tr = 1/2	
	

	pdf=lhapdf.mkAlphaS("NNFF10_PIp_nlo")
	

	def __init__(self,order):
		
		self.bmax = 1.
		self.scale= 10.58

		self.a = 0.05 ## GeV**2

		self.order = order		

		self.bmin = 2*e**(-euler_gamma)/10.58		
	
		self.bmax = .5
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


	def integrand(self,mu):
	
		if self.order == 1:
			integr = self.pdf.alphasQ(mu)/np.pi
			integr = integr*self.CF*(2*np.log(self.scale/mu) - 3/2)

		elif self.order == 2:
			integr1 = self.pdf.alphasQ(mu)/np.pi
			integr1 = integr1*self.CF*(2*np.log(self.scale/mu) - 3/2)

			integr2 = (self.pdf.alphasQ(mu)/np.pi)**2
			integr2 = integr2*self.CF*2*np.log(self.scale/mu)
			integr2 = integr2*(self.CA*(67/18 - np.pi**2/6) -self.Tr*self.nf*10/9 )
			
			integr = integr1 + integr2
		
		return integr		
		

	def soft_pert(self,bt):
	
		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			test = lambda xx: self.integrand(xx)/xx	
			
			mub_str = self.mu_b(bb)
			
			#print(mub_str)
			#print(self.scale)
			
			result = integrate.quad(test,mub_str,self.scale)
			out = np.append(out,np.exp(result[0]))
		
		
		return out

	def g_K(self,bt):
	
		g2 = 0.84
		Q_0 = np.sqrt(2.4) 
		
		esp =  g2*np.log(self.scale/Q_0)*np.log(bt/self.bt_str(bt))	
		
		return np.exp(-esp)

	def g_K_1h(self,bt):
	
		g2 = 0.84
		Q_0 = np.sqrt(2.4) 
		
		esp =  g2*np.log(self.scale/Q_0)*np.log(bt/self.bt_str(bt))/2	
		
		return np.exp(-esp)

		
'''		
	def soft_pert(self,bt):
	
		mub_str = self.mu_b(bt)
	
		xx = np.log(self.scale/mub_str)
	
		if self.order == 1:
			integr = self.pdf.alphasQ(mu)/np.pi
			integr = integr*self.CF*( xx**2 - 3/2*xx)

		elif self.order == 2:
			integr1 = self.pdf.alphasQ(mu)/np.pi
			integr1 = integr1*self.CF*( xx**2 - 3/2*xx)

			integr2 = (self.pdf.alphasQ(mu)/np.pi)**2
			integr2 = integr2*self.CF*xx**2
			integr2 = integr2*(self.CA*(67/18 - np.pi**2/6) -self.Tr*self.nf*10/9 )
			
			integr = integr1 + integr2
		
		return integr		
'''		
		
		
		
		
		
		
		
		
		
		
		
