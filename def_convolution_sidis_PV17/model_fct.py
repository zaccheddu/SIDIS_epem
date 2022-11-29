# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
from scipy import special
import numpy as np

class model_bt():
		
	def __init__(self):

		self.m = 1. # propagator mass GeV
		self.p= 2.

		self.unp_width = 0.2 ## ampiezza gaussiana non polarizzata

		self.mass = 0.

		self.frag2= 'dss'	
	
		self.qq = 10.58

		self.bmin = 2*e**(-euler_gamma)/self.qq 	### for mu_b_star
		self.bmax = 0.6



	def bt_str_mu(self,bt): 				 #define bt_star for mu_b_star with bmin
	
		bm= self.bmax	
		
		b_new = sqrt(bt**2 + self.bmin**2)

		bstar = b_new/sqrt(1+b_new**2/bm**2)			
		#print('in bt_star:' + str(bt))		
		return bstar

	def bt_c(self,bt): 
	
		b_new = sqrt(bt**2 + self.bmin**2)
		return b_new

	def bt_str(self,bt): 				 
	
		bm= self.bmax	
		
		bstar = bt/sqrt(1+bt**2/bm**2)			
		#print('in bt_star:' + str(bt))		
		return bstar


###############################################################################################
#
#		MODELS FOR G_J
#_______________________________________________________________________________


		#return mm

	def MD_bstar(self,bt,pwr,zz,xb2,mass):  ## power-law non polarizzata

		b_new = sqrt(bt**2 + self.bmin**2)
		#b_new_star = self.bt_str(b_new)
#
		#b_new = bt
		ml= 1.115


		eta_p= (1 - (ml**2/zz**2/self.qq**2)*xb2/(1-xb2))
		zp1 = zz*sqrt(eta_p)		# momentum fraction


		###############
		#pt0=mass1
		#mm2=(2*pwr*(1+pt0**2)/(1+5*pt0**2) -1)*pt0**2
		#mass=sqrt(mm2)
		#############


		b_new = b_new/zp1

		gmm = 2**(2-pwr)/special.gamma(pwr -1)
		
		bssl1 = special.kv(pwr-1,b_new*mass)
		mm1 = gmm*bssl1*(b_new*mass)**(pwr-1)		 # modified Bessel function of the second kind

		#bssl2 = special.kn(pwr-1,b_new_star*mass)
		#mm2 = gmm*bssl2*(b_new_star*mass)**(pwr-1)		 # modified Bessel function of the second kind
		

		
		return mm1 #/mm2




	def MD_gauss(self,bt,zz,xb2,unp_width):  ## gauss non polarizzata

		b_new = sqrt(bt**2 + self.bmin**2)
		#print(self.bmin)
		#b_new=bt
		ml= self.mass


		eta_p= (1 - (ml**2/zz**2/self.qq**2)*xb2/(1-xb2))
		zp1 = zz*sqrt(eta_p)		# momentum fraction


		wdt=   unp_width
		#print(wdt)
		esp= - wdt*b_new**2/4/zp1**2
		
#		esp= - wdt*b_new*self.bt_str(b_new)/4/zp1**2
#		esp= - wdt*self.bt_str(b_new)**2/4/zp1**2

		mm =   exp(esp)	
		#mm = mm/pi/sqrt(wdt*2)	

		return mm


		


	
	def MD_pv17(self,bt,zz):
	
		b_new = sqrt(bt**2 + self.bmin**2)
		#b_new = self.bt_str_mu(bt)
		#b_new = bt
		
		g3 = 0.21  #GeV**2
		g4 = 0.13  #GeV**2
		z_hat = 0.5
		
		beta = 1.65
		delta = 2.28
		gamma = 0.14
		lambda_f = 5.5	#GeV**-2	


		#ml= self.mass
		#eta_p= (1 - 4*ml**2/zz**2/self.qq**2)
		#zp1 = zz*sqrt(eta_p)		# momentum fraction

		
		
		gg3 = g3*(zz**beta + delta)*(1-zz)**gamma/(z_hat**beta + delta)/(1-z_hat)**gamma
		gg4 = g4*(zz**beta + delta)*(1-zz)**gamma/(z_hat**beta + delta)/(1-z_hat)**gamma

		num = gg3*exp(-b_new**2*gg3/4/zz**2) 
		num = num + (lambda_f/zz**2)*gg4**2*(1- b_new**2*gg4/4/zz**2)*exp(-b_new**2*gg4/4/zz**2)

		den = (gg3 +(lambda_f/zz**2)*gg4**2 )#*(2*pi*zz**2)
		
		return num/den

	def Mf_pv17(self,bt,xx):

		b_new = sqrt(bt**2 + self.bmin**2)

		x_hat = 0.1
		g1 = 0.28
		alpha = 2.95
		sigma = 0.17
		lambd = 0.86
		
		gg1 = g1*(1.-xx)**alpha*xx**sigma
		gg1 = gg1/(1.-x_hat)**alpha/x_hat**sigma
		
		fact = exp(-gg1*b_new**2/4)
		fact = fact*(1. - b_new**2*(lambd*gg1**2)/(1 + lambd*gg1)/4)
		
		return fact


















