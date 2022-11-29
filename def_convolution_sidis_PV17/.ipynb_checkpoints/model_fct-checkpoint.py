# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
from scipy import special


class model_bt():
		
	def __init__(self):

		self.m = 1. # propagator mass GeV
		self.p= 2.

		self.unp_width = 0.2 ## ampiezza gaussiana non polarizzata

		self.mass = 0.

		self.frag2= 'dss'	
	
		self.qq = 10.58

		self.bmin = 2*e**(-euler_gamma)/10.58	### for mu_b_star


	def MD(self,bt,pwr,mass):  ## power-law non polarizzata

		b_new = sqrt(bt**2 + self.bmin**2)
		
		gmm = 2**(2-pwr)/special.gamma(pwr -1)
		
		bssl = special.kn(pwr-1,b_new*mass)		 # modified Bessel function of the second kind

		mm = gmm*bssl*(b_new*mass)**(pwr-1)

		return mm


	def MD_gauss(self,bt,zz,unp_width):  ## gauss non polarizzata

		b_new = sqrt(bt**2 + self.bmin**2)

		ml= self.mass

		

		eta_p= (1 - 4*ml**2/zz**2/self.qq**2)
		zp1 = zz*sqrt(eta_p)		# momentum fraction


		wdt=   unp_width
		esp= - wdt*b_new**2/4/zp1**2

		mm =   exp(esp)	
		#mm = mm/pi/sqrt(wdt*2)	

		return mm

	def MD_log(self,bt,zz,width):
	
		b_new = sqrt(bt**2 + self.bmin**2)

		ml= self.mass

		

		eta_p= (1 - 4*ml**2/zz**2/self.qq**2)
		zp1 = zz*sqrt(eta_p)		# momentum fraction


		wdt=   unp_width
		esp= - b_new**2/4/zp1**2
		esp= esp*np.log(1/b_new/width)
		
		mm =   exp(esp)	
		#mm = mm/pi/sqrt(wdt*2)	

		return mm
		
		


	def MD_pol(self,bt,power,mass):  # power-law polarizzata

		#pt_0 = 0.329  # nel caso volessi mettere in relazione il parametro di massa con la posizione del picco

		#mass2 = (2*power -1)*pt_0**2
		#mass=sqrt(mass2)
			
		
		gmm = 2**(2-power)/special.gamma(power -1)
		
		bssl = special.kn(power-1,bt*mass) # modified Bessel function of the second kind

		mm = gmm*bssl*(bt*mass)**(power-1)

		return mm



	def MD_gauss_pol(self,bt,zz,wdt_pol):  # gauss polarizzata

		
		esp= - wdt_pol*bt**2/4/zz**2

		mm =   exp(esp) #/sqrt(2*wdt_pol)/pi		
		#mm=mm*bt/2/sqrt(2)

		return mm


