# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
import lhapdf
from andrea_model.Sudakov_and import* 
from def_crs import*
from model_fct import*

class polarization_1h_and :

	def __init__(self):

		self.mass = 1.115

		self.frag2= 'dss'

		self.scale = 10.58

		#self.coef = coefficient   


	def denominator(self,had1,z1,pt):

		fnt = cr_sec()
		fnt.mass = self.mass  

		mdl1 = model_bt()
		mdl1.mass = 1.115

		wdt1_unp = 0.2	

		scl = exp_fnc()
		gg = g_funct()
		#scl.scale = self.scale
		
		def fnc(btt):
			pwr=2
			mm=1
			res = btt*fnt.cross_sec1(had1,z1,gg.mu_b(btt))
			res = res*mdl1.MD(btt,pwr,mm)
			

			#res = res*scl.soft_pert_1h(btt)*scl.g_K(btt)
			res = res*scl.espn1(bt,QQ)*scl.espn2(bt,QQ)*scl.espn3(bt,z1,QQ)

			return res		


		test = lambda bt : fnc(bt)	
		
		ml = 0.
		eta_p= (1 - 4*ml**2/((z1**2)*self.scale**2))
		
		zp1 = z1*sqrt(eta_p)		# momentum fraction
	
		
		qts = pt/zp1

		N=10
		fbt = FBT(0)
		wfbt_unp = fbt.fbt(test,qts,N)
		#wfbt_unp = wfbt_unp*2*pi*qT_max	
		
		return wfbt_unp 


	def numerator(self,had1,z1,pt,param,pwr,mm):


		fnt = cr_sec()
		fnt.mass = self.mass   

		mdl1 = model_bt()
		mdl1.mass = 1.115

		scl = exp_fnc()
		scl.TT = 0.875
		gg = g_funct()

		def fnc(btt):

			res = btt**2*fnt.cross_sec1_polda_fm(had1,z1,gg.mu_b(btt),param)
			#print(fnt.cross_sec1_polda(had1,z1,scl.mu_b(btt),param))
			res = res*mdl1.MD(btt,pwr,mm)
			

			#res = res*scl.soft_pert_1h(btt)*scl.g_K(btt)
			res = res*scl.espn1(bt,QQ)*scl.espn2(bt,QQ)*scl.espn3(bt,z1,QQ)

			return res		

		test = lambda bt: fnc(bt)

		ml = 0.
		eta_p= (1 - 4*ml**2/((z1**2)*self.scale**2))
		#print('eta  ' + str(eta_p))
		zp1 = z1*sqrt(eta_p)		# momentum fraction
		#print('zp  ' + str(zp1))
		
		qts = pt/zp1


		N=10
		fbt1 = FBT(1)
		wfbt1 = fbt1.fbt(test,qts,N)
		#print(wfbt1)
		wfbt1 = wfbt1*self.mass
	
		return wfbt1 


	def ratio(self,had1,z1,pt,param,pwr,mm):
	
		
		fnt = cr_sec()
		fnt.mass = self.mass  
		#fact = fnt.fact_fst_mom(z1,0.2,wdt_pol)		
		
		num = self.numerator(had1,z1,pt,param,pwr,mm)
		den = self.denominator(had1,z1,pt)	
		
		result = num/den
		
		return result








