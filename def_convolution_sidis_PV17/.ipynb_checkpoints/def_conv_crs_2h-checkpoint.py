# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import*
from scipy import special
from def_crs import*
from model_fct import*
from g_funct import*
from FBT import FBT 
from model_fct import*
from Sudakov.sudakov_factor import*
from Sudakov.evolve import evolve, sng


class polarization :
    
	def __init__(self,coefficient):

		self.mass = 1.115

		self.frag2= 'dss'

		self.scale = 10.58

		self.coef = coefficient   
   	
		self.g_k = 'exp'
		self.bmax = 0.8
        	
	def denominator(self,had1,had2,z1,z2):

		fnt = cr_sec()
		fnt.mass = self.mass        
		fnt.frag2 = self.frag2

		mdl1 = model_bt()
		mdl1.mass = 1.115
		mdl2 = model_bt()
		mdl2.mass = 0.

		wdt1_unp = 0.2
		wdt2_unp = 0.2

		scl = Soft(1)
		scl.bmax = self.bmax 

		def fnc(btt):

			res = fnt.cross_sec2(had1,had2,z1,z2,scl.mu_b(btt))
			res = res*mdl1.MD_gauss(btt,z1,wdt1_unp)
			res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			

			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'old':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt)
	
			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			
			return res		
			
		test = lambda bt : fnc(bt)			


		qT_max = self.scale*self.coef

		N=10
		fbt = FBT(1)
		wfbt_unp = fbt.fbt(test,qT_max,N)
		wfbt_unp = wfbt_unp*2*pi*qT_max	

		return wfbt_unp


	def numerator(self,had1,had2,z1,z2,param,wdt_pol):

		fnt = cr_sec()
		fnt.mass = self.mass        
		fnt.frag2 = self.frag2

		mdl1 = model_bt()
		mdl1.mass = 1.115
		mdl2 = model_bt()
		mdl2.mass = 0.

		wdt1_unp = 0.2
		wdt2_unp = 0.2

		scl = Soft(1)
		scl.bmax = self.bmax 
		qT_max = self.scale*self.coef

		def fnc1(btt):

			res = btt*fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			res = res*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*special.struve(0,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'old':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt)

			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)

			return res		

		def fnc2(btt):

			res = btt*fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			res = res*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*special.struve(1,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'old':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt)

			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			
			return res		

			
		test1 = lambda bt : fnc1(bt)			
		test2 = lambda bt : fnc2(bt)

		

		N=10
		fbt1 = FBT(1)
		wfbt1 = fbt1.fbt(test1,qT_max,N)
		wfbt1 = wfbt1*pi**2*qT_max

		fbt2 = FBT(0)
		wfbt2 = fbt2.fbt(test2,qT_max,N)
		wfbt2 = wfbt2*pi**2*qT_max	

		wfbt_pol = wfbt1 - wfbt2


		return wfbt_pol	


	def ratio(self,had1,had2,z1,z2,param,wdt_pol):
	
		mass = 1.115	
		num = self.numerator(had1,had2,z1,z2,param,wdt_pol)
		den = self.denominator(had1,had2,z1,z2)
		
		fnt = cr_sec()
		fnt.mass = self.mass  
		fact = fnt.fact_fst_mom(z1,0.2,wdt_pol)		
		
		result = fact*mass*num/den
		
		return result










	
	
	
	
	
