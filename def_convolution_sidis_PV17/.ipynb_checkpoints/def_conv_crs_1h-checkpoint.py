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


class polarization_1h :

	def __init__(self):

		self.mass = 1.115

		self.frag2= 'dss'

		self.scale = 10.58

		#self.coef = coefficient   

		self.g_k = 'exp'
		self.bmax = 0.8

		self.mdl_den = 'gauss'
		self.mdl_num = 'gauss'

	def denominator(self,had1,z1,pt):

		fnt = cr_sec()
		fnt.mass = self.mass  

		mdl1 = model_bt()
		mdl1.mass = 1.115

		wdt1_unp = 0.2	

		scl = Soft(1)
		scl.scale = self.scale
		scl.bmax = self.bmax 
		
		def fnc(btt):

			res = btt*fnt.cross_sec1(had1,z1,scl.mu_b(btt))
			#res = res*mdl1.MD_gauss(btt,z1,wdt1_unp)

			if self.mdl_den == 'gauss' : res = res*mdl1.MD_gauss(btt,z1,wdt1_unp)
			elif self.mdl_den == 'pwr_lw' : res = res*mdl1.MD(btt,2.,1.)
			
			#res = res*scl.soft_pert_1h(btt)*scl.g_K(btt)
			if self.g_k == 'log_b': res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_1h(btt)
			elif self.g_k == 'blny':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_blny_1h(btt)			
			elif self.g_k == 'new':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_new_1h(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_1h(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_exp_1h(btt)
			###
			#
			elif self.g_k == 'log_b_lgm':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_1h_lgm(btt,z1)		
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_1h_lgm_n(btt,z1)		
			#
			###
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_1h_lgm(btt,z1)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_1h_lgm(btt,z1)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_1h_lgm(btt,z1)

			return res		


		test = lambda bt : fnc(bt)	
		
		ml = self.mass
		eta_p= (1 - 4*ml**2/((z1**2)*self.scale**2))
		
		zp1 = z1*sqrt(eta_p)		# momentum fraction
	
		
		qts = pt/zp1

		N=20
		fbt = FBT(0)
		wfbt_unp = fbt.fbt(test,qts,N)
		#wfbt_unp = wfbt_unp*2*pi*qT_max	
		
		#print('den = '+str(wfbt_unp))
		return wfbt_unp
		
		
	def numerator(self,had1,z1,pt,param,wdt_pol,mss):


		fnt = cr_sec()
		fnt.mass = self.mass   

		mdl1 = model_bt()
		mdl1.mass = 1.115

		scl = Soft(1)
		scl.scale = self.scale
		scl.bmax = self.bmax 


		def fnc(btt):

			res = btt**2*fnt.cross_sec1_polda(had1,z1,scl.mu_b(btt),param)
			#print(fnt.cross_sec1_polda(had1,z1,scl.mu_b(btt),param))
			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw' : res = res*mdl1.MD(btt,wdt_pol,mss)
			elif self.mdl_num == 'pwr_lw_pt' : res = res*mdl1.MD_pt(btt,wdt_pol,z1)
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,mss)


			
			#res = res*scl.soft_pert_1h(btt)*scl.g_K(btt)
			if self.g_k == 'log_b': res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_1h(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_blny_1h(btt)
			elif self.g_k == 'new': res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_new_1h(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_1h(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_exp_1h(btt)
			###
			#
			elif self.g_k == 'log_b_lgm':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_1h_lgm(btt,z1)		
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_1h_lgm_n(btt,z1)		
			#
			###
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_1h_lgm(btt,z1)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_1h_lgm(btt,z1)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_1h_lgm(btt,z1)


			return res		

		test = lambda bt: fnc(bt)

		ml = self.mass
		eta_p= (1 - 4*ml**2/((z1**2)*self.scale**2))
		#print('eta  ' + str(eta_p))
		zp1 = z1*sqrt(eta_p)		# momentum fraction
		#print('zp  ' + str(zp1))
		
		qts = pt/zp1


		N=20
		fbt1 = FBT(1)
		wfbt1 = fbt1.fbt(test,qts,N)
		#print(wfbt1)
		wfbt1 = wfbt1*self.mass
	
		#wfbt1=round(wfbt1,6)
		#print('--->> num = '+str(wfbt1))
		return wfbt1 
	
	def ratio(self,had1,z1,pt,param,wdt_pol,mss):
	
		
		fnt = cr_sec()
		fnt.mass = self.mass  
		fact =1# fnt.fact_fst_mom(z1,0.2,wdt_pol)		
		if self.mdl_num != 'gauss': fact=1.

		num = self.numerator(had1,z1,pt,param,wdt_pol,mss)
		den = self.denominator(had1,z1,pt)	
		
		result = fact*num/den
		#print('res = '+str(result))
		return result
		
		
		
		
		
		
			
