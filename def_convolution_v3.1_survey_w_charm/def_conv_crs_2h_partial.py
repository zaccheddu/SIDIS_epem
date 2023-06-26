# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import*
from scipy import special
from def_crs_partial import*
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
		self.bmax = 0.6
	

		self.mdl_den = 'gauss'
		self.mdl_num = 'gauss'

		self.su2 = 'no'
        	        	
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
			if self.mdl_den == 'gauss' : res = res*mdl1.MD_gauss(btt,z1,wdt1_unp)
			elif self.mdl_den == 'pwr_lw' : res = res*mdl1.MD(btt,2.,1.)
			elif self.mdl_den == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,2.,z1,1.)
			elif self.mdl_den == 'pwr_lw_ratio' : res = res*mdl1.MD_bstar_ratio(btt,2.,z1,1.)
			res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			#res = res*mdl2.MD(btt,2.,1.)
			

			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt,had2,z1,z2)
	
			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			####
			#
			elif self.g_k == 'log_b_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_lgm_n(btt,had2,z1,z2)		
			#
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)


			
			return res		
			
		test = lambda bt : fnc(bt)			


		qT_max = self.scale*self.coef

		N=10
		fbt = FBT(1)
		wfbt_unp = fbt.fbt(test,qT_max,N)
		wfbt_unp = wfbt_unp*2*pi*qT_max	

		return wfbt_unp


	def numerator(self,had1,had2,z1,z2,param,wdt_pol,mss):

		fnt = cr_sec()
		fnt.mass = self.mass        
		fnt.frag2 = self.frag2
		fnt.su2   = self.su2 

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

			#cross_tot,ff1,ff2,ff3,ff4 = fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			#res = btt*cross_tot
			
			res = btt*fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)[0]
			#print(res)
			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw' : res = res*mdl1.MD(btt,wdt_pol,mss)
			elif self.mdl_num == 'pwr_lw_pt' : res = res*mdl1.MD_pt(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_ratio' : res = res*mdl1.MD_bstar_ratio(btt,wdt_pol,z1,mss)
			
			res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			#res = res*mdl2.MD(btt,2.,1.)

			res = res*special.struve(0,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt,had2,z1,z2)

			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			####
			#
			elif self.g_k == 'log_b_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_lgm_n(btt,had2,z1,z2)		
			#
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)


			return res		

		def fnc2(btt):

			#cross_tot,ff1,ff2,ff3,ff4 = fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			#res = btt*cross_tot

			res = btt*fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)[0]

			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw' : res = res*mdl1.MD(btt,wdt_pol,mss)
			elif self.mdl_num == 'pwr_lw_pt' : res = res*mdl1.MD_pt(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_ratio' : res = res*mdl1.MD_bstar_ratio(btt,wdt_pol,z1,mss)
			
			res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			#res = res*mdl2.MD(btt,2.,1.)
		
			res = res*special.struve(1,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt,had2,z1,z2)

			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			####
			#
			elif self.g_k == 'log_b_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_lgm_n(btt,had2,z1,z2)		
			#
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)
			
			
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


	def ff1(self,had1,had2,z1,z2,param,wdt_pol,mss):

		fnt = cr_sec()
		fnt.mass = self.mass        
		fnt.frag2 = self.frag2
		fnt.su2   = self.su2 

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

			#cross_tot,ff1,ff2,ff3,ff4 = fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			#res = btt*cross_tot
			
			res = btt*fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)[1]
			#print(res)
			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw' : res = res*mdl1.MD(btt,wdt_pol,mss)
			elif self.mdl_num == 'pwr_lw_pt' : res = res*mdl1.MD_pt(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_ratio' : res = res*mdl1.MD_bstar_ratio(btt,wdt_pol,z1,mss)
			
			res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			#res = res*mdl2.MD(btt,2.,1.)

			res = res*special.struve(0,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt,had2,z1,z2)

			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			####
			#
			elif self.g_k == 'log_b_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_lgm_n(btt,had2,z1,z2)		
			#
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)


			return res		

		def fnc2(btt):

			#cross_tot,ff1,ff2,ff3,ff4 = fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			#res = btt*cross_tot

			res = btt*fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)[1]

			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw' : res = res*mdl1.MD(btt,wdt_pol,mss)
			elif self.mdl_num == 'pwr_lw_pt' : res = res*mdl1.MD_pt(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_ratio' : res = res*mdl1.MD_bstar_ratio(btt,wdt_pol,z1,mss)
			
			res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			#res = res*mdl2.MD(btt,2.,1.)
		
			res = res*special.struve(1,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt,had2,z1,z2)

			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			####
			#
			elif self.g_k == 'log_b_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_lgm_n(btt,had2,z1,z2)		
			#
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)
			
			
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
		
	def ff2(self,had1,had2,z1,z2,param,wdt_pol,mss):

		fnt = cr_sec()
		fnt.mass = self.mass        
		fnt.frag2 = self.frag2
		fnt.su2   = self.su2 

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

			#cross_tot,ff1,ff2,ff3,ff4 = fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			#res = btt*cross_tot
			
			res = btt*fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)[2]
			#print(res)
			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw' : res = res*mdl1.MD(btt,wdt_pol,mss)
			elif self.mdl_num == 'pwr_lw_pt' : res = res*mdl1.MD_pt(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_ratio' : res = res*mdl1.MD_bstar_ratio(btt,wdt_pol,z1,mss)
			
			res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			#res = res*mdl2.MD(btt,2.,1.)

			res = res*special.struve(0,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt,had2,z1,z2)

			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			####
			#
			elif self.g_k == 'log_b_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_lgm_n(btt,had2,z1,z2)		
			#
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)


			return res		

		def fnc2(btt):

			#cross_tot,ff1,ff2,ff3,ff4 = fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			#res = btt*cross_tot

			res = btt*fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)[2]

			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw' : res = res*mdl1.MD(btt,wdt_pol,mss)
			elif self.mdl_num == 'pwr_lw_pt' : res = res*mdl1.MD_pt(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_ratio' : res = res*mdl1.MD_bstar_ratio(btt,wdt_pol,z1,mss)
			
			res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			#res = res*mdl2.MD(btt,2.,1.)
		
			res = res*special.struve(1,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt,had2,z1,z2)

			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			####
			#
			elif self.g_k == 'log_b_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_lgm_n(btt,had2,z1,z2)		
			#
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)
			
			
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

	def ff3(self,had1,had2,z1,z2,param,wdt_pol,mss):

		fnt = cr_sec()
		fnt.mass = self.mass        
		fnt.frag2 = self.frag2
		fnt.su2   = self.su2 

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

			#cross_tot,ff1,ff2,ff3,ff4 = fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			#res = btt*cross_tot
			
			res = btt*fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)[3]
			#print(res)
			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw' : res = res*mdl1.MD(btt,wdt_pol,mss)
			elif self.mdl_num == 'pwr_lw_pt' : res = res*mdl1.MD_pt(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_ratio' : res = res*mdl1.MD_bstar_ratio(btt,wdt_pol,z1,mss)
			
			res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			#res = res*mdl2.MD(btt,2.,1.)

			res = res*special.struve(0,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt,had2,z1,z2)

			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			####
			#
			elif self.g_k == 'log_b_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_lgm_n(btt,had2,z1,z2)		
			#
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)


			return res		

		def fnc2(btt):

			#cross_tot,ff1,ff2,ff3,ff4 = fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			#res = btt*cross_tot

			res = btt*fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)[3]

			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw' : res = res*mdl1.MD(btt,wdt_pol,mss)
			elif self.mdl_num == 'pwr_lw_pt' : res = res*mdl1.MD_pt(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_ratio' : res = res*mdl1.MD_bstar_ratio(btt,wdt_pol,z1,mss)
			
			res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			#res = res*mdl2.MD(btt,2.,1.)
		
			res = res*special.struve(1,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt,had2,z1,z2)

			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			####
			#
			elif self.g_k == 'log_b_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_lgm_n(btt,had2,z1,z2)		
			#
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)
			
			
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
		
		
	def ff4(self,had1,had2,z1,z2,param,wdt_pol,mss):

		fnt = cr_sec()
		fnt.mass = self.mass        
		fnt.frag2 = self.frag2
		fnt.su2   = self.su2 

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

			#cross_tot,ff1,ff2,ff3,ff4 = fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			#res = btt*cross_tot
			
			res = btt*fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)[4]
			#print(res)
			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw' : res = res*mdl1.MD(btt,wdt_pol,mss)
			elif self.mdl_num == 'pwr_lw_pt' : res = res*mdl1.MD_pt(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_ratio' : res = res*mdl1.MD_bstar_ratio(btt,wdt_pol,z1,mss)
			
			res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			#res = res*mdl2.MD(btt,2.,1.)

			res = res*special.struve(0,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt,had2,z1,z2)

			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			####
			#
			elif self.g_k == 'log_b_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_lgm_n(btt,had2,z1,z2)		
			#
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)


			return res		

		def fnc2(btt):

			#cross_tot,ff1,ff2,ff3,ff4 = fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			#res = btt*cross_tot

			res = btt*fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)[4]

			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw' : res = res*mdl1.MD(btt,wdt_pol,mss)
			elif self.mdl_num == 'pwr_lw_pt' : res = res*mdl1.MD_pt(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,mss)
			elif self.mdl_num == 'pwr_lw_ratio' : res = res*mdl1.MD_bstar_ratio(btt,wdt_pol,z1,mss)
			
			res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			#res = res*mdl2.MD(btt,2.,1.)
		
			res = res*special.struve(1,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			elif self.g_k == 'blny': res = res*scl.analytic_sudakov(btt)*scl.g_K_blny(btt,had2,z1,z2)

			elif self.g_k == 'new': res = res*scl.analytic_sudakov(btt)*scl.g_K_new(btt)
			elif self.g_k == 'll': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll(btt)
			elif self.g_k == 'exp': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp(btt)
			####
			#
			elif self.g_k == 'log_b_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm_n':res = res*scl.analytic_sudakov_1h(btt)*scl.g_K_ll_lgm_n(btt,had2,z1,z2)		
			#
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)
			
			
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

	def ratio(self,had1,had2,z1,z2,param,wdt_pol,mss):
	
		mass = 1.115	
		num = self.numerator(had1,had2,z1,z2,param,wdt_pol,mss)
		den = self.denominator(had1,had2,z1,z2)

		ffs1 = self.ff1(had1,had2,z1,z2,param,wdt_pol,mss)
		ffs2 = self.ff2(had1,had2,z1,z2,param,wdt_pol,mss)
		ffs3 = self.ff3(had1,had2,z1,z2,param,wdt_pol,mss)
		ffs4 = self.ff4(had1,had2,z1,z2,param,wdt_pol,mss)

		#summ= ffs1+ ffs2+ ffs3+ ffs4
		#summ=around(summ/den,6)
		#num=around(num/den,6)
		#print('sums =' + str(summ) + '   numer =' + str(num))
		#print(num==summ)
		

		#fnt = cr_sec()
		#fnt.mass = self.mass  
		#fact = 1#fnt.fact_fst_mom(z1,0.2,wdt_pol)		
		#if self.mdl_num != 'gauss': fact=1.
		
#		result = fact*mass*num/den
		result = mass*num/den
		
		return result, mass*ffs1/den, mass*ffs2/den, mass*ffs3/den ,mass*ffs4/den










	
	
	
	
	
