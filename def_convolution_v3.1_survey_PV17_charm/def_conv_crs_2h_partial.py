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
   	
		self.g_k = 'PV17'
		self.bmax = 0.6
	

		self.mdl_den = 'gauss'
		self.mdl_num = 'gauss'
	
		self.charm ='no'
		self.nf=3	
        	
        	
	def denominator(self,had1,had2,z1,z2):

		fnt = cr_sec()
		fnt.mass = self.mass        
		fnt.frag2 = self.frag2
		fnt.qq=self.scale
		fnt.charm=self.charm


		mdl1 = model_bt()
		mdl1.mass = 1.115
		mdl1.qq = self.scale
		mdl2 = model_bt()
		mdl2.mass = 0.
		mdl2.qq = self.scale
		
		wdt1_unp = 0.2
		wdt2_unp = 0.2

		scl = Soft(1)
		scl.bmax = self.bmax 
		scl.scale = self.scale
		scl.nf=self.nf

		def fnc(btt):

			res = fnt.cross_sec2(had1,had2,z1,z2,scl.mu_b(btt))
			if self.mdl_den == 'gauss' : res = res*mdl1.MD_gauss(btt,z1,wdt1_unp)
			elif self.mdl_den == 'pwr_lw' : res = res*mdl1.MD(btt,2.,1.)
			elif self.mdl_den == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,2.,z1,1.)
			elif self.mdl_den == 'pwr_lw_ratio' : res = res*mdl1.MD_bstar_ratio(btt,2.,z1,1.)
			#res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*mdl2.MD_pv17(btt,z2) #PV17
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

			### PV17
			elif self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)

			
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
		fnt.qq=self.scale


		mdl1 = model_bt()
		mdl1.mass = 1.115
		mdl1.qq = self.scale
		mdl2 = model_bt()
		mdl2.mass = 0.
		mdl2.qq = self.scale
		
		wdt1_unp = 0.2
		wdt2_unp = 0.2

		scl = Soft(1)
		scl.bmax = self.bmax 
		scl.scale = self.scale
		qT_max = self.scale*self.coef
		scl.nf=self.nf

		def fnc1(btt):

			tot, UP, DO, ST, UPb, DOb, STb = fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			res=btt*tot

			UP=btt*UP
			DO=btt*DO
			ST=btt*ST
			UPb=btt*UPb
			DOb=btt*DOb
			STb=btt*STb

			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			res = res*mdl2.MD_pv17(btt,z2) #PV17
			res = res*special.struve(0,btt*qT_max)
			### PV17
			if self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)


			UP=btt*UP*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_pv17(btt,z2)*special.struve(0,btt*qT_max)*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)
			DO=btt*DO*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_pv17(btt,z2)*special.struve(0,btt*qT_max)*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)
			ST=btt*ST*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_pv17(btt,z2)*special.struve(0,btt*qT_max)*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)
			UPb=btt*UPb*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_pv17(btt,z2)*special.struve(0,btt*qT_max)*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)
			DOb=btt*DOb*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_pv17(btt,z2)*special.struve(0,btt*qT_max)*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)
			STb=btt*STb*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_pv17(btt,z2)*special.struve(0,btt*qT_max)*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)

			#print(res)
			return res, UP, DO, ST, UPb, DOb, STb 		

		def fnc2(btt):

			tot, UP, DO, ST, UPb, DOb, STb = fnt.cross_sec2_polda(had1,had2,z1,z2,scl.mu_b(btt),param)
			res=btt*tot
			
			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,wdt_pol)# h1 model
			res = res*mdl2.MD_pv17(btt,z2) #PV17
			res = res*special.struve(1,btt*qT_max)
			### PV17
			if self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)

			UP=btt*UP*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_pv17(btt,z2)*special.struve(1,btt*qT_max)*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)
			DO=btt*DO*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_pv17(btt,z2)*special.struve(1,btt*qT_max)*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)
			ST=btt*ST*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_pv17(btt,z2)*special.struve(1,btt*qT_max)*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)
			UPb=btt*UPb*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_pv17(btt,z2)*special.struve(1,btt*qT_max)*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)
			DOb=btt*DOb*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_pv17(btt,z2)*special.struve(1,btt*qT_max)*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)
			STb=btt*STb*mdl1.MD_gauss(btt,z1,wdt_pol)*mdl2.MD_pv17(btt,z2)*special.struve(1,btt*qT_max)*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,had2,z1,z2)

			
			return res, UP, DO, ST, UPb, DOb, STb 		

			
		test1 = lambda bt : fnc1(bt)[0]			
		test2 = lambda bt : fnc2(bt)[0]


		
		#print(qT_max)
		N=10
		fbt1 = FBT(1)
		wfbt1 = fbt1.fbt(test1,qT_max,N)
		wfbt1 = wfbt1*pi**2*qT_max



		fbt2 = FBT(0)
		wfbt2 = fbt2.fbt(test2,qT_max,N)
		wfbt2 = wfbt2*pi**2*qT_max	

		wfbt_pol = wfbt1 - wfbt2
###########################################################################
##### partial

		test1_UP = lambda bt : fnc1(bt)[1]			
		test2_UP = lambda bt : fnc2(bt)[1]

		test1_DO = lambda bt : fnc1(bt)[2]			
		test2_DO = lambda bt : fnc2(bt)[2]

		test1_ST = lambda bt : fnc1(bt)[3]			
		test2_ST = lambda bt : fnc2(bt)[3]

		test1_UPb = lambda bt : fnc1(bt)[4]			
		test2_UPb = lambda bt : fnc2(bt)[4]

		test1_DOb = lambda bt : fnc1(bt)[5]			
		test2_DOb = lambda bt : fnc2(bt)[5]

		test1_STb = lambda bt : fnc1(bt)[6]			
		test2_STb = lambda bt : fnc2(bt)[6]

		wfbt1_UP = fbt1.fbt(test1_UP,qT_max,N)
		wfbt1_UP = wfbt1_UP*pi**2*qT_max
		wfbt1_DO = fbt1.fbt(test1_DO,qT_max,N)
		wfbt1_DO = wfbt1_DO*pi**2*qT_max
		wfbt1_ST = fbt1.fbt(test1_ST,qT_max,N)
		wfbt1_ST = wfbt1_ST*pi**2*qT_max

		wfbt1_UPb = fbt1.fbt(test1_UPb,qT_max,N)
		wfbt1_UPb = wfbt1_UPb*pi**2*qT_max
		wfbt1_DOb = fbt1.fbt(test1_DOb,qT_max,N)
		wfbt1_DOb = wfbt1_DOb*pi**2*qT_max
		wfbt1_STb = fbt1.fbt(test1_STb,qT_max,N)
		wfbt1_STb = wfbt1_STb*pi**2*qT_max

		wfbt2_UP = fbt2.fbt(test2_UP,qT_max,N)
		wfbt2_UP = wfbt2_UP*pi**2*qT_max
		wfbt2_DO = fbt2.fbt(test2_DO,qT_max,N)
		wfbt2_DO = wfbt2_DO*pi**2*qT_max
		wfbt2_ST = fbt2.fbt(test2_ST,qT_max,N)
		wfbt2_ST = wfbt2_ST*pi**2*qT_max

		wfbt2_UPb = fbt2.fbt(test2_UPb,qT_max,N)
		wfbt2_UPb = wfbt2_UPb*pi**2*qT_max
		wfbt2_DOb = fbt2.fbt(test2_DOb,qT_max,N)
		wfbt2_DOb = wfbt2_DOb*pi**2*qT_max
		wfbt2_STb = fbt2.fbt(test2_STb,qT_max,N)
		wfbt2_STb = wfbt2_STb*pi**2*qT_max

		wfbt_pol_UP = wfbt1_UP - wfbt2_UP
		wfbt_pol_DO = wfbt1_DO - wfbt2_DO
		wfbt_pol_ST = wfbt1_ST - wfbt2_ST


		wfbt_pol_UPb = wfbt1_UPb - wfbt2_UPb
		wfbt_pol_DOb = wfbt1_DOb - wfbt2_DOb
		wfbt_pol_STb = wfbt1_STb - wfbt2_STb


		#wfbt_pol_UP =0# wfbt1_UP - wfbt2_UP
		#wfbt_pol_DO =0# wfbt1_DO - wfbt2_DO
		#wfbt_pol_ST =0# wfbt1_ST - wfbt2_ST
		#wfbt_pol_UPb =0# wfbt1_UPb - wfbt2_UPb
		#wfbt_pol_DOb =0# wfbt1_DOb - wfbt2_DOb
		#wfbt_pol_STb =0# wfbt1_STb - wfbt2_STb



		return wfbt_pol, wfbt_pol_UP, wfbt_pol_DO, wfbt_pol_ST, wfbt_pol_UPb, wfbt_pol_DOb, wfbt_pol_STb	


	def ratio(self,had1,had2,z1,z2,param,wdt_pol,mss):
	
		mass = 1.115	
		num, UP, DO, ST, UPb, DOb, STb  = self.numerator(had1,had2,z1,z2,param,wdt_pol,mss)
		den = self.denominator(had1,had2,z1,z2)

		#fnt = cr_sec()
		#fnt.mass = self.mass  
		#fact = 1#fnt.fact_fst_mom(z1,0.2,wdt_pol)		
		#if self.mdl_num != 'gauss': fact=1.
		
#		result = fact*mass*num/den
		result = mass*num/den
		
		return result, mass*UP/den, mass*DO/den, mass*ST/den, mass*UPb/den, mass*DOb/den, mass*STb/den 










	
	
	
	
	