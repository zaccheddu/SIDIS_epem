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
import lhapdf

class polarization :
    
	def __init__(self,coefficient):

		self.mass = 1.115

		self.frag2= 'dss'

		self.scale = 10.58

		self.coef = coefficient   
   	
		self.g_k = 'PV17'
		self.bmax = 0.6
	
		self.charm = 'no'

		self.mdl_den = 'gauss'
		self.mdl_num = 'gauss'
        	
	def denominator(self,had1,had2,z1,xb2):

		fnt = cr_sec()
		fnt.mass = self.mass        
		#fnt.frag2 = self.frag2
		fnt.qq=self.scale
		fnt.charm = self.charm

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

		def fnc(btt):

			res = fnt.cross_sec2(had1,had2,z1,xb2,scl.mu_b(btt))
			#print(scl.mu_b(btt))
			if self.mdl_den == 'gauss' : res = res*mdl1.MD_gauss(btt,z1,xb2,wdt1_unp)
			elif self.mdl_den == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,2.,z1,xb2,1.)
			#res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*mdl2.Mf_pv17(btt,xb2) #PV17
			#res = res*mdl2.MD(btt,2.,1.)
			

			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)

			### PV17
			elif self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,z1,xb2)

			
			return res		
			
		test = lambda bt : fnc(bt)			


		

		N=10
		fbt = FBT(1)
		wfbt_unp = fbt.fbt(test,qT_max,N)
		wfbt_unp = wfbt_unp*2*pi*qT_max	

		return wfbt_unp


	def numerator(self,had1,had2,z1,xb2,param,wdt_pol,mss):

		fnt = cr_sec()
		fnt.mass = self.mass        
		#fnt.frag2 = self.frag2
		fnt.qq=self.scale
		fnt.charm = self.charm

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

		def fnc1(btt):

			res = btt*fnt.cross_sec2_polda(had1,had2,z1,xb2,scl.mu_b(btt),param)
			#print(res)
			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,xb2,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,xb2,mss)
			#res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*mdl2.Mf_pv17(btt,xb2) #PV17
			#res = res*mdl2.MD(btt,2.,1.)

			res = res*special.struve(0,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)

			### PV17
			elif self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,z1,xb2)


			return res		

		def fnc2(btt):

			res = btt*fnt.cross_sec2_polda(had1,had2,z1,xb2,scl.mu_b(btt),param)

			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,xb2,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,xb2,mss)
			
			#res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*mdl2.Mf_pv17(btt,xb2)  #PV17
			#res = res*mdl2.MD(btt,2.,1.)
		
			res = res*special.struve(1,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)

			### PV17
			elif self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,z1,xb2)

			
			return res		

			
		test1 = lambda bt : fnc1(bt)			
		test2 = lambda bt : fnc2(bt)
		
		N=10
		
		fbt2 = FBT(0)
		wfbt2 = fbt2.fbt(test2,qT_max,N)
		wfbt2 = wfbt2*pi**2*qT_max	
		

		
		fbt1 = FBT(1)
		wfbt1 = fbt1.fbt(test1,qT_max,N)
		wfbt1 = wfbt1*pi**2*qT_max


		wfbt_pol = wfbt1 - wfbt2


		return wfbt_pol	




	def ratio(self,had1,had2,z1,xb,param,wdt_pol,mss):
	
		mass = 1.115	
		num = self.numerator(had1,had2,z1,xb,param,wdt_pol,mss)
		den = self.denominator(had1,had2,z1,xb)

		#fnt = cr_sec()
		#fnt.mass = self.mass  
		#fact = 1#fnt.fact_fst_mom(z1,0.2,wdt_pol)		
		#if self.mdl_num != 'gauss': fact=1.
		
#		result = fact*mass*num/den
		result = mass*num/den
		
		return result

###################################################################################


	def denominator_pt(self,had1,had2,z1,xb2,pt):

		fnt = cr_sec()
		fnt.mass = self.mass        
		#fnt.frag2 = self.frag2
		fnt.qq=self.scale
		fnt.charm = self.charm

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

		qT = pt/z1

		def fnc(btt):

			res = btt*fnt.cross_sec2(had1,had2,z1,xb2,scl.mu_b(btt))
			#print(scl.mu_b(btt))
			if self.mdl_den == 'gauss' : res = res*mdl1.MD_gauss(btt,z1,xb2,wdt1_unp)
			elif self.mdl_den == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,2.,z1,xb2,1.)
			#res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*mdl2.Mf_pv17(btt,xb2) #PV17
			#res = res*mdl2.MD(btt,2.,1.)
			

			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)

			### PV17
			elif self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,z1,xb2)

			
			return res		
			
		test = lambda bt : fnc(bt)			


		

		N=10
		fbt = FBT(0)
		wfbt_unp = fbt.fbt(test,qT,N)
		wfbt_unp = wfbt_unp*2*pi	

		return wfbt_unp


	def numerator_pt(self,had1,had2,z1,xb2,pt,param,wdt_pol,mss):

		fnt = cr_sec()
		fnt.mass = self.mass        
		#fnt.frag2 = self.frag2
		fnt.qq=self.scale
		fnt.charm = self.charm

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

		qT =  pt/z1

		def fnc1(btt):

			res = btt**2*fnt.cross_sec2_polda(had1,had2,z1,xb2,scl.mu_b(btt),param)
			#print(res)
			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,xb2,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,xb2,mss)
			#res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*mdl2.Mf_pv17(btt,xb2) #PV17
			#res = res*mdl2.MD(btt,2.,1.)

			#res = res*special.struve(0,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt)*scl.g_K(btt)
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt)*scl.g_K_exp_lgm(btt,had2,z1,z2)

			### PV17
			elif self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt)*scl.g_K_bac_lgm(btt,z1,xb2)


			return res		


			
		test1 = lambda bt : fnc1(bt)			
		#test2 = lambda bt : fnc2(bt)
		
		N=10
		
		fbt1 = FBT(1)
		wfbt1 = fbt1.fbt(test1,qT,N)
		wfbt1 = wfbt1*2*pi


		wfbt_pol = wfbt1 

		return wfbt_pol	
	
	def ratio_pt(self,had1,had2,z1,xb,pt,param,wdt_pol,mss):
	
		mass = 1.115	
		num = self.numerator_pt(had1,had2,z1,xb,pt,param,wdt_pol,mss)
		den = self.denominator_pt(had1,had2,z1,xb,pt)

		#fnt = cr_sec()
		#fnt.mass = self.mass  
		#fact = 1#fnt.fact_fst_mom(z1,0.2,wdt_pol)		
		#if self.mdl_num != 'gauss': fact=1.
		
#		result = fact*mass*num/den
		result = mass*num/den
		
		return result
	
	
	
