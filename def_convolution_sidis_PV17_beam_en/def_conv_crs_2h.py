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

import pycuba

class polarization :
    
	def __init__(self,coefficient,sep):

		self.mass = 1.115

		self.frag2= 'dss'

		#self.scale = 10.58

		self.coef = coefficient   
   	
		self.g_k = 'PV17'
		self.bmax = 0.6
	
		self.charm = 'no'

		self.mdl_den = 'gauss'
		self.mdl_num = 'gauss'

		self.sep=sep
		self.nf=3

		self.pdf_name='CT14IC'
		self.IC_num=0        	

		#self.nucleon='proton'# 'neutron', 'lead'

        	
	def denominator(self,had1,had2,z1,xb2,y):

		fnt = cr_sec(self.sep,self.IC_num,self.pdf_name)
		fnt.mass = self.mass        
		#fnt.frag2 = self.frag2
		#fnt.qq=self.scale
		fnt.charm = self.charm
		#fnt.nucleon=self.nucleon


		mdl1 = model_bt(self.sep)
		mdl1.mass = 1.115
		#mdl1.qq = self.scale
		mdl2 = model_bt(self.sep)
		mdl2.mass = 0.
		#mdl2.qq = self.scale
		
		
		wdt1_unp = 0.2
		wdt2_unp = 0.2

		scl = Soft(1,self.sep)
		scl.bmax = self.bmax 
		scl.nf=self.nf
		#scl.scale = self.scale	

		QQ = self.sep*xb2*y

		qT_max = QQ*self.coef

		#print('xb = '+str(xb2))
		#print('QQ = '+str(QQ))


		def fnc(btt):
			res = fnt.cross_sec2(had1,had2,z1,xb2,scl.mu_b(btt,xb2,y),y)
			#print(scl.mu_b(btt))
			#print('mu_b = '+str(scl.mu_b(btt,xb2,y)))

			if self.mdl_den == 'gauss' : res = res*mdl1.MD_gauss(btt,z1,xb2,y,wdt1_unp)
			elif self.mdl_den == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,2.,z1,xb2,y,1.)
			#res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*mdl2.Mf_pv17(btt,xb2,y) #PV17
			#res = res*mdl2.MD(btt,2.,1.)

			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K(btt)
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_exp_lgm(btt,had2,z1,z2)

			### PV17
			elif self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_bac_lgm(btt,z1,xb2,y)
			
			return res		
			
		test = lambda bt : fnc(bt)			


		

		N=20
		fbt = FBT(1)
		wfbt_unp = fbt.fbt(test,qT_max,N)
		wfbt_unp = wfbt_unp*2*pi*qT_max	

		return wfbt_unp


	def numerator(self,had1,had2,z1,xb2,y,param,wdt_pol,mss):

		fnt = cr_sec(self.sep,self.IC_num,self.pdf_name)
		fnt.mass = self.mass        
		#fnt.frag2 = self.frag2
		#fnt.qq=self.scale
		fnt.charm = self.charm
		#fnt.nucleon=self.nucleon

		mdl1 = model_bt(self.sep)
		mdl1.mass = 1.115
		#mdl1.qq = self.scale
		mdl2 = model_bt(self.sep)
		mdl2.mass = 0.
		#mdl2.qq = self.scale
		
		
		wdt1_unp = 0.2
		wdt2_unp = 0.2

		scl = Soft(1,self.sep)
		scl.bmax = self.bmax 
		#scl.scale = self.scale	
		scl.nf=self.nf
		QQ = self.sep*xb2*y

		qT_max = QQ*self.coef

		def fnc1(btt):

			res = btt*fnt.cross_sec2_polda(had1,had2,z1,xb2,scl.mu_b(btt,xb2,y),param,y)
			#print(res)
			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,xb2,y,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,xb2,mss)
			#res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*mdl2.Mf_pv17(btt,xb2,y) #PV17
			#res = res*mdl2.MD(btt,2.,1.)

			res = res*special.struve(0,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K(btt)
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_exp_lgm(btt,had2,z1,z2)

			### PV17
			elif self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_bac_lgm(btt,z1,xb2,y)

			return res		

		def fnc2(btt):

			res = btt*fnt.cross_sec2_polda(had1,had2,z1,xb2,scl.mu_b(btt,xb2,y),param,y)

			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,xb2,y,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,xb2,mss)
			
			#res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*mdl2.Mf_pv17(btt,xb2,y)  #PV17
			#res = res*mdl2.MD(btt,2.,1.)
		
			res = res*special.struve(1,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K(btt)
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_exp_lgm(btt,had2,z1,z2)

			### PV17
			elif self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_bac_lgm(btt,z1,xb2,y)

			return res		

			
		test1 = lambda bt : fnc1(bt)			
		test2 = lambda bt : fnc2(bt)
		
		N=20
		
		fbt2 = FBT(0)
		wfbt2 = fbt2.fbt(test2,qT_max,N)
		wfbt2 = wfbt2*pi**2*qT_max	
		

		
		fbt1 = FBT(1)
		wfbt1 = fbt1.fbt(test1,qT_max,N)
		wfbt1 = wfbt1*pi**2*qT_max


		wfbt_pol = wfbt1 - wfbt2


		return wfbt_pol	

	def ratio(self,had1,had2,z1,xb,param,wdt_pol,mss,y_dwn,y_up):
	
		mass = 1.115	
		#num = self.numerator(had1,had2,z1,xb,y,param,wdt_pol,mss)
		#den = self.denominator(had1,had2,z1,xb,y)

		def Integrand_den(ndim, xx, ncomp, ff, userdata):

			func = lambda yy : self.denominator(had1,had2,z1,xb,yy)

			x,y,z  = [xx[i] for i in range(ndim.contents.value)]
			x_v = x*(y_up-y_dwn) + y_dwn
			result = func(x_v)
			result = float(result)
			ff[0]  = result
			return 0

		def Integrand_num(ndim, xx, ncomp, ff, userdata):
		#

			func = lambda yy : self.numerator(had1,had2,z1,xb,yy,param,wdt_pol,mss)

			x,y,z  = [xx[i] for i in range(ndim.contents.value)]
			x_v =  x*(y_up-y_dwn) + y_dwn
			result = func(x_v)
			result = float(result)
			ff[0]  = result
			return 0
			
			
		NDIM = 3
		NCOMP = 1

		NNEW = 1000
		NMIN = 2
		FLATNESS = 50.

		KEY1 = 47
		KEY2 = 1
		KEY3 = 1
		MAXPASS = 5
		BORDER = 0.
		MAXCHISQ = 10.
		MINDEVIATION = .25
		NGIVEN = 0
		LDXGIVEN = NDIM
		NEXTRA = 0
		MINEVAL = 0
		MAXEVAL = 50000

		KEY=0

		#USERDATA=0.

		user=0.5

		itt=0		## seleziono l'integratore: 0 = VEGAS; 1 = CUHRE 

		dic_den = pycuba.Vegas(Integrand_den, NDIM, maxeval = 260, epsrel=1e-3 , epsabs=1e-11 ,verbose=0)

		dic_num = pycuba.Vegas(Integrand_num, NDIM, maxeval = 260, epsrel=1e-3 , epsabs=1e-11 ,verbose=0)

		cc_num=dic_num['fail']	
		a_num=dic_num['results']
		b_num=a_num[0] 

		cc_den=dic_den['fail']	
		a_den=dic_den['results']
		b_den=a_den[0] 

		result = mass*b_num['integral']/b_den['integral']
		return result


	def ratio_y(self,had1,had2,z1,xb,yy,param,wdt_pol,mss):
		#print(xb)
		mass = 1.115	
		den = self.denominator(had1,had2,z1,xb,yy)
		num =  self.numerator(had1,had2,z1,xb,yy,param,wdt_pol,mss)

		#fnt = cr_sec()
		#fnt.mass = self.mass  
		#fact = 1#fnt.fact_fst_mom(z1,0.2,wdt_pol)		
		#if self.mdl_num != 'gauss': fact=1.
		
#		result = fact*mass*num/den
		result = mass*num/den
		
		return result


###################################################################################


	def denominator_pt(self,had1,had2,z1,xb2,y,pt):
		fnt = cr_sec(self.sep,self.IC_num,self.pdf_name)

		#fnt = cr_sec(self.sep,self.IC_num)
		fnt.mass = self.mass        
		#fnt.frag2 = self.frag2
		#fnt.qq=self.scale
		fnt.charm = self.charm
		#fnt.nucleon=self.nucleon

		mdl1 = model_bt(self.sep)
		mdl1.mass = 1.115
		#mdl1.qq = self.scale
		mdl2 = model_bt(self.sep)
		mdl2.mass = 0.
		#mdl2.qq = self.scale
		
		
		wdt1_unp = 0.2
		wdt2_unp = 0.2

		scl = Soft(1,self.sep)
		scl.bmax = self.bmax 
		#scl.scale = self.scale	
		scl.nf=self.nf
		qT = pt/z1


		def fnc(btt):

			res = fnt.cross_sec2(had1,had2,z1,xb2,scl.mu_b(btt,xb2,y),y)
			#print(scl.mu_b(btt))
			if self.mdl_den == 'gauss' : res = res*mdl1.MD_gauss(btt,z1,xb2,y,wdt1_unp)
			elif self.mdl_den == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,2.,z1,xb2,y,1.)
			#res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*mdl2.Mf_pv17(btt,xb2,y) #PV17
			#res = res*mdl2.MD(btt,2.,1.)
			

			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K(btt)
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_exp_lgm(btt,had2,z1,z2)

			### PV17
			elif self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_bac_lgm(btt,z1,xb2,y)

			
			return res		


			
		test = lambda bt : fnc(bt)			


		

		N=10
		fbt = FBT(0)
		wfbt_unp = fbt.fbt(test,qT,N)
		wfbt_unp = wfbt_unp*2*pi	

		return wfbt_unp


	def numerator_pt(self,had1,had2,z1,xb2,y,pt,param,wdt_pol,mss):
		fnt = cr_sec(self.sep,self.IC_num,self.pdf_name)

		#fnt = cr_sec(self.sep,self.IC_num)
		fnt.mass = self.mass        
		#fnt.frag2 = self.frag2
		#fnt.qq=self.scale
		fnt.charm = self.charm
		#fnt.nucleon=self.nucleon

		mdl1 = model_bt(self.sep)
		mdl1.mass = 1.115
		#mdl1.qq = self.scale
		mdl2 = model_bt(self.sep)
		mdl2.mass = 0.
		#mdl2.qq = self.scale
		
		
		wdt1_unp = 0.2
		wdt2_unp = 0.2

		scl = Soft(1,self.sep)
		scl.bmax = self.bmax 
		#scl.scale = self.scale	
		scl.nf=self.nf
		qT =  pt/z1

		def fnc1(btt):

			res = btt**2*fnt.cross_sec2_polda(had1,had2,z1,xb2,scl.mu_b(btt,xb2,y),param,y)
			#print(res)
			if self.mdl_num == 'gauss' :res = res*mdl1.MD_gauss(btt,z1,xb2,y,wdt_pol)# h1 model
			elif self.mdl_num == 'pwr_lw_star' : res = res*mdl1.MD_bstar(btt,wdt_pol,z1,xb2,mss)
			#res = res*mdl2.MD_gauss(btt,z2,wdt2_unp)
			res = res*mdl2.Mf_pv17(btt,xb2,y) #PV17
			#res = res*mdl2.MD(btt,2.,1.)

			#res = res*special.struve(0,btt*qT_max)
			
			#res = res*scl.soft_pert(btt)*scl.g_K(btt)
			#res = res*scl.sudakov_integrated(btt)*scl.g_K(btt)
			if self.g_k == 'log_b':res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K(btt)
			####
			elif self.g_k == 'new_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_new_lgm(btt,had2,z1,z2)
			elif self.g_k == 'll_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_ll_lgm(btt,had2,z1,z2)
			elif self.g_k == 'exp_lgm': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_exp_lgm(btt,had2,z1,z2)

			### PV17
			elif self.g_k == 'PV17': res = res*scl.analytic_sudakov(btt,xb2,y)*scl.g_K_bac_lgm(btt,z1,xb2,y)


			return res		


			
		test1 = lambda bt : fnc1(bt)			
		#test2 = lambda bt : fnc2(bt)
		
		N=10
		
		fbt1 = FBT(1)
		wfbt1 = fbt1.fbt(test1,qT,N)
		wfbt1 = wfbt1*2*pi


		wfbt_pol = wfbt1 

		return wfbt_pol	
	

	def ratio_pt(self,had1,had2,z1,xb,pt,param,wdt_pol,mss,y_dwn,y_up):
	
		mass = 1.115	
		#num = self.numerator_pt(had1,had2,z1,xb,y,pt,param,wdt_pol,mss)
		#den = self.denominator_pt(had1,had2,z1,xb,y,pt)

		def Integrand_den(ndim, xx, ncomp, ff, userdata):
		
		
			func = lambda yy : self.denominator_pt(had1,had2,z1,xb,yy,pt)


			x,y,z  = [xx[i] for i in range(ndim.contents.value)]
			x_v = x*(y_up-y_dwn) + y_dwn
			result = func(x_v)
			result = float(result)
			ff[0]  = result
			return 0

		def Integrand_num(ndim, xx, ncomp, ff, userdata):
		#

			func = lambda yy : self.numerator_pt(had1,had2,z1,xb,yy,pt,param,wdt_pol,mss)

							
			x,y,z  = [xx[i] for i in range(ndim.contents.value)]
			x_v =  x*(y_up-y_dwn) + y_dwn
			result = func(x_v)
			result = float(result)
			ff[0]  = result
			return 0
		NDIM = 3
		NCOMP = 1

		NNEW = 1000
		NMIN = 2
		FLATNESS = 50.

		KEY1 = 47
		KEY2 = 1
		KEY3 = 1
		MAXPASS = 5
		BORDER = 0.
		MAXCHISQ = 10.
		MINDEVIATION = .25
		NGIVEN = 0
		LDXGIVEN = NDIM
		NEXTRA = 0
		MINEVAL = 0
		MAXEVAL = 50000

		KEY=0

		#USERDATA=0.

		user=0.5

		itt=0		## seleziono l'integratore: 0 = VEGAS; 1 = CUHRE 

		dic_den = pycuba.Vegas(Integrand_den, NDIM, maxeval = 260, epsrel=1e-3 , epsabs=1e-11 ,verbose=0)
		
		dic_num = pycuba.Vegas(Integrand_num, NDIM, maxeval = 260, epsrel=1e-3 , epsabs=1e-11 ,verbose=0)

		cc_num=dic_num['fail']	
		a_num=dic_num['results']
		b_num=a_num[0] 

		cc_den=dic_den['fail']	
		a_den=dic_den['results']
		b_den=a_den[0] 
		
#		result = fact*mass*num/den
		#result = mass*num/den
		result = mass*b_num['integral']/b_den['integral']		
		return result
	
	def ratio_pt_y(self,had1,had2,z1,xb,yy,pt,param,wdt_pol,mss):
	
		mass = 1.115	
		num = self.numerator_pt(had1,had2,z1,xb,yy,pt,param,wdt_pol,mss)
		den = self.denominator_pt(had1,had2,z1,xb,yy,pt)

		#fnt = cr_sec()
		#fnt.mass = self.mass  
		#fact = 1#fnt.fact_fst_mom(z1,0.2,wdt_pol)		
		#if self.mdl_num != 'gauss': fact=1.
		
#		result = fact*mass*num/den
		result = mass*num/den
		
		return result
	
	
