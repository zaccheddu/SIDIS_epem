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
from def_conv_crs_1h import*
from def_conv_crs_2h import*
import multiprocessing
import time 
from datetime import datetime
import pandas as pd
import sys 
warnings.filterwarnings("ignore")

class gauss_bands():

	def lp_band(self,num):

		model_numerator = 'gauss' 
		model_denominator =  'pwr_lw_star'
		gk_type =  'll_lgm'

		fnc = polarization(0.25)
		fnc.mdl_den = model_denominator
		fnc.mdl_num = model_numerator
		 
		fnc.mass = 1.115
		fnc.g_k = gk_type

		#dati_lp = pd.read_csv("fit_grid_plot/dati_lp_conv_0.25chi_1.192.csv")
		
		dati_lp1 = pd.read_csv('fit_parameters/4plot_/lines_lp_gauss.csv')
		#dati_lp1 = dati_lp.truncate(before=103)
		#dati_lk = pd.read_csv("fit_grid_plot/dati_lk_conv_0.25chi_1.192.csv")
		now = datetime.now()
		current_time = now.strftime("%H:%M:%S")
		print("Current Time =", current_time)

		if model_numerator == 'gauss': df_prm = pd.read_csv("fit_parameters/new_bands_/new_gauss_"+str(num)+".csv")
		chi_min = 104.9
		df_prm.loc[(df_prm['chi']<chi_min+15.79)]
		#if model_numerator == 'pwr_lw_star':df_prm = pd.read_csv("fit_parameters/new_pwrlw_tot.csv")


		mins_lp = zeros(len(dati_lp1.had1))
		maxx_lp = zeros(len(dati_lp1.had1))
		j=0


		for hads1,hads2,zs1,zs2 in zip(dati_lp1.had1,dati_lp1.had2,dati_lp1.z1,dati_lp1.z2):

			if hads1 == 300 : had1='lbd'
			elif hads1 == 310 : had1='lbd_b'

			if hads2 == 100 : had2='pi+'
			elif hads2 == 105 : had2='pi-'
			elif hads2 == 200 : had2='k+'
			elif hads2 == 205 : had2='k-'

			num_tmp=zeros(len(df_prm.NUP))
			i=0
			for nup, ndo, nst, nsea, ast, bup, bsea, pp in zip(df_prm.NUP,df_prm.NDO,df_prm.NST,df_prm.NSEA,df_prm.AST,df_prm.BUP,df_prm.BSEA,df_prm.PP2):

				f_prm=arange(0.,18.,1.)
				f_prm[0] = nup
				f_prm[1] = ndo#par[0]  #do
				f_prm[2] = nst  #par[0]      #st
				f_prm[3] = nsea  #par[0]      #upb
				f_prm[4] = f_prm[3] #par[0]      #dob
				f_prm[5] = f_prm[3] # par[0]      #stb
				#¯¯¯¯¯¯
				f_prm[6] = 0 #par[0]      #aup
				f_prm[7] = 0 #par[0]      #ado
				f_prm[8] = ast  # par[0]      #ast
				f_prm[9] = 0# par[0]      #aupb
				f_prm[10] = f_prm[9] #par[0]      #adob
				f_prm[11] = f_prm[9] #par[0]      #astb
				#¯#¯¯¯¯¯¯¯
				f_prm[12] = bup #par[0]      #bup
				f_prm[13] = 0 #par[0]      #bdo
				f_prm[14] = 0# par[0]      #bst
				f_prm[15] = bsea #par[0]      #bupb
				f_prm[16] = f_prm[15] #par[0]      #bdob
				f_prm[17] = f_prm[15] # par[0]      #bstb
				pt_pp=pp
				#print(pt_pp)
				mss = 0# float(df['MSS'])
				
				num_tmp[i] = fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,mss)
				i+=1
			mins_lp[j] = min(num_tmp)
			maxx_lp[j] = max(num_tmp)

			
			j+=1	
		dati_lp1['mins']= mins_lp
		dati_lp1['maxx']= maxx_lp
		dati_lp1.to_csv(r'fit_parameters/new_bands_/lp_point_gauss_'+str(num)+'.2.csv')

		now = datetime.now()
		current_time = now.strftime("%H:%M:%S")
		print("Current Time =", current_time)
					
			
	def lk_band(self,num):

		model_numerator = 'gauss' 
		model_denominator =  'pwr_lw_star'
		gk_type =  'll_lgm'

		fnc = polarization(0.25)
		fnc.mdl_den = model_denominator
		fnc.mdl_num = model_numerator
		 
		fnc.mass = 1.115
		fnc.g_k = gk_type

		#dati_lp = pd.read_csv("fit_grid_plot/dati_lp_conv_0.25chi_1.192.csv")
		dati_lk1 =  pd.read_csv('fit_parameters/4plot_/lines_lk_gauss.csv')
		#dati_lk1 = dati_lk.truncate(before=103)
		if model_numerator == 'gauss': df_prm =pd.read_csv("fit_parameters/new_bands_/new_gauss_"+str(num)+".csv")
		#if model_numerator == 'pwr_lw_star':df_prm = pd.read_csv("fit_parameters/new_pwrlw_tot.csv")
		chi_min = 104.9
		df_prm.loc[(df_prm['chi']<chi_min+15.79)]


		j=0
		mins_lk = zeros(len(dati_lk1.had1))
		maxx_lk = zeros(len(dati_lk1.had1))

		now = datetime.now()
		current_time = now.strftime("%H:%M:%S")
		print("Current Time =", current_time)
		for hads1,hads2,zs1,zs2 in zip(dati_lk1.had1,dati_lk1.had2,dati_lk1.z1,dati_lk1.z2):

			if hads1 == 300 : had1='lbd'
			elif hads1 == 310 : had1='lbd_b'

			if hads2 == 100 : had2='pi+'
			elif hads2 == 105 : had2='pi-'
			elif hads2 == 200 : had2='k+'
			elif hads2 == 205 : had2='k-'

			num_tmp=zeros(len(df_prm.NUP))
			i=0
			for nup, ndo, nst, nsea, ast, bup, bsea, pp in zip(df_prm.NUP,df_prm.NDO,df_prm.NST,df_prm.NSEA,df_prm.AST,df_prm.BUP,df_prm.BSEA,df_prm.PP2):

				f_prm=arange(0.,18.,1.)
				f_prm[0] = nup
				f_prm[1] = ndo#par[0]  #do
				f_prm[2] = nst  #par[0]      #st
				f_prm[3] = nsea  #par[0]      #upb
				f_prm[4] = f_prm[3] #par[0]      #dob
				f_prm[5] = f_prm[3] # par[0]      #stb
				#¯¯¯¯¯¯
				f_prm[6] = 0 #par[0]      #aup
				f_prm[7] = 0 #par[0]      #ado
				f_prm[8] = ast  # par[0]      #ast
				f_prm[9] = 0# par[0]      #aupb
				f_prm[10] = f_prm[9] #par[0]      #adob
				f_prm[11] = f_prm[9] #par[0]      #astb
				#¯#¯¯¯¯¯¯¯
				f_prm[12] = bup #par[0]      #bup
				f_prm[13] = 0 #par[0]      #bdo
				f_prm[14] = 0# par[0]      #bst
				f_prm[15] = bsea #par[0]      #bupb
				f_prm[16] = f_prm[15] #par[0]      #bdob
				f_prm[17] = f_prm[15] # par[0]      #bstb
				pt_pp=pp
				#print(pt_pp)
				mss = 0# float(df['MSS'])

				num_tmp[i] = fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,mss)
				i+=1

			mins_lk[j] = min(num_tmp)
			maxx_lk[j] = max(num_tmp)
			j+=1
			
		dati_lk1['mins'] = mins_lk
		dati_lk1['maxx'] = maxx_lk

		dati_lk1.to_csv(r'fit_parameters/new_bands_/lk_point_gauss_'+str(num)+'.2.csv')
		
		now = datetime.now()
		current_time = now.strftime("%H:%M:%S")
		print("Current Time =", current_time)
	
	
	
class pwrlw_bands():	
	
	def lp_band(self,num):

		model_numerator = 'pwr_lw_star' 
		model_denominator =  'pwr_lw_star'
		gk_type =  'll_lgm'

		fnc = polarization(0.25)
		fnc.mdl_den = model_denominator
		fnc.mdl_num = model_numerator
		 
		fnc.mass = 1.115
		fnc.g_k = gk_type

		dati_lp = pd.read_csv('fit_parameters/4plot_/dati_lp_conv_0.25_chi_1.21_mdl_pwr_lw_star.csv')
		#dati_lk = pd.read_csv("fit_grid_plot/dati_lk_conv_0.25chi_1.192.csv")

		#if model_numerator == 'gauss': df_prm = pd.read_csv("fit_parameters/new_gauss_tot.csv")
		if model_numerator == 'pwr_lw_star':df_prm =pd.read_csv("fit_parameters/new_bands_/new_pwrlw_70k_set.csv")
		chi_min = 105.3
		df_prm.loc[(df_prm['chi']<chi_min+17.21)]


		mins_lp = zeros(len(dati_lp.had1))
		maxx_lp = zeros(len(dati_lp.had1))
		j=0


		for hads1,hads2,zs1,zs2 in zip(dati_lp.had1,dati_lp.had2,dati_lp.z1,dati_lp.z2):

			if hads1 == 300 : had1='lbd'
			elif hads1 == 310 : had1='lbd_b'

			if hads2 == 100 : had2='pi+'
			elif hads2 == 105 : had2='pi-'
			elif hads2 == 200 : had2='k+'
			elif hads2 == 205 : had2='k-'

			num_tmp=zeros(len(df_prm.NUP))
			i=0
			for nup, ndo, nst, nsea, ast, bup, bsea, pp, mss in zip(df_prm.NUP,df_prm.NDO,df_prm.NST,df_prm.NSEA,df_prm.AST,df_prm.BUP,df_prm.BSEA,df_prm.PP2,df_prm.MSS2):

				f_prm=arange(0.,18.,1.)
				f_prm[0] = nup
				f_prm[1] = ndo#par[0]  #do
				f_prm[2] = nst  #par[0]      #st
				f_prm[3] = nsea  #par[0]      #upb
				f_prm[4] = f_prm[3] #par[0]      #dob
				f_prm[5] = f_prm[3] # par[0]      #stb
				#¯¯¯¯¯¯
				f_prm[6] = 0 #par[0]      #aup
				f_prm[7] = 0 #par[0]      #ado
				f_prm[8] = ast  # par[0]      #ast
				f_prm[9] = 0# par[0]      #aupb
				f_prm[10] = f_prm[9] #par[0]      #adob
				f_prm[11] = f_prm[9] #par[0]      #astb
				#¯#¯¯¯¯¯¯¯
				f_prm[12] = bup #par[0]      #bup
				f_prm[13] = 0 #par[0]      #bdo
				f_prm[14] = 0# par[0]      #bst
				f_prm[15] = bsea #par[0]      #bupb
				f_prm[16] = f_prm[15] #par[0]      #bdob
				f_prm[17] = f_prm[15] # par[0]      #bstb
				pt_pp=pp
				#print(pt_pp)
				#mss = 0# float(df['MSS'])
				
				num_tmp[i] = fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,mss)
				i+=1
			mins_lp[j] = min(num_tmp)
			maxx_lp[j] = max(num_tmp)
			j+=1
		dati_lp['mins']	= mins_lp
		dati_lp['maxx']	= maxx_lp

		dati_lp.to_csv(r'fit_parameters/4plot_/lines_lp_pwrlw.csv')
				
				
		now = datetime.now()
		current_time = now.strftime("%H:%M:%S")
		print("Current Time =", current_time)
			
	def lk_band(self,num):

		model_numerator = 'pwr_lw_star' 
		model_denominator =  'pwr_lw_star'
		gk_type =  'll_lgm'

		fnc = polarization(0.25)
		fnc.mdl_den = model_denominator
		fnc.mdl_num = model_numerator
		 
		fnc.mass = 1.115
		fnc.g_k = gk_type

		#dati_lp = pd.read_csv("fit_grid_plot/dati_lp_conv_0.25chi_1.192.csv")
		dati_lk = pd.read_csv('fit_parameters/4plot_/dati_lk_conv_0.25_chi_1.21_mdl_pwr_lw_star.csv')
		#print(dati_lk)
		#if model_numerator == 'gauss': df_prm = pd.read_csv("fit_parameters/new_gauss_tot.csv")
		if model_numerator == 'pwr_lw_star':df_prm =pd.read_csv("fit_parameters/new_bands_/new_pwrlw_70k_set.csv")
		chi_min = 105.3
		df_prm.loc[(df_prm['chi']<chi_min+17.21)]


		j=0
		mins_lk = zeros(len(dati_lk.had1))
		maxx_lk = zeros(len(dati_lk.had1))

		now = datetime.now()
		current_time = now.strftime("%H:%M:%S")
		print("Current Time =", current_time)
		for hads1,hads2,zs1,zs2 in zip(dati_lk.had1,dati_lk.had2,dati_lk.z1,dati_lk.z2):

			if hads1 == 300 : had1='lbd'
			elif hads1 == 310 : had1='lbd_b'

			if hads2 == 100 : had2='pi+'
			elif hads2 == 105 : had2='pi-'
			elif hads2 == 200 : had2='k+'
			elif hads2 == 205 : had2='k-'

			num_tmp=zeros(len(df_prm.NUP))
			i=0
			for nup, ndo, nst, nsea, ast, bup, bsea, pp, mss in zip(df_prm.NUP,df_prm.NDO,df_prm.NST,df_prm.NSEA,df_prm.AST,df_prm.BUP,df_prm.BSEA,df_prm.PP2,df_prm.MSS2):

				f_prm=arange(0.,18.,1.)
				f_prm[0] = nup
				f_prm[1] = ndo#par[0]  #do
				f_prm[2] = nst  #par[0]      #st
				f_prm[3] = nsea  #par[0]      #upb
				f_prm[4] = f_prm[3] #par[0]      #dob
				f_prm[5] = f_prm[3] # par[0]      #stb
				#¯¯¯¯¯¯
				f_prm[6] = 0 #par[0]      #aup
				f_prm[7] = 0 #par[0]      #ado
				f_prm[8] = ast  # par[0]      #ast
				f_prm[9] = 0# par[0]      #aupb
				f_prm[10] = f_prm[9] #par[0]      #adob
				f_prm[11] = f_prm[9] #par[0]      #astb
				#¯#¯¯¯¯¯¯¯
				f_prm[12] = bup #par[0]      #bup
				f_prm[13] = 0 #par[0]      #bdo
				f_prm[14] = 0# par[0]      #bst
				f_prm[15] = bsea #par[0]      #bupb
				f_prm[16] = f_prm[15] #par[0]      #bdob
				f_prm[17] = f_prm[15] # par[0]      #bstb
				pt_pp=pp
				#print(pt_pp)
				#mss = 0# float(df['MSS'])

				num_tmp[i] = fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,mss)
				i+=1

			mins_lk[j] = min(num_tmp)
			maxx_lk[j] = max(num_tmp)
			j+=1
			
		dati_lk['mins'] = mins_lk
		dati_lk['maxx'] = maxx_lk
		print(dati_lk)
		dati_lk.to_csv(r'fit_parameters/4plot_/lines_lp_pwrlw.csv')
		
				
		now = datetime.now()
		current_time = now.strftime("%H:%M:%S")
		print("Current Time =", current_time)
		
'''	
bd_gss = gauss_bands()
bd_pwrlw = pwrlw_bands()

p1 = multiprocessing.Process(target=bd_gss.lp_band,args=(211,))
p2 = multiprocessing.Process(target=bd_gss.lp_band,args=(345,))
p3 = multiprocessing.Process(target=bd_gss.lp_band,args=(3673,))


p4 = multiprocessing.Process(target=bd_gss.lk_band,args=(211,))
p5 = multiprocessing.Process(target=bd_gss.lk_band,args=(345,))
p6 = multiprocessing.Process(target=bd_gss.lk_band,args=(3673,))

#p3 = multiprocessing.Process(target=bd_pwrlw.lp_band)
#p4 = multiprocessing.Process(target=bd_pwrlw.lk_band)


p1.start()
p2.start()
p3.start()
p4.start()
p5.start()
p6.start()


p1.join()
p2.join()
p3.join()
p4.join()
p5.join()
p6.join()
'''
	
	
		
