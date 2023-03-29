# -*- coding: utf-8 -*-


from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
import pandas as pd
from iminuit import Minuit, minuit
from least_squares_lh import*
import time 
from datetime import datetime
import multiprocessing

class fitter:
	"""   """

	def __init__(self,type_fit):
		self.type = type_fit
		
		self.nup = 0.19
		self.ndo = -0.22
		self.nst = -0.32
		self.nsea = -0.11

		self.aup = 0.
		self.ado = 0.
		self.ast = 2.1
		self.asea = 0.
		
		self.bup = 3.2
		self.bdo = 0.
		self.bst = 0.
		self.bsea = 2.3

		self.pp=1.2
		self.mss= .6
##_____________6
##--------------
		self.nup_fix = False
		self.ndo_fix = False
		self.nst_fix = False
		self.nsea_fix = False

		self.aup_fix = True

		self.ado_fix = True
		self.ast_fix = False
		self.asea_fix = True
		
		self.bup_fix = False
		self.bdo_fix = True
		self.bst_fix = True
		self.bsea_fix = False

		self.pp_fix = False
		self.mss_fix = False

		self.nst_lt_lw = -1
		self.nst_lt_up = 1

		self.pp_down_lm = 1.
		self.pp_up_lm  = None
		
##################################################à
		'''
		self.nup =0.197
		self.ndo = -0.222
		self.nst = -0.34
		self.nsea = -0.113

		self.aup = 0.
		self.ado = 0.
		self.ast =  2.5
		self.asea = 0.
		
		self.bup = 3.294
		self.bdo = 0.
		self.bst = 0.
		self.bsea = 2.1

		self.pp=2.1
		self.mss= .3
##______________
##--------------
		self.nup_fix = True
		self.ndo_fix = True
		self.nst_fix = True
		self.nsea_fix = True

		self.aup_fix = True

		self.ado_fix = True
		self.ast_fix = True
		self.asea_fix = True
		
		self.bup_fix = True
		self.bdo_fix = True
		self.bst_fix = True
		self.bsea_fix = True

		self.pp_fix = False
		self.mss_fix = False

		self.nst_lt_lw = -1
		self.nst_lt_up = 1

		self.pp_down_lm = 2.
		self.pp_up_lm  =None
		'''
################################################	




		self.mass = 1.115

		self.wd = 0.2

		self.mns = 'n'
		
		self.cut_h2= 5   # 0 : global, 1: cut umb; 2 : my cut; 5: cut_articolo
		self.cut_h1 ='all' #'cut_lwz' #'all' #'cut_lwz' #'cut_hgz' #'cut_lwz' #'all' #'cut_lwz' # 'all' #'cut_hgz' # 'all' #'cut_lwz' #  # 'all' #

		self.ff2= 'dss'

		self.coef = 0.25 # # 0.3 #

		self.tollerance = 2.	## standard 20	
		
		self.g_k_1h = 'PV17'
		self.g_k_2h = 'PV17'


		#self.g_k_2h = 'll_lgm'
		#self.g_k_1h = self.g_k_2h
		
		self.mdl_num = 'pwr_lw_star' #'gauss' #	
		self.mdl_den = 'pwr_lw_star'	

		self.su2='no'
		self.charm= 'no'
		self.nf=3
		self.scale = 10.58

		self.correct='no'
	def fit(self): 


		now = datetime.now()

		current_time = now.strftime("%H:%M:%S")
		print("Current Time =", current_time)


		ct_h2 = self.cut_h2
		ct_h1 = self.cut_h1
		
		charm = self.charm
		mdl_den=self.mdl_den
		
		
		lst = least_sq(ct_h2,ct_h1,charm,mdl_den)
		lst.mm = self.mass
		#lst.data_cut = self.cut
		lst.unp_wd = self.wd
		lst.f2 = self.ff2
	
		lst.coef = self.coef

		lst.g_k_1h = self.g_k_1h 
		lst.g_k_2h = self.g_k_2h 

		lst.mdl_num = self.mdl_num
		lst.mdl_den = self.mdl_den

		lst.su2 = self.su2
		lst.nf=self.nf
		#lst.charm = self.charm
		lst.scale = self.scale 
		lst.correct = self.correct
		#print(lst.unp_cr)
		#print(self.charm)
		#print(lst.dataframe)
		

		start = time.time()

		if self.type=='hadron' :
			print('lambda_had')
			print(self.mdl_num)
			fit6 = Minuit(lst.least_squares_lh, NUP=self.nup,NDO= self.ndo,NST=self.nst,NSEA=self.nsea,\
			AUP=self.aup,ADO=self.ado,AST=self.ast,ASEA=self.asea,BUP=self.bup,BDO=self.bdo,BST=self.bst,\
			BSEA=self.bsea,PP=self.pp,MSS=self.mss)
			fit6.tol = self.tollerance
			print(fit6.tol)
###___________________________________________________________________________________________________
			fit6.errordef=1  # set = 1 for least-square minimization
			nn_call = 4000

			##########################
			#    limits

			fit6.limits['NUP'] = (None,None)
			fit6.limits['NDO'] = (None,None)
			fit6.limits['NST'] = (None,None)
			fit6.limits['NSEA'] = (None,None)

			fit6.limits['AUP'] = (0, None)
			fit6.limits['ADO'] = (0, None)
			fit6.limits['AST'] = (0, None)
			fit6.limits['ASEA'] = (0, None)

			fit6.limits['BUP'] = (0, None)
			fit6.limits['BDO'] = (0, None)
			fit6.limits['BST'] = (0, None)
			fit6.limits['BSEA'] = (0, None)

			fit6.limits['PP'] = (self.pp_down_lm,self.pp_up_lm)
			fit6.limits['MSS'] = (0, None)
			
			##########################
			#  fix params

			fit6.fixed['NUP'] = self.nup_fix
			fit6.fixed['NDO'] = self.ndo_fix
			fit6.fixed['NST'] = self.nst_fix 
			fit6.fixed['NSEA'] = self.nsea_fix

			fit6.fixed['AUP'] = self.aup_fix
			fit6.fixed['ADO'] = self.ado_fix
			fit6.fixed['AST'] = self.ast_fix
			fit6.fixed['ASEA'] = self.asea_fix

			fit6.fixed['BUP'] = self.bup_fix 
			fit6.fixed['BDO'] = self.bdo_fix
			fit6.fixed['BST'] = self.bst_fix
			fit6.fixed['BSEA'] = self.bsea_fix

			fit6.fixed['PP'] = self.pp_fix
			fit6.fixed['MSS'] = self.mss_fix

			##########################
			#  step 


			fit6.errors['NUP'] = 0.01
			fit6.errors['NDO'] = 0.01
			fit6.errors['NST'] = 0.01
			fit6.errors['NSEA'] = 0.01

			fit6.errors['AUP'] = 0.01
			fit6.errors['ADO'] = 0.01
			fit6.errors['AST'] = 0.01
			fit6.errors['ASEA'] = 0.01

			fit6.errors['BUP'] = 0.01
			fit6.errors['BDO'] = 0.01
			fit6.errors['BST'] = 0.01
			fit6.errors['BSEA'] = 0.01

			fit6.errors['PP'] = 0.01
			fit6.errors['MSS'] = 0.01


			#m.migrad()
			#print(m.migrad())
			print(fit6.params)
			fit6.migrad(nn_call)
			#fit6.migrad()

			print(fit6.params)
			print(fit6.fmin)
			#fit6.get_param_states()
			mean=np.array(fit6.values)
			print('lambda_had')
			print('the chi square d.o.f value is:')
			print(fit6.fval/(len(lst.z1) -fit6.nfit ))
			print('parameters:')
			#print(fit6.np_values())
			
			print( 'lambda mass = ' + str(self.mass) )
			print( 'fragmentation set pion/k = ' + str(self.ff2) )
			print( 'coef = ' + str(self.coef) )
			print('g_k type = ' + str(self.g_k_2h))

			print( 'unpolarized gauss. width = ' + str(lst.unp_wd) )
			if self.mns == 'y' :fit6.minos()
			fit6.params
			print(fit6.params)

			end = time.time()

			mins=(end -start)/60
			print('¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯')
			print('time passed:')

			print(str(end - start) +'   ' + 'sec')
			print(str((end - start)/60) +'   ' + 'min')

			print('data point: 2-h = ' + str(len(lst.z1))  )
			cov=fit6.covariance.correlation()

			chi_min=fit6.fval

		elif self.type=='both' :
			print('lambda_had + lambda_thrust')
			fit6 = Minuit(lst.least_squares_ljh, NUP=self.nup,NDO= self.ndo,NST=self.nst,NSEA=self.nsea,\
			AUP=self.aup,ADO=self.ado,AST=self.ast,ASEA=self.asea,BUP=self.bup,BDO=self.bdo,BST=self.bst,\
			BSEA=self.bsea,PP=self.pp,MSS=self.mss)
			fit6.tol = self.tollerance
			print(fit6.tol)

###___________________________________________________________________________________________________
			fit6.errordef=1  # set = 1 for least-square minimization
			nn_call = 4000

			##########################
			#    limits

			fit6.limits['NUP'] = (None, None)
			fit6.limits['NDO'] = (None, None)
			fit6.limits['NST'] = (None, None)
			fit6.limits['NSEA'] = (None, None)

			fit6.limits['AUP'] = (0, None)
			fit6.limits['ADO'] = (0, None)
			fit6.limits['AST'] = (0, None)
			fit6.limits['ASEA'] = (0, None)

			fit6.limits['BUP'] = (0, None)
			fit6.limits['BDO'] = (0, None)
			fit6.limits['BST'] = (0, None)
			fit6.limits['BSEA'] = (0, None)

			fit6.limits['PP'] = (self.pp_down_lm, self.pp_up_lm)
			fit6.limits['MSS'] = (0., None)
			
			##########################
			#  fix params

			fit6.fixed['NUP'] = self.nup_fix
			fit6.fixed['NDO'] = self.ndo_fix
			fit6.fixed['NST'] = self.nst_fix 
			fit6.fixed['NSEA'] = self.nsea_fix

			fit6.fixed['AUP'] = self.aup_fix
			fit6.fixed['ADO'] = self.ado_fix
			fit6.fixed['AST'] = self.ast_fix
			fit6.fixed['ASEA'] = self.asea_fix

			fit6.fixed['BUP'] = self.bup_fix 
			fit6.fixed['BDO'] = self.bdo_fix
			fit6.fixed['BST'] = self.bst_fix
			fit6.fixed['BSEA'] = self.bsea_fix

			fit6.fixed['PP'] = self.pp_fix
			fit6.fixed['MSS'] = self.mss_fix

			##########################
			#  step 


			fit6.errors['NUP'] = 0.01
			fit6.errors['NDO'] = 0.01
			fit6.errors['NST'] = 0.01
			fit6.errors['NSEA'] = 0.01

			fit6.errors['AUP'] = 0.01
			fit6.errors['ADO'] = 0.01
			fit6.errors['AST'] = 0.01
			fit6.errors['ASEA'] = 0.01

			fit6.errors['BUP'] = 0.01
			fit6.errors['BDO'] = 0.01
			fit6.errors['BST'] = 0.01
			fit6.errors['BSEA'] = 0.01

			fit6.errors['PP'] = 0.01
			fit6.errors['MSS'] = 0.01
			
			
			
			#m.migrad()
			#print(m.migrad())
			print(fit6.params)
			fit6.migrad(nn_call)
			print(fit6.params)
			print(fit6.fmin)
			#fit6.get_param_states()
			mean=np.array(fit6.values)
			print('lambda_had + lambda_thrust')
			print('the chi square d.o.f value is:')
			print(fit6.fval/(len(lst.z1) + len(lst.zh) -fit6.nfit ))
			print('parameters:')
			#print(fit6.np_values())
			
			print( 'lambda mass = ' + str(self.mass) )
			#print( 'fragmentation set pion/k = ' + str(self.ff2) )
			print( 'coef = ' + str(self.coef) )
			print('g_k type = ' + str(self.g_k_2h))

			print( 'unpolarized gauss. width = ' + str(lst.unp_wd) )
			if self.mns == 'y' :fit6.minos()
			fit6.params
			print(fit6.params)

			end = time.time()

			mins=(end -start)/60
			print('¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯')
			print('time passed:')

			print(str(end - start) +'   ' + 'sec')
			print(str((end - start)/60) +'   ' + 'min')

			print('data point: 2-h = ' + str(len(lst.z1))  )
			cov=fit6.covariance.correlation()

			chi_min=fit6.fval
			


		elif self.type=='single' :
			print(' lambda_thrust')
			fit6 = Minuit(lst.least_squares_lj, NUP=self.nup,NDO= self.ndo,NST=self.nst,NSEA=self.nsea,\
			AUP=self.aup,ADO=self.ado,AST=self.ast,ASEA=self.asea,BUP=self.bup,BDO=self.bdo,BST=self.bst,\
			BSEA=self.bsea,PP=self.pp,MSS=self.mss)
###___________________________________________________________________________________________________
			fit6.errordef=1  # set = 1 for least-square minimization
			nn_call = 2000

			##########################
			#    limits

			fit6.limits['NUP'] = (None, None)
			fit6.limits['NDO'] = (None, None)
			fit6.limits['NST'] = (None, None)
			fit6.limits['NSEA'] = (None, None)

			fit6.limits['AUP'] = (0, None)
			fit6.limits['ADO'] = (0, None)
			fit6.limits['AST'] = (0, None)
			fit6.limits['ASEA'] = (0, None)

			fit6.limits['BUP'] = (0, None)
			fit6.limits['BDO'] = (0, None)
			fit6.limits['BST'] = (0, None)
			fit6.limits['BSEA'] = (0, None)

			fit6.limits['PP'] = (self.pp_down_lm, self.pp_up_lm)
			fit6.limits['MSS'] = (0., None)
			
			##########################
			#  fix params

			fit6.fixed['NUP'] = self.nup_fix
			fit6.fixed['NDO'] = self.ndo_fix
			fit6.fixed['NST'] = self.nst_fix 
			fit6.fixed['NSEA'] = self.nsea_fix

			fit6.fixed['AUP'] = self.aup_fix
			fit6.fixed['ADO'] = self.ado_fix
			fit6.fixed['AST'] = self.ast_fix
			fit6.fixed['ASEA'] = self.asea_fix

			fit6.fixed['BUP'] = self.bup_fix 
			fit6.fixed['BDO'] = self.bdo_fix
			fit6.fixed['BST'] = self.bst_fix
			fit6.fixed['BSEA'] = self.bsea_fix

			fit6.fixed['PP'] = self.pp_fix
			fit6.fixed['MSS'] = self.mss_fix

			##########################
			#  step 


			fit6.errors['NUP'] = 0.01
			fit6.errors['NDO'] = 0.01
			fit6.errors['NST'] = 0.01
			fit6.errors['NSEA'] = 0.01

			fit6.errors['AUP'] = 0.01
			fit6.errors['ADO'] = 0.01
			fit6.errors['AST'] = 0.01
			fit6.errors['ASEA'] = 0.01

			fit6.errors['BUP'] = 0.01
			fit6.errors['BDO'] = 0.01
			fit6.errors['BST'] = 0.01
			fit6.errors['BSEA'] = 0.01

			fit6.errors['PP'] = 0.01
			fit6.errors['MSS'] = 0.02
			
			
			
			#m.migrad()
			#print(m.migrad())
			print(fit6.params)
			fit6.migrad(nn_call)
			print(fit6.params)
			print(fit6.fmin)
			#fit6.get_param_states()
			mean=np.array(fit6.values)
			print('lambda_had + lambda_thrust')
			print('the chi square d.o.f value is:')
			print(fit6.fval/( len(lst.zh) -fit6.nfit ))
			print('parameters:')
			#print(fit6.np_values())
			
			print( 'lambda mass = ' + str(self.mass) )
			#print( 'fragmentation set pion/k = ' + str(self.ff2) )multiprocessing
			print( 'coef = ' + str(self.coef) )
			print('g_k type = ' + str(self.g_k_2h))

			print( 'unpolarized gauss. width = ' + str(lst.unp_wd) )
			if self.mns == 'y' :fit6.minos()
			fit6.params
			print(fit6.params)

			end = time.time()

			mins=(end -start)/60
			print('¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯')
			print('time passed:')

			print(str(end - start) +'   ' + 'sec')
			print(str((end - start)/60) +'   ' + 'min')

			print('data point: 2-h = ' + str(len(lst.z1))  )
			cov=fit6.covariance.correlation()

			chi_min=fit6.fval
			#print(chi_min)


		if self.type=='hadron' :chi_dof = fit6.fval/(len(lst.z1) -fit6.nfit )
		elif self.type=='both' : chi_dof = fit6.fval/(len(lst.z1) + len(lst.zh) -fit6.nfit )
		elif self.type=='single' : chi_dof = fit6.fval/( len(lst.zh) -fit6.nfit )


		clm = fit6.parameters
		clm = list(clm)
		clm.append('coef')
		clm.append('chi_sq')
		df = pd.DataFrame(columns = clm)
		
		values = np.array(fit6.values)
		values = np.append(values,self.coef)
		print(chi_dof)
		chi_dof = round(chi_dof,3)
		values = np.append(values,chi_dof)

		df.loc[len(df), :] = values
		df.to_csv(r'fit_parameters/fit_'+str(self.type)+'_coef_'+ str(self.coef)+'_chi_'+str(chi_dof)+'__'+str(fit6.valid)+'_gk_'+str(self.g_k_2h)+'_su_'+str(self.su2)+'_charm'+str(self.charm)\
			+'_correction_'+str(self.correct) +'.csv',index=False)

		if self.cut_h2==0:
			ct_h2 = 'pions'
		elif self.cut_h2==5:
			ct_h2 = 'pions_cut'
		
		lst1 = least_sq(ct_h2,ct_h1,charm,mdl_den)
		lst1.mm = self.mass
		#lst.data_cut = self.cut
		lst1.unp_wd = self.wd
		lst1.f2 = self.ff2
	
		lst1.coef = self.coef

		lst1.g_k_1h = self.g_k_1h 
		lst1.g_k_2h = self.g_k_2h 

		lst1.mdl_num = self.mdl_num
		lst1.mdl_den = self.mdl_den

		lst1.su2 = self.su2
		lst1.nf=self.nf
		#lst.charm = self.charm
		lst1.scale = self.scale 
		lst1.correct = self.correct
		#prm=pd.read_csv('fit_parameters/fit_hadron_coef_0.27_chi_1.259__True_gk_PV17_su_no_charmyes_correction_no.csv')

		prm=pd.read_csv(r'fit_parameters/fit_'+str(self.type)+'_coef_'+ str(self.coef)+'_chi_'+str(chi_dof)+'__'+str(fit6.valid)+'_gk_'+str(self.g_k_2h)+'_su_'+str(self.su2)+'_charm'+str(self.charm)\
			+'_correction_'+str(self.correct) +'.csv')


		chi_pion=lst1.least_squares_lh(prm['NUP'].to_numpy(),prm['NDO'].to_numpy(),prm['NST'].to_numpy(),prm['NSEA'].to_numpy(),\
                         prm['AUP'].to_numpy(),prm['ADO'].to_numpy(),prm['AST'].to_numpy(),prm['ASEA'].to_numpy(),\
                         prm['BUP'].to_numpy(),prm['BDO'].to_numpy(),prm['BST'].to_numpy(),prm['BSEA'].to_numpy(),prm['PP'].to_numpy(),prm['MSS'].to_numpy())
                         
		#chi_pion=lst1.least_squares_lh(df['NUP'].to_numpy(),df['NDO'].to_numpy(),df['NST'].to_numpy(),df['NSEA'].to_numpy(),\
                #         df['AUP'].to_numpy(),df['ADO'].to_numpy(),df['AST'].to_numpy(),df['ASEA'].to_numpy(),\
                #         df['BUP'].to_numpy(),df['BDO'].to_numpy(),df['BST'].to_numpy(),df['BSEA'].to_numpy(),df['PP'].to_numpy(),df['MSS'].to_numpy())

		if self.cut_h2==0:
			ct_h2 = 'kaons'
		elif self.cut_h2==5:
			ct_h2 = 'kaons_cut'
		
		lst2 = least_sq(ct_h2,ct_h1,charm,mdl_den)
		lst2.mm = self.mass
		#lst.data_cut = self.cut
		lst2.unp_wd = self.wd
		lst2.f2 = self.ff2
	
		lst2.coef = self.coef

		lst2.g_k_1h = self.g_k_1h 
		lst2.g_k_2h = self.g_k_2h 

		lst2.mdl_num = self.mdl_num
		lst2.mdl_den = self.mdl_den

		lst2.su2 = self.su2
		lst2.nf=self.nf
		#lst.charm = self.charm
		lst2.scale = self.scale 
		lst2.correct = self.correct

		chi_kaon=lst2.least_squares_lh(prm['NUP'].to_numpy(),prm['NDO'].to_numpy(),prm['NST'].to_numpy(),prm['NSEA'].to_numpy(),\
                         prm['AUP'].to_numpy(),prm['ADO'].to_numpy(),prm['AST'].to_numpy(),prm['ASEA'].to_numpy(),\
                         prm['BUP'].to_numpy(),prm['BDO'].to_numpy(),prm['BST'].to_numpy(),prm['BSEA'].to_numpy(),prm['PP'].to_numpy(),prm['MSS'].to_numpy())


		
		sourceFile = open(r'fit_parameters/fit_'+str(self.type)+'_coef_'+ str(self.coef)+'_chi_'+str(chi_dof)+'__'+str(fit6.valid)+'_gk_'+str(self.g_k_2h)+'_su_'+str(self.su2)+'_charm'+str(self.charm) \
			+'_correction_'+str(self.correct)+'.txt', 'w')
		print(fit6.init_params, file = sourceFile)
		print(fit6, file = sourceFile)
		#
		if self.type=='hadron' :print('lambda_had',file = sourceFile)
		elif self.type=='both' : print('lambda_had + thrust',file = sourceFile)
		elif self.type=='single' :print('single',file = sourceFile)
		
		print('the chi square d.o.f value is:',file = sourceFile)
		print(chi_dof,file = sourceFile)
		#print('parameters:')
		#print(fit6.np_values())
		
		print( 'lambda mass = ' + str(self.mass),file = sourceFile )
		#print( 'fragmentation set pion/k = ' + str(self.ff2) )
		print( 'coef = ' + str(self.coef),file = sourceFile )
		print('g_k type = ' + str(self.g_k_2h),file = sourceFile)
		print('model numerator = '+str(self.mdl_num),file = sourceFile)
		print('model denominator = '+str(self.mdl_den),file = sourceFile)
		print('data cut had2 = '+str(self.cut_h2),file = sourceFile)
		print('data cut had1 = '+ str(self.cut_h1),file = sourceFile)
		
		print('________________________',file = sourceFile)
		print('1h points = '+ str(len(lst.zh)),file = sourceFile)
		print('2h points = '+ str(len(lst.z1)),file = sourceFile)
		print('chi_square = '+ str(fit6.fval),file = sourceFile)
		print('number parameters = '+ str(fit6.nfit),file = sourceFile)
		print('________________________',file = sourceFile)
		print('chi_square pions = '+ str(chi_pion)+ ' ;  points = '+ str(len(lst1.z1)),file = sourceFile)
		print('chi_square kaons = '+ str(chi_kaon)+ ' ;  points = '+ str(len(lst2.z1)),file = sourceFile)

		print('________________________',file = sourceFile)

		print('SU2 simmetry = '+ str(self.su2),file = sourceFile)
		print('Charm = '+ str(self.charm),file = sourceFile)
		print('Correction_alpha_decay = '+ str(self.correct),file = sourceFile)
		print('________________________',file = sourceFile)
		
		
		print('________________________',file = sourceFile)
		print('covariance',file = sourceFile)
		print('         ',file = sourceFile)
		print(cov,file = sourceFile)

		print('________________________',file = sourceFile)
		print('mean',file = sourceFile)
		print('         ',file = sourceFile)
		print(mean,file = sourceFile)

		sourceFile.close()



		return mean, cov, chi_min, fit6


### how to run
ft,ft2,gt,gt2,ht,ht2 = fitter('hadron'),fitter('hadron'),fitter('hadron'),fitter('hadron'),fitter('hadron'),fitter('hadron')

#print(ft.g_k_2h)
ft.mdl_den = 'pwr_lw_star'
ft.mdl_num = 'gauss'
ft.pp = 0.1
ft.pp_down_lm = 0.
ft.pp_up_lm  = 0.2
ft.mss=0.
ft.mss_fix = True
ft.correct = 'no'
ft.coef=0.27

ft.su2='no'
ft.charm= 'no'
ft.nf=3
ft.scale = 10.58
ft.cut_h2= 0
ft.correct = 'no'

#ft.fit()

###
ft2.mdl_den = 'pwr_lw_star'
ft2.mdl_num = 'gauss'
ft2.pp = .1
ft2.pp_down_lm = 0.
ft2.pp_up_lm  = 0.2
ft2.mss=0.
ft2.mss_fix = True

ft2.su2='no'
ft2.charm= 'yes'
ft2.nf=4
ft2.scale = 10.58
ft2.cut_h2= 0
ft2.coef = 0.27
ft2.ado=1.2
ft2.ado_fix = False
ft2.correct = 'no'
#ft2.bsea = 0.
#ft2.bsea_fix = True



ht.mdl_den = 'pwr_lw_star'
ht.mdl_num = 'gauss'
ht.pp =.1
ht.pp_down_lm = 0.
ht.pp_up_lm  = .2
ht.mss=0.
ht.mss_fix = True

ht.su2='yes'
ht.charm= 'yes'
ht.nf=4
ht.scale = 10.58
ht.cut_h2= 0
ht.correct = 'no'

ht.aup=0.
ht.ado=0.

ht.bup = 1.2
ht.bdo = 1.
ht.bst = 0.
ht.bsea = 0.

ht.aup_fix = True
ht.ado_fix = True

ht.bup_fix = False
ht.bdo_fix = False
ht.bst_fix = True
ht.bsea_fix = True
ht.coef=0.27


p1 = multiprocessing.Process(target=ft.fit)
p2 = multiprocessing.Process(target=ft2.fit)
#p3 = multiprocessing.Process(target=gt.fit)
#p4 = multiprocessing.Process(target=gt2.fit)
p5 = multiprocessing.Process(target=ht.fit)
#p6 = multiprocessing.Process(target=ht2.fit)


#p1.start()
#p2.start()
#p3.start()
#p4.start()
p5.start()
#p6.start()



#p1.join()
#p2.join()
#p3.join()
#p4.join()
p5.join()
#p6.join()



'''
ht.mdl_den = 'pwr_lw_star'
ht.mdl_num = 'gauss'
ht.pp =.1
ht.pp_down_lm = 0.
ht.pp_up_lm  = .2
ht.mss=0.
ht.mss_fix = True

ht.su2='yes'
ht.charm= 'yes'
ht.scale = 10.58
ht.cut_h2= 5

ht.aup=0.
ht.ado=0.

ht.bup = 1.2
ht.bdo = 1.
ht.bst = 0.
ht.bsea = 0.

ht.aup_fix = True
ht.ado_fix = True

ht.bup_fix = False
ht.bdo_fix = False
ht.bst_fix = True
ht.bsea_fix = True
ht.coef=0.27


#ft.fit()

###
ht2.mdl_den = 'pwr_lw_star'
ht2.mdl_num = 'gauss'
ht2.pp =.1
ht2.pp_down_lm = 0.
ht2.pp_up_lm  = .2
ht2.mss=0.
ht2.mss_fix = True

ht2.su2='yes'
ht2.charm= 'yes'
ht2.scale = 10.58
ht2.cut_h2= 0

ht2.aup=0.
ht2.ado=0.

ht2.bup = 1.2
ht2.bdo = 1.
ht2.bst = 0.
ht2.bsea = 0.

ht2.aup_fix = True
ht2.ado_fix = True

ht2.bup_fix = False
ht2.bdo_fix = False
ht2.bst_fix = True
ht2.bsea_fix = True

ht2.coef=0.27


ht.mdl_den = 'pwr_lw_star'
ht.mdl_num = 'gauss'
ht.pp =.1
ht.pp_down_lm = 0.
ht.pp_up_lm  = .2
ht.mss=0.
ht.mss_fix = True

ht.su2='yes'
ht.charm= 'yes'
ht.scale = 10.58
ht.cut_h2= 5

ht.aup=0.
ht.ado=0.

ht.bup = 1.2
ht.bdo = 1.
ht.bst = 0.
ht.bsea = 0.

ht.aup_fix = True
ht.ado_fix = True

ht.bup_fix = False
ht.bdo_fix = False
ht.bst_fix = True
ht.bsea_fix = True
ht.coef=0.27


#ft.fit()

###
ht2.mdl_den = 'pwr_lw_star'
ht2.mdl_num = 'gauss'
ht2.pp =.1
ht2.pp_down_lm = 0.
ht2.pp_up_lm  = .2
ht2.mss=0.
ht2.mss_fix = True

ht2.su2='yes'
ht2.charm= 'yes'
ht2.scale = 10.58
ht2.cut_h2= 0

ht2.aup=0.
ht2.ado=0.

ht2.bup = 1.2
ht2.bdo = 1.
ht2.bst = 0.
ht2.bsea = 0.

ht2.aup_fix = True
ht2.ado_fix = True

ht2.bup_fix = False
ht2.bdo_fix = False
ht2.bst_fix = True
ht2.bsea_fix = True

ht2.coef=0.27

'''













