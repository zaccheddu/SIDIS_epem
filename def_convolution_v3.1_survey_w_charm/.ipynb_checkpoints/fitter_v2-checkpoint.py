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

class fitter:
	"""   """

	def __init__(self,type_fit):
		self.type = type_fit

		self.nup = 0.1
		self.ndo = -0.1
		self.nst = -0.1
		self.nsea = -0.2

		self.aup = 0.
		self.ado = 0.
		self.ast = 1.
		self.asea = 0.
		
		self.bup = 1.
		self.bdo = 0.
		self.bst = 0.
		self.bsea = 1.

		self.pp=0.05
		self.mss= 0.
##______________
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

		self.nst_lt_lw = -1
		self.nst_lt_up = 1
	


		self.mass = 1.115

		self.wd = 0.2

		self.mns = 'n'
		
		self.cut_h2= 5   # 0 : global, 1: cut umb; 2 : my cut; 5: cut_articolo
		self.cut_h1 = 'all'

		self.ff2= 'dss'

		self.coef = 0.27

		self.g_k_1h = 'log_b'
		self.g_k_2h = 'log_b'
		
		self.mdl_num = 'pwr_lw'	
		self.mdl_den = 'pwr_lw'	
	def fit(self): 


		now = datetime.now()

		current_time = now.strftime("%H:%M:%S")
		print("Current Time =", current_time)


		ct_h2 = self.cut_h2
		ct_h1 = self.cut_h1
		
		lst = least_sq(ct_h2,ct_h1)
		lst.mm = self.mass
		#lst.data_cut = self.cut
		lst.unp_wd = self.wd
		lst.f2 = self.ff2
	
		lst.coef = self.coef

		lst.g_k_1h = self.g_k_1h 
		lst.g_k_2h = self.g_k_2h 

		lst.mdl_num = self.mdl_num
		lst.mdl_den = self.mdl_den

		start = time.time()

		if self.type=='hadron' :
			print('lambda_had')
			fit6 = Minuit(lst.least_squares_lh, NUP=self.nup,NDO= self.ndo,NST=self.nst,NSEA=self.nsea,\
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

			fit6.limits['PP'] = (2.,None)
			fit6.limits['MSS'] = (1, None)
			
			##########################
			#  fix params

			fit6.fixed['NUP'] = False
			fit6.fixed['NDO'] = False
			fit6.fixed['NST'] = False
			fit6.fixed['NSEA'] = False

			fit6.fixed['AUP'] = True
			fit6.fixed['ADO'] = True
			fit6.fixed['AST'] = False
			fit6.fixed['ASEA'] = True

			fit6.fixed['BUP'] = False
			fit6.fixed['BDO'] = True
			fit6.fixed['BST'] = True
			fit6.fixed['BSEA'] = False

			fit6.fixed['PP'] = False
			fit6.fixed['MSS'] = False

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

			fit6.limits['PP'] = (2, None)
			fit6.limits['MSS'] = (1, None)
			
			##########################
			#  fix params

			fit6.fixed['NUP'] = False
			fit6.fixed['NDO'] = False
			fit6.fixed['NST'] = False
			fit6.fixed['NSEA'] = False

			fit6.fixed['AUP'] = True
			fit6.fixed['ADO'] = True
			fit6.fixed['AST'] = False
			fit6.fixed['ASEA'] = True

			fit6.fixed['BUP'] = False
			fit6.fixed['BDO'] = True
			fit6.fixed['BST'] = True
			fit6.fixed['BSEA'] = False

			fit6.fixed['PP'] = False
			fit6.fixed['MSS'] = False

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

		if self.type=='hadron' :chi_dof = fit6.fval/(len(lst.z1) -fit6.nfit )
		elif self.type=='both' : chi_dof = fit6.fval/(len(lst.z1) + len(lst.zh) -fit6.nfit )
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
		df.to_csv(r'fit_parameters/fit_'+str(self.type)+'_coef_'+ str(self.coef)+'_chi_'+str(chi_dof)+'.csv',index=False)






		return mean, cov, chi_min, fit6







