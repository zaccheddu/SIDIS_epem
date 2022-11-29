# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
import pandas as pd
from iminuit import Minuit
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

		self.pp=2.1
		self.mss = 1.
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
		self.mss_fix = True



		self.nst_lt_lw = -1
		self.nst_lt_up = 1

		self.pp_down_lm = 2.
		self.pp_up_lm  =None

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
		self.mdl_den = 'gauss'

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
			fit6 = Minuit(lst.least_squares_lh, NUP= self.nup ,NDO= self.ndo, NST=self.nst, NSEA=self.nsea,AUP=self.aup,ADO=self.ado,AST=self.ast,ASEA=self.asea,\
BUP=self.bup,BDO=self.bdo,BST=self.bst,BSEA=self.bsea,PP=self.pp,MSS=self.mss, error_NUP=0.001,error_NDO=0.001, error_NST=0.001,\
error_NSEA=0.001,error_AUP=0.001,error_ADO=0.001,error_AST=0.001,error_BUP=0.001,error_BDO=0.001,error_BSEA=0.001,error_PP=0.005,error_MSS=0.005,\
fix_NUP=self.nup_fix,fix_NDO=self.ndo_fix  ,fix_NST=self.nst_fix  ,fix_NSEA=self.nsea_fix ,fix_AUP=self.aup_fix ,fix_ADO=self.ado_fix ,fix_AST=self.ast_fix ,fix_ASEA=self.asea_fix,\
fix_BUP=self.bup_fix,fix_BDO=self.bdo_fix, fix_BST=self.bst_fix,fix_BSEA=self.bsea_fix, fix_PP=self.pp_fix,fix_MSS=self.mss_fix, limit_NUP=(None,None),limit_NDO=(None,None),limit_NST=(None,None ), limit_NSEA=(None,None), limit_BUP=(0., 20.),limit_AST=(0., None),limit_BSEA=(0., None), limit_PP=(self.pp_down_lm, self.pp_up_lm ),limit_MSS=(0., None),  errordef=1)
			#m.migrad()
			#print(m.migrad())
			print(fit6.get_param_states())
			fit6.migrad()
			print(fit6.get_param_states())
			#fit6.get_param_states()
			mean=fit6.np_values()
			print('lambda_had')
			print('the chi square d.o.f value is:')
			print(fit6.fval/(len(lst.z1) -count_nonzero(mean)))
			print('parameters:')
			print(fit6.np_values())
			
			print( 'lambda mass = ' + str(self.mass) )
			print( 'fragmentation set pion/k = ' + str(self.ff2) )
			print( 'coef = ' + str(self.coef) )
			print('g_k type = ' + str(self.g_k_2h))

			print( 'unpolarized gauss. width = ' + str(lst.unp_wd) )
			if self.mns == 'y' :fit6.minos()
			fit6.get_param_states()
			print(fit6.get_param_states())

			end = time.time()

			mins=(end -start)/60
			print('¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯')
			print('time passed:')

			print(str(end - start) +'   ' + 'sec')
			print(str((end - start)/60) +'   ' + 'min')

			print('data point: 2-h = ' + str(len(lst.z1))  )
			cov=fit6.np_matrix()

			chi_min=fit6.fval

		elif self.type=='both' :
			print('lambda_had + lambda_thrust')
			fit6 = Minuit(lst.least_squares_ljh, NUP= self.nup ,NDO= self.ndo, NST=self.nst, NSEA=self.nsea,AUP=self.aup,ADO=self.ado,AST=self.ast,ASEA=self.asea,\
BUP=self.bup,BDO=self.bdo,BST=self.bst,BSEA=self.bsea,PP=self.pp,MSS=self.mss, error_NUP=0.001,error_NDO=0.001, error_NST=0.001,\
error_NSEA=0.001,error_AUP=0.001,error_ADO=0.001,error_AST=0.001,error_BUP=0.001,error_BDO=0.001,error_BSEA=0.001,error_PP=0.005,error_MSS=0.005,\
fix_NUP=self.nup_fix,fix_NDO=self.ndo_fix  ,fix_NST=self.nst_fix  ,fix_NSEA=self.nsea_fix ,fix_AUP=self.aup_fix ,fix_ADO=self.ado_fix ,fix_AST=self.ast_fix ,fix_ASEA=self.asea_fix,\
fix_BUP=self.bup_fix,fix_BDO=self.bdo_fix, fix_BST=self.bst_fix,fix_BSEA=self.bsea_fix, fix_PP=self.pp_fix,fix_MSS=self.mss_fix, limit_NUP=(None,None),limit_NDO=(None,None),limit_NST=(None,None ), limit_NSEA=(None,None), limit_BUP=(0., 20.),limit_AST=(0., None),limit_BSEA=(0., None), limit_PP=(self.pp_down_lm, self.pp_up_lm ),limit_MSS=(0., None),  errordef=1)
			#m.migrad()
			#print(m.migrad())
			print(fit6.get_param_states())
			fit6.migrad()
			#print(fit6.get_param_states())
			#fit6.get_param_states()
			mean=fit6.np_values()
			print('lambda_had + lambda_thrust')
			print('the chi square d.o.f value is:')
			print(fit6.fval/(len(lst.z1) + len(lst.zh) -count_nonzero(mean)))
			print('parameters:')
			print(fit6.np_values())
			
			print( 'lambda masss= ' + str(self.mass) )
			print( 'fragmentation set pion/k = ' + str(self.ff2) )
			print( 'coef = ' + str(self.coef) )
			print( 'unpolarized gauss. width = ' + str(lst.unp_wd) )
			if self.mns == 'y' :fit6.minos()
			#fit6.get_param_states()
			#print(fit6.get_param_states())

			end = time.time()

			mins=(end -start)/60
			print('¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯')
			print('time :')

			print(str(end - start) +'   ' + 'sec')
			print(str((end - start)/60) +'   ' + 'min')
			
			print('data point: 2-h = :' + str(len(lst.z1)) +'  ' +' ; 1-h = :' + str(len(lst.zh)) )
			
			cov=fit6.np_matrix()

			chi_min=fit6.fval

		if self.type=='hadron' :chi_dof = fit6.fval/(len(lst.z1) -count_nonzero(mean))
		elif self.type=='both' : chi_dof = fit6.fval/(len(lst.z1) + len(lst.zh) -count_nonzero(mean))
		clm = fit6.parameters
		clm.append('coef')
		clm.append('chi_sq')
		df = pd.DataFrame(columns = clm)
		
		values = fit6.np_values()
		values = np.append(values,self.coef)
		print(chi_dof)
		chi_dof = round(chi_dof,3)
		values = np.append(values,chi_dof)

		df.loc[len(df), :] = values
		df.to_csv(r'fit_parameters/fit_'+str(self.type)+'_coef_'+ str(self.coef)+'_chi_'+str(chi_dof)+'.csv',index=False)






		return mean, cov, chi_min, fit6

### how to run

'''
ft = fitter('hadron')  # 'both'
ft.mdl_den = 'gauss'
ft.mdl_num = 'gauss'

ft.g_k_2h = 'log_b'

prms, cov_m, chi_min, fit_info = ft.fit()
'''






