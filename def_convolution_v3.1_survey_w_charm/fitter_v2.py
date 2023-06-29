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

		self.tollerance = 20.		
		
		self.g_k_1h = 'll_lgm'
		self.g_k_2h = 'll_lgm'


		#self.g_k_2h = 'll_lgm'
		#self.g_k_1h = self.g_k_2h
		
		self.mdl_num = 'pwr_lw_star' #'gauss' #	
		self.mdl_den = 'pwr_lw_star'	

		self.su2='no'

		self.charm = 'yes'

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

		lst.su2 = self.su2
		lst.charm = self.charm 
		start = time.time()

		if self.type=='hadron' :
			print('lambda_had')
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
		df.to_csv(r'fit_parameters/fit_'+str(self.type)+'_coef_'+ str(self.coef)+'_chi_'+str(chi_dof)+'__'+str(fit6.valid)+'_gk_'+str(self.g_k_2h)+ '.csv',index=False)

		
		sourceFile = open(r'fit_parameters/fit_'+str(self.type)+'_coef_'+ str(self.coef)+'_chi_'+str(chi_dof)+'__'+str(fit6.valid)+'_gk_'+str(self.g_k_2h)+'.txt', 'w')
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
		print('SU2 simmetry = '+ str(self.su2),file = sourceFile)
		print('Charm contribution = '+ str(self.charm),file = sourceFile)
		
		print('unpolarized FF2 = '+ str(self.ff2),file = sourceFile)
		
		
		
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
#ft = fitter('hadron')
#ft.fit()
'''
ft,ft2,gt,gt2 = fitter('hadron'),fitter('hadron'),fitter('hadron'),fitter('hadron')


ft.mdl_den = 'pwr_lw_star'
ft.mdl_num = 'gauss'
ft.g_k_1h = 'll_lgm'
ft.g_k_2h = 'll_lgm'

ft.pp = .1
ft.pp_down_lm = 0.
ft.pp_up_lm  = 0.2

ft.mss=0
ft.mss_fix = True

ft.su2='no'
ft.coef = 0.25
ft.ff2= 'dss'

ft.fit()
'''
ft,ft2,gt,gt2 = fitter('hadron'),fitter('hadron'),fitter('hadron'),fitter('hadron')


ft.mdl_den = 'pwr_lw_star'
ft.mdl_num = 'pwr_lw_star' # 'gauss' # 
ft.g_k_1h = 'll_lgm'
ft.g_k_2h = 'll_lgm'

ft.pp =2.1 # .1
ft.pp_down_lm = 1.#0.
ft.pp_up_lm  = None #0.2

ft.mss=0.6#0
ft.mss_fix =False #True

ft.su2='no'
ft.charm = 'no'
ft.coef = 0.25
ft.ff2= 'dss'




ft2.mdl_den = 'pwr_lw_star'
ft2.mdl_num = 'pwr_lw_star' # 'gauss'
ft2.g_k_1h = 'll_lgm'
ft2.g_k_2h = 'll_lgm'

ft2.pp =2.1 # .1
ft2.pp_down_lm = 1.#0.
ft2.pp_up_lm  = None #0.2

ft2.mss=0.6#0
ft2.mss_fix =False #True

ft2.su2='no'
ft2.charm = 'no'
ft2.coef = 0.25
ft2.cut_h2= 0
ft2.ff2= 'dss'


######
gt.mdl_den = 'pwr_lw_star'
gt.mdl_num = 'pwr_lw_star' # 'gauss'
gt.g_k_1h = 'll_lgm'
gt.g_k_2h = 'll_lgm'

gt.pp =2.1 # .1
gt.pp_down_lm = 1.#0.
gt.pp_up_lm  = None #0.2

gt.mss=0.6#0
gt.mss_fix =False #True

gt.su2='yes'
gt.charm = 'no'
gt.coef = 0.25
gt.cut_h2= 5

gt.bup = 1.2
gt.bdo = 1.
gt.bst = 1.
gt.bsea = 0.

gt.bup_fix = False
gt.bdo_fix = False
gt.bst_fix = False
gt.bsea_fix = True

gt.ff2= 'dss'


####
gt2.mdl_den = 'pwr_lw_star'
gt2.mdl_num = 'pwr_lw_star' # 'gauss'
gt2.g_k_1h = 'll_lgm'
gt2.g_k_2h = 'll_lgm'

gt2.pp =2.1 # .1
gt2.pp_down_lm = 1.#0.
gt2.pp_up_lm  = None #0.2

gt2.mss=0.6#0
gt2.mss_fix =False #True

gt2.su2='yes'
gt2.charm = 'no'
gt2.coef = 0.25
gt2.cut_h2= 0

gt2.bup = 1.2
gt2.bdo = 1.
gt2.bst = 1.
gt2.bsea = 0.

gt2.bup_fix = False
gt2.bdo_fix = False
gt2.bst_fix = False
gt2.bsea_fix = True

gt2.ff2= 'dss'





p1 = multiprocessing.Process(target=ft.fit)
p2 = multiprocessing.Process(target=ft2.fit)
p3 = multiprocessing.Process(target=gt.fit)
p4 = multiprocessing.Process(target=gt2.fit)
#p5 = multiprocessing.Process(target=jt.fit)
#p6 = multiprocessing.Process(target=jt2.fit)
#p7 = multiprocessing.Process(target=ht.fit)
#p8 = multiprocessing.Process(target=ht2.fit)
#p9 = multiprocessing.Process(target=kt.fit)
#p10 = multiprocessing.Process(target=kt2.fit)
#p11 = multiprocessing.Process(target=lt.fit)
#p12 = multiprocessing.Process(target=lt2.fit)



p1.start()
p2.start()
p3.start()
p4.start()
#p5.start()
#p6.start()
#p7.start()
#p8.start()
#p9.start()
#p10.start()
#p11.start()
#p12.start()

p1.join()
p2.join()
p3.join()
p4.join()
#p5.join()
#p6.join()
#p7.join()
#p8.join()
#p9.join()
#p10.join()
#p11.join()
#p12.join()









