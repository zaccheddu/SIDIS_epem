from pylab import*
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import*
from scipy import special
from scipy.stats import multivariate_normal
import pandas as pd
from least_squares_lh import*
from datetime import datetime
from def_conv_crs_1h import*
from def_conv_crs_2h import*
import time 
from datetime import datetime
#warnings.filterwarnings("ignore")
import multiprocessing

def bands_gauss(su2_charm):
	#print(num)
	
	mean = array([0.168, -0.138 , -0.161, -0.104,  2.195,   4.023, 2.905, 0.103  ])

	if su2_charm=='no_no':   ## 8 parametri chi_1.174 chi_min 103,3 ncall 767
	
		mean =array([ 0.15108095, -0.13516462, -0.14654559, -0.095736, 2.06062647, 3.83497502,  2.88930635, 0.10388819 ])
	
		cov =[[ 0.000755, -0.000192, -0.000159, -0.000368,   0.00218,   0.00582,  -0.00043,  0.000345],
		      [-0.000192,  0.000214,  -1.9e-05,   0.00012, -8.18e-05,  0.000746, -0.000141, -0.000157],
		      [-0.000159,  -1.9e-05,  0.000838,  9.97e-05,  -0.00795,  -0.00143,  0.000234, -4.26e-05],
		      [-0.000368,   0.00012,  9.97e-05,  0.000264,  -0.00126, -0.000236,  8.95e-05, -0.000244],
		      [  0.00218, -8.18e-05,  -0.00795,  -0.00126,    0.0832,    0.0188,  -0.00287,  0.000643],
		      [  0.00582,  0.000746,  -0.00143, -0.000236,    0.0188,     0.152,  -0.00488, -6.26e-05],
		      [ -0.00043, -0.000141,  0.000234,  8.95e-05,  -0.00287,  -0.00488,   0.00615, -0.000118],
		      [ 0.000345, -0.000157, -4.26e-05, -0.000244,  0.000643, -6.26e-05, -0.000118,  0.000263]]


	elif su2_charm=='no_yes':   ## 9 parametri chi_1.262 chi_min 109.8 ncall 837
	
		mean =array([ 0.16940709, -0.36707107, -0.11713259, -0.12222033, 0.79289446,  0.86638731, 2.69376086,1.59579872, 0.09534966 ])
	
		cov= [[  0.00872,  -0.00625, -0.000475,  -0.00719,    0.0146,  0.000515,     0.105,     0.111,   0.00272],
		      [ -0.00625,    0.0408,  0.000423,   0.00764,   -0.0997,  -0.00692,    -0.081,    -0.138,   -0.0016],
		      [-0.000475,  0.000423,  0.000448,  0.000425,  -0.00196,  -0.00451,  -0.00624,  -0.00589, -0.000102],
		      [ -0.00719,   0.00764,  0.000425,   0.00624,   -0.0182, -0.000981,   -0.0852,   -0.0977,  -0.00228],
		      [   0.0146,   -0.0997,  -0.00196,   -0.0182,     0.268,    0.0307,     0.193,     0.302,   0.00329],
		      [ 0.000515,  -0.00692,  -0.00451, -0.000981,    0.0307,    0.0554,    0.0129,   0.00224, -0.000566],
		      [    0.105,    -0.081,  -0.00624,   -0.0852,     0.193,    0.0129,      1.33,      1.38,    0.0308],
		      [    0.111,    -0.138,  -0.00589,   -0.0977,     0.302,   0.00224,      1.38,      1.64,    0.0342],
		      [  0.00272,   -0.0016, -0.000102,  -0.00228,   0.00329, -0.000566,    0.0308,    0.0342,  0.000941]]	


	elif su2_charm=='yes_yes':   ## 8 parametri chi_1.447 chi_min 127.4 ncall 422

		mean =array([ 0.20754979, -0.12811842, -0.25272118, -0.16499243, 2.11244742, 2.43826013,  1.15411948, 0.11715731 ])

		cov= [[0.000473,  -0.00039, -0.000328, -0.000261,   0.00222,   0.00196,   0.00539,  0.000182],
		      [-0.00039,  0.000902,  0.000387,  0.000285,  -0.00206,  0.000829,   -0.0136, -0.000209],
		      [-0.000328,  0.000387,   0.00344,  0.000231,   -0.0206, -0.000865,  -0.00657, -9.59e-05],
		      [-0.000261,  0.000285,  0.000231,  0.000274,  -0.00153,   0.00198,  -0.00291, -0.000201],
		      [ 0.00222,  -0.00206,   -0.0206,  -0.00153,     0.134,   0.00734,    0.0305,  0.000571],
		      [ 0.00196,  0.000829, -0.000865,   0.00198,   0.00734,    0.0887,     0.018,  -0.00155],
		      [ 0.00539,   -0.0136,  -0.00657,  -0.00291,    0.0305,     0.018,     0.259,   0.00207],
		      [0.000182, -0.000209, -9.59e-05, -0.000201,  0.000571,  -0.00155,   0.00207,  0.000175]]



	
	n_set = 2000
	     
#	np.random.seed(428)
	np.random.seed(np.random.random_integers(0,500))
	
	y1 = np.random.multivariate_normal(mean, cov, n_set).T
	
	df1 = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],\
			   'AST':y1[4],'BUP':y1[5],'BSEA':y1[6],'PP2':y1[7]})

#	np.random.seed(818)
	np.random.seed(np.random.random_integers(0,500))

	y1 = np.random.multivariate_normal(mean, cov, n_set).T
	
	df2 = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],\
			   'AST':y1[4],'BUP':y1[5],'BSEA':y1[6],'PP2':y1[7]})

#	np.random.seed(324)
	np.random.seed(np.random.random_integers(0,500))

	y1 = np.random.multivariate_normal(mean, cov, n_set).T
	
	df3 = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],\
			   'AST':y1[4],'BUP':y1[5],'BSEA':y1[6],'PP2':y1[7]})

#	np.random.seed(215)
	np.random.seed(np.random.random_integers(0,500))

	y1 = np.random.multivariate_normal(mean, cov, n_set).T
	
	df4 = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],\
			   'AST':y1[4],'BUP':y1[5],'BSEA':y1[6],'PP2':y1[7]})


#	np.random.seed(21)
	np.random.seed(np.random.random_integers(0,500))

	y1 = np.random.multivariate_normal(mean, cov, n_set).T
	
	df5 = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],\
			   'AST':y1[4],'BUP':y1[5],'BSEA':y1[6],'PP2':y1[7]})

#	np.random.seed(222)
	np.random.seed(np.random.random_integers(0,500))

	y1 = np.random.multivariate_normal(mean, cov, n_set).T
	
	df6 = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],\
			   'AST':y1[4],'BUP':y1[5],'BSEA':y1[6],'PP2':y1[7]})

#	np.random.seed(43)
	np.random.seed(np.random.random_integers(0,500))

	y1 = np.random.multivariate_normal(mean, cov, n_set).T
	
	df7 = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],\
			   'AST':y1[4],'BUP':y1[5],'BSEA':y1[6],'PP2':y1[7]})


	df = pd.concat([df1,df2,df3,df4,df5,df6,df7])

			   
	lst =least_sq(5,'all')
	lst.g_k_2h = 'PV17'
	lst.g_k_1h = 'PV17'

	lst.mdl_num = 'gauss'
	lst.mdl_den = 'pwr_lw_star'
	lst.coef=0.25
	
	
	now = datetime.now()
	current_time = now.strftime("%H:%M:%S")
	print("Current Time =", current_time)

	a=zeros(len(df.NUP))
	i=0
	for nup,ndo,nst,nsea,ast,bup,bsea,pp2 in zip(df.NUP, df.NDO, df.NST, df.NSEA, df.AST,df.BUP, df.BSEA, df.PP2):
		a[i]=lst.least_squares_lh(nup,ndo,nst,nsea,0,0,ast,0,bup,0,0,bsea,pp2,0)
		df['chi']=a
		df.to_csv(r'fit_parameters/new_bands_pv17/new_gauss_iter_set.csv',index=False)

		#print('gauss '+str(a[i])+' serie =' +str(num))
		i+=1
#	df['chi']=a    
#	df.to_csv(r'fit_parameters/new_bands_pv17/new_gauss_70k_set.csv',index=False)

	now = datetime.now()
	current_time = now.strftime("%H:%M:%S")
	print("Current Time =", current_time)	    

	
def bands_gauss_2(su2_charm):
	#print(num)
	
	mean = array([0.168, -0.138 , -0.161, -0.104,  2.195,   4.023, 2.905, 0.103  ])

	if su2_charm=='no_no':   ## 8 parametri chi_1.174 chi_min 103,3 ncall 767
	
		mean =array([ 0.15108095, -0.13516462, -0.14654559, -0.095736, 2.06062647, 3.83497502,  2.88930635, 0.10388819 ])
	
		cov =[[ 0.000755, -0.000192, -0.000159, -0.000368,   0.00218,   0.00582,  -0.00043,  0.000345],
		      [-0.000192,  0.000214,  -1.9e-05,   0.00012, -8.18e-05,  0.000746, -0.000141, -0.000157],
		      [-0.000159,  -1.9e-05,  0.000838,  9.97e-05,  -0.00795,  -0.00143,  0.000234, -4.26e-05],
		      [-0.000368,   0.00012,  9.97e-05,  0.000264,  -0.00126, -0.000236,  8.95e-05, -0.000244],
		      [  0.00218, -8.18e-05,  -0.00795,  -0.00126,    0.0832,    0.0188,  -0.00287,  0.000643],
		      [  0.00582,  0.000746,  -0.00143, -0.000236,    0.0188,     0.152,  -0.00488, -6.26e-05],
		      [ -0.00043, -0.000141,  0.000234,  8.95e-05,  -0.00287,  -0.00488,   0.00615, -0.000118],
		      [ 0.000345, -0.000157, -4.26e-05, -0.000244,  0.000643, -6.26e-05, -0.000118,  0.000263]]


	elif su2_charm=='no_yes':   ## 9 parametri chi_1.262 chi_min 109.8 ncall 837
	
		mean =array([ 0.16940709, -0.36707107, -0.11713259, -0.12222033, 0.79289446,  0.86638731, 2.69376086,1.59579872, 0.09534966 ])
	
		cov= [[  0.00872,  -0.00625, -0.000475,  -0.00719,    0.0146,  0.000515,     0.105,     0.111,   0.00272],
		      [ -0.00625,    0.0408,  0.000423,   0.00764,   -0.0997,  -0.00692,    -0.081,    -0.138,   -0.0016],
		      [-0.000475,  0.000423,  0.000448,  0.000425,  -0.00196,  -0.00451,  -0.00624,  -0.00589, -0.000102],
		      [ -0.00719,   0.00764,  0.000425,   0.00624,   -0.0182, -0.000981,   -0.0852,   -0.0977,  -0.00228],
		      [   0.0146,   -0.0997,  -0.00196,   -0.0182,     0.268,    0.0307,     0.193,     0.302,   0.00329],
		      [ 0.000515,  -0.00692,  -0.00451, -0.000981,    0.0307,    0.0554,    0.0129,   0.00224, -0.000566],
		      [    0.105,    -0.081,  -0.00624,   -0.0852,     0.193,    0.0129,      1.33,      1.38,    0.0308],
		      [    0.111,    -0.138,  -0.00589,   -0.0977,     0.302,   0.00224,      1.38,      1.64,    0.0342],
		      [  0.00272,   -0.0016, -0.000102,  -0.00228,   0.00329, -0.000566,    0.0308,    0.0342,  0.000941]]	


	elif su2_charm=='yes_yes':   ## 8 parametri chi_1.447 chi_min 127.4 ncall 422

		mean =array([ 0.20754979, -0.12811842, -0.25272118, -0.16499243, 2.11244742, 2.43826013,  1.15411948, 0.11715731 ])

		cov= [[0.000473,  -0.00039, -0.000328, -0.000261,   0.00222,   0.00196,   0.00539,  0.000182],
		      [-0.00039,  0.000902,  0.000387,  0.000285,  -0.00206,  0.000829,   -0.0136, -0.000209],
		      [-0.000328,  0.000387,   0.00344,  0.000231,   -0.0206, -0.000865,  -0.00657, -9.59e-05],
		      [-0.000261,  0.000285,  0.000231,  0.000274,  -0.00153,   0.00198,  -0.00291, -0.000201],
		      [ 0.00222,  -0.00206,   -0.0206,  -0.00153,     0.134,   0.00734,    0.0305,  0.000571],
		      [ 0.00196,  0.000829, -0.000865,   0.00198,   0.00734,    0.0887,     0.018,  -0.00155],
		      [ 0.00539,   -0.0136,  -0.00657,  -0.00291,    0.0305,     0.018,     0.259,   0.00207],
		      [0.000182, -0.000209, -9.59e-05, -0.000201,  0.000571,  -0.00155,   0.00207,  0.000175]]



	
	n_set = 1000
	     
#	np.random.seed(428)
	np.random.seed(np.random.random_integers(0,200))
	
	y1 = np.random.multivariate_normal(mean, cov, n_set).T
	
	if su2_charm=='no_no': df = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],'AST':y1[4],'BUP':y1[5],'BSEA':y1[6],'PP2':y1[7]})

	elif su2_charm=='no_yes': df = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],'ADO':y1[4],'AST':y1[5],'BUP':y1[6],'BSEA':y1[7],'PP2':y1[8]})
	
	elif su2_charm=='yes_yes': df = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],'AST':y1[4],'BUP':y1[5],'BDO':y1[6],'PP2':y1[7]})
#	
	data_cut=5
	data_cut_thr='all'
	charm = 'no'
	if su2_charm=='no_yes': charm = 'yes' 
	elif su2_charm=='yes_yes': charm = 'yes' 
	
	mdl_den='pwr_lw_star'
		   
	lst =least_sq(data_cut, data_cut_thr,charm,mdl_den)
	lst.g_k_2h = 'PV17'
	lst.g_k_1h = 'PV17'

	if su2_charm=='no_yes': lst.su2 = 'no'
	elif su2_charm=='yes_yes': lst.su2 = 'yes' 


	lst.mdl_num = 'gauss'
	lst.mdl_den = 'pwr_lw_star'
	lst.coef=0.27
	
	
	now = datetime.now()
	current_time = now.strftime("%H:%M:%S")
	print("Current Time =", current_time)


	start = time.time()

	if su2_charm=='no_no': 
		a=zeros(len(df.NUP))
		i=0
		for nup,ndo,nst,nsea,ast,bup,bsea,pp2 in zip(df.NUP, df.NDO, df.NST, df.NSEA, df.AST,df.BUP, df.BSEA, df.PP2):
			a[i]=lst.least_squares_lh(nup,ndo,nst,nsea,0,0,ast,0,bup,0,0,bsea,pp2,0)
			df['chi']=a
			df.to_csv(r'fit_parameters/new_gauss_su_charm_'+str(su2_charm)+'.csv',index=False)
			i+=1
	elif su2_charm=='no_yes': 
		a=zeros(len(df.NUP))
		i=0
		for nup,ndo,nst,nsea,ado,ast,bup,bsea,pp2 in zip(df.NUP, df.NDO, df.NST, df.NSEA, df.ADO,df.AST,df.BUP, df.BSEA, df.PP2):
			a[i]=lst.least_squares_lh(nup,ndo,nst,nsea,0,ado,ast,0,bup,0,0,bsea,pp2,0)
			df['chi']=a
			df.to_csv(r'fit_parameters/new_gauss_su_charm_'+str(su2_charm)+'.csv',index=False)
			i+=1
	elif su2_charm=='yes_yes':
		a=zeros(len(df.NUP))
		i=0
		for nup,ndo,nst,nsea,ast,bup,bdo,pp2 in zip(df.NUP, df.NDO, df.NST, df.NSEA, df.AST,df.BUP, df.BDO , df.PP2):
			a[i]=lst.least_squares_lh(nup,ndo,nst,nsea,0,0,ast,0,bup,bdo,0,0,pp2,0)
			df['chi']=a
			df.to_csv(r'fit_parameters/new_gauss_su_charm_'+str(su2_charm)+'.csv',index=False)
			i+=1

		#print('gauss '+str(a[i])+' serie =' +str(num))
		
#	df['chi']=a    
#	df.to_csv(r'fit_parameters/new_bands_pv17/new_gauss_70k_set.csv',index=False)


	now = datetime.now()
	current_time = now.strftime("%H:%M:%S")
	print("Current Time =", current_time)	    

	end = time.time()
	mins=(end -start)/60
	print('¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯')
	print('time passed:')

	print(str(end - start) +'   ' + 'sec')
	print(str((end - start)/60) +'   ' + 'min')
		
	
	#return df		
	
#df1 =bands_gauss_2('no_no')
#df2 =bands_gauss_2('no_yes')
#df3 =bands_gauss_2('yes_yes')
#bands_gauss()
#bands_pwrlw()


		
p1 = multiprocessing.Process(target=bands_gauss_2,args=('no_no',))
p2 = multiprocessing.Process(target=bands_gauss_2,args=('no_yes',))
p3 = multiprocessing.Process(target=bands_gauss_2,args=('yes_yes',))


p1.start()
p2.start()
p3.start()


p1.join()
p2.join()
p3.join()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
