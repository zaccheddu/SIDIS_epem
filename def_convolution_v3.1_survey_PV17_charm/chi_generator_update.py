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
	
		mean=array([ 0.19191584, -0.16685054, -0.16756225, -0.1225011,   0., 0.,  1.99427729,  0.,   3.93156437,  0.,0.,2.99347891,  0.09985722,  0. ])
	
		cov=[[   0.00406, -0.000527, -0.000276 , -0.00283 ,        0   ,      0,  0.000676 ,        0,    0.0508 ,        0  ,       0 ,   0.0538,  0.000934,         0], 
		 [-0.000527 , 0.000337, -8.96e-05 , 0.000452  ,       0    ,     0 ,  0.00103  ,       0 , -0.00487  ,       0   ,      0 ,  -0.0101 ,-0.000189 ,        0 ],
		 [-0.000276 ,-8.96e-05,   0.00113 , 0.000167  ,       0    ,     0 , -0.00932  ,       0 , -0.00288  ,       0   ,      0 ,-0.000203 , -3.6e-05 ,        0 ],
		 [ -0.00283 , 0.000452,  0.000167 ,  0.00218  ,       0    ,     0 ,  1.3e-05  ,       0 ,  -0.0333  ,       0   ,      0  , -0.0417 ,-0.000733 ,        0 ],
		  [       0 ,        0,         0  ,       0  ,       0    ,     0 ,        0  ,       0 ,        0  ,       0   ,      0  ,       0  ,       0 ,        0 ],
		   [      0 ,        0,         0  ,       0  ,       0    ,     0 ,        0  ,       0 ,        0  ,       0   ,      0  ,       0  ,       0 ,        0 ],
		  [0.000676 ,  0.00103,  -0.00932  , 1.3e-05  ,       0    ,     0 ,   0.0854  ,       0 ,  0.00138  ,       0   ,      0  , -0.0332 ,-0.000141 ,        0 ],
		   [      0 ,        0,         0  ,       0  ,       0    ,     0 ,        0  ,       0 ,        0  ,       0   ,      0  ,       0 ,        0 ,        0 ],
		   [ 0.0508 , -0.00487,  -0.00288  , -0.0333  ,       0     ,    0 ,  0.00138  ,       0 ,    0.739  ,       0   ,      0  ,   0.713 ,  0.00935 ,        0 ],
		    [     0 ,        0,         0  ,       0  ,       0     ,    0 ,        0  ,       0 ,        0  ,       0   ,      0  ,       0 ,        0 ,        0 ],
		    [    0  ,       0 ,        0   ,      0   ,      0      ,   0  ,       0   ,      0  ,       0   ,      0    ,     0   ,      0  ,       0  ,       0 ],
		    [.0538  , -0.0101, -0.000203  , -0.0417   ,      0     ,    0  , -0.0332   ,      0  ,   0.713   ,      0    ,     0   ,  0.923  ,  0.0123  ,       0 ],
		  [0.000934 ,-0.000189,  -3.6e-05, -0.000733  ,       0   ,      0 ,-0.000141  ,       0 ,  0.00935  ,       0   ,      0  ,  0.0123 , 0.000306 ,        0 ],
		[0        , 0        , 0,         0,         0 ,        0 ,        0,         0,         0,         0,         0 ,        0,         0,         0 ]]


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

	
def bands_gauss_2(su2_charm,coef):
	#print(num)
	
	mean = array([0.168, -0.138 , -0.161, -0.104,  2.195,   4.023, 2.905, 0.103  ])

	if su2_charm=='no_no':   ## 8 parametri chi_1.174 chi_min 103,3 ncall 767
	
		mean=array([ 0.19191584, -0.16685054, -0.16756225, -0.1225011,   0., 0.,  1.99427729,  0.,   3.93156437,  0.,0.,2.99347891,  0.09985722,  0. ])
	
		cov=[[   0.00406, -0.000527, -0.000276 , -0.00283 ,        0   ,      0,  0.000676 ,        0,    0.0508 ,        0  ,       0 ,   0.0538,  0.000934,         0], 
		 [-0.000527 , 0.000337, -8.96e-05 , 0.000452  ,       0    ,     0 ,  0.00103  ,       0 , -0.00487  ,       0   ,      0 ,  -0.0101 ,-0.000189 ,        0 ],
		 [-0.000276 ,-8.96e-05,   0.00113 , 0.000167  ,       0    ,     0 , -0.00932  ,       0 , -0.00288  ,       0   ,      0 ,-0.000203 , -3.6e-05 ,        0 ],
		 [ -0.00283 , 0.000452,  0.000167 ,  0.00218  ,       0    ,     0 ,  1.3e-05  ,       0 ,  -0.0333  ,       0   ,      0  , -0.0417 ,-0.000733 ,        0 ],
		  [       0 ,        0,         0  ,       0  ,       0    ,     0 ,        0  ,       0 ,        0  ,       0   ,      0  ,       0  ,       0 ,        0 ],
		   [      0 ,        0,         0  ,       0  ,       0    ,     0 ,        0  ,       0 ,        0  ,       0   ,      0  ,       0  ,       0 ,        0 ],
		  [0.000676 ,  0.00103,  -0.00932  , 1.3e-05  ,       0    ,     0 ,   0.0854  ,       0 ,  0.00138  ,       0   ,      0  , -0.0332 ,-0.000141 ,        0 ],
		   [      0 ,        0,         0  ,       0  ,       0    ,     0 ,        0  ,       0 ,        0  ,       0   ,      0  ,       0 ,        0 ,        0 ],
		   [ 0.0508 , -0.00487,  -0.00288  , -0.0333  ,       0     ,    0 ,  0.00138  ,       0 ,    0.739  ,       0   ,      0  ,   0.713 ,  0.00935 ,        0 ],
		    [     0 ,        0,         0  ,       0  ,       0     ,    0 ,        0  ,       0 ,        0  ,       0   ,      0  ,       0 ,        0 ,        0 ],
		    [    0  ,       0 ,        0   ,      0   ,      0      ,   0  ,       0   ,      0  ,       0   ,      0    ,     0   ,      0  ,       0  ,       0 ],
		    [.0538  , -0.0101, -0.000203  , -0.0417   ,      0     ,    0  , -0.0332   ,      0  ,   0.713   ,      0    ,     0   ,  0.923  ,  0.0123  ,       0 ],
		  [0.000934 ,-0.000189,  -3.6e-05, -0.000733  ,       0   ,      0 ,-0.000141  ,       0 ,  0.00935  ,       0   ,      0  ,  0.0123 , 0.000306 ,        0 ],
		[0        , 0        , 0,         0,         0 ,        0 ,        0,         0,         0,         0,         0 ,        0,         0,         0 ]]


	elif su2_charm=='no_yes':   ## 9 parametri chi_1.262 chi_min 109.8 ncall 837

		mean=array([ 0.20147387, -0.2367614,  -0.1202825,  -0.12543309,  0., 0.32660506, 0.78635502,  0.,3.02174084,  0., 0.,1.49666629,0.09772576,  0.])
	
		#mean =array([ 0.16940709, -0.36707107, -0.11713259, -0.12222033, 0.79289446,  0.86638731, 2.69376086,1.59579872, 0.09534966 ])

		cov= [[2.83e-05, -3.06e-07,  3.15e-05, -5.94e-08,         0, -8.89e-07, -0.000448,         0,  0.000538,         0,         0,  -2.2e-06,   8.9e-06,         0 ],
		 [-3.06e-07 , 5.75e-07, -3.77e-07 , 5.23e-10 ,        0 , 7.18e-09,   4.9e-06 ,        0, -5.83e-06 ,        0 ,        0 , 1.98e-08, -1.79e-07 ,        0 ],
		 [ 3.15e-05 ,-3.77e-07,   0.00041 ,-4.76e-08 ,        0 ,-1.04e-06,  -0.00421 ,        0,    0.0006 ,        0 ,        0 ,-2.11e-06,  4.35e-05 ,        0 ],
		 [-5.94e-08 , 5.23e-10, -4.76e-08 , 2.23e-08 ,        0 , 1.46e-09,  6.85e-07 ,        0, -1.13e-06 ,        0 ,        0 , 3.85e-09, -3.41e-08 ,        0 ],
		 [        0 ,        0,        0  ,       0  ,       0  ,       0 ,        0  ,       0 ,        0  ,       0  ,       0  ,       0 ,        0  ,       0 ],
		 [-8.89e-07 , 7.18e-09, -1.04e-06 , 1.46e-09 ,        0 , 5.82e-06,  1.42e-05 ,        0, -1.69e-05 ,        0 ,        0 , 5.71e-08, -5.89e-07 ,        0 ],
		 [-0.000448 ,  4.9e-06,  -0.00421 , 6.85e-07 ,        0 , 1.42e-05,    0.049  ,       0 , -0.00854  ,       0  ,       0  ,2.79e-05 ,-0.000615  ,       0 ],
		 [       0  ,       0 ,        0  ,       0 ,        0  ,       0 ,        0  ,       0 ,        0  ,       0  ,       0  ,       0 ,        0  ,       0 ],
		 [ 0.000538 ,-5.83e-06,    0.0006 ,-1.13e-06,         0 ,-1.69e-05,  -0.00854 ,        0,    0.0103 ,        0 ,        0 ,-4.19e-05 ,  0.00017 ,        0 ],
		 [       0  ,       0 ,       0   ,      0  ,       0   ,      0  ,      0    ,     0   ,      0    ,     0    ,     0    ,     0    ,     0    ,     0 ],
		 [       0  ,       0 ,        0  ,       0 ,        0  ,       0 ,        0  ,       0 ,        0  ,       0  ,       0  ,       0 ,        0  ,       0 ],
		 [ -2.2e-06 , 1.98e-08, -2.11e-06 , 3.85e-09,         0 , 5.71e-08,  2.79e-05 ,        0, -4.19e-05 ,        0 ,        0 , 1.83e-05, -1.05e-06 ,        0 ],
		 [  8.9e-06 ,-1.79e-07,  4.35e-05 ,-3.41e-08,         0 ,-5.89e-07, -0.000615 ,        0,   0.00017 ,        0 ,        0 ,-1.05e-06,  3.01e-05 ,        0 ],
		 [       0 ,        0,         0 ,        0,         0 ,        0,         0 ,        0,         0 ,        0 ,        0 ,        0,         0 ,        0 ]]
		'''
		cov= [[  0.00872,  -0.00625, -0.000475,  -0.00719,    0.0146,  0.000515,     0.105,     0.111,   0.00272],
		      [ -0.00625,    0.0408,  0.000423,   0.00764,   -0.0997,  -0.00692,    -0.081,    -0.138,   -0.0016],
		      [-0.000475,  0.000423,  0.000448,  0.000425,  -0.00196,  -0.00451,  -0.00624,  -0.00589, -0.000102],
		      [ -0.00719,   0.00764,  0.000425,   0.00624,   -0.0182, -0.000981,   -0.0852,   -0.0977,  -0.00228],
		      [   0.0146,   -0.0997,  -0.00196,   -0.0182,     0.268,    0.0307,     0.193,     0.302,   0.00329],
		      [ 0.000515,  -0.00692,  -0.00451, -0.000981,    0.0307,    0.0554,    0.0129,   0.00224, -0.000566],
		      [    0.105,    -0.081,  -0.00624,   -0.0852,     0.193,    0.0129,      1.33,      1.38,    0.0308],
		      [    0.111,    -0.138,  -0.00589,   -0.0977,     0.302,   0.00224,      1.38,      1.64,    0.0342],
		      [  0.00272,   -0.0016, -0.000102,  -0.00228,   0.00329, -0.000566,    0.0308,    0.0342,  0.000941]]	
		'''

	elif su2_charm=='yes_yes':   ## 8 parametri chi_1.447 chi_min 127.4 ncall 422

		mean =array([ 0.22219531, -0.13806872, -0.26662248, -0.17717025,  0.,0.,  2.10203422,  0., 2.42733213,  1.15821458,  0., 0.,  0.11528321,  0.        ])

		#mean =array([ 0.20754979, -0.12811842, -0.25272118, -0.16499243, 2.11244742, 2.43826013,  1.15411948, 0.11715731 ])

		cov= [[  0.000557, -0.000487 ,-0.000397, -0.000309 ,        0  ,       0,   0.00251 ,        0 ,  0.00203 ,  0.00589 ,        0 ,        0,  0.000204 ,        0 ],
		      [-0.000487 ,  0.00119 , 0.000487 , 0.000359  ,       0   ,      0 , -0.00243  ,       0 ,  0.00105 ,  -0.0157  ,       0  ,       0 ,-0.000248  ,       0 ],
		      [-0.000397 , 0.000487 ,  0.00401 , 0.000279  ,       0   ,      0 ,  -0.0226  ,       0 ,-0.000992 , -0.00713  ,       0  ,       0 ,-0.000112  ,       0 ],
		      [ -0.000309,  0.000359,  0.000279,  0.000322 ,        0 ,        0,  -0.00172 ,        0 ,   0.0021,  -0.00319 ,        0 ,        0, -0.000225 ,        0], 
		      [        0 ,        0 ,        0,         0  ,       0  ,       0 ,        0 ,        0  ,       0 ,        0  ,       0  ,       0 ,        0 ,        0], 
		      [       0  ,       0  ,       0 ,        0   ,      0   ,      0  ,       0  ,       0   ,      0  ,       0   ,      0   ,      0  ,       0  ,       0], 
		      [  0.00251 , -0.00243 ,  -0.0226,  -0.00172  ,       0   ,      0 ,    0.138 ,        0 ,  0.00773 ,   0.0313  ,       0  ,       0 , 0.000626 ,        0], 
		      [       0  ,       0  ,       0,         0   ,      0   ,      0  ,       0  ,       0   ,      0 ,        0   ,      0   ,      0  ,       0  ,       0], 
		      [ 0.00203 ,  0.00105 ,-0.000992,    0.0021   ,      0  ,       0  , 0.00773  ,       0 ,   0.0851 ,   0.0166   ,      0   ,      0 , -0.00156  ,       0], 
		      [ 0.00589 ,  -0.0157 , -0.00713,  -0.00319   ,      0 ,        0  ,  0.0313  ,       0 ,   0.0166 ,    0.261   ,      0  ,       0 ,  0.00215  ,       0], 
		      [       0 ,        0 ,        0,         0   ,      0,         0  ,       0  ,       0 ,        0 ,        0   ,      0 ,        0 ,        0  ,       0 ],
		      [     0   ,      0   ,      0  ,       0     ,    0  ,       0   ,      0    ,     0   ,      0   ,      0     ,    0   ,      0   ,      0   ,      0 ],
		      [0.000204, -0.000248 ,-0.000112, -0.000225   ,      0,         0 , 0.000626 ,        0 , -0.00156 ,  0.00215   ,      0 ,        0 , 0.000185,         0 ],
		      [       0 ,        0  ,       0,         0    ,     0 ,        0  ,       0 ,        0  ,       0 ,        0    ,     0  ,       0   ,      0 ,        0 ]]

		'''
		cov= [[0.000473,  -0.00039, -0.000328, -0.000261,   0.00222,   0.00196,   0.00539,  0.000182],
		      [-0.00039,  0.000902,  0.000387,  0.000285,  -0.00206,  0.000829,   -0.0136, -0.000209],
		      [-0.000328,  0.000387,   0.00344,  0.000231,   -0.0206, -0.000865,  -0.00657, -9.59e-05],
		      [-0.000261,  0.000285,  0.000231,  0.000274,  -0.00153,   0.00198,  -0.00291, -0.000201],
		      [ 0.00222,  -0.00206,   -0.0206,  -0.00153,     0.134,   0.00734,    0.0305,  0.000571],
		      [ 0.00196,  0.000829, -0.000865,   0.00198,   0.00734,    0.0887,     0.018,  -0.00155],
		      [ 0.00539,   -0.0136,  -0.00657,  -0.00291,    0.0305,     0.018,     0.259,   0.00207],
		      [0.000182, -0.000209, -9.59e-05, -0.000201,  0.000571,  -0.00155,   0.00207,  0.000175]]
		'''


	
	n_set = 2000
	     
#	np.random.seed(428)
	np.random.seed(np.random.random_integers(0,200))
	
	y1 = np.random.multivariate_normal(mean, cov, n_set).T
	
	if su2_charm=='no_no': df = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],'AST':y1[6],'BUP':y1[8],'BSEA':y1[11],'PP2':y1[12]})

	elif su2_charm=='no_yes': df = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],'ADO':y1[5],'AST':y1[6],'BUP':y1[8],'BSEA':y1[11],'PP2':y1[12]})
	
	elif su2_charm=='yes_yes': df = pd.DataFrame({'NUP':y1[0],'NDO':y1[1],'NST':y1[2],'NSEA':y1[3],'AST':y1[6],'BUP':y1[8],'BDO':y1[9],'PP2':y1[12]})
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
	lst.coef=coef
	
	
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
			df.to_csv(r'fit_parameters/dfs_2/new_gauss_su_charm_'+str(su2_charm)+'.csv',index=False)
			i+=1
	elif su2_charm=='no_yes': 
		a=zeros(len(df.NUP))
		i=0
		for nup,ndo,nst,nsea,ado,ast,bup,bsea,pp2 in zip(df.NUP, df.NDO, df.NST, df.NSEA, df.ADO,df.AST,df.BUP, df.BSEA, df.PP2):
			a[i]=lst.least_squares_lh(nup,ndo,nst,nsea,0,ado,ast,0,bup,0,0,bsea,pp2,0)
			df['chi']=a
			df.to_csv(r'fit_parameters/dfs_2/new_gauss_su_charm_'+str(su2_charm)+'.csv',index=False)
			i+=1
	elif su2_charm=='yes_yes':
		a=zeros(len(df.NUP))
		i=0
		for nup,ndo,nst,nsea,ast,bup,bdo,pp2 in zip(df.NUP, df.NDO, df.NST, df.NSEA, df.AST,df.BUP, df.BDO , df.PP2):
			a[i]=lst.least_squares_lh(nup,ndo,nst,nsea,0,0,ast,0,bup,bdo,0,0,pp2,0)
			df['chi']=a
			df.to_csv(r'fit_parameters/dfs_2/new_gauss_su_charm_'+str(su2_charm)+'.csv',index=False)
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


		
p1 = multiprocessing.Process(target=bands_gauss_2,args=('no_no',0.25))
p2 = multiprocessing.Process(target=bands_gauss_2,args=('no_yes',0.27))
p3 = multiprocessing.Process(target=bands_gauss_2,args=('yes_yes',0.27))


p1.start()
p2.start()
p3.start()


p1.join()
p2.join()
p3.join()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
