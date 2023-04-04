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


		mean=array([ 0.151, -0.135, -0.147, -0.096,  2.061, 3.835, 2.889, 0.104 ])

		cov= [[ 0.000755, -0.000192, -0.000159, -0.000368,     0.00218,   0.00582,  -0.00043,  0.000345], 
		 [-0.000192 , 0.000214,  -1.9e-05  , 0.00012,   -8.18e-05 , 0.000746, -0.000141, -0.000157], 
		 [-0.000159 , -1.9e-05,  0.000838 , 9.97e-05,    -0.00795 , -0.00143,  0.000234, -4.26e-05],
		 [-0.000368 ,  0.00012,  9.97e-05 , 0.000264,    -0.00126 ,-0.000236,  8.95e-05, -0.000244],
		 [  0.00218 ,-8.18e-05,  -0.00795 , -0.00126,      0.0832 ,   0.0188,  -0.00287,  0.000643],
		 [  0.00582 , 0.000746,  -0.00143 ,-0.000236,      0.0188 ,    0.152,  -0.00488, -6.26e-05],
		 [ -0.00043 ,-0.000141,  0.000234 , 8.95e-05,    -0.00287 , -0.00488,   0.00615, -0.000118],
		 [ 0.000345 ,-0.000157, -4.26e-05 ,-0.000244,    0.000643 ,-6.26e-05, -0.000118,  0.000263]]
	

	elif su2_charm=='no_yes':   ## 9 parametri chi_1.262 chi_min 109.8 ncall 837

		mean =array([ 0.178, -0.376, -0.121,  -0.127,  0.774,   0.848, 2.709, 1.589, 0.093 ])

		cov=  [[ 0.00572,  -0.00323, -0.000304,   -0.0046,   0.00732,  0.000235,    0.0658,     0.068,    0.0017],
		  [-0.00323 ,   0.0385,  0.000271 ,  0.00516,   -0.0919 , -0.00677,    -0.043 ,  -0.0961, -0.000606],
		 [-0.000304 , 0.000271,  0.000462 , 0.000284,  -0.00162 ,  -0.0046,  -0.00405 , -0.00351, -4.21e-05],
		 [  -0.0046 ,  0.00516,  0.000284 ,  0.00403,   -0.0122 ,-0.000801,   -0.0518 ,  -0.0606,   -0.0014],
		 [  0.00732 ,  -0.0919,  -0.00162 ,  -0.0122,     0.244 ,   0.0299,     0.101 ,    0.198,  0.000933],
		 [ 0.000235 , -0.00677,   -0.0046 ,-0.000801,    0.0299 ,   0.0549,   0.00918 ,-0.000948, -0.000658],
		 [ 0.0658   , -0.043  ,-0.00405   ,-0.0518  ,     0.101 ,  0.00918,      0.83 ,    0.824,    0.0176],
		 [    0.068 ,  -0.0961,  -0.00351 ,  -0.0606,     0.198 ,-0.000948,     0.824 ,     1.03,    0.0196],
		 [   0.0017 ,-0.000606, -4.21e-05 ,  -0.0014,  0.000933 ,-0.000658,    0.0176 ,   0.0196,  0.000594]]

	elif su2_charm=='yes_yes':   ## 8 parametri chi_1.447 chi_min 127.4 ncall 422


		mean =array([ 0.208, -0.128, -0.253, -0.165,  2.112, 2.438,  1.154, 0.117 ])

		cov= [[ 0.000473,  -0.00039, -0.000328, -0.000261,     0.00222,   0.00196,   0.00539,  0.000182],
		  [-0.00039,  0.000902 , 0.000387,  0.000285 ,   -0.00206,  0.000829 ,  -0.0136, -0.000209],
		 [-0.000328,  0.000387 ,  0.00344,  0.000231 ,    -0.0206, -0.000865 , -0.00657, -9.59e-05],
		 [-0.000261,  0.000285 , 0.000231,  0.000274 ,   -0.00153,   0.00198 , -0.00291, -0.000201],
		 [  0.00222,  -0.00206 ,  -0.0206,  -0.00153 ,      0.134,   0.00734 ,   0.0305,  0.000571],
		 [  0.00196,  0.000829 ,-0.000865,   0.00198 ,    0.00734,    0.0887 ,    0.018,  -0.00155],
		 [  0.00539,   -0.0136 , -0.00657,  -0.00291 ,     0.0305,     0.018 ,    0.259,   0.00207],
		 [ 0.000182, -0.000209 ,-9.59e-05, -0.000201 ,   0.000571,  -0.00155 ,  0.00207,  0.000175]]

	
	n_set = 2000
	     
#	np.random.seed(428)
	np.random.seed(np.random.randint(0,500))
	
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
	lst.coef=0.27
	
	
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

	
def bands_gauss_2(su2_charm,coef,num):
	#print(num)
	
	mean = array([0.168, -0.138 , -0.161, -0.104,  2.195,   4.023, 2.905, 0.103  ])

	if su2_charm=='no_no':   ## 8 parametri chi_1.174 chi_min 103,3 ncall 767


		mean=array([ 0.151, -0.135, -0.147, -0.096,  2.061, 3.835, 2.889, 0.104 ])

		cov= [[ 0.000755, -0.000192, -0.000159, -0.000368,     0.00218,   0.00582,  -0.00043,  0.000345], 
		 [-0.000192 , 0.000214,  -1.9e-05  , 0.00012,   -8.18e-05 , 0.000746, -0.000141, -0.000157], 
		 [-0.000159 , -1.9e-05,  0.000838 , 9.97e-05,    -0.00795 , -0.00143,  0.000234, -4.26e-05],
		 [-0.000368 ,  0.00012,  9.97e-05 , 0.000264,    -0.00126 ,-0.000236,  8.95e-05, -0.000244],
		 [  0.00218 ,-8.18e-05,  -0.00795 , -0.00126,      0.0832 ,   0.0188,  -0.00287,  0.000643],
		 [  0.00582 , 0.000746,  -0.00143 ,-0.000236,      0.0188 ,    0.152,  -0.00488, -6.26e-05],
		 [ -0.00043 ,-0.000141,  0.000234 , 8.95e-05,    -0.00287 , -0.00488,   0.00615, -0.000118],
		 [ 0.000345 ,-0.000157, -4.26e-05 ,-0.000244,    0.000643 ,-6.26e-05, -0.000118,  0.000263]]
	

	elif su2_charm=='no_yes':   ## 9 parametri chi_1.262 chi_min 109.8 ncall 837

		mean =array([ 0.178, -0.376, -0.121,  -0.127,  0.774,   0.848, 2.709, 1.589, 0.093 ])

		cov=  [[ 0.00572,  -0.00323, -0.000304,   -0.0046,   0.00732,  0.000235,    0.0658,     0.068,    0.0017],
		  [-0.00323 ,   0.0385,  0.000271 ,  0.00516,   -0.0919 , -0.00677,    -0.043 ,  -0.0961, -0.000606],
		 [-0.000304 , 0.000271,  0.000462 , 0.000284,  -0.00162 ,  -0.0046,  -0.00405 , -0.00351, -4.21e-05],
		 [  -0.0046 ,  0.00516,  0.000284 ,  0.00403,   -0.0122 ,-0.000801,   -0.0518 ,  -0.0606,   -0.0014],
		 [  0.00732 ,  -0.0919,  -0.00162 ,  -0.0122,     0.244 ,   0.0299,     0.101 ,    0.198,  0.000933],
		 [ 0.000235 , -0.00677,   -0.0046 ,-0.000801,    0.0299 ,   0.0549,   0.00918 ,-0.000948, -0.000658],
		 [ 0.0658   , -0.043  ,-0.00405   ,-0.0518  ,     0.101 ,  0.00918,      0.83 ,    0.824,    0.0176],
		 [    0.068 ,  -0.0961,  -0.00351 ,  -0.0606,     0.198 ,-0.000948,     0.824 ,     1.03,    0.0196],
		 [   0.0017 ,-0.000606, -4.21e-05 ,  -0.0014,  0.000933 ,-0.000658,    0.0176 ,   0.0196,  0.000594]]

	elif su2_charm=='yes_yes':   ## 8 parametri chi_1.447 chi_min 127.4 ncall 422


		mean =array([ 0.208, -0.128, -0.253, -0.165,  2.112, 2.438,  1.154, 0.117 ])

		cov= [[ 0.000473,  -0.00039, -0.000328, -0.000261,     0.00222,   0.00196,   0.00539,  0.000182],
		  [-0.00039,  0.000902 , 0.000387,  0.000285 ,   -0.00206,  0.000829 ,  -0.0136, -0.000209],
		 [-0.000328,  0.000387 ,  0.00344,  0.000231 ,    -0.0206, -0.000865 , -0.00657, -9.59e-05],
		 [-0.000261,  0.000285 , 0.000231,  0.000274 ,   -0.00153,   0.00198 , -0.00291, -0.000201],
		 [  0.00222,  -0.00206 ,  -0.0206,  -0.00153 ,      0.134,   0.00734 ,   0.0305,  0.000571],
		 [  0.00196,  0.000829 ,-0.000865,   0.00198 ,    0.00734,    0.0887 ,    0.018,  -0.00155],
		 [  0.00539,   -0.0136 , -0.00657,  -0.00291 ,     0.0305,     0.018 ,    0.259,   0.00207],
		 [ 0.000182, -0.000209 ,-9.59e-05, -0.000201 ,   0.000571,  -0.00155 ,  0.00207,  0.000175]]

	

	
	n_set = 100
	     
#	np.random.seed(428)
	np.random.seed(np.random.randint(0,200))
	y = np.random.multivariate_normal(mean, cov, n_set).T
	
	np.random.seed(np.random.randint(0,500))
	y2 = np.random.multivariate_normal(mean, cov, n_set).T

	np.random.seed(np.random.randint(0,500))
	y3 = np.random.multivariate_normal(mean, cov, n_set).T

	np.random.seed(np.random.randint(0,500))
	y4 = np.random.multivariate_normal(mean, cov, n_set).T

	np.random.seed(np.random.randint(0,500))
	y5 = np.random.multivariate_normal(mean, cov, n_set).T

	y1=np.concatenate((y,y2,y3,y4,y5),axis=1)
	
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

	lst.su2='no'
	lst.nf=3
	
	if su2_charm=='no_yes': 
		lst.su2 = 'no'
		lst.nf=4
	elif su2_charm=='yes_yes': 
		lst.su2 = 'yes' 
		lst.nf=4

	lst.mdl_num = 'gauss'
	lst.mdl_den = 'pwr_lw_star'
	lst.coef=coef
	lst.correct='no'
	
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
			i+=1
		df.to_csv(r'fit_parameters/dfs_to_use/new_gauss_su_charm_'+str(su2_charm)+'_'+str(num)+'.csv',index=False)
				
	elif su2_charm=='no_yes': 
		a=zeros(len(df.NUP))
		i=0
		for nup,ndo,nst,nsea,ado,ast,bup,bsea,pp2 in zip(df.NUP, df.NDO, df.NST, df.NSEA, df.ADO,df.AST,df.BUP, df.BSEA, df.PP2):
			a[i]=lst.least_squares_lh(nup,ndo,nst,nsea,0,ado,ast,0,bup,0,0,bsea,pp2,0)
			df['chi']=a
			i+=1
		df.to_csv(r'fit_parameters/dfs_to_use/new_gauss_su_charm_'+str(su2_charm)+'_'+str(num)+'.csv',index=False)
				
	elif su2_charm=='yes_yes':
		a=zeros(len(df.NUP))
		i=0
		for nup,ndo,nst,nsea,ast,bup,bdo,pp2 in zip(df.NUP, df.NDO, df.NST, df.NSEA, df.AST,df.BUP, df.BDO , df.PP2):
			a[i]=lst.least_squares_lh(nup,ndo,nst,nsea,0,0,ast,0,bup,bdo,0,0,pp2,0)
			df['chi']=a
			i+=1
		df.to_csv(r'fit_parameters/dfs_to_use/new_gauss_su_charm_'+str(su2_charm)+'_'+str(num)+'.csv',index=False)
				

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


num=7		
p1 = multiprocessing.Process(target=bands_gauss_2,args=('no_no',0.27,num))
p2 = multiprocessing.Process(target=bands_gauss_2,args=('no_yes',0.27,num))
p3 = multiprocessing.Process(target=bands_gauss_2,args=('yes_yes',0.27,num))


p1.start()
p2.start()
p3.start()


p1.join()
p2.join()
p3.join()
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
