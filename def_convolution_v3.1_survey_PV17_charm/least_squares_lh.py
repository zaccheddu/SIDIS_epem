# -*- coding: utf-8 -*-


from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
import pandas as pd
from iminuit import Minuit
from def_conv_crs_2h import*
from def_conv_crs_1h import*



class least_sq:

	def __init__(self,data_cut, data_cut_thr,charm,mdl_den):

		self.data_cut = data_cut
		self.data_cut_thr = data_cut_thr
		

		self.coef = 0.27


		self.unp_wd = 0.2
		self.mm = 1.115   #massa lambda fare richiamo

		self.f2 = 'dss'   # scelta della frammentazione

### double hadron data


		if data_cut == 0:
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			unpol = pd.read_csv('unpol_dfs/full.csv')  # tutti i punti
			#dati_lh2['P_exp']= dati_lh2['P_exp']*0.877
		elif data_cut == 1: dati_lh2=pd.read_csv("/exp_data/lambda_had_global_zkzj.dat", delimiter=r"\s+", header=0, engine='python') # taglio i punti di zk e zl >0.5

		#dati_lh2=dati_lh2.loc[(dati_lh2['z1']>0.3)]
		

		elif data_cut==2:  # tiene i punti per cui zk < 0.5
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[~(dati_lh2['z2']>0.5) | ~(dati_lh2['h2']==200)]
			dati_lh2=dati_lh2.loc[~(dati_lh2['z2']>0.5) | ~(dati_lh2['h2']==205)]
		elif data_cut==3: # tiene i punti per cui zl < 0.5
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[~(dati_lh2['z1']>0.5) | ~(dati_lh2['h2']==200)]
			dati_lh2=dati_lh2.loc[~(dati_lh2['z1']>0.5) | ~(dati_lh2['h2']==205)]
		elif data_cut==4:
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[~(dati_lh2['z1']>0.5) ]
		elif data_cut==5:
			unpol = pd.read_csv('unpol_dfs/full.csv')
			unpol = unpol.loc[~(unpol['z2']>0.5)]
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[~(dati_lh2['z2']>0.5) ]
			
			#dati_lh2['P_exp']= dati_lh2['P_exp']*0.877
		elif data_cut==6:
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[~(dati_lh2['z2']>0.5) ]
			dati_lh2=dati_lh2.loc[~(dati_lh2['z1']>0.5) ]
		#
		#PIONS
		elif data_cut=='pions':
			unpol = pd.read_csv('unpol_dfs/full.csv')
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[(dati_lh2['h2']==100)|(unpol['h2']==105)]
			unpol=unpol.loc[(unpol['h2']==100)|(unpol['h2']==105)]

		elif data_cut=='pions_cut':
			unpol = pd.read_csv('unpol_dfs/full.csv')
			unpol = unpol.loc[~(unpol['z2']>0.5)]
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[~(dati_lh2['z2']>0.5) ]
			dati_lh2=dati_lh2.loc[(dati_lh2['h2']==100)|(unpol['h2']==105)]
			unpol=unpol.loc[(unpol['h2']==100)|(unpol['h2']==105)]

		elif data_cut=='pions+':
			unpol = pd.read_csv('unpol_dfs/full.csv')
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[(dati_lh2['h2']==100)]
			unpol=unpol.loc[(unpol['h2']==100)]

		elif data_cut=='pions+_cut':
			unpol = pd.read_csv('unpol_dfs/full.csv')
			unpol = unpol.loc[~(unpol['z2']>0.5)]
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[~(dati_lh2['z2']>0.5) ]
			dati_lh2=dati_lh2.loc[(dati_lh2['h2']==100)]
			unpol=unpol.loc[(unpol['h2']==100)]
		elif data_cut=='pions-':
			unpol = pd.read_csv('unpol_dfs/full.csv')
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[(dati_lh2['h2']==105)]
			unpol=unpol.loc[(unpol['h2']==105)]
		elif data_cut=='pions-_cut':
			unpol = pd.read_csv('unpol_dfs/full.csv')
			unpol = unpol.loc[~(unpol['z2']>0.5)]
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[~(dati_lh2['z2']>0.5) ]
			dati_lh2=dati_lh2.loc[(dati_lh2['h2']==105)]
			unpol=unpol.loc[(unpol['h2']==105)]
#			
#
#
		#KAONS
		elif data_cut=='kaons':
			unpol = pd.read_csv('unpol_dfs/full.csv')
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[(dati_lh2['h2']==200)|(unpol['h2']==205)]
			unpol=unpol.loc[(unpol['h2']==200)|(unpol['h2']==205)]
		elif data_cut=='kaons_cut':
			unpol = pd.read_csv('unpol_dfs/full.csv')
			unpol = unpol.loc[~(unpol['z2']>0.5)]
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[~(dati_lh2['z2']>0.5) ]
			dati_lh2=dati_lh2.loc[(dati_lh2['h2']==200)|(unpol['h2']==205)]
			unpol=unpol.loc[(unpol['h2']==200)|(unpol['h2']==205)]
		
		elif data_cut=='kaons+':
			unpol = pd.read_csv('unpol_dfs/full.csv')
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[(dati_lh2['h2']==200)]
			unpol=unpol.loc[(unpol['h2']==200)]
		elif data_cut=='kaons+_cut':
			unpol = pd.read_csv('unpol_dfs/full.csv')
			unpol = unpol.loc[~(unpol['z2']>0.5)]
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[~(dati_lh2['z2']>0.5) ]
			dati_lh2=dati_lh2.loc[(dati_lh2['h2']==200)]
			unpol=unpol.loc[(unpol['h2']==200)]
		elif data_cut=='kaons-':
			unpol = pd.read_csv('unpol_dfs/full.csv')
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[(dati_lh2['h2']==205)]
			unpol=unpol.loc[(unpol['h2']==205)]
		elif data_cut=='kaons-_cut':
			unpol = pd.read_csv('unpol_dfs/full.csv')
			unpol = unpol.loc[~(unpol['z2']>0.5)]
			dati_lh2=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')
			dati_lh2=dati_lh2.loc[~(dati_lh2['z2']>0.5) ]
			dati_lh2=dati_lh2.loc[(dati_lh2['h2']==205)]
			unpol=unpol.loc[(unpol['h2']==205)]




		self.dati_h2=dati_lh2	
		dati_lh2.astype(float)

#		self.hh1=dati_lh2["h1"].get_values()
#		self.hh2=dati_lh2["h2"].get_values()
#		self.z1=dati_lh2["z1"].get_values()
#		self.z2=dati_lh2["z2"].get_values()
#		self.p_exp=dati_lh2["P_exp"].get_values()
#		self.p_the=dati_lh2["P(theo)"].get_values()
#		self.err=dati_lh2["err"].get_values()

		self.hh1=dati_lh2["h1"].to_numpy()
		self.hh2=dati_lh2["h2"].to_numpy()
		self.z1=dati_lh2["z1"].to_numpy()
		self.z2=dati_lh2["z2"].to_numpy()
		self.p_exp=dati_lh2["P_exp"].to_numpy()
		self.p_the=dati_lh2["P(theo)"].to_numpy()
		self.err=dati_lh2["err"].to_numpy()


## thrust (1 hadron) data

		if data_cut_thr=='all' :
			ml = 1.115
			data = pd.read_csv("exp_data/lambda_jet_global.dat", delimiter=r"\s+", header=0, engine='python') 
			#data = data.loc[~(data['z1']<0.3 )] 
			data['zp'] = data['z1']*sqrt(1- 4*ml**2/10.58**2/data['z1']**2)
			data['qt'] =data['pt']/data['zp']
			data = data.loc[~(data['P_exp']>0.001 )]
			#data['pt/zq'] = data['pt']/data['z1']/10.58 
			#data = data.loc[~(data['pt/zq']> 0.28 )] 
			#data = data.loc[~(data['qt']>10.58*self.coef )]

		if data_cut_thr=='cut_lwz' :
			ml = 1.115
			data = pd.read_csv("exp_data/lambda_jet_global.dat", delimiter=r"\s+", header=0, engine='python') 
			#data = data.loc[~(data['z1']<0.3 )] 
			data['zp'] = data['z1']#*sqrt(1- 4*ml**2/10.58**2/data['z1']**2)
			data['qt'] =data['pt']/data['zp']
			data = data.loc[~(data['P_exp']>0.05 )]  #0.001
			#data['pt/zq'] = data['pt']/data['z1']/10.58 
			#data = data.loc[~(data['pt/zq']> 0.28 )] 
			#data = data.loc[~(data['qt']>10.58*self.coef )]
			data = data.loc[~(data['z1']< 0.299)]
		if data_cut_thr=='cut_hgz' :
			ml = 1.115
			data = pd.read_csv("exp_data/lambda_jet_global.dat", delimiter=r"\s+", header=0, engine='python') 
			#data = data.loc[~(data['z1']<0.3 )] 
			data['zp'] = data['z1']#*sqrt(1- 4*ml**2/10.58**2/data['z1']**2)
			data['qt'] =data['pt']/data['zp']
			data = data.loc[~(data['P_exp']>0.05 )]  #0.001
			#data['pt/zq'] = data['pt']/data['z1']/10.58 
			#data = data.loc[~(data['pt/zq']> 0.28 )] 
			#data = data.loc[~(data['qt']>10.58*self.coef )]
			data = data.loc[~(data['z1']> 0.5)]

		if data_cut_thr=='1bin' :
			ml = 1.115
			data = pd.read_csv("exp_data/lambda_jet_global.dat", delimiter=r"\s+", header=0, engine='python') 
			#data = data.loc[~(data['z1']<0.3 )] 
			data['zp'] = data['z1']#*sqrt(1- 4*ml**2/10.58**2/data['z1']**2)
			data['qt'] =data['pt']/data['zp']
			data = data.loc[~(data['P_exp']>0.05 )]  #0.001
			#data['pt/zq'] = data['pt']/data['z1']/10.58 
			#data = data.loc[~(data['pt/zq']> 0.28 )] 
			#data = data.loc[~(data['qt']>10.58*self.coef )]
			data = data.loc[~(data['z1']> 0.299)]
		if data_cut_thr=='2bin' :
			ml = 1.115
			data = pd.read_csv("exp_data/lambda_jet_global.dat", delimiter=r"\s+", header=0, engine='python') 
			#data = data.loc[~(data['z1']<0.3 )] 
			data['zp'] = data['z1']#*sqrt(1- 4*ml**2/10.58**2/data['z1']**2)
			data['qt'] =data['pt']/data['zp']
			data = data.loc[~(data['P_exp']>0.05 )]  #0.001
			#data['pt/zq'] = data['pt']/data['z1']/10.58 
			#data = data.loc[~(data['pt/zq']> 0.28 )] 
			#data = data.loc[~(data['qt']>10.58*self.coef )]
			data = data.loc[(data['z1']> 0.299) & (data['z1']< 0.399)]
		if data_cut_thr=='3bin' :
			ml = 1.115
			data = pd.read_csv("exp_data/lambda_jet_global.dat", delimiter=r"\s+", header=0, engine='python') 
			#data = data.loc[~(data['z1']<0.3 )] 
			data['zp'] = data['z1']#*sqrt(1- 4*ml**2/10.58**2/data['z1']**2)
			data['qt'] =data['pt']/data['zp']
			data = data.loc[~(data['P_exp']>0.05 )]  #0.001
			#data['pt/zq'] = data['pt']/data['z1']/10.58 
			#data = data.loc[~(data['pt/zq']> 0.28 )] 
			#data = data.loc[~(data['qt']>10.58*self.coef )]
			data = data.loc[~(data['z1']> 0.5)]
			data = data.loc[(data['z1']> 0.399) & (data['z1']< 0.499)]
		if data_cut_thr=='4bin' :
			ml = 1.115
			data = pd.read_csv("exp_data/lambda_jet_global.dat", delimiter=r"\s+", header=0, engine='python') 
			#data = data.loc[~(data['z1']<0.3 )] 
			data['zp'] = data['z1']#*sqrt(1- 4*ml**2/10.58**2/data['z1']**2)
			data['qt'] =data['pt']/data['zp']
			data = data.loc[~(data['P_exp']>0.05 )]  #0.001
			#data['pt/zq'] = data['pt']/data['z1']/10.58 
			#data = data.loc[~(data['pt/zq']> 0.28 )] 
			#data = data.loc[~(data['qt']>10.58*self.coef )]
			data = data.loc[~(data['z1']< 0.5)]


		self.had1   =  data['h1'].to_numpy()
		self.zh    =  data['z1'].to_numpy()
		self.pt    =  data['pt'].to_numpy()
		self.p_exp_h1 =  data['P_exp'].to_numpy()
		self.err_h1   =  data['err'].to_numpy()


		self.g_k_2h = 'PV17'
		self.g_k_1h = 'PV17'

		self.mdl_num = 'gauss'	#pwr_lw_star	
		self.mdl_den = 'pwr_lw_star'		
		self.su2='no'	

		self.charm =None
		self.scale = 10.58
		#('unpol_dfs/dataframe_pwr_nocharm.csv')

			
		if mdl_den == 'pwr_lw_star': 
			if charm == 'no':  self.unp_cr = unpol['unpol1'].to_numpy()
			elif charm == 'yes': self.unp_cr = unpol['unpol2'].to_numpy()
		elif mdl_den == 'gauss': 
			if charm == 'no':  self.unp_cr = unpol['unpol3'].to_numpy()
			elif charm == 'yes': self.unp_cr = unpol['unpol4'].to_numpy()
		
		
		self.dataframe = unpol
		self.nf=3 
		self.correct='no'
	def least_squares_lh(self,NUP,NDO,NST,NSEA,AUP,ADO,AST,ASEA,BUP,BDO,BST,BSEA,PP,MSS):  #chi-square lambda-had
		#print(NUP)
		#print(self.unp_cr)

		fnc = polarization(self.coef)
		fnc.scale = self.scale
		fnc.mass = self.mm 
		fnc.frag2 = self.f2
		fnc.mdl_num = self.mdl_num 
		fnc.mdl_den = self.mdl_den 

		fnc.g_k = self.g_k_2h 
		fnc.charm=self.charm
		fnc.nf=self.nf


		#NUP.astype(float)

		f_prm=arange(0.,18.,1.)
		if self.su2=='yes':
			f_prm[0] = NUP # par[0]   #up
			f_prm[1] = NUP #par[0]  #do
			f_prm[2] = NST #par[0]      #st
			f_prm[3] = NSEA #par[0]      #upb
			f_prm[4] = NSEA #par[0]      #dob
			f_prm[5] = NDO # par[0]      #stb
		#¯¯¯¯¯¯
			f_prm[6] = AUP #par[0]      #aup
			f_prm[7] = AUP #par[0]      #ado
			f_prm[8] = AST # par[0]      #ast
			f_prm[9] = ASEA# par[0]      #aupb
			f_prm[10] =ASEA #par[0]      #adob
			f_prm[11] = ADO #par[0]      #astb
		#¯#¯¯¯¯¯¯¯
			f_prm[12] = BUP #par[0]      #bup
			f_prm[13] = BUP #par[0]      #bdo
			f_prm[14] = BST # par[0]      #bst
			f_prm[15] = BSEA #par[0]      #bupb
			f_prm[16] = BSEA #par[0]      #bdob
			f_prm[17] = BDO # par[0]      #bstb
			pt_pp=PP	

		elif self.su2=='no':

			f_prm[0] = NUP # par[0]   #up
			f_prm[1] = NDO #par[0]  #do
			f_prm[2] = NST #par[0]      #st
			f_prm[3] = NSEA #par[0]      #upb
			f_prm[4] = f_prm[3] #par[0]      #dob
			f_prm[5] = f_prm[3] # par[0]      #stb
		#¯¯¯¯¯¯
			f_prm[6] = AUP #par[0]      #aup
			f_prm[7] = ADO #par[0]      #ado
			f_prm[8] = AST # par[0]      #ast
			f_prm[9] = ASEA# par[0]      #aupb
			f_prm[10] = f_prm[9] #par[0]      #adob
			f_prm[11] = f_prm[9] #par[0]      #astb
		#¯#¯¯¯¯¯¯¯
			f_prm[12] = BUP #par[0]      #bup
			f_prm[13] = BDO #par[0]      #bdo
			f_prm[14] = BST # par[0]      #bst
			f_prm[15] = BSEA #par[0]      #bupb
			f_prm[16] = f_prm[15] #par[0]      #bdob
			f_prm[17] = f_prm[15] # par[0]      #bstb
			pt_pp=PP	
			
			
		elif self.su2=='mix':
			f_prm[0] = NUP # par[0]   #up
			f_prm[1] = NUP #par[0]  #do
			f_prm[2] = NST #par[0]      #st
			f_prm[3] = NSEA #par[0]      #upb
			f_prm[4] = NSEA #par[0]      #dob
			f_prm[5] = NSEA # par[0]      #stb
		#¯¯¯¯¯¯
			f_prm[6] = AUP #par[0]      #aup
			f_prm[7] = AUP #par[0]      #ado
			f_prm[8] = AST # par[0]      #ast
			f_prm[9] = ASEA# par[0]      #aupb
			f_prm[10] =ASEA #par[0]      #adob
			f_prm[11] = ASEA #par[0]      #astb
		#¯#¯¯¯¯¯¯¯
			f_prm[12] = BUP #par[0]      #bup
			f_prm[13] = BUP #par[0]      #bdo
			f_prm[14] = BST # par[0]      #bst
			f_prm[15] = BSEA #par[0]      #bupb
			f_prm[16] = BSEA #par[0]      #bdob
			f_prm[17] = BSEA # par[0]      #bstb
			pt_pp=PP	

			
		pt_pl1=self.unp_wd
		pt_pl2=0.2

		zz1=self.z1
		zz2=self.z2

		num=zeros(len(zz1))
		den=zeros(len(zz1))
		#a=-.8
		#b=0.68

		ratio=zeros(len(zz1))
		faxx=zeros(len(zz1))
		pol2=zeros(len(zz1))
		
		polar=zeros(len(zz1))
		for i in range(len(zz1)):

			if self.hh1[i]== 300: had1='lbd'
			elif self.hh1[i]== 310: had1='lbd_b'
			if self.hh2[i]== 100: had2='pi+'
			elif self.hh2[i]== 105: had2='pi-'
			elif self.hh2[i]== 200: had2='k+'
			elif self.hh2[i]==  205: had2='k-'


			#num[i]=fnc.cross_sec2_polda(had1,had2,self.z1[i],self.z2[i],10.58,f_prm)
			#den[i]=fnc.cross_sec2(had1,had2,self.z1[i],self.z2[i],10.58)	

			
			#pol2[i] = fnc.ratio(had1,had2,self.z1[i],self.z2[i],f_prm,pt_pp,MSS)

			pol2[i] = fnc.numerator(had1,had2,self.z1[i],self.z2[i],f_prm,pt_pp,MSS)
			pol2[i] = pol2[i]/self.unp_cr[i] 



	#		ratio[i]=ratio_pol(had1,had2,zp1[i],z2[i],10.58,f_prm)
			#faxx[i]=fnc.fact_pol(self.z1[i],self.z2[i],pt_pl1,pt_pl2,pt_pp,10.58)
	#		polar[i]=polarisation(had1,had2,zp1[i],z2[i],10.58,pt_pl1,pt_pl2,pt_pp,f_prm)
			#print(pol2[i])
		#print(pol2)	
		#ratio=num/den	
		#pol2=ratio*faxx
		#print(pol2)
		fact=1.
		if self.correct=='yes':	fact=642/732
		exp_pnts=self.p_exp
		exp_pnts=exp_pnts*fact
		
		errors=self.err
		errors=errors*fact


		chi_sq= sum ((exp_pnts - pol2)**2/errors**2)
		#chi_sq=round(chi_sq,4)
		#print(chi_sq)
		#print(((exp_pnts - pol2)**2/errors**2))
		#print(self.dati_h2)
		#tr=((exp_pnts - pol2)**2/errors**2)
		return chi_sq 



	def least_squares_ljh(self,NUP,NDO,NST,NSEA,AUP,ADO,AST,ASEA,BUP,BDO,BST,BSEA,PP,MSS):  #chi-square lambda-had and lambda-thrust

		fnc = polarization(self.coef)
		fnc.mass = self.mm 
		fnc.frag2 = self.f2

		fnc.mdl_num = self.mdl_num 
		fnc.mdl_den = self.mdl_den 

		fnc.g_k = self.g_k_2h 

		f_prm=arange(0.,18.,1.)
		f_prm[0] = NUP # par[0]   #up
		f_prm[1] = NDO #par[0]  #do
		f_prm[2] = NST #par[0]      #st
		f_prm[3] = NSEA #par[0]      #upb
		f_prm[4] = f_prm[3] #par[0]      #dob
		f_prm[5] = f_prm[3] # par[0]      #stb
	#¯¯¯¯¯¯
		f_prm[6] = AUP #par[0]      #aup
		f_prm[7] = ADO #par[0]      #ado
		f_prm[8] = AST # par[0]      #ast
		f_prm[9] = ASEA# par[0]      #aupb
		f_prm[10] = f_prm[9] #par[0]      #adob
		f_prm[11] = f_prm[9] #par[0]      #astb
	#¯#¯¯¯¯¯¯¯
		f_prm[12] = BUP #par[0]      #bup
		f_prm[13] = BDO #par[0]      #bdo
		f_prm[14] = BST # par[0]      #bst
		f_prm[15] = BSEA #par[0]      #bupb
		f_prm[16] = f_prm[15] #par[0]      #bdob
		f_prm[17] = f_prm[15] # par[0]      #bstb
		pt_pp=PP	

		pt_pl1=self.unp_wd
		pt_pl2=0.2

		zz1=self.z1
		zz2=self.z2

		#num=zeros(len(zz1))
		#den=zeros(len(zz1))


		#ratio=zeros(len(zz1))
		#faxx=zeros(len(zz1))
		pol2=zeros(len(zz1))
		
		polar=zeros(len(zz1))
		for i in range(len(zz1)):

			if self.hh1[i]== 300: had1='lbd'
			elif self.hh1[i]== 310: had1='lbd_b'
			if self.hh2[i]== 100: had2='pi+'
			elif self.hh2[i]== 105: had2='pi-'
			elif self.hh2[i]== 200: had2='k+'
			elif self.hh2[i]==  205: had2='k-'


			#num[i]=fnc.cross_sec2_polda(had1,had2,self.z1[i],self.z2[i],10.58,f_prm)
			#den[i]=fnc.cross_sec2(had1,had2,self.z1[i],self.z2[i],10.58)	

			
			pol2[i] = fnc.ratio(had1,had2,self.z1[i],self.z2[i],f_prm,pt_pp,MSS)

		
## thrust
		
		fnc_1h = polarization_1h()
		fnc_1h.mass = self.mm
		fnc_1h.mdl_num = self.mdl_num 
		fnc_1h.mdl_den = self.mdl_den 


		fnc_1h.g_k = self.g_k_1h 
		
		zh1 = self.zh
		pol1 = zeros(len(zh1))
		
		for j in range(len(zh1)):
		
			if self.had1[j]== 300: had1='lbd'
			elif self.had1[j]== 310: had1='lbd_b'
			
			pol1[j] = fnc_1h.ratio(had1,self.zh[j],self.pt[j],f_prm,pt_pp,MSS)
			

		
		pol_exp = np.concatenate((self.p_exp,self.p_exp_h1))
		errors = np.concatenate((self.err,self.err_h1))
		
		pol_theo = np.concatenate((pol2,pol1))
		
		chi_sq = sum ((pol_exp - pol_theo)**2/errors**2)

		return chi_sq

	def least_squares_lj(self,NUP,NDO,NST,NSEA,AUP,ADO,AST,ASEA,BUP,BDO,BST,BSEA,PP,MSS):  #chi-square lambda-had and lambda-thrust


		f_prm=arange(0.,18.,1.)
		f_prm[0] = NUP # par[0]   #up
		f_prm[1] = NDO #par[0]  #do
		f_prm[2] = NST #par[0]      #st
		f_prm[3] = NSEA #par[0]      #upb
		f_prm[4] = f_prm[3] #par[0]      #dob
		f_prm[5] = f_prm[3] # par[0]      #stb
	#¯¯¯¯¯¯
		f_prm[6] = AUP #par[0]      #aup
		f_prm[7] = ADO #par[0]      #ado
		f_prm[8] = AST # par[0]      #ast
		f_prm[9] = ASEA# par[0]      #aupb
		f_prm[10] = f_prm[9] #par[0]      #adob
		f_prm[11] = f_prm[9] #par[0]      #astb
	#¯#¯¯¯¯¯¯¯
		f_prm[12] = BUP #par[0]      #bup
		f_prm[13] = BDO #par[0]      #bdo
		f_prm[14] = BST # par[0]      #bst
		f_prm[15] = BSEA #par[0]      #bupb
		f_prm[16] = f_prm[15] #par[0]      #bdob
		f_prm[17] = f_prm[15] # par[0]      #bstb
		pt_pp=PP	

		pt_pl1=self.unp_wd
		pt_pl2=0.2

		
## thrust
		
		fnc_1h = polarization_1h()
		fnc_1h.mass = self.mm
		fnc_1h.mdl_num = self.mdl_num 
		fnc_1h.mdl_den = self.mdl_den 


		fnc_1h.g_k = self.g_k_1h 
		
		zh1 = self.zh
		pol1 = zeros(len(zh1))
		
		for j in range(len(zh1)):
		
			if self.had1[j]== 300: had1='lbd'
			elif self.had1[j]== 310: had1='lbd_b'
			
			pol1[j] = fnc_1h.ratio(had1,self.zh[j],self.pt[j],f_prm,pt_pp,MSS)
			

		
		pol_exp = self.p_exp_h1
		errors = self.err_h1
		
		pol_theo = pol1
		
		chi_sq = sum ((pol_exp - pol_theo)**2/errors**2)
		#print(chi_sq)
		return chi_sq





