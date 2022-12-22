# -*- coding: utf-8 -*-
from pylab import*
import numpy as np
import lhapdf
import scipy.integrate as integrate
from Sudakov.evolve import evolve, sng
#from evolve import evolve, sng
import pycuba


class Soft:

	#CA= 3.
	#CF= 4/3
	#nf= 3	#active flavors
	#Tr = 1/2	
	

#	pdf=lhapdf.mkAlphaS("NNFF10_PIp_nlo")
	#pdf=lhapdf.mkPDF("CT10nlo")
	#pdf=lhapdf.mkPDF("CT14nnloIC")	

	def __init__(self,order,sep):


		self.CA= 3.
		self.CF= 4/3
		self.nf= 4 #3	#active flavors
		self.Tf = 1/2	

		self.sep= sep   # GeV
		
		#self.scale= 10.58  # GeV

		self.a = 0.05 ## GeV**2

		self.order = order		

	
		self.bmax = 0.6
		self.mu= 10.58

		self.a = 0.05 ## GeV**2
		

		#self.bmin = 2*e**(-euler_gamma)/self.scale	### for mu_b_star


	def bt_str_mu(self,bt,x,y): 				 #define bt_star for mu_b_star with bmin
	
		QQ = self.sep*x*y
		bmin = 2*e**(-euler_gamma)/QQ
		
		
		bm= self.bmax	
		
		b_new = sqrt(bt**2 + bmin**2)

		bstar = b_new/sqrt(1+b_new**2/bm**2)			
		#print('in bt_star:' + str(bt))		
		return bstar

	def bt_str_mu_pv17(self,bt,x,y): 				 #define bt_star for mu_b_star with bmin

		QQ = self.sep*x*y
		bmin = 2*e**(-euler_gamma)/QQ
	
		#print(bt)
		bm= self.bmax	
		bstar = bm

		#print('bstar '+str(bstar))
		if type(bt)==float:

			bstar = bstar*((1-np.exp(-bt**4/bm**4))**(1/4))/(1-np.exp(-bt**4/bmin**4))**(1/4)
			if bt==0: bstar=bmin

		#elif type(bt)==np.array:
		#	print('cacca')

		#print(type(bt))
		#else:
		#	if bt==0: bstar=bmin
		#print('bstar '+str((1-np.exp(-bt**4/bm**4))**(1/4)))

		#if bstar !=0. :bstar = bstar/(1-np.exp(-bt**4/self.bmin**4))**(1/4)
		#elif bstar ==0. : bstar = self.bmin
		#print('bstar '+str(bstar))
		
		return bstar



	def mu_b(self,bt,x,y): 				 #define mu_b

		#print(bt)
		#bt_str_mu
		bstr=self.bt_str_mu_pv17(bt,x,y)
		#print(bstr)
		#bstr=self.bt_str_mu(bt,x,y)

		
		mub=2*e**(-euler_gamma)/bstr
		#print(mub)
		#print('in mu_b:' + str(bt))
		#print('mub  ' + str(mub))
		#print(mub)
		return mub

###################### 		NO BMIN

	def bt_str(self,bt): 				 #define bt_star --- no bmin
	
		bm= self.bmax	
		
		#b_new =  sqrt(bt**2 + self.bmin**2)

		bstar = bt/sqrt(1+bt**2/bm**2)			
		#print('in bt_star:' + str(bt))		
		return bstar

	def mu_b_no(self,bt): 				 #define mu_b --- no bmin

		bstr=self.bt_str(bt)

		mub=2*e**(-euler_gamma)/bstr
		#print('in mu_b:' + str(bt))
		return mub



#
#	sSUDAKOV TERMS


######################################
#  ANALYTIC SUDAKOV

	def analytic_sudakov(self,bt,x,y):

		if self.nf==3:lbd_qcd = 0.2123 ## nf= 3
		elif self.nf==4:lbd_qcd = 0.1737 ## nf=4

		beta0 = (11*self.CA - 4*self.Tf*self.nf)/3  ## 4 pi beta0
		gammad1 = 6*self.CF 
		gammak1 = 8*self.CF
		gammak2 = self.CA*self.CF*(536/9 - 8*np.pi**2/3) - 80*self.CF*self.nf/9

		QQ = self.sep*x*y
		mu_bstr = self.mu_b(bt,x,y)
		#bmin = 2*e**(-euler_gamma)/QQ
	
		#lbd_qcd = 0.2123 
		#QQ = self.scale
		LL1 = np.log(QQ/lbd_qcd)
		LL2 = np.log(mu_bstr/lbd_qcd)
		LL3 = np.log(QQ/mu_bstr)
		
		fact = gammad1*np.log(LL1/LL2)/beta0 +	gammak1*(LL3 - LL1*np.log(LL1/LL2))/beta0 + gammak2*(-LL3/LL2 + np.log(LL1/LL2))/2/beta0**2		
		
		return np.exp(fact)	

	def analytic_sudakov_1h(self,bt):
		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 


			#b_new = sqrt(bb**2 + self.bmin**2)
		
			lbd_qcd = 0.2123 
			QQ = self.scale
			AA = 12*np.pi/(33 - 2*self.nf)
			
			mu_bstr = self.mu_b(bb)
			
			LL1 = np.log(QQ/lbd_qcd)
			LL2 = np.log(mu_bstr/lbd_qcd)
			
			fact = 2*AA/np.pi*( np.log(LL1/LL2) -4/3*LL1*np.log(LL1/LL2) + 4/3*np.log(QQ/mu_bstr) )
			#out = np.append(out,np.exp(fact)*sng(self.scale,b_new))
			#out = np.append(out,np.exp(fact)*self.non_glob(self.scale,bb))
			
			out = np.append(out,np.exp(fact)*self.sng_imp(self.scale,bb))

			

			#out = np.append(out,np.exp(fact))
			#out = np.append(out,1.)
		
		return out	
####################################################################################
#
#	g_K FUNCTIONS   ( with exp(-gk) )
#_____________________________________________________________________________

	def g_K(self,bt):	# LOG G_K
	
	
		b_new = sqrt(bt**2 + self.bmin**2)

		g2 = 0.84
		Q_0 = np.sqrt(2.4) # GeV
		
		esp =  g2*np.log(self.scale/Q_0)*np.log(b_new/self.bt_str_mu(bt))	
		out =  np.exp(-esp)
		
		return out #esp #

	def g_K_1h(self,bt):

		b_new = sqrt(bt**2 + self.bmin**2)
	
		g2 = 0.84
		Q_0 = np.sqrt(2.4) # GeV
		
		esp =  g2*np.log(self.scale/Q_0)*np.log(b_new/self.bt_str_mu(bt))	
		#print('fact =' + str(g2) + '; log = '+str(np.log(b_new/self.bt_str_mu(bt))) )
		#print('   ')
		#print('lgM_ =' +str(np.log(self.scale/Q_0)))
		#print('________________________')
		
		out = np.exp(-esp/2)
		
		return  out #esp #

#################################	gk log_b con lgm


	def g_K_lgm(self,bt,had2,z1,z2):	# LOG G_K
	
		ml1= 1.115 # lambda mass
		
		if had2 == 'pi+' or had2 == 'pi-' : ml2 = 0.13957  # pion mass
		elif had2 == 'k+' or had2 == 'k-' : ml2 = 0.4937  # kaon mass
		 
		
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale**2*zlc1*z2
		term = term/ml1/ml2

		logs = np.log(term)

	
		b_new = sqrt(bt**2 + self.bmin**2)

		g2 = 0.84
		Q_0 = np.sqrt(2.4) # GeV
		
		esp =  g2*logs*np.log(b_new/self.bt_str_mu(bt))	
		out =  np.exp(-esp)
		
		return out #esp #

	def g_K_1h_lgm(self,bt,z1):

		ml1= 1.115 # lambda mass
		
	
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale*zlc1
		term = term/ml1

		logs = np.log(term)

		b_new = sqrt(bt**2 + self.bmin**2)
	
		g2 = 0.84
		#Q_0 = np.sqrt(2.4) # GeV
		
		esp =  g2*logs*np.log(b_new/self.bt_str_mu(bt))	
		out = np.exp(-esp)
		
		return  out #esp #






############


	def g_K_blny(self,bt,had2,z1,z2):	# LOG G_K

		ml1= 1.115 # lambda mass
		
		if had2 == 'pi+' or had2 == 'pi-' : ml2 = 0.13957  # pion mass
		elif had2 == 'k+' or had2 == 'k-' : ml2 = 0.4937  # kaon mass
		 
		
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale**2*zlc1*z2
		term = term/ml1/ml2

		logs = np.log(term)

		b_new = sqrt(bt**2 + self.bmin**2)

		g2 = 0.68
		#Q_0 = 3.2 # GeV
		
		esp =  g2*b_new**2	
		esp =  esp*logs
		out =  np.exp(-esp)
		
		return out #esp #


	def g_K_blny_1h(self,bt,z1):	# LOG G_K

		ml1= 1.115 # lambda mass
		
	
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale*zlc1
		term = term/ml1

		logs = np.log(term)
	
	
		b_new = sqrt(bt**2 + self.bmin**2)

		g2 = 0.68
		#Q_0 = 3.2 # GeV
		
		esp =  g2*b_new**2
		esp =  esp*logs	
		out =  np.exp(-esp/2)
		
		return out #esp #


############# Bacchetta et al


	def g_K_bac_lgm(self,bt,z1,xb2,y):	# LOG G_K


		QQ = self.sep*xb2*y
		bmin = 2*e**(-euler_gamma)/QQ

		ml1= 1.115 # lambda mass
		
		ml2 = 0.938
		#if had2 == 'pi+' or had2 == 'pi-' : ml2 = 0.13957  # pion mass
		#elif had2 == 'k+' or had2 == 'k-' : ml2 = 0.4937  # kaon mass
		 
		eta_p= (1 - (ml1**2/z1**2/QQ**2)*xb2/(1-xb2))
		zp1 = z1*sqrt(eta_p)		# momentum fraction
		
		#eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		#zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = QQ**2*zlc1#*xb2
		term = term/ml1/ml2/xb2

		logs = np.log(term)	
	
		b_new = sqrt(bt**2 + bmin**2)

		g2 = 0.13
		
		esp =  g2*logs*b_new**2	
		out =  np.exp(-esp)
		
		return out #esp #


	def g_K_bac_1h_lgm(self,bt,z1):	# LOG G_K
	
		ml1= 1.115 # lambda mass
		
	
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale*zlc1
		term = term/ml1

		logs = np.log(term)
	
		b_new = sqrt(bt**2 + self.bmin**2)

		g2 = 0.13
		
		esp =  g2*logs*b_new**2	
		out =  np.exp(-esp/2)
		
		return out #esp #












###### missing log(Q/Q0)

	def g_K_new(self,bt): # QUADRATIC G_K

		out = array([])	
		res1 = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 
	
			b_new = sqrt(bb**2 + self.bmin**2)
	
			alpha = self.pdf.alphasQ(self.mu_b(bb))	
			
			res = self.CF*b_new**2*alpha/np.pi/self.bmax**2
			#res1 = np.append(res1,res)
			out = np.append(out,np.exp(-res))
		
		return out	

	def g_K_new_1h(self,bt):

		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			b_new = sqrt(bb**2 + self.bmin**2)
	
			alpha = self.pdf.alphasQ(self.mu_b(bb))	
			
			res = self.CF*b_new**2*alpha/np.pi/self.bmax**2/2
			out = np.append(out,np.exp(-res))
		
		return out	

	def g_K_new_lgm(self,bt,had2,z1,z2): # QUADRATIC G_K

		ml11= 1.115 # lambda mass
		ml1=0.
		
		if had2 == 'pi+' or had2 == 'pi-' : ml2 = 0.13957  # pion mass
		elif had2 == 'k+' or had2 == 'k-' : ml2 = 0.4937  # kaon mass
		 
		
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale**2*zlc1*z2
		term = term/ml11/ml2

		logs = np.log(term)


		out = array([])	
		res1 = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 
	
			b_new = sqrt(bb**2 + self.bmin**2)
	
			alpha = self.pdf.alphasQ(self.mu_b(bb))	
			
			res = self.CF*b_new**2*alpha/np.pi/self.bmax**2
			#res1 = np.append(res1,res)
			out = np.append(out,np.exp(-res*logs))
		
		return out	

	def g_K_new_1h_lgm(self,bt,z1):

		ml1= 1.115 # lambda mass
		
	
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale*zlc1
		term = term/ml1

		logs = np.log(term)


		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			b_new = sqrt(bb**2 + self.bmin**2)
	
			alpha = self.pdf.alphasQ(self.mu_b(bb))	
			
			res = self.CF*b_new**2*alpha/np.pi/self.bmax**2
			out = np.append(out,np.exp(-res*logs))
		
		return out	




########## ## logs collins rogers
	def g_K_ll(self,bt):
	
		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			b_new = sqrt(bb**2 + self.bmin**2)

			#alpha = self.pdf.alphasQ(2*e**(-euler_gamma)/self.bmax)
			alpha = self.pdf.alphasQ(self.mu_b(bb))
			
			res = self.CF*alpha/np.pi*np.log(1+ b_new**2/self.bmax**2)
			out = np.append(out,np.exp(-res))
	
		return out

	def g_K_ll_1h(self,bt):

		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			b_new = sqrt(bb**2 + self.bmin**2)
		
			alpha = self.pdf.alphasQ(self.mu_b(bb))

			#g2 = 0.84
			#Q_0 = np.sqrt(2.4) # GeV

			res = self.CF*alpha/np.pi*np.log(1+ b_new**2/self.bmax**2)/2
			#res = g2*np.log(self.scale/Q_0)**np.log(1+ b_new**2/self.bmax**2)/4
			out = np.append(out,np.exp(-res))
	
		return out


	def g_K_ll_lgm(self,bt,had2,z1,z2):

		ml1= 1.115 # lambda mass
		
		if had2 == 'pi+' or had2 == 'pi-' : ml2 = 0.13957  # pion mass
		elif had2 == 'k+' or had2 == 'k-' : ml2 = 0.4937  # kaon mass
		 
		
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale**2*zlc1*z2
		term = term/ml1/ml2

		logs = np.log(term)

	
		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			b_new = sqrt(bb**2 + self.bmin**2)

			#alpha = self.pdf.alphasQ(2*e**(-euler_gamma)/self.bmax)
			alpha = self.pdf.alphasQ(self.mu_b(bb))
			
			res = self.CF*alpha/np.pi*np.log(1+ b_new**2/self.bmax**2)
			out = np.append(out,np.exp(-res*logs))
	
		return out

	def g_K_ll_1h_lgm(self,bt,z1):

		ml1= 1.115 # lambda mass
		
	
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale*(zlc1)
		term = term/ml1

		logs = np.log(term)


		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			b_new = sqrt(bb**2 + self.bmin**2)
		
			alpha = self.pdf.alphasQ(self.mu_b(bb))

			#g2 = 0.84
			#Q_0 = np.sqrt(2.4) # GeV

			res = self.CF*alpha/np.pi*np.log(1+ b_new**2/self.bmax**2)
			#print('fact =' + str(self.CF*alpha/np.pi) + '; log = '+str(np.log(1+ b_new**2/self.bmax**2)) )
			#print('   ')
			#print('lgM_ =' +str(logs))
			#print('________________________')
			#res = g2*np.log(self.scale/Q_0)**np.log(1+ b_new**2/self.bmax**2)/4
			out = np.append(out,np.exp(-res*logs))
	
		return out
##################################################################################
	def g_K_ll_lgm_n(self,bt,had2,z1,z2):

		ml1= 1.115 # lambda mass
		
		if had2 == 'pi+' or had2 == 'pi-' : ml2 = 0.13957  # pion mass
		elif had2 == 'k+' or had2 == 'k-' : ml2 = 0.4937  # kaon mass
		 
		
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale**2*zlc1*z2
		term = term/ml1/ml2

		logs = np.log(term)

	
		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			b_new = sqrt(bb**2 + self.bmin**2)

			#alpha = self.pdf.alphasQ(2*e**(-euler_gamma)/self.bmax)
			alpha = self.pdf.alphasQ(self.mu_b(bb))
			res = self.CF*alpha/np.pi/self.bmax**2
			res = res*np.log(1+ b_new**2/self.bmax**2)
			out = np.append(out,np.exp(-res*logs))
	
		return out

	def g_K_ll_1h_lgm_n(self,bt,z1):

		ml1= 1.115 # lambda mass
		
	
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale*zlc1
		term = term/ml1

		logs = np.log(term)


		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			b_new = sqrt(bb**2 + self.bmin**2)
		
			alpha = self.pdf.alphasQ(self.mu_b(bb))

			#g2 = 0.84
			#Q_0 = np.sqrt(2.4) # GeV
			res = self.CF*alpha/np.pi/self.bmax**2
			#print('fact =' + str(res) + '; log = '+str(np.log(1+ b_new**2/self.bmax**2)) )
			#print('   ')
			res = res*np.log(1+ b_new**2/self.bmax**2)

			#print('lgM_ =' +str(logs))
			#print('________________________')
			#res = g2*np.log(self.scale/Q_0)**np.log(1+ b_new**2/self.bmax**2)/4
			out = np.append(out,np.exp(-res*logs))
	
		return out






####################___ : exp collins rogers
#

	def g_K_exp(self,bt):
	
		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			b_new = sqrt(bb**2 + self.bmin**2)
		
			alpha = self.pdf.alphasQ(self.mu_b(bb))
			
			g_0 = 0.55
			
			res = g_0*(1 - np.exp(- self.CF*alpha*b_new**2/np.pi/g_0/self.bmax**2)  )	
			#type(res)
			out = np.append(out,np.exp(-res))
			
		return out	
			
	def g_K_exp_1h(self,bt):

	
		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			b_new = sqrt(bb**2 + self.bmin**2)
		
			alpha = self.pdf.alphasQ(self.mu_b(bb))
			
			g_0 = 0.55
			
			res = g_0*(1 - np.exp(- self.CF*alpha*b_new**2/np.pi/g_0/self.bmax**2)  )/2	
			
			out = np.append(out,np.exp(-res))
			
		return out




####################___ : exp collins rogers
	def g_K_exp_lgm(self,bt,had2,z1,z2):

		ml1= 1.115 # lambda mass
		
		if had2 == 'pi+' or had2 == 'pi-' : ml2 = 0.13957  # pion mass
		elif had2 == 'k+' or had2 == 'k-' : ml2 = 0.4937  # kaon mass
		 
		
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale**2*zlc1*z2
		term = term/ml1/ml2

		logs = np.log(term)

	
		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			b_new = sqrt(bb**2 + self.bmin**2)
		
			alpha = self.pdf.alphasQ(self.mu_b(bb))
			
			g_0 = 0.55
			
			res = g_0*(1 - np.exp(- self.CF*alpha*b_new**2/np.pi/g_0/self.bmax**2)  )	
			#type(res)
			out = np.append(out,np.exp(-res*logs))
			
		return out	
			
	def g_K_exp_1h_lgm(self,bt,z1):

		ml1= 1.115 # lambda mass
		
	
		eta_p= (1 - 4*ml1**2/z1**2/self.scale**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction

		term = self.scale*zlc1
		term = term/ml1

		logs = np.log(term)

	
		out = array([])	
		if type(bt)==float or type(bt)== np.float64 :bt=array([bt])
		for bb in bt: 

			b_new = sqrt(bb**2 + self.bmin**2)
		
			alpha = self.pdf.alphasQ(self.mu_b(bb))
			
			g_0 = .7
			
			res = g_0*(1 - np.exp(- self.CF*alpha*b_new**2/np.pi/g_0/self.bmax**2)  )	
			
			out = np.append(out,np.exp(-res*logs))
			
		return out



		
		
		
		
		
		
		
		
		
