# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
from set_FF.FF_akk import akk_ff
from set_FF.FF_dass import dass_ff
from set_FF.FF_dass07 import dass07_ff
from set_FF.FF_dass14 import dass14_ff
import lhapdf
######¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
#	definizioni della sezione d'urto non polarizzata e polarizzatta
#	per la produzione associata lambda adrone e lambda-jet
#____________________________________________________________________
pip_ff = lhapdf.mkPDF("NNFF10_PIp_lo")
pim_ff = lhapdf.mkPDF("NNFF10_PIm_lo")

kap_ff = lhapdf.mkPDF("NNFF10_KAp_lo")
kam_ff = lhapdf.mkPDF("NNFF10_KAm_lo")



class cr_sec:


	def __init__(self):

		self.mass = 0.

		self.frag2= 'dss'	
	
		self.qq = 10.58

		self.charm = 'yes'

	def cross_sec2(self,had1,had2,z1,z2,q):	#lambda-had non pol

		ml= self.mass
		
		eta_p= (1 - 4*ml**2/z1**2/self.qq**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction


		if had1=='lbd': 

			fact=1/(2-zp1) ##1/(1+(1-z)**2)	
			fact_bar=(1-zp1)/(2-zp1)  ## ((1-z)**2)/(1+(1-z)**2)	

			#### high suppression
			#fact=1/(1+(1-zp1)**2)	
			#fact_bar=((1-zp1)**2)/(1+(1-zp1)**2)	


			hadron1='lbd'
		#
			dss1=akk_ff()
			dup1=dss1.D1(hadron1,"u",zp1,q)
			dup1=dup1*fact	
		#	
			dupb1=dss1.D1(hadron1,"ub",zp1,q)
			dupb1=dupb1*fact_bar		
		#
			ddo1=dss1.D1(hadron1,"d",zp1,q)
			ddo1=ddo1*fact		

			ddob1=dss1.D1(hadron1,"db",zp1,q)
			ddob1=ddob1*fact_bar
		#
			dst1=dss1.D1(hadron1,"s",zp1,q)
			dst1=dst1*fact
		#		
			dstb1=dss1.D1(hadron1,"sb",zp1,q)	
			dstb1=dstb1*fact_bar

			dch1=dss1.D1(hadron1,"c",zp1,q)
			dch1=dch1/2#*fact

			dchb1=dss1.D1(hadron1,"cb",zp1,q)
			dchb1=dchb1/2#*fact_bar
			
			
		
		elif had1=='lbd_b':

			fact_bar=1/(2-zp1)
			fact=(1-zp1)/(2-zp1)

			#### high suppression
			#fact_bar=1/(1+(1-zp1)**2)	
			#fact=((1-zp1)**2)/(1+(1-zp1)**2)	


			hadron1='lbd'
		#
			dss1=akk_ff()
			dup1=dss1.D1(hadron1,"u",zp1,q)
			dup1=dup1*fact	
		#	
			dupb1=dss1.D1(hadron1,"ub",zp1,q)
			dupb1=dupb1*fact_bar		
		#
			ddo1=dss1.D1(hadron1,"d",zp1,q)
			ddo1=ddo1*fact		

			ddob1=dss1.D1(hadron1,"db",zp1,q)
			ddob1=ddob1*fact_bar
		#
			dst1=dss1.D1(hadron1,"s",zp1,q)
			dst1=dst1*fact
		#		
			dstb1=dss1.D1(hadron1,"sb",zp1,q)	
			dstb1=dstb1*fact_bar

			dch1=dss1.D1(hadron1,"c",zp1,q)
			dch1=dch1/2#*fact
		#
			dchb1=dss1.D1(hadron1,"cb",zp1,q)
			dchb1=dchb1/2#*fact_bar


	
		fr=self.frag2

		if fr=='dss':
			dss2=dass_ff()

			dup2=dss2.D3(had2,"u",z2,q,0)/z2
			dupb2=dss2.D3(had2,"ub",z2,q,0)/z2
			ddo2=dss2.D3(had2,"d",z2,q,0)/z2
			ddob2=dss2.D3(had2,"db",z2,q,0)/z2
			dst2=dss2.D3(had2,"s",z2,q,0)/z2
			dstb2=dss2.D3(had2,"sb",z2,q,0)/z2	

			dch2=dss2.D3(had2,"c",z2,q,0)/z2
			dchb2=dss2.D3(had2,"c",z2,q,0)/z2	
			
			
		elif fr=='dss07':
			dss2=dass07_ff()

			dup2=dss2.D07(had2,"u",z2,q)/z2
			dupb2=dss2.D07(had2,"ub",z2,q)/z2
			ddo2=dss2.D07(had2,"d",z2,q)/z2
			ddob2=dss2.D07(had2,"db",z2,q)/z2
			dst2=dss2.D07(had2,"s",z2,q)/z2
			dstb2=dss2.D07(had2,"sb",z2,q)/z2	

		elif fr=='dss14':
			dss2=dass14_ff()

			dup2=dss2.D14(had2,"u",z2,q)/z2
			dupb2=dss2.D14(had2,"ub",z2,q)/z2
			ddo2=dss2.D14(had2,"d",z2,q)/z2
			ddob2=dss2.D14(had2,"db",z2,q)/z2
			dst2=dss2.D14(had2,"s",z2,q)/z2
			dstb2=dss2.D14(had2,"sb",z2,q)/z2	


			dch2=dss2.D14(had2,"c",z2,q)/z2
			dchb2=dss2.D14(had2,"c",z2,q)/z2


		elif fr =='nnpdf':
		
			if had2=='pi+':
			
				dup2=self.pip_ff.xfxQ(2,z2,q)/z2
				dupb2=self.pip_ff.xfxQ(-2,z2,q)/z2
				
				ddo2=self.pip_ff.xfxQ(1,z2,q)/z2
				ddob2=self.pip_ff.xfxQ(-1,z2,q)/z2
				
				dst2=self.pip_ff.xfxQ(3,z2,q)/z2
				dstb2=self.pip_ff.xfxQ(-3,z2,q)/z2	
				#dgl = self.pip_ff.xfxQ(21,z2,q)/z2

				dch2=self.pip_ff.xfxQ(4,z2,q)/z2
				dchb2=self.pip_ff.xfxQ(-4,z2,q)/z2	
			
			
			elif had2=='pi-':
				dup2=self.pim_ff.xfxQ(2,z2,q)/z2
				dupb2=self.pim_ff.xfxQ(-2,z2,q)/z2
				
				ddo2=self.pim_ff.xfxQ(1,z2,q)/z2
				ddob2=self.pim_ff.xfxQ(-1,z2,q)/z2
				
				dst2=self.pim_ff.xfxQ(3,z2,q)/z2
				dstb2=self.pim_ff.xfxQ(-3,z2,q)/z2	
				#dgl = self.pip_ff.xfxQ(21,z2,q)/z2

				dch2=self.pim_ff.xfxQ(4,z2,q)/z2
				dchb2=self.pim_ff.xfxQ(-4,z2,q)/z2	
			

			elif had2=='k+':
				dup2=self.kap_ff.xfxQ(2,z2,q)/z2
				dupb2=self.kap_ff.xfxQ(-2,z2,q)/z2
				
				ddo2=self.kap_ff.xfxQ(1,z2,q)/z2
				ddob2=self.kap_ff.xfxQ(-1,z2,q)/z2
				
				dst2=self.kap_ff.xfxQ(3,z2,q)/z2
				dstb2=self.kap_ff.xfxQ(-3,z2,q)/z2	
				#dgl = self.pip_ff.xfxQ(21,z2,q)/z2

				dch2=self.kap_ff.xfxQ(4,z2,q)/z2
				dchb2=self.kap_ff.xfxQ(-4,z2,q)/z2	


			elif had2=='k-':
				dup2=self.kam_ff.xfxQ(2,z2,q)/z2
				dupb2=self.kam_ff.xfxQ(-2,z2,q)/z2
				
				ddo2=self.kam_ff.xfxQ(1,z2,q)/z2
				ddob2=self.kam_ff.xfxQ(-1,z2,q)/z2
				
				dst2=self.kam_ff.xfxQ(3,z2,q)/z2
				dstb2=self.kam_ff.xfxQ(-3,z2,q)/z2	
				#dgl = self.pip_ff.xfxQ(21,z2,q)/z2

				dch2=self.kam_ff.xfxQ(4,z2,q)/z2
				dchb2=self.kam_ff.xfxQ(-4,z2,q)/z2	



		cs2= 4/9*(dup1*dupb2 + dupb1*dup2) + 1/9*(ddo1*ddob2 + ddob1*ddo2) + 1/9*(dst1*dstb2 + dstb1*dst2) 
		if self.charm == 'yes': 
			cs2= cs2 + 4/9*(dch1*dchb2 + dchb1*dch2)
			#print('yes ss')
	#
		return cs2


	def cross_sec2_polda(self,had1,had2,z1,z2,q,param): # lambda-had polarizzata   ,z1 z2 energy fractions

		ml= self.mass

		eta_p= (1 - 4*ml**2/((z1**2)*self.qq**2))
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 =  (z1 + zp1)/2	#light-cone momentum fraction



		nup=param[0]
		ndo=param[1]
		nst=param[2]
		nupb=param[3]
		ndob=param[4]
		nstb=param[5]
	#¯¯¯¯¯¯
		aup=param[6]
		ado=param[7]
		ast=param[8]
		aupb=param[9]
		adob=param[10]
		astb=param[11]
	#¯¯¯¯¯¯¯
		bup=param[12]
		bdo=param[13]
		bst=param[14]
		bupb=param[15]
		bdob=param[16]
		bstb=param[17]#### metto i parametri



		fact=1/(2-zp1)
		fact_bar=(1-zp1)/(2-zp1)

		#### high suppression
		#fact=1/(1+(1-zp1)**2)	
		#fact_bar=((1-zp1)**2)/(1+(1-zp1)**2)	


		hadron1='lbd'
	#
		dss1=akk_ff()
		dup1=dss1.D1(hadron1,"u",zp1,q)
		dup1=dup1*fact	
	#	
		dupb1=dss1.D1(hadron1,"ub",zp1,q)
		dupb1=dupb1*fact_bar		
	#
		ddo1=dss1.D1(hadron1,"d",zp1,q)
		ddo1=ddo1*fact		

		ddob1=dss1.D1(hadron1,"db",zp1,q)
		ddob1=ddob1*fact_bar
	#
		dst1=dss1.D1(hadron1,"s",zp1,q)
		dst1=dst1*fact
	#		
		dstb1=dss1.D1(hadron1,"sb",zp1,q)	
		dstb1=dstb1*fact_bar

		dup1=dup1*nup*zlc1**aup*(1.-zlc1)**bup/ ( aup**aup * bup**bup / (aup+bup)**(aup+bup) )
		ddo1=ddo1*ndo*zlc1**ado*(1.-zlc1)**bdo/ ( ado**ado * bdo**bdo / (ado+bdo)**(ado+bdo) )
		dst1=dst1*nst*zlc1**ast*(1.-zlc1)**bst/ ( ast**ast * bst**bst / (ast+bst)**(ast+bst) )

		dupb1=dupb1*nupb*zlc1**aupb*(1.-zlc1)**bupb/ ( aupb**aupb * bupb**bupb / (aupb+bupb)**(aupb+bupb) )
		ddob1=ddob1*ndob*zlc1**adob*(1.-zlc1)**bdob/ ( adob**adob * bdob**bdob / (adob+bdob)**(adob+bdob) )
		dstb1=dstb1*nstb*zlc1**astb*(1.-zlc1)**bstb/ ( astb**astb * bstb**bstb / (astb+bstb)**(astb+bstb) )


		
		if had1=='lbd_b':
			tmp   = dupb1
			dupb1 = dup1
			dup1  = tmp
			tmp   = ddob1
			ddob1 = ddo1
			ddo1  = tmp
			tmp   = dstb1
			dstb1 = dst1
			dst1  = tmp

		fr=self.frag2
		if fr=='dss':
			dss2=dass_ff()

			dup2=dss2.D3(had2,"u",z2,q,0)/z2
			dupb2=dss2.D3(had2,"ub",z2,q,0)/z2
			ddo2=dss2.D3(had2,"d",z2,q,0)/z2
			ddob2=dss2.D3(had2,"db",z2,q,0)/z2
			dst2=dss2.D3(had2,"s",z2,q,0)/z2
			dstb2=dss2.D3(had2,"sb",z2,q,0)/z2	
		elif fr=='dss07':
			dss2=dass07_ff()

			dup2=dss2.D07(had2,"u",z2,q)/z2
			dupb2=dss2.D07(had2,"ub",z2,q)/z2
			ddo2=dss2.D07(had2,"d",z2,q)/z2
			ddob2=dss2.D07(had2,"db",z2,q)/z2
			dst2=dss2.D07(had2,"s",z2,q)/z2
			dstb2=dss2.D07(had2,"sb",z2,q)/z2	

		elif fr=='dss14':
			dss2=dass14_ff()

			dup2=dss2.D14(had2,"u",z2,q)/z2
			dupb2=dss2.D14(had2,"ub",z2,q)/z2
			ddo2=dss2.D14(had2,"d",z2,q)/z2
			ddob2=dss2.D14(had2,"db",z2,q)/z2
			dst2=dss2.D14(had2,"s",z2,q)/z2
			dstb2=dss2.D14(had2,"sb",z2,q)/z2	

		elif fr =='nnpdf':
		
			if had2=='pi+':
			
				dup2=self.pip_ff.xfxQ(2,z2,q)/z2
				dupb2=self.pip_ff.xfxQ(-2,z2,q)/z2
				
				ddo2=self.pip_ff.xfxQ(1,z2,q)/z2
				ddob2=self.pip_ff.xfxQ(-1,z2,q)/z2
				
				dst2=self.pip_ff.xfxQ(3,z2,q)/z2
				dstb2=self.pip_ff.xfxQ(-3,z2,q)/z2	
				#dgl = self.pip_ff.xfxQ(21,z2,q)/z2

				dch2=self.pip_ff.xfxQ(4,z2,q)/z2
				dchb2=self.pip_ff.xfxQ(-4,z2,q)/z2	
			
			
			elif had2=='pi-':
				dup2=self.pim_ff.xfxQ(2,z2,q)/z2
				dupb2=self.pim_ff.xfxQ(-2,z2,q)/z2
				
				ddo2=self.pim_ff.xfxQ(1,z2,q)/z2
				ddob2=self.pim_ff.xfxQ(-1,z2,q)/z2
				
				dst2=self.pim_ff.xfxQ(3,z2,q)/z2
				dstb2=self.pim_ff.xfxQ(-3,z2,q)/z2	
				#dgl = self.pip_ff.xfxQ(21,z2,q)/z2

				dch2=self.pim_ff.xfxQ(4,z2,q)/z2
				dchb2=self.pim_ff.xfxQ(-4,z2,q)/z2	
			

			elif had2=='k+':
				dup2=self.kap_ff.xfxQ(2,z2,q)/z2
				dupb2=self.kap_ff.xfxQ(-2,z2,q)/z2
				
				ddo2=self.kap_ff.xfxQ(1,z2,q)/z2
				ddob2=self.kap_ff.xfxQ(-1,z2,q)/z2
				
				dst2=self.kap_ff.xfxQ(3,z2,q)/z2
				dstb2=self.kap_ff.xfxQ(-3,z2,q)/z2	
				#dgl = self.pip_ff.xfxQ(21,z2,q)/z2

				dch2=self.kap_ff.xfxQ(4,z2,q)/z2
				dchb2=self.kap_ff.xfxQ(-4,z2,q)/z2	


			elif had2=='k-':
				dup2=self.kam_ff.xfxQ(2,z2,q)/z2
				dupb2=self.kam_ff.xfxQ(-2,z2,q)/z2
				
				ddo2=self.kam_ff.xfxQ(1,z2,q)/z2
				ddob2=self.kam_ff.xfxQ(-1,z2,q)/z2
				
				dst2=self.kam_ff.xfxQ(3,z2,q)/z2
				dstb2=self.kam_ff.xfxQ(-3,z2,q)/z2	
				#dgl = self.pip_ff.xfxQ(21,z2,q)/z2

				dch2=self.kam_ff.xfxQ(4,z2,q)/z2
				dchb2=self.kam_ff.xfxQ(-4,z2,q)/z2	

		

		cs2= 4/9*(dup1*dupb2 + dupb1*dup2) + 1/9*(ddo1*ddob2 + ddob1*ddo2) + 1/9*(dst1*dstb2 + dstb1*dst2)
		
		
		UP = 4/9*(dup1*dupb2)
		DO = 1/9*(ddo1*ddob2)
		ST = 1/9*(dst1*dstb2)
		
		UPb = 4/9*(dupb1*dup2)
		DOb = 1/9*(ddob1*ddo2)
		STb = 1/9*(dstb1*dst2)
	#
		return cs2, UP, DO, ST, UPb, DOb, STb






	def ratio_pol(self,had1,had2,z1,z2,q,param):  #rapporto fra numeratore e denominatore

	
		num, UP, DO, ST, UPb, DOb, STb=self.cross_sec2_polda(had1,had2,z1,z2,q,param)
		den=self.cross_sec2(had1,had2,z1,z2,q)	

		ratio=num/den
		return ratio

	def ratio_pol_partial(self,had1,had2,z1,z2,q,param):  #rapporto fra numeratore e denominatore e rapporti parziali

	
		num, UP, DO, ST, UPb, DOb, STb=self.cross_sec2_polda(had1,had2,z1,z2,q,param)
		den=self.cross_sec2(had1,had2,z1,z2,q)	

		ratio=num/den
		return ratio, UP/den, DO/den, ST/den, UPb/den, DOb/den, STb/den


	def fact_pol(self,z1,z2,pt_pl1,pt_pl2,pt_pp,q):  #fattore per le rotazioni

		ml= self.mass

		#ml= 1.115

		eta_p= (1 - 4*ml**2/((z1**2)*q**2))
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction


##### per fittare la delta D, con il bound integrato 
#		print(type(pt_pl1)

		mp=sqrt( (pt_pl1*pt_pp)/(pt_pl1-pt_pp))
		pre_fact =sqrt(e*pi/2)/mp
		pre_fact =pre_fact*pt_pp**2/pt_pl1
	#
		bt2=(z2**2)/((zp1**2)*pt_pl2 + (z2**2)*pt_pp )
		bt=sqrt(bt2)

		factt=pre_fact*bt


##### per fittare la delta D 


		#mp=sqrt( (pt_pl1*pt_pp)/(pt_pl1-pt_pp))
#		pre_fact=sqrt(pi)/2
#		pre_fact=pre_fact*pt_pp
	#
#		bt2=(z2**2)/((zp1**2)*pt_pl2 + (z2**2)*pt_pp )
#		bt=sqrt(bt2)

#		factt=pre_fact*bt

#### per fittare il primo momento uso come prefattore invece

#		pre_fact=sqrt(pi)/2
#		pre_fact=pre_fact*2*zlc1*ml

#		bt2=(z2**2)/((zp1**2)*pt_pl2 + (z2**2)*pt_pp )
#		bt=sqrt(bt2)

#		factt=pre_fact*bt

		return factt

	def fact_fst_mom(self,z1,pt_pl1,pt_pp):

		q = 10.58

		mp=sqrt( (pt_pl1*pt_pp)/(pt_pl1-pt_pp))

		ml = self.mass

		eta_p= (1 - 4*ml**2/((z1**2)*q**2))
		zp1 = z1*sqrt(eta_p)		# momentum fraction
		
		factt = sqrt(2*e)/2/mp
		factt = factt/zp1/ml
		factt = factt*pt_pp**2/pt_pl1

		return factt

	def polarisation(self,had1,had2,z1,z2,q,pt_pl1,pt_pl2,pt_pp,param):  # polarizzazione

		ratio=self.ratio_pol(had1,had2,z1,z2,q,param)
		faxx=self.fact_pol(z1,z2,pt_pl1,pt_pl2,pt_pp,q)
		pols=ratio*faxx

		return pols
		

########################################################################



	def cross_sec1(self,had1,z1,q): # lambda-jet non polarizzata

		ml= self.mass

		qq = 10.58

		eta_p= (1 - 4*ml**2/((z1**2)*qq**2))
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = 1/2*(z1 + zp1)	#light-cone momentum fraction


		if had1=='lbd': 

			fact=1/(2-zp1)
			fact_bar=(1-zp1)/(2-zp1)

			hadron1='lbd'
		#
			dss1=akk_ff()
			dup1=dss1.D1(hadron1,"u",zp1,q)
			dup1=dup1*fact	
		#	
			dupb1=dss1.D1(hadron1,"ub",zp1,q)
			dupb1=dupb1*fact_bar		
		#
			ddo1=dss1.D1(hadron1,"d",zp1,q)
			ddo1=ddo1*fact		

			ddob1=dss1.D1(hadron1,"db",zp1,q)
			ddob1=ddob1*fact_bar
		#
			dst1=dss1.D1(hadron1,"s",zp1,q)
			dst1=dst1*fact
		#		
			dstb1=dss1.D1(hadron1,"sb",zp1,q)	
			dstb1=dstb1*fact_bar
		
		elif had1=='lbd_b':

			fact_bar=1/(2-zp1)
			fact=(1-zp1)/(2-zp1)

			hadron1='lbd'
		#
			dss1=akk_ff()
			dup1=dss1.D1(hadron1,"u",zp1,q)
			dup1=dup1*fact	
		#	
			dupb1=dss1.D1(hadron1,"ub",zp1,q)
			dupb1=dupb1*fact_bar		
		#
			ddo1=dss1.D1(hadron1,"d",zp1,q)
			ddo1=ddo1*fact		

			ddob1=dss1.D1(hadron1,"db",zp1,q)
			ddob1=ddob1*fact_bar
		#
			dst1=dss1.D1(hadron1,"s",zp1,q)
			dst1=dst1*fact
		#		
			dstb1=dss1.D1(hadron1,"sb",zp1,q)	
			dstb1=dstb1*fact_bar

		cs2= 4/9*(dup1 + dupb1) + 1/9*(ddo1 + ddob1) + 1/9*(dst1 + dstb1)
	#
		return cs2

	def cross_sec1_polda(self,had1,z1,q,param):  # lambda-jet polarizzata

		ml= self.mass
		qq=10.58

		eta_p= (1 - 4*ml**2/z1**2/qq**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = 1/2*(z1 + zp1)	#light-cone momentum fraction


		nup=param[0]
		ndo=param[1]
		nst=param[2]
		nupb=param[3]
		ndob=param[4]
		nstb=param[5]
	#¯¯¯¯¯¯
		aup=param[6]
		ado=param[7]
		ast=param[8]
		aupb=param[9]
		adob=param[10]
		astb=param[11]
	#¯¯¯¯¯¯¯
		bup=param[12]
		bdo=param[13]
		bst=param[14]
		bupb=param[15]
		bdob=param[16]
		bstb=param[17]#### metto i parametri



		fact=1/(2-zp1)
		fact_bar=(1-zp1)/(2-zp1)

		hadron1='lbd'
	#
		dss1=akk_ff()
		dup1=dss1.D1(hadron1,"u",zp1,q)
		dup1=dup1*fact	
	#	
		dupb1=dss1.D1(hadron1,"ub",zp1,q)
		dupb1=dupb1*fact_bar		
	#
		ddo1=dss1.D1(hadron1,"d",zp1,q)
		ddo1=ddo1*fact		

		ddob1=dss1.D1(hadron1,"db",zp1,q)
		ddob1=ddob1*fact_bar
	#
		dst1=dss1.D1(hadron1,"s",zp1,q)
		dst1=dst1*fact
	#		
		dstb1=dss1.D1(hadron1,"sb",zp1,q)	
		dstb1=dstb1*fact_bar

		dup1=dup1*nup*zlc1**aup*(1.-zlc1)**bup*(aup+bup)**(aup+bup)/aup**aup/bup**bup 
		ddo1=ddo1*ndo*zlc1**ado*(1.-zlc1)**bdo*(ado+bdo)**(ado+bdo)/ado**ado/bdo**bdo 
		dst1=dst1*nst*zlc1**ast*(1.-zlc1)**bst*(ast+bst)**(ast+bst)/ast**ast/bst**bst 

		dupb1=dupb1*nupb*zlc1**aupb*(1.-zlc1)**bupb*(aupb+bupb)**(aupb+bupb)/aupb**aupb/bupb**bupb 
		ddob1=ddob1*ndob*zlc1**adob*(1.-zlc1)**bdob*(adob+bdob)**(adob+bdob)/adob**adob/bdob**bdob 
		dstb1=dstb1*nstb*zlc1**astb*(1.-zlc1)**bstb*(astb+bstb)**(astb+bstb)/astb**astb/bstb**bstb 


		
		if had1=='lbd_b':
			tmp   = dupb1
			dupb1 = dup1
			dup1  = tmp
			tmp   = ddob1
			ddob1 = ddo1
			ddo1  = tmp
			tmp   = dstb1
			dstb1 = dst1
			dst1  = tmp


		cs2= 4/9*(dup1 + dupb1) + 1/9*(ddo1 + ddob1) + 1/9*(dst1 + dstb1)
		
		
	#
		return cs2


	def cross_sec1_polda_fm(self,had1,z1,q,param):  # lambda-jet polarizzata

		ml= self.mass
		qq=10.58

		eta_p= (1 - 4*ml**2/z1**2/qq**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = 1/2*(z1 + zp1)	#light-cone momentum fraction


		nup=param[0]
		ndo=param[1]
		nst=param[2]
		nupb=param[3]
		ndob=param[4]
		nstb=param[5]
	#¯¯¯¯¯¯
		aup=param[6]
		ado=param[7]
		ast=param[8]
		aupb=param[9]
		adob=param[10]
		astb=param[11]
	#¯¯¯¯¯¯¯
		bup=param[12]
		bdo=param[13]
		bst=param[14]
		bupb=param[15]
		bdob=param[16]
		bstb=param[17]#### metto i parametri



		fact=1/(2-zp1)
		fact_bar=(1-zp1)/(2-zp1)

		hadron1='lbd'
	#
		dss1=akk_ff()
		dup1=dss1.D1(hadron1,"u",zp1,q)
		dup1=dup1*fact	
	#	
		dupb1=dss1.D1(hadron1,"ub",zp1,q)
		dupb1=dupb1*fact_bar		
	#
		ddo1=dss1.D1(hadron1,"d",zp1,q)
		ddo1=ddo1*fact		

		ddob1=dss1.D1(hadron1,"db",zp1,q)
		ddob1=ddob1*fact_bar
	#
		dst1=dss1.D1(hadron1,"s",zp1,q)
		dst1=dst1*fact
	#		
		dstb1=dss1.D1(hadron1,"sb",zp1,q)	
		dstb1=dstb1*fact_bar

		dup1=dup1*nup*zlc1**aup*(1.-zlc1)**bup#*(aup+bup)**(aup+bup)/aup**aup/bup**bup 
		ddo1=ddo1*ndo*zlc1**ado*(1.-zlc1)**bdo#*(ado+bdo)**(ado+bdo)/ado**ado/bdo**bdo 
		dst1=dst1*nst*zlc1**ast*(1.-zlc1)**bst#*(ast+bst)**(ast+bst)/ast**ast/bst**bst 

		dupb1=dupb1*nupb*zlc1**aupb*(1.-zlc1)**bupb#*(aupb+bupb)**(aupb+bupb)/aupb**aupb/bupb**bupb 
		ddob1=ddob1*ndob*zlc1**adob*(1.-zlc1)**bdob#*(adob+bdob)**(adob+bdob)/adob**adob/bdob**bdob 
		dstb1=dstb1*nstb*zlc1**astb*(1.-zlc1)**bstb#*(astb+bstb)**(astb+bstb)/astb**astb/bstb**bstb 


		
		if had1=='lbd_b':
			tmp   = dupb1
			dupb1 = dup1
			dup1  = tmp
			tmp   = ddob1
			ddob1 = ddo1
			ddo1  = tmp
			tmp   = dstb1
			dstb1 = dst1
			dst1  = tmp


		cs2= 4/9*(dup1 + dupb1) + 1/9*(ddo1 + ddob1) + 1/9*(dst1 + dstb1)
	#
		return cs2


	def ratio_pol2(self,had1,z1,q,param,pt,pt_pl1,pt_pp):  #rapporto numeratore denominatore per l-jet con   

								# con fattori già inclusi
#		mp=sqrt( (pt_pl1*pt_pp)/(pt_pl1-pt_pp))		
#restituisce il valore della polarizzazione
		ml= self.mass

		eta_p= (1 - 4*ml**2/z1**2/q**2)
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = 1/2*(z1 + zp1)	#light-cone momentum fraction

		

		fact=(e**(-pt**2/pt_pl1))/(pi*pt_pl1)

		cost=1/pt_pp**2/pi
		cost = cost*2*zlc1*ml
		fact_delt=cost*pt*e**(-pt**2/pt_pp)

		den=self.cross_sec1(had1,z1,q)*fact
		num=self.cross_sec1_polda(had1,z1,q,param)*fact_delt

		ratio=num/den
		return ratio


	def denominator(self,had1,z1,q,pt,pt_pl1):  #rapporto numeratore denominatore per l-jet con   

								# con fattori già inclusi
#		mp=sqrt( (pt_pl1*pt_pp)/(pt_pl1-pt_pp))		
#restituisce il valore della polarizzazione

		

		fact=(e**(-pt**2/pt_pl1))/(pi*pt_pl1)


		den=self.cross_sec1(had1,z1,q)*fact

		return den

	def numerator(self,had1,z1,q,param,pt,pt_pp):  #rapporto numeratore denominatore per l-jet con   

								# con fattori già inclusi
#		mp=sqrt( (pt_pl1*pt_pp)/(pt_pl1-pt_pp))		
#restituisce il valore della polarizzazione




		cost=1/pt_pp/pi
		fact_delt=cost*pt*e**(-pt**2/pt_pp)

#		cost=sqrt(2*e)/mp/pt_pl1/pi
#		fact_delt=cost*pt*e**(-pt**2/pt_pp)


		num=self.cross_sec1_polda(had1,z1,q,param)*fact_delt
#		num=self.cross_sec1_polda(had1,z1,q,param)*pt/pt

		return num



















































