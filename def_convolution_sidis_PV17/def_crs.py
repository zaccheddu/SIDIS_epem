# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
import numpy as np
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


#pdf=lhapdf.mkPDF("cteq6l1")
pdf=lhapdf.mkPDF("CT10nlo")



class cr_sec:



	def __init__(self):

		self.mass = 0.

		#self.frag2= 'dss'	
	
		self.qq = 10.58

		self.charm = 'yes'

	def cross_sec2(self,had1,had2,z1,xb2,q):	#lambda-had non pol

		ml= self.mass
		
		eta_p= (1 - (ml**2/z1**2/self.qq**2)*xb2/(1-xb2))
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		zlc1 = (z1 + zp1)/2	#light-cone momentum fraction


		if had1=='lbd': 

			fact=1/(2-zp1)
			fact_bar=(1-zp1)/(2-zp1)

			hadron1='lbd'
		#
			dss1=akk_ff()
			dup1=dss1.D1(hadron1,"u",zp1,q)
			dup1=dup1*fact	
			#print('UP ' + str(dup1))
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
			dch1=dch1*fact

			dchb1=dss1.D1(hadron1,"cb",zp1,q)
			dchb1=dchb1*fact_bar

		
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

			dch1=dss1.D1(hadron1,"c",zp1,q)
			dch1=dch1*fact

			dchb1=dss1.D1(hadron1,"cb",zp1,q)
			dchb1=dchb1*fact_bar


		#print( dup1)
	
#		dup2=self.pdf.xfxQ(2,xb2,q)/xb2
#		dupb2=self.pdf.xfxQ(-2,xb2,q)/xb2
#		ddo2=self.pdf.xfxQ(1,xb2,q)/xb2
#		ddob2=self.pdf.xfxQ(-1,xb2,q)/xb2
#		dst2=self.pdf.xfxQ(3,xb2,q)/xb2
#		dstb2=self.pdf.xfxQ(-3,xb2,q)/xb2	
	
#		dch2=self.pdf.xfxQ(4,xb2,q)/xb2
#		dchb2=self.pdf.xfxQ(-4,xb2,q)/xb2	
#		print(q)
#		print(type(q))
		#q=np.array([q])
		if type(q)!=ndarray:
			dup2=pdf.xfxQ(2,xb2,q)/xb2
			dupb2=pdf.xfxQ(-2,xb2,q)/xb2
			ddo2=pdf.xfxQ(1,xb2,q)/xb2
			ddob2=pdf.xfxQ(-1,xb2,q)/xb2
			dst2=pdf.xfxQ(3,xb2,q)/xb2
			dstb2=pdf.xfxQ(-3,xb2,q)/xb2	
		
			dch2=pdf.xfxQ(4,xb2,q)/xb2
			dchb2=pdf.xfxQ(-4,xb2,q)/xb2	

		if type(q)==ndarray:
			#print(q)

			dup2= np.array([])
			dupb2= np.array([])
			ddo2= np.array([])
			ddob2= np.array([])
			dst2= np.array([])
			dstb2= np.array([])	
		
			dch2= np.array([])
			dchb2= np.array([])	
			
			for qq in q:#
				dup2= np.append(dup2,pdf.xfxQ(2,xb2,qq)/xb2)
				dupb2=np.append(dupb2,pdf.xfxQ(-2,xb2,qq)/xb2)
				ddo2=np.append(ddo2,pdf.xfxQ(1,xb2,qq)/xb2)
				ddob2=np.append(ddob2,pdf.xfxQ(-1,xb2,qq)/xb2)
				dst2=np.append(dst2,pdf.xfxQ(3,xb2,qq)/xb2)
				dstb2=np.append(dstb2,pdf.xfxQ(-3,xb2,qq)/xb2)	
			
				dch2=np.append(dch2,pdf.xfxQ(4,xb2,qq)/xb2)
				dchb2=np.append(dchb2,pdf.xfxQ(-4,xb2,qq)/xb2)	
				#print(type(dup2))



		cs2= 4/9*(dup1*dup2 + dupb1*dupb2) + 1/9*(ddo1*ddo2 + ddob1*ddob2) + 1/9*(dst1*dst2 + dstb1*dstb2) 
		if self.charm == 'yes': cs2= cs2 + 4/9*(dch1*dch2 + dchb1*dchb2)
	#
		return cs2


	def cross_sec2_polda(self,had1,had2,z1,xb2,q,param): # lambda-had polarizzata   ,z1 z2 energy fractions

		ml= self.mass

		eta_p= (1 - (ml**2/z1**2/self.qq**2)*xb2/(1-xb2))
		zp1 = z1*sqrt(eta_p)		# momentum fraction

		#zlc1 =  (z1 + zp1)/2	#light-cone momentum fraction
		zlc1 = z1


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

		if type(q)!=ndarray:
			dup2=pdf.xfxQ(2,xb2,q)/xb2
			dupb2=pdf.xfxQ(-2,xb2,q)/xb2
			ddo2=pdf.xfxQ(1,xb2,q)/xb2
			ddob2=pdf.xfxQ(-1,xb2,q)/xb2
			dst2=pdf.xfxQ(3,xb2,q)/xb2
			dstb2=pdf.xfxQ(-3,xb2,q)/xb2	
		
			dch2=pdf.xfxQ(4,xb2,q)/xb2
			dchb2=pdf.xfxQ(-4,xb2,q)/xb2	

		if type(q)==ndarray:
			#print(q)

			dup2= np.array([])
			dupb2= np.array([])
			ddo2= np.array([])
			ddob2= np.array([])
			dst2= np.array([])
			dstb2= np.array([])	
		
			dch2= np.array([])
			dchb2= np.array([])	
			
			for qq in q:#
				dup2= np.append(dup2,pdf.xfxQ(2,xb2,qq)/xb2)
				dupb2=np.append(dupb2,pdf.xfxQ(-2,xb2,qq)/xb2)
				ddo2=np.append(ddo2,pdf.xfxQ(1,xb2,qq)/xb2)
				ddob2=np.append(ddob2,pdf.xfxQ(-1,xb2,qq)/xb2)
				dst2=np.append(dst2,pdf.xfxQ(3,xb2,qq)/xb2)
				dstb2=np.append(dstb2,pdf.xfxQ(-3,xb2,qq)/xb2)	
			
				dch2=np.append(dch2,pdf.xfxQ(4,xb2,qq)/xb2)
				dchb2=np.append(dchb2,pdf.xfxQ(-4,xb2,qq)/xb2)	
				#print(type(dup2))
	

		cs2= 4/9*(dup1*dup2 + dupb1*dupb2) + 1/9*(ddo1*ddo2 + ddob1*ddob2) + 1/9*(dst1*dst2 + dstb1*dstb2) 
	#
		return cs2


















































