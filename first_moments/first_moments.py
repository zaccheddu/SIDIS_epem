# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
from set_FF.FF_akk import akk_ff
from set_FF.FF_dass import dass_ff
from set_FF.FF_dass07 import dass07_ff
from set_FF.FF_dass14 import dass14_ff

######¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯
#	definizioni della sezione d'urto non polarizzata e polarizzatta
#	per la produzione associata lambda adrone e lambda-jet
#____________________________________________________________________



class cr_sec:


	def __init__(self):

		self.mass = 1.115

		self.frag2= 'dss'	
	
		self.qq = 10.58


	def jack(self,z1):

		ml=self.mass
		q=self.qq

		eta_p= (1 - 4*ml**2/z1**2/q**2)
		zp1 = z1*sqrt(eta_p)	# momentum fraction

		zlc1 = 1/2*(z1 + zp1)	#light-cone momentum fraction

		jack= 2/(1 + sqrt(eta_p))


		return jack, zlc1



	def cross_sec2_polda(self,had1,z1,q,param): # lambda-had polarizzata   ,z1 z2 energy fractions

		ml= self.mass

		eta_p= (1 - 4*ml**2/((z1**2)*self.qq**2))
		#eta_p = 1 - ml**2/z1**2/self.qq**2
		zp1 = z1*eta_p		# momentum fraction

		#zlc1 = z1 # (z1 + zp1)/2	#light-cone momentum fraction
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


	#
		return dup1, ddo1, dst1, dupb1, ddob1, dstb1


	def cross_sec2_polda_fst(self,had1,z1,q,param): # lambda-had polarizzata   ,z1 z2 energy fractions

		ml= self.mass

		#eta_p= (1 - 4*ml**2/((z1**2)*self.qq**2))
		eta_p = 1 - ml**2/2/z1**2/self.qq**2
		zp1 = z1*eta_p		# momentum fraction

		zlc1 = z1 # (z1 + zp1)/2	#light-cone momentum fraction
		#zlc1 =  (z1 + zp1)/2	#light-cone momentum fraction



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


	#
		return dup1, ddo1, dst1, dupb1, ddob1, dstb1



	def cross_sec2_polda_jk(self,had1,z1,q,param): # lambda-had polarizzata   ,z1 z2 energy fractions

		ml= self.mass

		eta_p= (1 - 4*ml**2/((z1**2)*self.qq**2))
		#eta_p = 1 - ml**2/z1**2/self.qq**2
		zp1 = z1*eta_p		# momentum fraction

		#zlc1 = z1 # (z1 + zp1)/2	#light-cone momentum fraction
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

		jack, zlc=self.jack(z1)
	#
		return dup1*jack*zlc1, ddo1*jack*zlc1, dst1*jack*zlc1, dupb1*jack*zlc1, ddob1*jack*zlc1, dstb1*jack*zlc1




















