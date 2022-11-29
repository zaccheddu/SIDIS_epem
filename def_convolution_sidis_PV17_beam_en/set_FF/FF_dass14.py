# -*- coding: utf-8 -*-
#!/home/zack/Documenti/AKK_f2py/__main__.py
import sys
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
import os
percorso = os.getcwd()
path = os.path.join(percorso, "set_FF")
#print(percorso)
sys.path.insert(1, path)

import dss_14_k
import dss_14_pi

class dass14_ff():


	def D14(self,hadron,parton,x,q):


		io=1
		if hadron=='k+':
			self.ddh=np.array(dss_14_k.fdss(2,2,1,io,x,q**2))
			if parton == "gl": return self.ddh[8]
			elif parton == "u": return self.ddh[0]
			elif parton == "ub": return self.ddh[1]
			elif parton == "d": return self.ddh[2]
			elif parton == "db": return self.ddh[3]
			elif parton == "s": return self.ddh[4]
			elif parton == "sb": return self.ddh[5]
			elif parton == "c": return self.ddh[6]
			elif parton == "b": return self.ddh[7]
		elif hadron=='k-':
			self.ddh=np.array(dss_14_k.fdss(2,2,-1,io,x,q**2))
			if parton == "gl": return self.ddh[8]
			elif parton == "u": return self.ddh[0]
			elif parton == "ub": return self.ddh[1]
			elif parton == "d": return self.ddh[2]
			elif parton == "db": return self.ddh[3]
			elif parton == "s": return self.ddh[4]
			elif parton == "sb": return self.ddh[5]
			elif parton == "c": return self.ddh[6]
			elif parton == "b": return self.ddh[7]
	

		elif hadron=='pi+':
			self.ddh=np.array(dss_14_pi.fdss(2,1,1,io,x,q**2))
			if parton == "gl": return self.ddh[8]
			elif parton == "u": return self.ddh[0]
			elif parton == "ub": return self.ddh[1]
			elif parton == "d": return self.ddh[2]
			elif parton == "db": return self.ddh[3]
			elif parton == "s": return self.ddh[4]
			elif parton == "sb": return self.ddh[5]
			elif parton == "c": return self.ddh[6]
			elif parton == "b": return self.ddh[7]
###______________________________________________________________________
		elif hadron=='pi-':
			self.ddh=np.array(dss_14_pi.fdss(2,1,-1,io,x,q**2))
			if parton == "gl": return self.ddh[8]
			elif parton == "u": return self.ddh[0]
			elif parton == "ub": return self.ddh[1]
			elif parton == "d": return self.ddh[2]
			elif parton == "db": return self.ddh[3]
			elif parton == "s": return self.ddh[4]
			elif parton == "sb": return self.ddh[5]
			elif parton == "c": return self.ddh[6]
			elif parton == "b": return self.ddh[7]



#----- esempio
#	dss=dass14_ff()
#	dup1=dss.D14('pi+',"u",0.5,91.2)



