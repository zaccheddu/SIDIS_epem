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

import dss_07_k
import dss_07_pi


class dass07_ff():


	def D07(self,hadron,parton,x,q):


		io=0
		if hadron=='k+':
			self.ddh=np.array(dss_07_k.fdss(1,2,1,io,x,q**2))
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
			self.ddh=np.array(dss_07_k.fdss(1,2,-1,io,x,q**2))
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
			self.ddh=np.array(dss_07_pi.fdss(1,1,1,io,x,q**2))
			if parton == "gl": return self.ddh[8]
			elif parton == "u": return self.ddh[0]
			elif parton == "ub": return self.ddh[1]
			elif parton == "d": return self.ddh[2]
			elif parton == "db": return self.ddh[3]
			elif parton == "s": return self.ddh[4]
			elif parton == "sb": return self.ddh[5]
			elif parton == "c": return self.ddh[6]
			elif parton == "b": return self.ddh[7]
####______________________________________________________________________
		elif hadron=='pi-':
			self.ddh=np.array(dss_07_pi.fdss(1,1,-1,io,x,q**2))
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
#	dss1=dass07_ff()
#	dup=dss1.D07('pi+',"u",0.5,91.2)




