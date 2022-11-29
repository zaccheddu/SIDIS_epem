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

import dss_pi
import dss_k




class dass_ff():

	def D3(self,hadron,parton,x,q,io):

		if hadron=='pi+':
			self.ddh=np.array(dss_pi.fdssda(1,1,io,x,q**2))
			if parton == "gl": return self.ddh[0]
			elif parton == "u": return self.ddh[1]
			elif parton == "ub": return self.ddh[2]
			elif parton == "d": return self.ddh[3]
			elif parton == "db": return self.ddh[4]
			elif parton == "s": return self.ddh[5]
			elif parton == "sb": return self.ddh[6]
			elif parton == "c": return self.ddh[7]
			elif parton == "cb": return self.ddh[8]
			elif parton == "b": return self.ddh[9]
			elif parton == "bb": return self.ddh[10]
####______________________________________________________________________
		elif hadron=='pi-':
			self.ddh=np.array(dss_pi.fdssda(1,-1,io,x,q**2))
			if parton == "gl": return self.ddh[0]
			elif parton == "u": return self.ddh[1]
			elif parton == "ub": return self.ddh[2]
			elif parton == "d": return self.ddh[3]
			elif parton == "db": return self.ddh[4]
			elif parton == "s": return self.ddh[5]
			elif parton == "sb": return self.ddh[6]
			elif parton == "c": return self.ddh[7]
			elif parton == "cb": return self.ddh[8]
			elif parton == "b": return self.ddh[9]
			elif parton == "bb": return self.ddh[10]
####______________________________________________________________________
		elif hadron=='pi':
			self.ddh=np.array(dss_pi.fdssda(1,0,io,x,q**2))
			if parton == "gl": return self.ddh[0]
			elif parton == "u": return self.ddh[1]
			elif parton == "ub": return self.ddh[2]
			elif parton == "d": return self.ddh[3]
			elif parton == "db": return self.ddh[4]
			elif parton == "s": return self.ddh[5]
			elif parton == "sb": return self.ddh[6]
			elif parton == "c": return self.ddh[7]
			elif parton == "cb": return self.ddh[8]
			elif parton == "b": return self.ddh[9]
			elif parton == "bb": return self.ddh[10]
###___________________________________________________________________
		elif hadron=='k+':
			self.ddh=np.array(dss_k.fdssda(2,1,io,x,q**2))
			if parton == "gl": return self.ddh[0]
			elif parton == "u": return self.ddh[1]
			elif parton == "ub": return self.ddh[2]
			elif parton == "d": return self.ddh[3]
			elif parton == "db": return self.ddh[4]
			elif parton == "s": return self.ddh[5]
			elif parton == "sb": return self.ddh[6]
			elif parton == "c": return self.ddh[7]
			elif parton == "cb": return self.ddh[8]
			elif parton == "b": return self.ddh[9]
			elif parton == "bb": return self.ddh[10]
####______________________________________________________________________
		elif hadron=='k-':
			self.ddh=np.array(dss_k.fdssda(2,-1,io,x,q**2))
			if parton == "gl": return self.ddh[0]
			elif parton == "u": return self.ddh[1]
			elif parton == "ub": return self.ddh[2]
			elif parton == "d": return self.ddh[3]
			elif parton == "db": return self.ddh[4]
			elif parton == "s": return self.ddh[5]
			elif parton == "sb": return self.ddh[6]
			elif parton == "c": return self.ddh[7]
			elif parton == "cb": return self.ddh[8]
			elif parton == "b": return self.ddh[9]
			elif parton == "bb": return self.ddh[10]
###___________________________________________________________________
		elif hadron=='k':
			self.ddh=np.array(dss_k.fdssda(2,1,io,x,q**2))
			if parton == "gl": return self.ddh[0]
			elif parton == "u": return self.ddh[1]
			elif parton == "ub": return self.ddh[2]
			elif parton == "d": return self.ddh[3]
			elif parton == "db": return self.ddh[4]
			elif parton == "s": return self.ddh[5]
			elif parton == "sb": return self.ddh[6]
			elif parton == "c": return self.ddh[7]
			elif parton == "cb": return self.ddh[8]
			elif parton == "b": return self.ddh[9]
			elif parton == "bb": return self.ddh[10]
	

	
#----- esempio
#	dss=dass_ff()
#	dup1=dss.D3('pi+',"u",0.5,91.2,0)
#	dss1=dass_ff()
#	dup2=dss1.D14('k+',"u",0.5,91.2,0)








