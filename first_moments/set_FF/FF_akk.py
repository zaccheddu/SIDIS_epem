# -*- coding: utf-8 -*-
#!/home/zack/Documenti/AKK_f2py/__main__.py
import sys
import os
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*

percorso = os.getcwd()
path = os.path.join(percorso, "set_FF")
#print(percorso)
sys.path.insert(1, path)
#sys.path.insert(1, '/home/zackmrc/Documenti/git_project/lambda_fit/def_convolution/set_FF')
import akkff
##__________________________________________________________________________
#  DEFINISCO LE CLASSI PER LE FRAMMENTAZIONI

class akk_ff():


	def D1(self,hadron,parton,z,q):
		if hadron=='lbd':######  LAMBDA FF
			self.ddh=np.array(akkff.akk(1,z,q))
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
###_______________________________________________________________________
		elif hadron == 'pi'  :###### 	PION FF
			self.ddh=np.array(akkff.akk(2,z,q))
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
###_________________________________________________________________________
		elif hadron == 'k'   :###### KAON FF
			self.ddh=np.array(akkff.akk(3,z,q))
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
###_________________________________________________________________________
		elif hadron == 'piCSA' :###### PION CSA
			self.ddh=np.array(akkff.akk(4,z,q))
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
###_________________________________________________________________________
		elif hadron == 'kCSA'   :###### KAON CSA
			self.ddh=np.array(akkff.akk(5,z,q))
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
###_________________________________________________________________________
		elif hadron == 'pi+'   :###### pion +
			self.ddh1=np.array(akkff.akk(2,z,q))
			self.ddh2=np.array(akkff.akk(4,z,q))
			if parton == "gl": return (self.ddh1[0] + self.ddh2[0])/2
			elif parton == "u": return (self.ddh1[1] + self.ddh2[1])/2
			elif parton == "ub": return (self.ddh1[2] + self.ddh2[2])/2
			elif parton == "d": return (self.ddh1[3] + self.ddh2[3])/2
			elif parton == "db": return (self.ddh1[4] + self.ddh2[4])/2
			elif parton == "s": return (self.ddh1[5] + self.ddh2[5])/2
			elif parton == "sb": return (self.ddh1[6] + self.ddh2[6])/2
			elif parton == "c": return (self.ddh1[7] + self.ddh2[7])/2
			elif parton == "cb": return (self.ddh1[8] + self.ddh2[8])/2
			elif parton == "b": return (self.ddh1[9] + self.ddh2[9])/2
			elif parton == "bb": return (self.ddh1[10] + self.ddh2[10])/2
###_________________________________________________________________________
		elif hadron == 'pi-'   :###### pion +
			self.ddh1=np.array(akkff.akk(2,z,q))
			self.ddh2=np.array(akkff.akk(4,z,q))
			if parton == "gl": return (self.ddh1[0] - self.ddh2[0])/2
			elif parton == "u": return (self.ddh1[1] - self.ddh2[1])/2
			elif parton == "ub": return (self.ddh1[2] - self.ddh2[2])/2
			elif parton == "d": return (self.ddh1[3] - self.ddh2[3])/2
			elif parton == "db": return (self.ddh1[4] - self.ddh2[4])/2
			elif parton == "s": return (self.ddh1[5] - self.ddh2[5])/2
			elif parton == "sb": return (self.ddh1[6] - self.ddh2[6])/2
			elif parton == "c": return (self.ddh1[7] - self.ddh2[7])/2
			elif parton == "cb": return (self.ddh1[8] - self.ddh2[8])/2
			elif parton == "b": return (self.ddh1[9] - self.ddh2[9])/2
			elif parton == "bb": return (self.ddh1[10] - self.ddh2[10])/2
###_________________________________________________________________________
		elif hadron == 'k+'   :###### kaon +
			self.ddh1=np.array(akkff.akk(3,z,q))
			self.ddh2=np.array(akkff.akk(5,z,q))
			if parton == "gl": return (self.ddh1[0] + self.ddh2[0])/2
			elif parton == "u": return (self.ddh1[1] + self.ddh2[1])/2
			elif parton == "ub": return (self.ddh1[2] + self.ddh2[2])/2
			elif parton == "d": return (self.ddh1[3] + self.ddh2[3])/2
			elif parton == "db": return (self.ddh1[4] + self.ddh2[4])/2
			elif parton == "s": return (self.ddh1[5] + self.ddh2[5])/2
			elif parton == "sb": return (self.ddh1[6] + self.ddh2[6])/2
			elif parton == "c": return (self.ddh1[7] + self.ddh2[7])/2
			elif parton == "cb": return (self.ddh1[8] + self.ddh2[8])/2
			elif parton == "b": return (self.ddh1[9] + self.ddh2[9])/2
			elif parton == "bb": return (self.ddh1[10] + self.ddh2[10])/2
###_________________________________________________________________________
		elif hadron == 'k-'   :###### kaon+
			self.ddh1=np.array(akkff.akk(3,z,q))
			self.ddh2=np.array(akkff.akk(5,z,q))
			if parton == "gl": return (self.ddh1[0] - self.ddh2[0])/2
			elif parton == "u": return (self.ddh1[1] - self.ddh2[1])/2
			elif parton == "ub": return (self.ddh1[2] - self.ddh2[2])/2
			elif parton == "d": return (self.ddh1[3] - self.ddh2[3])/2
			elif parton == "db": return (self.ddh1[4] - self.ddh2[4])/2
			elif parton == "s": return (self.ddh1[5] - self.ddh2[5])/2
			elif parton == "sb": return (self.ddh1[6] - self.ddh2[6])/2
			elif parton == "c": return (self.ddh1[7] - self.ddh2[7])/2
			elif parton == "cb": return (self.ddh1[8] - self.ddh2[8])/2
			elif parton == "b": return (self.ddh1[9] - self.ddh2[9])/2
			elif parton == "bb": return (self.ddh1[10] - self.ddh2[10])/2
###_________________________________________________________________________

#----- esempio
#	dss=akk_ff()
#	dup1=dss.D1('lbd',"u",0.5,91.2)

  




























