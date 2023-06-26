# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
import pandas as pd
#from iminuit import Minuit
#from def_crs import*
#from least_squares import*
import time 
#from fitter import*

y_lbl=48
x_lbl =50


ers_bar = 1.2  # dimensione barre errore
pnt_dat = 6 # dimensione punti

title_s=35 #dimensione titolo bin

fit_l=3.  # dimensione linea fit

coefs = [0.3,0.4,0.5]

for cf in coefs:

	dati_lp = pd.read_csv("grid_conv/dati_lp_conv_"+ str(cf)+".csv")
	dati_lk = pd.read_csv("grid_conv/dati_lk_conv_"+ str(cf)+".csv")
	dati_exp=pd.read_csv("lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')


	#################################

	#___ PLOT PER PIONI  +



	fig, axes = plt.subplots(1,4)
	z1=[0.25,.35,.45,.6]
	lim=[-0.15,0.055]
	ct=1

	z_err =[[0.05,0.05,0.05,0.1],[0.05,0.05,0.05,0.3]]

	for zs,ax in zip(z1,axes):
		
		fig.suptitle('$\Lambda$ - $\pi^+$ -- coef = '+str(cf),fontsize=30)	
		dt = dati_lp.loc[(dati_lp['had1']==300) & (dati_lp['z1']==zs)& (dati_lp['had2']==100)]

		pnt = dati_exp.loc[(dati_exp['h1']==300) & (dati_exp['z1']==zs)& (dati_exp['h2']==100)]	

		ax=plt.subplot(1,4,ct)
		ax.plot(dt.z2,dt.fit_line,label='fit_line',linewidth=fit_l)

		ax.plot(dt.z2,dt.conv,label='conv',linewidth=fit_l-0.3,linestyle='--',color='red')
	#	ax.plot(dt.z2,dt.down,label='down',linewidth=fit_l-0.3,linestyle='-.',color='blue')
	#	ax.plot(dt.z2,dt.strange,label='strange',linewidth=fit_l-0.3,linestyle=':',color='purple')
	#	ax.plot(dt.z2,dt.sea,label='sea',linewidth=fit_l-0.3,linestyle=(0, (3, 5, 1, 5, 1, 5)),color='green')
		axhline(linewidth=1.5, ls=':', color='g')

		ax.errorbar(pnt.z2, pnt.P_exp, pnt.err, z_err, fmt='o', markersize=4, color='blue',elinewidth=0.7, label= '$\Lambda$ - $\pi^+$')



		if zs == 0.25: title("0.2<$z_{\Lambda}$<0.3 ",fontsize=title_s,x=0.5, y=1.)
		if zs == 0.35: title("0.3<$z_{\Lambda}$<0.4 ",fontsize=title_s,x=0.5, y=1.)
		if zs == 0.45: title("0.4<$z_{\Lambda}$<0.5 ",fontsize=title_s,x=0.5, y=1.)
		if zs == 0.6: title("0.5<$z_{\Lambda}$<0.9 ",fontsize=title_s,x=0.5, y=1.)

		ax.tick_params(axis='both', which='major', labelsize=28)
		if ct >1: ax.set_yticklabels([])
		
		if ct ==1 :legend(loc='lower left', fontsize=15,frameon=True), ylabel('Polarization',size=y_lbl) 

		ax.set_ylim(lim)
		ct+=1


	fig.subplots_adjust(top=0.8,bottom=0.2,left=0.105,right=0.99,hspace=0.2,wspace=0.0)
	fig.set_size_inches(19.5, 10.5, forward=True)
	fig.savefig('grid_conv/Lb_pi_p_partial_'+ str(cf)+'.pdf')
	#fig.show()

	############################

	##_  PLOT PER KAONI +


	fig1, axes = plt.subplots(1,4)
	z1=[0.25,.35,.45,.6]
	#lim=[-0.10,0.055]
	ct=1
	for zs,ax in zip(z1,axes):
		fig1.suptitle('$\Lambda$ - $K^+$',fontsize=30)	

		dt = dati_lk.loc[(dati_lk['had1']==300) & (dati_lk['z1']==zs)& (dati_lk['had2']==200)]

		pnt = dati_exp.loc[(dati_exp['h1']==300) & (dati_exp['z1']==zs)& (dati_exp['h2']==200)]	

		
		ax=plt.subplot(1,4,ct)
		ax.plot(dt.z2,dt.fit_line,label='fit_line',linewidth=fit_l)

		ax.plot(dt.z2,dt.conv,label='conv',linewidth=fit_l-0.3,linestyle='--',color='red')
	#	ax.plot(dt.z2,dt.down,label='down',linewidth=fit_l-0.3,linestyle='-.',color='blue')
	#	ax.plot(dt.z2,dt.strange,label='strange',linewidth=fit_l-0.3,linestyle=':',color='purple')
	#	ax.plot(dt.z2,dt.sea,label='sea',linewidth=fit_l-0.3,linestyle=(0, (3, 5, 1, 5, 1, 5)),color='green')
		axhline(linewidth=1.5, ls=':', color='g')

		ax.errorbar(pnt.z2, pnt.P_exp, pnt.err, z_err, fmt='o', markersize=4, color='blue',elinewidth=0.7, label= '$\Lambda$ - $K^+$')


		if zs == 0.25: title("0.2<$z_{\Lambda}$<0.3 ",fontsize=title_s,x=0.5, y=1.)
		if zs == 0.35: title("0.3<$z_{\Lambda}$<0.4 ",fontsize=title_s,x=0.5, y=1.)
		if zs == 0.45: title("0.4<$z_{\Lambda}$<0.5 ",fontsize=title_s,x=0.5, y=1.)
		if zs == 0.6: title("0.5<$z_{\Lambda}$<0.9 ",fontsize=title_s,x=0.5, y=1.)

		ax.tick_params(axis='both', which='major', labelsize=28)
		if ct >1: ax.set_yticklabels([])
		
		if ct ==1 :legend(loc='lower left', fontsize=15,frameon=True), ylabel('Polarization',size=y_lbl) 

		ax.set_ylim(lim)
		ct+=1

		
	fig1.subplots_adjust(top=0.8,bottom=0.2,left=0.105,right=0.99,hspace=0.2,wspace=0.0)
	fig1.set_size_inches(19.5, 10.5, forward=True)
	#fig1.savefig('Lb_k_p_partial.pdf')
	#fig1.show()




	#################################

	#___ PLOT PER PIONI  -



	fig2, axes = plt.subplots(1,4)
	z1=[0.25,.35,.45,.6]
	lim=[-0.06,0.15]
	ct=1

	z_err =[[0.05,0.05,0.05,0.1],[0.05,0.05,0.05,0.3]]

	for zs,ax in zip(z1,axes):
		
		fig2.suptitle('$\Lambda$ - $\pi^-$ -- coef = '+str(cf),fontsize=30)	
		dt = dati_lp.loc[(dati_lp['had1']==300) & (dati_lp['z1']==zs)& (dati_lp['had2']==105)]

		pnt = dati_exp.loc[(dati_exp['h1']==300) & (dati_exp['z1']==zs)& (dati_exp['h2']==105)]	

		ax=plt.subplot(1,4,ct)
		ax.plot(dt.z2,dt.fit_line,label='fit_line',linewidth=fit_l)

		ax.plot(dt.z2,dt.conv,label='conv',linewidth=fit_l-0.3,linestyle='--',color='red')
	#	ax.plot(dt.z2,dt.down,label='down',linewidth=fit_l-0.3,linestyle='-.',color='blue')
	#	ax.plot(dt.z2,dt.strange,label='strange',linewidth=fit_l-0.3,linestyle=':',color='purple')
	#	ax.plot(dt.z2,dt.sea,label='sea',linewidth=fit_l-0.3,linestyle=(0, (3, 5, 1, 5, 1, 5)),color='green')
		axhline(linewidth=1.5, ls=':', color='g')

		ax.errorbar(pnt.z2, pnt.P_exp, pnt.err, z_err, fmt='o', markersize=4, color='blue',elinewidth=0.7, label= '$\Lambda$ - $\pi^-$')



		if zs == 0.25: title("0.2<$z_{\Lambda}$<0.3 ",fontsize=title_s,x=0.5, y=1.)
		if zs == 0.35: title("0.3<$z_{\Lambda}$<0.4 ",fontsize=title_s,x=0.5, y=1.)
		if zs == 0.45: title("0.4<$z_{\Lambda}$<0.5 ",fontsize=title_s,x=0.5, y=1.)
		if zs == 0.6: title("0.5<$z_{\Lambda}$<0.9 ",fontsize=title_s,x=0.5, y=1.)

		ax.tick_params(axis='both', which='major', labelsize=28)
		if ct >1: ax.set_yticklabels([])
		
		if ct ==1 :ylabel('Polarization',size=y_lbl) 
		if ct ==4 :legend(loc='upper right', fontsize=15,frameon=True) 

		ax.set_ylim(lim)
		ct+=1


	fig2.subplots_adjust(top=0.8,bottom=0.2,left=0.105,right=0.99,hspace=0.2,wspace=0.0)
	fig2.set_size_inches(19.5, 10.5, forward=True)
	fig2.savefig('grid_conv/Lb_pi_m_partial_'+ str(cf)+'.pdf')
	#fig2.show()


	###############################################à
	#___ PLOT PER KAONI  -



	fig3, axes = plt.subplots(1,4)
	z1=[0.25,.35,.45,.6]
	lim=[-0.06,0.15]
	ct=1

	z_err =[[0.05,0.05,0.05,0.1],[0.05,0.05,0.05,0.3]]

	for zs,ax in zip(z1,axes):
		
		fig3.suptitle('$\Lambda$ - $K^-$',fontsize=30)	
		dt = dati_lk.loc[(dati_lk['had1']==300) & (dati_lk['z1']==zs)& (dati_lk['had2']==205)]

		pnt = dati_exp.loc[(dati_exp['h1']==300) & (dati_exp['z1']==zs)& (dati_exp['h2']==205)]	

		ax=plt.subplot(1,4,ct)
		ax.plot(dt.z2,dt.fit_line,label='fit_line',linewidth=fit_l)

		ax.plot(dt.z2,dt.conv,label='conv',linewidth=fit_l-0.3,linestyle='--',color='red')
	#	ax.plot(dt.z2,dt.down,label='down',linewidth=fit_l-0.3,linestyle='-.',color='blue')
	#	ax.plot(dt.z2,dt.strange,label='strange',linewidth=fit_l-0.3,linestyle=':',color='purple')
	#	ax.plot(dt.z2,dt.sea,label='sea',linewidth=fit_l-0.3,linestyle=(0, (3, 5, 1, 5, 1, 5)),color='green')
		axhline(linewidth=1.5, ls=':', color='g')

		ax.errorbar(pnt.z2, pnt.P_exp, pnt.err, z_err, fmt='o', markersize=4, color='blue',elinewidth=0.7, label= '$\Lambda$ - $K^-$')



		if zs == 0.25: title("0.2<$z_{\Lambda}$<0.3 ",fontsize=title_s,x=0.5, y=1.)
		if zs == 0.35: title("0.3<$z_{\Lambda}$<0.4 ",fontsize=title_s,x=0.5, y=1.)
		if zs == 0.45: title("0.4<$z_{\Lambda}$<0.5 ",fontsize=title_s,x=0.5, y=1.)
		if zs == 0.6: title("0.5<$z_{\Lambda}$<0.9 ",fontsize=title_s,x=0.5, y=1.)

		ax.tick_params(axis='both', which='major', labelsize=28)
		if ct >1: ax.set_yticklabels([])
		
		if ct ==1 :ylabel('Polarization',size=y_lbl) 
		if ct ==4 :legend(loc='upper right', fontsize=15,frameon=True) 

		ax.set_ylim(lim)
		ct+=1


	fig3.subplots_adjust(top=0.8,bottom=0.2,left=0.105,right=0.99,hspace=0.2,wspace=0.0)
	fig3.set_size_inches(19.5, 10.5, forward=True)
	#fig3.savefig('Lb_k_m_partial.pdf')
	#fig3.show()


























