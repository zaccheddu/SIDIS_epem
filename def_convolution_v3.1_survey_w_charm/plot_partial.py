# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
import pandas as pd
from iminuit import Minuit
from def_crs import*
from least_squares import*
import time 
from fitter import*

y_lbl=48
x_lbl =50


ers_bar = 1.3  # dimensione barre errore
pnt_dat = 6 # dimensione punti

title_s=35 #dimensione titolo bin

fit_l=7  # dimensione linea fit
leg_sz=27

dati_lp = pd.read_csv("dati_lp_partial.csv")
dati_lk = pd.read_csv("dati_lk_partial.csv")
dati_exp=pd.read_csv("lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')


#################################

#___ PLOT PER PIONI  +



fig, axes = plt.subplots(1,4)
z1=[0.25,.35,.45,.6]
lim=[-0.15,0.055]
ct=1

z_err =[[0.05,0.05,0.05,0.1],[0.05,0.05,0.05,0.3]]

for zs,ax in zip(z1,axes):
	
	#fig.suptitle('$\Lambda$ - $\pi^+$',fontsize=30)	
	dt = dati_lp.loc[(dati_lp['had1']==300) & (dati_lp['z1']==zs)& (dati_lp['had2']==100)]

	pnt = dati_exp.loc[(dati_exp['h1']==300) & (dati_exp['z1']==zs)& (dati_exp['h2']==100)]	

	ax=plt.subplot(1,4,ct)
	ax.plot(dt.z2,dt.fit_line,label='_nolegend_',linewidth=fit_l-2.)

	ax.plot(dt.z2,dt.up,label='_nolegend_',linewidth=fit_l-0.3,linestyle='-',color='red')
	ax.plot(dt.z2,dt.down,label='_nolegend_',linewidth=fit_l-0.3,linestyle='--',color='blue')
	ax.plot(dt.z2,dt.strange,label='_nolegend_',linewidth=fit_l-0.3,linestyle='-.',color='purple')
	ax.plot(dt.z2,dt.sea,label='_nolegend_',linewidth=fit_l-0.3,linestyle=':',color='green')  # (0, (3, 5, 1, 5, 1, 5))
	axhline(linewidth=1.5, ls=':', color='g')

	ax.errorbar(pnt.z2, pnt.P_exp, pnt.err, z_err, fmt='o', markersize=6, color='blue',elinewidth=ers_bar, label= '$\Lambda$ - $\pi^+$')



	if zs == 0.25: title("0.2<$z_{\Lambda}$<0.3 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.35: title("0.3<$z_{\Lambda}$<0.4 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.45: title("0.4<$z_{\Lambda}$<0.5 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.6: title("0.5<$z_{\Lambda}$<0.9 ",fontsize=title_s,x=0.5, y=1.)

	ax.tick_params(axis='both', which='major', labelsize=28)
	if ct >1: ax.set_yticklabels([])
	
	if ct ==1 :legend(loc='lower left', fontsize=leg_sz,frameon=True), ylabel('Polarization',size=y_lbl) 
	#if ct ==1 :ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
	ax.set_ylim(lim)
	ct+=1
	fig.text(0.5, 0.1, r'$z_{\pi}$', ha='center',size=x_lbl)
	fig.text(0.6, 0.7, '(c)', ha='center',size=x_lbl)
fig.subplots_adjust(top=0.8,bottom=0.2,left=0.105,right=0.99,hspace=0.2,wspace=0.0)
fig.set_size_inches(19.5, 10.5, forward=True)

############################

##_  PLOT PER KAONI +


fig1, axes = plt.subplots(1,4)
z1=[0.25,.35,.45,.6]
#lim=[-0.10,0.055]
ct=1
for zs,ax in zip(z1,axes):
	#fig1.suptitle('$\Lambda$ - $K^+$',fontsize=30)	

	dt = dati_lk.loc[(dati_lk['had1']==300) & (dati_lk['z1']==zs)& (dati_lk['had2']==200)]

	pnt = dati_exp.loc[(dati_exp['h1']==300) & (dati_exp['z1']==zs)& (dati_exp['h2']==200)]	

	
	ax=plt.subplot(1,4,ct)
	ax.plot(dt.z2,dt.fit_line,label='_nolegend_',linewidth=fit_l-2.)

	ax.plot(dt.z2,dt.up,label='_nolegend_',linewidth=fit_l-0.3,linestyle='-',color='red')
	ax.plot(dt.z2,dt.down,label='_nolegend_',linewidth=fit_l-0.3,linestyle='--',color='blue')
	ax.plot(dt.z2,dt.strange,label='_nolegend_',linewidth=fit_l-0.3,linestyle='-.',color='purple')
	ax.plot(dt.z2,dt.sea,label='_nolegend_',linewidth=fit_l-0.3,linestyle=':',color='green')
	axhline(linewidth=1.5, ls=':', color='g')

	ax.errorbar(pnt.z2, pnt.P_exp, pnt.err, z_err, fmt='o', markersize=6, color='blue',elinewidth=ers_bar, label= '$\Lambda$ - $K^+$')


	if zs == 0.25: title("0.2<$z_{\Lambda}$<0.3 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.35: title("0.3<$z_{\Lambda}$<0.4 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.45: title("0.4<$z_{\Lambda}$<0.5 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.6: title("0.5<$z_{\Lambda}$<0.9 ",fontsize=title_s,x=0.5, y=1.)

	ax.tick_params(axis='both', which='major', labelsize=28)
	if ct >1: ax.set_yticklabels([])
	
	if ct ==1 :legend(loc='lower left', fontsize=leg_sz,frameon=True), ylabel('Polarization',size=y_lbl) 

	ax.set_ylim(lim)
	ct+=1
	fig1.text(0.5, 0.1, r'$z_{K}$', ha='center',size=x_lbl)
	fig1.text(0.6, 0.7, '(d)', ha='center',size=x_lbl)
fig1.subplots_adjust(top=0.8,bottom=0.2,left=0.105,right=0.99,hspace=0.2,wspace=0.0)
fig1.set_size_inches(19.5, 10.5, forward=True)




#################################

#___ PLOT PER PIONI  -



fig2, axes = plt.subplots(1,4)
z1=[0.25,.35,.45,.6]
lim=[-0.06,0.15]
ct=1

z_err =[[0.05,0.05,0.05,0.1],[0.05,0.05,0.05,0.3]]

for zs,ax in zip(z1,axes):
	
	#fig2.suptitle('$\Lambda$ - $\pi^-$',fontsize=30)	
	dt = dati_lp.loc[(dati_lp['had1']==300) & (dati_lp['z1']==zs)& (dati_lp['had2']==105)]

	pnt = dati_exp.loc[(dati_exp['h1']==300) & (dati_exp['z1']==zs)& (dati_exp['h2']==105)]	

	ax=plt.subplot(1,4,ct)
	ax.plot(dt.z2,dt.fit_line,label='_nolegend_',linewidth=fit_l-2.)

	ax.plot(dt.z2,dt.up,label='up',linewidth=fit_l-0.3,linestyle='-',color='red')
	ax.plot(dt.z2,dt.down,label='down',linewidth=fit_l-0.3,linestyle='--',color='blue')
	ax.plot(dt.z2,dt.strange,label='strange',linewidth=fit_l-0.3,linestyle='-.',color='purple')
	ax.plot(dt.z2,dt.sea,label='sea',linewidth=fit_l-0.3,linestyle=':',color='green')
	axhline(linewidth=1.5, ls=':', color='g')

	ax.errorbar(pnt.z2, pnt.P_exp, pnt.err, z_err, fmt='o', markersize=6, color='blue',elinewidth=ers_bar, label= '$\Lambda$ - $\pi^-$')



	if zs == 0.25: title("0.2<$z_{\Lambda}$<0.3 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.35: title("0.3<$z_{\Lambda}$<0.4 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.45: title("0.4<$z_{\Lambda}$<0.5 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.6: title("0.5<$z_{\Lambda}$<0.9 ",fontsize=title_s,x=0.5, y=1.)

	ax.tick_params(axis='both', which='major', labelsize=28)
	if ct >1: ax.set_yticklabels([])
	
	if ct ==1 :ylabel('Polarization',size=y_lbl) 
	if ct ==4 :legend(loc='upper right', fontsize=leg_sz,frameon=True) 

	ax.set_ylim(lim)
	ct+=1
	fig2.text(0.5, 0.1, r'$z_{\pi}$', ha='center',size=x_lbl)
	fig2.text(0.6, 0.7, '(a)', ha='center',size=x_lbl)
fig2.subplots_adjust(top=0.8,bottom=0.2,left=0.105,right=0.99,hspace=0.2,wspace=0.0)
fig2.set_size_inches(19.5, 10.5, forward=True)



###############################################à
#___ PLOT PER KAONI  -



fig3, axes = plt.subplots(1,4)
z1=[0.25,.35,.45,.6]
lim=[-0.06,0.15]
ct=1

z_err =[[0.05,0.05,0.05,0.1],[0.05,0.05,0.05,0.3]]

for zs,ax in zip(z1,axes):
	
	#fig3.suptitle('$\Lambda$ - $K^-$',fontsize=30)	
	dt = dati_lk.loc[(dati_lk['had1']==300) & (dati_lk['z1']==zs)& (dati_lk['had2']==205)]

	pnt = dati_exp.loc[(dati_exp['h1']==300) & (dati_exp['z1']==zs)& (dati_exp['h2']==205)]	

	ax=plt.subplot(1,4,ct)
	ax.plot(dt.z2,dt.fit_line,label='_nolegend_',linewidth=fit_l-2.)

	ax.plot(dt.z2,dt.up,label='_nolegend_',linewidth=fit_l-0.3,linestyle='-',color='red')
	ax.plot(dt.z2,dt.down,label='_nolegend_',linewidth=fit_l-0.3,linestyle='--',color='blue')
	ax.plot(dt.z2,dt.strange,label='_nolegend_',linewidth=fit_l-0.3,linestyle='-.',color='purple')
	ax.plot(dt.z2,dt.sea,label='_nolegend_',linewidth=fit_l-0.3,linestyle=':',color='green')
	axhline(linewidth=1.5, ls=':', color='g')

	ax.errorbar(pnt.z2, pnt.P_exp, pnt.err, z_err, fmt='o', markersize=6, color='blue',elinewidth=ers_bar, label= '$\Lambda$ - $K^-$')



	if zs == 0.25: title("0.2<$z_{\Lambda}$<0.3 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.35: title("0.3<$z_{\Lambda}$<0.4 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.45: title("0.4<$z_{\Lambda}$<0.5 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.6: title("0.5<$z_{\Lambda}$<0.9 ",fontsize=title_s,x=0.5, y=1.)

	ax.tick_params(axis='both', which='major', labelsize=28)
	if ct >1: ax.set_yticklabels([])
	
	if ct ==1 :ylabel('Polarization',size=y_lbl) 
	if ct ==4 :legend(loc='upper right', fontsize=leg_sz,frameon=True) 

	ax.set_ylim(lim)
	ct+=1
	fig3.text(0.5, 0.1, r'$z_{K}$', ha='center',size=x_lbl)
	fig3.text(0.6, 0.7, '(b)', ha='center',size=x_lbl)
fig3.subplots_adjust(top=0.8,bottom=0.2,left=0.105,right=0.99,hspace=0.2,wspace=0.0)
fig3.set_size_inches(19.5, 10.5, forward=True)


fig.savefig('Lb_pi_p_partial.pdf')
fig.show()


fig1.savefig('Lb_k_p_partial.pdf')
fig1.show()


fig2.savefig('Lb_pi_m_partial.pdf')
fig2.show()

fig3.savefig('Lb_k_m_partial.pdf')
fig3.show()


























