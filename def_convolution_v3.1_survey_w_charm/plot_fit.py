# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
import pandas as pd
import time
import sys
import warnings
warnings.filterwarnings('ignore')

y_lbl=48
x_lbl =50


ers_bar = 1.2  # dimensione barre errore
pnt_dat = 6 # dimensione punti

title_s=35 #dimensione titolo bin

fit_l=3.  # dimensione linea fit

#coef=0.27
#chi=1.337

coef=  0.25 # np.float(sys.argv[1])#0.27
print(type(coef))
chi= 1.214  #sys.argv[2]# 1.337

dati_lp = pd.read_csv('fit_grid_plot/dati_lp_conv_'+ str(coef)+'chi_'+str(chi)+'.csv')
dati_lk = pd.read_csv('fit_grid_plot/dati_lk_conv_'+ str(coef)+'chi_'+str(chi)+'.csv')
dati_exp=pd.read_csv("exp_data/lambda_had_global.dat", delimiter=r"\s+", header=0, engine='python')


#undati_exp = dati_exp.loc[(dati_exp['z2']>0.5 )]
#dati_exp = dati_exp.loc[(dati_exp['z2']<0.5 )]

fig, axes = plt.subplots(1,4)
z1=[0.25,.35,.45,.6]
lim=[-0.15,0.055]
ct=1

z_err =[[0.05,0.05,0.05,0.1],[0.05,0.05,0.05,0.3]]

for zs,ax in zip(z1,axes):
	
	fig.suptitle('$\Lambda$ - $\pi^+$ -- coef = '+str(coef)+' $\chi^2_{dof}$ = '+str(chi),fontsize=30)	
	dt = dati_lp.loc[(dati_lp['had1']==300) & (dati_lp['z1']==zs)& (dati_lp['had2']==100)]

	pnt = dati_exp.loc[(dati_exp['h1']==300) & (dati_exp['z1']==zs)& (dati_exp['h2']==100)]	
	#pnt_white = undati_exp.loc[(undati_exp['h1']==300) & (undati_exp['z1']==zs)& (undati_exp['h2']==100)]	


	ax=plt.subplot(1,4,ct)
	ax.plot(dt.z2,dt.fit_line,label='reference fit',linewidth=fit_l)

	ax.plot(dt.z2,dt.conv,label='new fit',linewidth=fit_l-0.3,linestyle='--',color='red')
#	ax.plot(dt.z2,dt.down,label='down',linewidth=fit_l-0.3,linestyle='-.',color='blue')
#	ax.plot(dt.z2,dt.strange,label='strange',linewidth=fit_l-0.3,linestyle=':',color='purple')
#	ax.plot(dt.z2,dt.sea,label='sea',linewidth=fit_l-0.3,linestyle=(0, (3, 5, 1, 5, 1, 5)),color='green')
	axhline(linewidth=1.5, ls=':', color='g')

	ax.errorbar(pnt.z2, pnt.P_exp, pnt.err, z_err, fmt='o', markersize=4, color='blue',elinewidth=0.7, label= '$\Lambda$ - $\pi^+$')
	#xlabel("$z_{\pi}$",size=12)


	if zs == 0.25: title("0.2<$z_{\Lambda}$<0.3 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.35: title("0.3<$z_{\Lambda}$<0.4 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.45: title("0.4<$z_{\Lambda}$<0.5 ",fontsize=title_s,x=0.5, y=1.)
	if zs == 0.6: title("0.5<$z_{\Lambda}$<0.9 ",fontsize=title_s,x=0.5, y=1.)

	ax.tick_params(axis='both', which='major', labelsize=28)
	if ct >1: ax.set_yticklabels([])
	
	if ct ==1 :legend(loc='lower left', fontsize=15,frameon=True), ylabel('Polarization',size=y_lbl) 

	ax.set_ylim(lim)
	ct+=1
	fig.text(0.5, 0.1, r'$z_{\pi}$', ha='center',size=x_lbl)

fig.subplots_adjust(top=0.8,bottom=0.2,left=0.105,right=0.99,hspace=0.2,wspace=0.0)
fig.set_size_inches(19.5, 10.5, forward=True)
#fig.savefig('grid_conv/Lb_pi_p_partial_'+ str(cf)+'.pdf')
#fig.show()



############################

##_  PLOT  KAON +


fig1, axes = plt.subplots(1,4)
z1=[0.25,.35,.45,.6]
#lim=[-0.10,0.055]
ct=1
for zs,ax in zip(z1,axes):
	fig1.suptitle('$\Lambda$ - $K^+$-- coef = '+str(coef)+' $\chi^2_{dof}$ = '+str(chi),fontsize=30)	

	dt = dati_lk.loc[(dati_lk['had1']==300) & (dati_lk['z1']==zs)& (dati_lk['had2']==200)]

	pnt = dati_exp.loc[(dati_exp['h1']==300) & (dati_exp['z1']==zs)& (dati_exp['h2']==200)]	
	#pnt_white = undati_exp.loc[(undati_exp['h1']==300) & (undati_exp['z1']==zs)& (undati_exp['h2']==200)]	

	
	ax=plt.subplot(1,4,ct)
	ax.plot(dt.z2,dt.fit_line,label='reference fit',linewidth=fit_l)

	ax.plot(dt.z2,dt.conv,label='new fit',linewidth=fit_l-0.3,linestyle='--',color='red')
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
	fig1.text(0.5, 0.1, r'$z_K$', ha='center',size=x_lbl)
	
fig1.subplots_adjust(top=0.8,bottom=0.2,left=0.105,right=0.99,hspace=0.2,wspace=0.0)
fig1.set_size_inches(19.5, 10.5, forward=True)
#fig1.savefig('Lb_k_p_partial.pdf')
#fig1.show()

#################################

#___ PLOT  PION  -



fig2, axes = plt.subplots(1,4)
z1=[0.25,.35,.45,.6]
lim=[-0.06,0.15]
ct=1

z_err =[[0.05,0.05,0.05,0.1],[0.05,0.05,0.05,0.3]]

for zs,ax in zip(z1,axes):
	
	fig2.suptitle('$\Lambda$ - $\pi^-$ coef = '+str(coef)+' $\chi^2_{dof}$ = '+str(chi),fontsize=30)	
	dt = dati_lp.loc[(dati_lp['had1']==300) & (dati_lp['z1']==zs)& (dati_lp['had2']==105)]

	pnt = dati_exp.loc[(dati_exp['h1']==300) & (dati_exp['z1']==zs)& (dati_exp['h2']==105)]	
	#pnt_white = undati_exp.loc[(undati_exp['h1']==300) & (undati_exp['z1']==zs)& (undati_exp['h2']==105)]	

	ax=plt.subplot(1,4,ct)
	ax.plot(dt.z2,dt.fit_line,label='reference fit',linewidth=fit_l)

	ax.plot(dt.z2,dt.conv,label='new fit',linewidth=fit_l-0.3,linestyle='--',color='red')
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
	fig2.text(0.5, 0.1, r'$z_{\pi}$', ha='center',size=x_lbl)

fig2.subplots_adjust(top=0.8,bottom=0.2,left=0.105,right=0.99,hspace=0.2,wspace=0.0)
fig2.set_size_inches(19.5, 10.5, forward=True)
#fig2.savefig('grid_conv/Lb_pi_m_partial_'+ str(cf)+'.pdf')
#fig2.show()



###############################################à
#___ PLOT  KAON  -



fig3, axes = plt.subplots(1,4)
z1=[0.25,.35,.45,.6]
lim=[-0.06,0.15]
ct=1

z_err =[[0.05,0.05,0.05,0.1],[0.05,0.05,0.05,0.3]]

for zs,ax in zip(z1,axes):
	
	fig3.suptitle('$\Lambda$ - $K^-$ coef = '+str(coef)+' $\chi^2_{dof}$ = '+str(chi),fontsize=30)	
	dt = dati_lk.loc[(dati_lk['had1']==300) & (dati_lk['z1']==zs)& (dati_lk['had2']==205)]

	pnt = dati_exp.loc[(dati_exp['h1']==300) & (dati_exp['z1']==zs)& (dati_exp['h2']==205)]	
	#pnt_white = undati_exp.loc[(undati_exp['h1']==300) & (undati_exp['z1']==zs)& (undati_exp['h2']==205)]	

	ax=plt.subplot(1,4,ct)
	ax.plot(dt.z2,dt.fit_line,label='reference fit',linewidth=fit_l)

	ax.plot(dt.z2,dt.conv,label='new fit',linewidth=fit_l-0.3,linestyle='--',color='red')
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
	fig3.text(0.5, 0.1, r'$z_K$', ha='center',size=x_lbl)

fig3.subplots_adjust(top=0.8,bottom=0.2,left=0.105,right=0.99,hspace=0.2,wspace=0.0)
fig3.set_size_inches(19.5, 10.5, forward=True)
#fig3.savefig('Lb_k_m_partial.pdf')
#fig3.show()



##################################################

#dati_exp_lj=pd.read_csv("exp_data/lambda_jet_global.dat", delimiter=r"\s+", header=0, engine='python')

dati_exp_lj=pd.read_csv("exp_data/lambda_jet_global_err.dat")


ml = 1.115
dati_exp_lj['zp'] = dati_exp_lj['z1']#*sqrt(1- 4*ml**2/10.58**2/dati_exp_lj['z1']**2)
dati_exp_lj['qt'] =dati_exp_lj['pt']/dati_exp_lj['zp']
coefs=0.5

undati_exp_lj = dati_exp_lj.loc[(dati_exp_lj['qt']>10.58*coefs ) | (dati_exp_lj['P_exp']>0.05 )]

dati_exp_lj = dati_exp_lj.loc[~(dati_exp_lj['P_exp']>0.05 )]
dati_exp_lj = dati_exp_lj.loc[~(dati_exp_lj['qt']>10.58*coefs )]



dati_lj = pd.read_csv('fit_grid_plot/dati_lj_conv_'+ str(coef)+'chi_'+str(chi)+'.csv')
fig4, axes = plt.subplots(1,4)

cl_r='red'
cl_b='mediumblue'

lim=[-0.15,0.15]
z1=[0.25,.35,.45,.6]
#z1=[.35,.45,.6]

ct=1
#pt1_err =[[0.13,0.14,0.1,0.08],[0.07,0.16,0.2,0.72]]
#pt2_err =[[0.13,0.14,0.12,0.16],[0.07,0.16,0.18,0.64]]
#pt3_err =[[0.13,0.15,0.13,0.19],[0.07,0.15,0.17,0.61]]
#pt4_err =[[0.13,0.15,0.13,0.24],[0.07,0.15,0.17,0.56]]

for zs, ax in zip(z1,axes):
	#print(zs)
	#zmin=zs-0.05
	#zmin=zs-0.25
#	if ct==1 : pt_err =[[0.13,0.14,0.1,0.08],[0.07,0.16,0.2,0.72]]
#	elif ct==2 :pt_err =[[0.13,0.14,0.12,0.16],[0.07,0.16,0.18,0.64]]
#	elif ct==3 :pt_err =[[0.13,0.15,0.13,0.19],[0.07,0.15,0.17,0.61]]
#	elif ct==4 :pt_err =[[0.13,0.15,0.13,0.24],[0.07,0.15,0.17,0.56]]	
	#print(ax)
	lim=[-0.15,0.15]
	
	fig4.suptitle('$\Lambda$ - thrust coef = '+str(coef)+' $\chi^2_{dof}$ = '+str(chi),fontsize=30)
	dt = dati_lj.loc[(dati_lj['had1']==300) & (dati_lj['z1']==zs)] 
	
	pnt1 = dati_exp_lj.loc[(dati_exp_lj['h1']==300) & (dati_exp_lj['z1']< zs+0.05) & (dati_exp_lj['z1']> zs-0.05)]  
	pnt2 = dati_exp_lj.loc[(dati_exp_lj['h1']==310) & (dati_exp_lj['z1']< zs+0.05) & (dati_exp_lj['z1']> zs-0.05)] 

	unpnt1 = undati_exp_lj.loc[(undati_exp_lj['h1']==300) & (undati_exp_lj['z1']< zs+0.05) & (undati_exp_lj['z1']> zs-0.05)]  
	unpnt2 = undati_exp_lj.loc[(undati_exp_lj['h1']==310) & (undati_exp_lj['z1']< zs+0.05) & (undati_exp_lj['z1']> zs-0.05)] 


	pnt1.astype(float)
	pnt2.astype(float)
	#print(pnt1)
	#print(pnt2)
	ax=plt.subplot(1,4,ct)
	#ax.plot(dt.pt,dt.fit_line,label='fit_line',linewidth=fit_l)
	#if ct !=1 :
	ax.plot(dt.pt,dt.conv,label='new fit',linewidth=fit_l-0.3,linestyle='--',color='red')
	ax.set_ylim(lim)

	pt1_err=(pnt1['err_x_l'].tolist(),pnt1['err_x_r'].tolist())
	pt2_err=(pnt2['err_x_l'].tolist(),pnt2['err_x_r'].tolist())

	ax.errorbar(pnt1.pt, pnt1.P_exp, pnt1.err,pt1_err,  fmt='o', markersize=4, color=cl_r,elinewidth=0.7,label= '$\Lambda$ $')	
	ax.errorbar(pnt2.pt, pnt2.P_exp, pnt2.err, pt2_err, fmt='o', markersize=4, color=cl_b,elinewidth=0.7, label= '$\Lambda$ $')	

	unpt1_err=(unpnt1['err_x_l'].tolist(),unpnt1['err_x_r'].tolist())
	unpt2_err=(unpnt2['err_x_l'].tolist(),unpnt2['err_x_r'].tolist())

	ax.set_ylim(lim)
	ax.errorbar(unpnt1.pt, unpnt1.P_exp, unpnt1.err, unpt1_err,  fmt='o', markersize=4, color=cl_r,elinewidth=0.7,mfc='White',label= '$\Lambda$ $')	
	ax.errorbar(unpnt2.pt, unpnt2.P_exp, unpnt2.err, unpt2_err, fmt='o', markersize=4, color=cl_b,elinewidth=0.7,mfc='White', label= '$\Lambda$ $')	
	tick_params(axis='both', which='major', labelsize=25)

	if zs==0.25 :	title("0.2<$z_{\Lambda}$<0.3  ",fontsize=title_s,x=0.5, y=1.)
	elif zs==0.35 :	title("0.3<$z_{\Lambda}$<0.4",fontsize=title_s,x=0.5, y=1.)
	elif zs==0.45 :	title("0.4<$z_{\Lambda}$<0.5",fontsize=title_s,x=0.5, y=1.)
	elif zs==0.6 :	title("0.5<$z_{\Lambda}$<0.9",fontsize=title_s,x=0.5, y=1.)

	if ct==1 : ylabel('Polarization',size=y_lbl)
	if ct !=1 :ax.set_yticklabels([])
	axhline(linewidth=1.5, ls=':', color='g')
	ct+=1
	fig4.text(0.5, 0.1, r'$p_{\perp}$', ha='center',size=x_lbl)
fig.savefig('fit_grid_plot/pres_/Lb_pi_p_partial_'+ str(coef)+'_chi_'+str(chi)+'.pdf')
fig1.savefig('fit_grid_plot/pres_/Lb_k_p_partial_'+ str(coef)+'_chi_'+str(chi)+'.pdf')
fig2.savefig('fit_grid_plot/pres_/Lb_pi_m_partial_'+ str(coef)+'_chi_'+str(chi)+'.pdf')
fig3.savefig('fit_grid_plot/pres_/Lb_k_m_partial_'+ str(coef)+'_chi_'+str(chi)+'.pdf')

fig.show()
fig1.show()
fig2.show()
fig3.show()


fig4.subplots_adjust(top=0.8,bottom=0.2,left=0.1,right=0.99,hspace=0.2,wspace=0.0)

fig4.set_size_inches(19.5, 10.5, forward=True)
fig4.savefig('fit_grid_plot/pres_/LB_thrust_'+ str(coef)+'_chi_'+str(chi)+'_4.pdf')
fig4.show()




