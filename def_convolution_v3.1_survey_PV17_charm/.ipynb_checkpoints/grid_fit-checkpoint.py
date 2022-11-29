#!/usr/bin/env python
# coding: utf-8

# In[1]:


# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.pyplot import*
from scipy import special
from def_crs import*
from model_fct import*
from g_funct import*
from FBT import FBT 
from model_fct import*
from def_conv_crs_1h import*
from def_conv_crs_2h import*

import pandas as pd
import sys 

# In[2]:
coef=np.float(sys.argv[1])#0.27
print(type(coef))

chi=sys.argv[2] #1.337


fit_type= sys.argv[3] #'hadron'
df = pd.read_csv("fit_parameters/fit_"+str(fit_type)+"_coef_"+str(coef)+"_chi_"+str(chi)+".csv") 
gk_type = sys.argv[4] # 'blny'

#print (sys.argv[1])
#print (sys.argv[1])

#coef=0.5
fnc = polarization(coef)
fnc.mass = 1.115
fnc.g_k = gk_type

print(fit_type)
print(gk_type)
# In[3]:


f_prm=arange(0.,18.,1.)
f_prm[0] = df['NUP']
f_prm[1] = df['NDO']#par[0]  #do
f_prm[2] = df['NST']  #par[0]      #st
f_prm[3] = df['NSEA']  #par[0]      #upb
f_prm[4] = f_prm[3] #par[0]      #dob
f_prm[5] = f_prm[3] # par[0]      #stb
#¯¯¯¯¯¯
f_prm[6] = 0 #par[0]      #aup
f_prm[7] = 0 #par[0]      #ado
f_prm[8] = df['AST']  # par[0]      #ast
f_prm[9] = 0# par[0]      #aupb
f_prm[10] = f_prm[9] #par[0]      #adob
f_prm[11] = f_prm[9] #par[0]      #astb
#¯#¯¯¯¯¯¯¯
f_prm[12] = df['BUP'] #par[0]      #bup
f_prm[13] = 0 #par[0]      #bdo
f_prm[14] = 0# par[0]      #bst
f_prm[15] =  df['BSEA'] #par[0]      #bupb
f_prm[16] = f_prm[15] #par[0]      #bdob
f_prm[17] = f_prm[15] # par[0]      #bstb
pt_pp=float(df['PP'])


# In[4]:


dati_lp = pd.read_csv("set_bande_pols/dati_lp_flh.csv")
dati_lk = pd.read_csv("set_bande_pols/dati_lk_flh.csv")


# In[5]:


num=zeros(len(dati_lp.had1))
i=0


# In[6]:


for hads1,hads2,zs1,zs2 in zip(dati_lp.had1,dati_lp.had2,dati_lp.z1,dati_lp.z2):

    if hads1 == 300 : had1='lbd'
    elif hads1 == 310 : had1='lbd_b'

    if hads2 == 100 : had2='pi+'
    elif hads2 == 105 : had2='pi-'
    elif hads2 == 200 : had2='k-'
    elif hads2 == 205 : had2='k-'

    num[i]= fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp)
    i+=1


# In[7]:


dati_lp['conv'] = num


# In[8]:


dati_lp
dati_lp.to_csv(r'fit_grid_plot/dati_lp_conv_'+ str(coef)+'chi_'+str(chi)+'.csv',index=False)


# In[9]:


num=zeros(len(dati_lk.had1))
i=0
for hads1,hads2,zs1,zs2 in zip(dati_lk.had1,dati_lk.had2,dati_lk.z1,dati_lk.z2):
    if hads1 == 300 : had1='lbd'
    elif hads1 == 310 : had1='lbd_b'

    if hads2 == 100 : had2='pi+'
    elif hads2 == 105 : had2='pi-'
    elif hads2 == 200 : had2='k+'
    elif hads2 == 205 : had2='k-'

    num[i]= fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp)
    i+=1


# In[10]:


dati_lk['conv'] = num
dati_lk.to_csv(r'fit_grid_plot/dati_lk_conv_'+ str(coef)+'chi_'+str(chi)+'.csv',index=False)

######

dati_lj = pd.read_csv("set_bande_pols/dati_lj_flh.csv")
dati_lj  = dati_lj.loc[~(dati_lj['pt']==0.0 )]	
fnc = polarization_1h()
fnc.mass = 1.115
fnc.g_k= gk_type


num=zeros(len(dati_lj.had1))
i=0

for hads1,zs1,pt in zip(dati_lj.had1,dati_lj.z1,dati_lj.pt):

    if hads1 == 300 : had1='lbd'
    elif hads1 == 310 : had1='lbd_b'

    #print(i)
    num[i] = fnc.ratio(had1,zs1,pt,f_prm,pt_pp)  
    i+=1
	
dati_lj['conv'] = num
dati_lj.to_csv(r'fit_grid_plot/dati_lj_conv_'+ str(coef)+'chi_'+str(chi)+'.csv',index=False)


fnc.ratio('lbd',0.35,0.2,f_prm,pt_pp) 


