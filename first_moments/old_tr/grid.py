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
from def_conv_crs import*
import pandas as pd


# In[2]:

coef=0.5
fnc = polarization(coef)
fnc.mass = 1.115


# In[3]:


f_prm=arange(0.,18.,1.)
f_prm[0] = 0.67486074
f_prm[1] = -0.65605671#par[0]  #do
f_prm[2] = -0.99999975  #par[0]      #st
f_prm[3] = -0.39721677  #par[0]      #upb
f_prm[4] = f_prm[3] #par[0]      #dob
f_prm[5] = f_prm[3] # par[0]      #stb
#¯¯¯¯¯¯
f_prm[6] = 0 #par[0]      #aup
f_prm[7] = 0 #par[0]      #ado
f_prm[8] = 2.1669248  # par[0]      #ast
f_prm[9] = 0# par[0]      #aupb
f_prm[10] = f_prm[9] #par[0]      #adob
f_prm[11] = f_prm[9] #par[0]      #astb
#¯#¯¯¯¯¯¯¯
f_prm[12] = 3.35974926 #par[0]      #bup
f_prm[13] = 0 #par[0]      #bdo
f_prm[14] = 0# par[0]      #bst
f_prm[15] =  2.1412581 #par[0]      #bupb
f_prm[16] = f_prm[15] #par[0]      #bdob
f_prm[17] = f_prm[15] # par[0]      #bstb
pt_pp=0.05622766


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
dati_lp.to_csv(r'dati_lp_conv_'+ str(coef)+'.csv',index=False)


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
dati_lk.to_csv(r'dati_lk_conv_'+ str(coef)+'.csv',index=False)







