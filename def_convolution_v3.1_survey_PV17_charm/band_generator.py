#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
#warnings.filterwarnings("ignore")
import multiprocessing

# In[2]:


df1 = pd.read_csv('fit_parameters/to_use/fit_hadron_coef_0.27_chi_1.174__True_gk_PV17_su_no_charmno_correction_no.csv')
dfs1= pd.read_csv('fit_parameters/dfs_to_use/new_gauss_su_charm_no_no_def.csv')

df2 = pd.read_csv('fit_parameters/to_use/fit_hadron_coef_0.27_chi_1.259__True_gk_PV17_su_no_charmyes_correction_no.csv')
dfs2= pd.read_csv('fit_parameters/dfs_to_use/new_gauss_su_charm_no_yes_def.csv')

df3 = pd.read_csv('fit_parameters/to_use/fit_hadron_coef_0.27_chi_1.361__True_gk_PV17_su_yes_charmyes_correction_no.csv')
dfs3= pd.read_csv('fit_parameters/dfs_to_use/new_gauss_su_charm_yes_yes_def.csv')


# In[3]:


dati_lp = pd.read_csv("fit_parameters/lp_point.csv")
dati_lk = pd.read_csv("fit_parameters/lk_point.csv")


# In[ ]:





# In[4]:


def grids_lp(df,su2,charm,scale,coef,nf):
    mdl_den = 'pwr_lw_star'
    mdl_num = 'gauss'
    g_k_2h = 'PV17'
    coef = coef #0.27

    dati_lp = pd.read_csv("fit_parameters/lp_point.csv")


    fnc = polarization(coef)
    fnc.scale = scale
    fnc.mass = 1.115 
    fnc.frag2 = 'dss'
    fnc.mdl_num = mdl_num 
    fnc.mdl_den = mdl_den 
    fnc.charm =charm
    fnc.nf= nf
    
    fnc.g_k = g_k_2h 

    if su2=='no':

        f_prm=arange(0.,18.,1.)
        f_prm[0] = df['NUP']
        f_prm[1] = df['NDO']#par[0]  #do
        f_prm[2] = df['NST']  #par[0]      #st
        f_prm[3] = df['NSEA']  #par[0]      #upb
        f_prm[4] = f_prm[3] #par[0]      #dob
        f_prm[5] = f_prm[3] # par[0]      #stb
        #¯¯¯¯¯¯
        f_prm[6] = df['AUP']      #aup
        f_prm[7] = df['ADO']      #ado
        f_prm[8] = df['AST']  # par[0]      #ast
        f_prm[9] = df['ASEA']      #aupb
        f_prm[10] = f_prm[9] #par[0]      #adob
        f_prm[11] = f_prm[9] #par[0]      #astb
        #¯#¯¯¯¯¯¯¯
        f_prm[12] = df['BUP'] #par[0]      #bup
        f_prm[13] = df['BDO']      #bdo
        f_prm[14] = df['BST']      #bst
        f_prm[15] =  df['BSEA'] #par[0]      #bupb
        f_prm[16] = f_prm[15] #par[0]      #bdob
        f_prm[17] = f_prm[15] # par[0]      #bstb
        pt_pp=float(df['PP'])

    elif su2=='yes':

        f_prm=arange(0.,18.,1.)
        f_prm[0] = df['NUP']
        f_prm[1] = df['NUP']#par[0]  #do
        f_prm[2] = df['NST']  #par[0]      #st
        f_prm[3] = df['NSEA']  #par[0]      #upb
        f_prm[4] = df['NSEA'] #par[0]      #dob
        f_prm[5] = df['NDO']# par[0]      #stb
        #¯¯¯¯¯¯
        f_prm[6] = df['AUP']     #aup
        f_prm[7] = df['AUP']      #ado
        f_prm[8] = df['AST']  # par[0]      #ast
        f_prm[9] = df['ASEA']      #aupb
        f_prm[10] = df['ASEA'] #par[0]      #adob
        f_prm[11] = df['ADO'] #par[0]      #astb
        #¯#¯¯¯¯¯¯¯
        f_prm[12] = df['BUP'] #par[0]      #bup
        f_prm[13] = df['BUP']      #bdo
        f_prm[14] = df['BST']      #bst
        f_prm[15] =  df['BSEA'] #par[0]      #bupb
        f_prm[16] = df['BSEA'] #par[0]      #bdob
        f_prm[17] = df['BDO'] # par[0]      #bstb
        pt_pp=float(df['PP'])

    #print(f_prm)
    num=zeros(len(dati_lp.had1))
    i=0


    # In[6]:


    for hads1,hads2,zs1,zs2 in zip(dati_lp.had1,dati_lp.had2,dati_lp.z1,dati_lp.z2):

        if hads1 == 300 : had1='lbd'
        elif hads1 == 310 : had1='lbd_b'

        if hads2 == 100 : had2='pi+'
        elif hads2 == 105 : had2='pi-'
        elif hads2 == 200 : had2='k+'
        elif hads2 == 205 : had2='k-'

        num[i]= fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,0)
        i+=1

    dati_lp['conv'] = num

    dati_lp=dati_lp.drop(columns=['mins','maxx','fit_line'])
    
    return dati_lp


# In[5]:


def grids_lk(df,su2,charm,scale,coef,nf):
    mdl_den = 'pwr_lw_star'
    mdl_num = 'gauss'
    g_k_2h = 'PV17'
    coef = coef #0.27

    dati_lk = pd.read_csv("fit_parameters/lk_point.csv")


    fnc = polarization(coef)
    fnc.scale = scale
    fnc.mass = 1.115 
    fnc.frag2 = 'dss'
    fnc.mdl_num = mdl_num 
    fnc.mdl_den = mdl_den 
    fnc.charm =charm
    fnc.nf= nf
    
    fnc.g_k = g_k_2h 

    if su2=='no':

        f_prm=arange(0.,18.,1.)
        f_prm[0] = df['NUP']
        f_prm[1] = df['NDO']#par[0]  #do
        f_prm[2] = df['NST']  #par[0]      #st
        f_prm[3] = df['NSEA']  #par[0]      #upb
        f_prm[4] = f_prm[3] #par[0]      #dob
        f_prm[5] = f_prm[3] # par[0]      #stb
        #¯¯¯¯¯¯
        f_prm[6] = df['AUP']      #aup
        f_prm[7] = df['ADO']      #ado
        f_prm[8] = df['AST']  # par[0]      #ast
        f_prm[9] = df['ASEA']      #aupb
        f_prm[10] = f_prm[9] #par[0]      #adob
        f_prm[11] = f_prm[9] #par[0]      #astb
        #¯#¯¯¯¯¯¯¯
        f_prm[12] = df['BUP'] #par[0]      #bup
        f_prm[13] = df['BDO']      #bdo
        f_prm[14] = df['BST']      #bst
        f_prm[15] =  df['BSEA'] #par[0]      #bupb
        f_prm[16] = f_prm[15] #par[0]      #bdob
        f_prm[17] = f_prm[15] # par[0]      #bstb
        pt_pp=float(df['PP'])

    elif su2=='yes':

        f_prm=arange(0.,18.,1.)
        f_prm[0] = df['NUP']
        f_prm[1] = df['NUP']#par[0]  #do
        f_prm[2] = df['NST']  #par[0]      #st
        f_prm[3] = df['NSEA']  #par[0]      #upb
        f_prm[4] = df['NSEA'] #par[0]      #dob
        f_prm[5] = df['NDO']# par[0]      #stb
        #¯¯¯¯¯¯
        f_prm[6] = df['AUP']     #aup
        f_prm[7] = df['AUP']      #ado
        f_prm[8] = df['AST']  # par[0]      #ast
        f_prm[9] = df['ASEA']      #aupb
        f_prm[10] = df['ASEA'] #par[0]      #adob
        f_prm[11] = df['ADO'] #par[0]      #astb
        #¯#¯¯¯¯¯¯¯
        f_prm[12] = df['BUP'] #par[0]      #bup
        f_prm[13] = df['BUP']      #bdo
        f_prm[14] = df['BST']      #bst
        f_prm[15] =  df['BSEA'] #par[0]      #bupb
        f_prm[16] = df['BSEA'] #par[0]      #bdob
        f_prm[17] = df['BDO'] # par[0]      #bstb
        pt_pp=float(df['PP'])


    num=zeros(len(dati_lk.had1))
    i=0


    # In[6]:
    num=zeros(len(dati_lk.had1))
    i=0
    for hads1,hads2,zs1,zs2 in zip(dati_lk.had1,dati_lk.had2,dati_lk.z1,dati_lk.z2):
        if hads1 == 300 : had1='lbd'
        elif hads1 == 310 : had1='lbd_b'

        if hads2 == 100 : had2='pi+'
        elif hads2 == 105 : had2='pi-'
        elif hads2 == 200 : had2='k+'
        elif hads2 == 205 : had2='k-'

        num[i]= fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,0)
        i+=1


    dati_lk['conv'] = num

    dati_lk=dati_lk.drop(columns=['mins','maxx','fit_line'])
    
    return  dati_lk


# In[6]:


scale = 12.58
print(scale)
su2_1= 'no'
charm1='no'
su2_2= 'no'
charm2='yes'
su2_3= 'yes'
charm3='yes'

cut_dfs=500

dati_lp1=grids_lp(df1,su2_1,charm1,scale,0.27,3)
dati_lk1=grids_lk(df1,su2_1,charm1,scale,0.27,3)
dati_lp2=grids_lp(df2,su2_2,charm2,scale,0.27,4)
dati_lk2=grids_lk(df2,su2_2,charm2,scale,0.27,4)
dati_lp3=grids_lp(df3,su2_3,charm3,scale,0.27,4)
dati_lk3=grids_lk(df3,su2_3,charm3,scale,0.27,4)



# In[39]:


def grids_lp_bands(dati_lp,df_prm,su2,charm,scale,coef,nf):
	
    df_prm=df_prm.loc[(df_prm.index<cut_dfs)]
    mdl_den = 'pwr_lw_star'
    mdl_num = 'gauss'
    g_k_2h = 'PV17'
    coef = coef #0.27


    fnc = polarization(coef)
    fnc.scale = scale
    fnc.mass = 1.115 
    fnc.frag2 = 'dss'
    fnc.mdl_num = mdl_num 
    fnc.mdl_den = mdl_den 
    fnc.charm =charm
    fnc.nf= nf
    
    fnc.g_k = g_k_2h 

    num_tmp=zeros(len(df_prm.NUP))
    mins_lp = zeros(len(dati_lp.had1))
    maxx_lp = zeros(len(dati_lp.had1))
    i=0 
    for hads1,hads2,zs1,zs2 in zip(dati_lp.had1,dati_lp.had2,dati_lp.z1,dati_lp.z2):
        
        if hads1 == 300 : had1='lbd'
        elif hads1 == 310 : had1='lbd_b'

        if hads2 == 100 : had2='pi+'
        elif hads2 == 105 : had2='pi-'
        elif hads2 == 200 : had2='k+'
        elif hads2 == 205 : had2='k-'

    
        if su2=='no' and charm=='no':
            j=0
            for nup, ndo, nst, nsea, ast, bup, bsea, pp in zip(df_prm.NUP,df_prm.NDO,df_prm.NST,df_prm.NSEA,df_prm.AST, df_prm.BUP,df_prm.BSEA,df_prm.PP2):

                f_prm=arange(0.,18.,1.)
                f_prm[0] = nup
                f_prm[1] = ndo#par[0]  #do
                f_prm[2] = nst  #par[0]      #st
                f_prm[3] = nsea  #par[0]      #upb
                f_prm[4] = f_prm[3] #par[0]      #dob
                f_prm[5] = f_prm[3] # par[0]      #stb
                #¯¯¯¯¯¯
                f_prm[6] = 0 #par[0]      #aup
                f_prm[7] = 0 #par[0]      #ado
                f_prm[8] = ast  # par[0]      #ast
                f_prm[9] = 0# par[0]      #aupb
                f_prm[10] = f_prm[9] #par[0]      #adob
                f_prm[11] = f_prm[9] #par[0]      #astb
                #¯#¯¯¯¯¯¯¯
                f_prm[12] = bup #par[0]      #bup
                f_prm[13] = 0 #par[0]      #bdo
                f_prm[14] = 0# par[0]      #bst
                f_prm[15] = bsea #par[0]      #bupb
                f_prm[16] = f_prm[15] #par[0]      #bdob
                f_prm[17] = f_prm[15] # par[0]      #bstb
                pt_pp=pp
                #print(pt_pp)
                mss = 0# float(df['MSS'])

                num_tmp[j] = fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,mss)
                j+=1

        elif su2=='no' and charm=='yes':
            
            j=0
            for nup, ndo, nst, nsea, ado, ast, bup, bsea, pp in zip(df_prm.NUP,df_prm.NDO,df_prm.NST,df_prm.NSEA,df_prm.ADO,df_prm.AST, df_prm.BUP,df_prm.BSEA,df_prm.PP2):

                f_prm=arange(0.,18.,1.)
                f_prm[0] = nup
                f_prm[1] = ndo#par[0]  #do
                f_prm[2] = nst  #par[0]      #st
                f_prm[3] = nsea  #par[0]      #upb
                f_prm[4] = f_prm[3] #par[0]      #dob
                f_prm[5] = f_prm[3] # par[0]      #stb
                #¯¯¯¯¯¯
                f_prm[6] = 0 #par[0]      #aup
                f_prm[7] = ado #par[0]      #ado
                f_prm[8] = ast  # par[0]      #ast
                f_prm[9] = 0# par[0]      #aupb
                f_prm[10] = f_prm[9] #par[0]      #adob
                f_prm[11] = f_prm[9] #par[0]      #astb
                #¯#¯¯¯¯¯¯¯
                f_prm[12] = bup #par[0]      #bup
                f_prm[13] = 0 #par[0]      #bdo
                f_prm[14] = 0# par[0]      #bst
                f_prm[15] = bsea #par[0]      #bupb
                f_prm[16] = f_prm[15] #par[0]      #bdob
                f_prm[17] = f_prm[15] # par[0]      #bstb
                pt_pp=pp
                #print(pt_pp)
                mss = 0# float(df['MSS'])

                num_tmp[j] = fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,mss)
                j+=1

        elif su2=='yes' and charm=='yes':
            
            j=0
            for nup, ndo, nst, nsea,aup, ast, bup,bdo, pp in zip(df_prm.NUP,df_prm.NDO,df_prm.NST,df_prm.NSEA,df_prm.AUP,df_prm.AST,df_prm.BUP,df_prm.BDO,df_prm.PP2):

                f_prm=arange(0.,18.,1.)
                f_prm[0] = nup
                f_prm[1] = nup#par[0]  #do
                f_prm[2] = nst  #par[0]      #st
                f_prm[3] = nsea  #par[0]      #upb
                f_prm[4] = nsea #par[0]      #dob
                f_prm[5] = ndo # par[0]      #stb
                #¯¯¯¯¯¯
                f_prm[6] = aup #par[0]      #aup
                f_prm[7] = aup #par[0]      #ado
                f_prm[8] = ast  # par[0]      #ast
                f_prm[9] = 0# par[0]      #aupb
                f_prm[10] = 0 #par[0]      #adob
                f_prm[11] = 0 #par[0]      #astb
                #¯#¯¯¯¯¯¯¯
                f_prm[12] = bup #par[0]      #bup
                f_prm[13] = bup #par[0]      #bdo
                f_prm[14] = 0# par[0]      #bst
                f_prm[15] = 0 #par[0]      #bupb
                f_prm[16] = 0 #par[0]      #bdob
                f_prm[17] = bdo # par[0]      #bstb
                pt_pp=pp
                #print(pt_pp)
                mss = 0# float(df['MSS'])

                num_tmp[j] = fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,mss)
                j+=1

                
            
        mins_lp[i]=min(num_tmp)# fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,0)
        maxx_lp[i]=max(num_tmp)# fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,0)

        i+=1

    dati_lp['mins'] = mins_lp
    dati_lp['maxx'] = maxx_lp
    
    dati_lp.to_csv(r'fit_parameters/bands_/dati_lp_su_charm_'+str(su2)+'_'+str(charm)+'_scale_'+str(scale)+'_def.csv',index=False)
    
    #return dati_lp




# In[42]:


def grids_lk_bands(dati_lk,df_prm,su2,charm,scale,coef,nf):

    df_prm=df_prm.loc[(df_prm.index<cut_dfs)]
    mdl_den = 'pwr_lw_star'
    mdl_num = 'gauss'
    g_k_2h = 'PV17'
    coef = coef #0.27


    fnc = polarization(coef)
    fnc.scale = scale
    fnc.mass = 1.115 
    fnc.frag2 = 'dss'
    fnc.mdl_num = mdl_num 
    fnc.mdl_den = mdl_den 
    fnc.charm =charm
    fnc.nf= nf
    
    fnc.g_k = g_k_2h 

    num_tmp=zeros(len(df_prm.NUP))
    mins_lk = zeros(len(dati_lk.had1))
    maxx_lk = zeros(len(dati_lk.had1))
    i=0 
    for hads1,hads2,zs1,zs2 in zip(dati_lk.had1,dati_lk.had2,dati_lk.z1,dati_lk.z2):
        
        if hads1 == 300 : had1='lbd'
        elif hads1 == 310 : had1='lbd_b'

        if hads2 == 100 : had2='pi+'
        elif hads2 == 105 : had2='pi-'
        elif hads2 == 200 : had2='k+'
        elif hads2 == 205 : had2='k-'

    
        if su2=='no' and charm=='no':
            j=0
            for nup, ndo, nst, nsea, ast, bup, bsea, pp in zip(df_prm.NUP,df_prm.NDO,df_prm.NST,df_prm.NSEA,df_prm.AST,   df_prm.BUP,df_prm.BSEA,df_prm.PP2):

                f_prm=arange(0.,18.,1.)
                f_prm[0] = nup
                f_prm[1] = ndo#par[0]  #do
                f_prm[2] = nst  #par[0]      #st
                f_prm[3] = nsea  #par[0]      #upb
                f_prm[4] = f_prm[3] #par[0]      #dob
                f_prm[5] = f_prm[3] # par[0]      #stb
                #¯¯¯¯¯¯
                f_prm[6] = 0 #par[0]      #aup
                f_prm[7] = 0 #par[0]      #ado
                f_prm[8] = ast  # par[0]      #ast
                f_prm[9] = 0# par[0]      #aupb
                f_prm[10] = f_prm[9] #par[0]      #adob
                f_prm[11] = f_prm[9] #par[0]      #astb
                #¯#¯¯¯¯¯¯¯
                f_prm[12] = bup #par[0]      #bup
                f_prm[13] = 0 #par[0]      #bdo
                f_prm[14] = 0# par[0]      #bst
                f_prm[15] = bsea #par[0]      #bupb
                f_prm[16] = f_prm[15] #par[0]      #bdob
                f_prm[17] = f_prm[15] # par[0]      #bstb
                pt_pp=pp
                #print(pt_pp)
                mss = 0# float(df['MSS'])

                num_tmp[j] = fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,mss)
                j+=1

        elif su2=='no' and charm=='yes':
            
            j=0
            for nup, ndo, nst, nsea, ado, ast, bup, bsea, pp in zip(df_prm.NUP,df_prm.NDO,df_prm.NST,df_prm.NSEA,df_prm.ADO,df_prm.AST,  df_prm.BUP,df_prm.BSEA,df_prm.PP2):

                f_prm=arange(0.,18.,1.)
                f_prm[0] = nup
                f_prm[1] = ndo#par[0]  #do
                f_prm[2] = nst  #par[0]      #st
                f_prm[3] = nsea  #par[0]      #upb
                f_prm[4] = f_prm[3] #par[0]      #dob
                f_prm[5] = f_prm[3] # par[0]      #stb
                #¯¯¯¯¯¯
                f_prm[6] = 0 #par[0]      #aup
                f_prm[7] = ado #par[0]      #ado
                f_prm[8] = ast  # par[0]      #ast
                f_prm[9] = 0# par[0]      #aupb
                f_prm[10] = f_prm[9] #par[0]      #adob
                f_prm[11] = f_prm[9] #par[0]      #astb
                #¯#¯¯¯¯¯¯¯
                f_prm[12] = bup #par[0]      #bup
                f_prm[13] = 0 #par[0]      #bdo
                f_prm[14] = 0# par[0]      #bst
                f_prm[15] = bsea #par[0]      #bupb
                f_prm[16] = f_prm[15] #par[0]      #bdob
                f_prm[17] = f_prm[15] # par[0]      #bstb
                pt_pp=pp
                #print(pt_pp)
                mss = 0# float(df['MSS'])

                num_tmp[j] = fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,mss)
                j+=1

        elif su2=='yes' and charm=='yes':
            
            j=0
            for nup, ndo, nst, nsea, aup, ast, bup,bdo, pp in zip(df_prm.NUP,df_prm.NDO,df_prm.NST,df_prm.NSEA,df_prm.AUP,df_prm.AST,  df_prm.BUP,df_prm.BDO,df_prm.PP2):

                f_prm=arange(0.,18.,1.)
                f_prm[0] = nup
                f_prm[1] = nup#par[0]  #do
                f_prm[2] = nst  #par[0]      #st
                f_prm[3] = nsea  #par[0]      #upb
                f_prm[4] = nsea #par[0]      #dob
                f_prm[5] = ndo # par[0]      #stb
                #¯¯¯¯¯¯
                f_prm[6] = aup #par[0]      #aup
                f_prm[7] = aup #par[0]      #ado
                f_prm[8] = ast  # par[0]      #ast
                f_prm[9] = 0# par[0]      #aupb
                f_prm[10] = 0 #par[0]      #adob
                f_prm[11] = 0 #par[0]      #astb
                #¯#¯¯¯¯¯¯¯
                f_prm[12] = bup #par[0]      #bup
                f_prm[13] = bup #par[0]      #bdo
                f_prm[14] = 0# par[0]      #bst
                f_prm[15] = 0 #par[0]      #bupb
                f_prm[16] = 0 #par[0]      #bdob
                f_prm[17] = bdo # par[0]      #bstb
                pt_pp=pp
                #print(pt_pp)
                mss = 0# float(df['MSS'])

                num_tmp[j] = fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,mss)
                j+=1

                
            
        mins_lk[i]=min(num_tmp)# fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,0)
        maxx_lk[i]=max(num_tmp)# fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,0)

        i+=1

    dati_lk['mins'] = mins_lk
    dati_lk['maxx'] = maxx_lk

    dati_lk.to_csv(r'fit_parameters/bands_/dati_lk_su_charm_'+str(su2)+'_'+str(charm)+'_scale_'+str(scale)+'_def.csv',index=False)

    #return dati_lk

# In[40]:


#dfs_lp1=grids_lp_bands(dati_lp1,dfs1,su2_1,charm1,scale,0.27,3)
#dfs_lp2=grids_lp_bands(dati_lp2,dfs2,su2_2,charm2,scale,0.27,4)
#dfs_lp3=grids_lp_bands(dati_lp3,dfs3,su2_3,charm3,scale,0.27,4)

#p1 = multiprocessing.Process(target=grids_lp_bands,args=(dati_lp1,dfs1,su2_1,charm1,scale,0.27,3))
#p2 = multiprocessing.Process(target=grids_lp_bands,args=(dati_lp2,dfs2,su2_2,charm2,scale,0.27,4))
#p3 = multiprocessing.Process(target=grids_lp_bands,args=(dati_lp3,dfs3,su2_3,charm3,scale,0.27,4))


p4 = multiprocessing.Process(target=grids_lk_bands,args=(dati_lk1,dfs1,su2_1,charm1,scale,0.27,3))
p5 = multiprocessing.Process(target=grids_lk_bands,args=(dati_lk2,dfs2,su2_2,charm2,scale,0.27,4))
p6 = multiprocessing.Process(target=grids_lk_bands,args=(dati_lk3,dfs3,su2_3,charm3,scale,0.27,4))


# In[43]:


#dfs_lk1=grids_lk_bands(dati_lk1,dfs1,su2_1,charm1,scale,0.27,3)
#dfs_lk2=grids_lk_bands(dati_lk2,dfs2,su2_2,charm2,scale,0.27,4)
#dfs_lk3=grids_lk_bands(dati_lk3,dfs3,su2_3,charm3,scale,0.27,4)


# In[ ]:


#p1.start()
#p2.start()
#p3.start()
p4.start()
p5.start()
p6.start()


#p1.join()
#p2.join()
#p3.join()
p4.join()
p5.join()
p6.join()

