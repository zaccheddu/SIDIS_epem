#!/usr/bin/env python
# coding: utf-8

# In[3]:


# -*- coding: utf-8 -*-
from pylab import*
import matplotlib.pyplot as plt
from numpy import*
from matplotlib.pyplot import*
import pandas as pd
import time
import sys
from def_conv_crs_2h import*
from numpy import*
import time 
from datetime import datetime
from matplotlib.legend_handler import HandlerTuple
#scale=10.58
lhapdf.setVerbosity(0)
import multiprocessing


# In[4]:


#sep=44.7
#IC=1
#nucleon='proton'# 'proton','neutron','helium','deuterium','lead'
#pdf_name='CT10'# 'NNPDF40'#, 'CT14IC', 'CT10'


# In[5]:


dfs1 = pd.read_csv('fit_parameters/dfs_/new_gauss_su_charm_no_no_def.csv') 
dfs2 = pd.read_csv('fit_parameters/dfs_/new_gauss_su_charm_no_yes_def.csv')
dfs3 = pd.read_csv('fit_parameters/dfs_/new_gauss_su_charm_yes_yes_def.csv')

df1 = pd.read_csv('fit_parameters/dfs_/fit_hadron_coef_0.27_chi_1.174__True_gk_PV17_su_no_charmno.csv')
df2 = pd.read_csv('fit_parameters/dfs_/fit_hadron_coef_0.27_chi_1.259__True_gk_PV17_su_no_charmyes_correction_no.csv')
df3 = pd.read_csv('fit_parameters/dfs_/fit_hadron_coef_0.27_chi_1.447__True_gk_PV17_su_yes_charmyes.csv')


# In[6]:


def grids_lp(df,su2,charm,sep,IC,pdf_name):

    yy=0.4
    dati_lp = pd.read_csv("fit_parameters/lprot_point.csv")
    if sep == 28.6:dati_lp=dati_lp.loc[(dati_lp['xb']>0.05)]

    fnc = polarization(0.27,sep)
    fnc.mass = 1.115 
    fnc.frag2 = 'dss'
    fnc.g_k = 'PV17'
    fnc.mdl_den = 'pwr_lw_star'
    fnc.mdl_num = 'gauss'
    fnc.charm =charm
    if charm=='yes':nf=4
    elif charm=='no':nf=3
    fnc.nf= nf
    fnc.IC_num=IC
    fnc.pdf_name=pdf_name

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
    num=zeros(len(dati_lp.hads1))
    i=0

    for hads1,hads2,zz,xbs in zip(dati_lp.hads1,dati_lp.hads2,dati_lp.z1,dati_lp.xb):

        if hads1 == 300 : had1='lbd'
        elif hads1 == 310 : had1='lbd_b'

        if hads2 == 100 : had2='proton'
        elif hads2 == 105 : had2='neutron'
        elif hads2 == 200 : had2='k+'
        elif hads2 == 205 : had2='k-'
        #print()
        num[i]= fnc.ratio_y(had1,had2,zz,xbs,yy,f_prm,pt_pp,0.)
        i+=1

    dati_lp['conv'] = num

    #dati_lp=dati_lp.drop(columns=['mins','maxx','fit_line'])
    
    return dati_lp


# In[7]:


'''
df1 = pd.read_csv('fit_parameters/dfs_/fit_hadron_coef_0.27_chi_1.174__True_gk_PV17_su_no_charmno.csv')
su2='no'
charm='no'
sep=40.7
IC=0
pdf_name='CT10'
a=grids_lp(df1,su2,charm,sep,IC,pdf_name)
a
'''


# In[8]:


def grids_ln(df,su2,charm,sep,IC,pdf_name):

    yy=0.4
    dati_lp = pd.read_csv("fit_parameters/lneutr_point.csv")
    if sep == 28.6:dati_lp=dati_lp.loc[(dati_lp['xb']>0.05)]

    fnc = polarization(0.27,sep)
    fnc.mass = 1.115 
    fnc.frag2 = 'dss'
    fnc.g_k = 'PV17'
    fnc.mdl_den = 'pwr_lw_star'
    fnc.mdl_num = 'gauss'
    fnc.charm =charm
    if charm=='yes':nf=4
    elif charm=='no':nf=3
    fnc.nf= nf
    fnc.IC_num=IC
    fnc.pdf_name=pdf_name

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
    num=zeros(len(dati_lp.hads1))
    i=0

    for hads1,hads2,zz,xbs in zip(dati_lp.hads1,dati_lp.hads2,dati_lp.z1,dati_lp.xb):

        if hads1 == 300 : had1='lbd'
        elif hads1 == 310 : had1='lbd_b'

        if hads2 == 100 : had2='proton'
        elif hads2 == 105 : had2='deuterium' # 'proton','neutron','helium','deuterium','lead'
        #print()
        num[i]= fnc.ratio_y(had1,had2,zz,xbs,yy,f_prm,pt_pp,0.)
        i+=1

    dati_lp['conv'] = num

    #dati_lp=dati_lp.drop(columns=['mins','maxx','fit_line'])
    
    return dati_lp


# In[9]:


'''
df1 = pd.read_csv('fit_parameters/dfs_/fit_hadron_coef_0.27_chi_1.174__True_gk_PV17_su_no_charmno.csv')
su2='no'
charm='no'
sep=40.7
IC=0
pdf_name='CT10'
a=grids_ln(df1,su2,charm,sep,IC,pdf_name)
a
'''


# In[10]:


'''
sep = 44.7 # 28.6, 44.7, 63.2, 104.9, 140.7
pdf_name='CT10'
IC=0

su2_1= 'no'
charm1='no'
su2_2= 'no'
charm2='yes'
su2_3= 'yes'
charm3='yes'
cut_dfs=20
dati_lp1=grids_lp(df1,su2_1,charm1,sep,IC,pdf_name)
dati_ln1=grids_ln(df1,su2_1,charm1,sep,IC,pdf_name)

dati_lp2=grids_lp(df2,su2_2,charm2,sep,IC,pdf_name)
dati_ln2=grids_ln(df2,su2_2,charm2,sep,IC,pdf_name)

dati_lp3=grids_lp(df3,su2_3,charm3,sep,IC,pdf_name)
dati_ln3=grids_ln(df3,su2_3,charm3,sep,IC,pdf_name)
'''


# In[15]:


cut_dfs=600


# In[16]:


def grids_lp_bands(df_prm,df,su2,charm,sep,IC,pdf_name):
    df_prm=df_prm.loc[(df_prm.index<cut_dfs)]
    
    yy=0.4
    dati_lp = pd.read_csv("fit_parameters/lprot_point.csv")
    if sep == 28.6:dati_lp=dati_lp.loc[(dati_lp['xb']>0.05)]

    fnc = polarization(0.27,sep)
    fnc.mass = 1.115 
    fnc.frag2 = 'dss'
    fnc.g_k = 'PV17'
    fnc.mdl_den = 'pwr_lw_star'
    fnc.mdl_num = 'gauss'
    fnc.charm =charm
    if charm=='yes':nf=4
    elif charm=='no':nf=3
    fnc.nf= nf
    fnc.IC_num=IC
    fnc.pdf_name=pdf_name

    dati_lp=grids_lp(df,su2,charm,sep,IC,pdf_name)
    num_tmp=zeros(len(df_prm.NUP))
    mins_lp = zeros(len(dati_lp.hads1))
    maxx_lp = zeros(len(dati_lp.hads1))
    i=0 

    for hads1,hads2,zz,xbs in zip(dati_lp.hads1,dati_lp.hads2,dati_lp.z1,dati_lp.xb):

        if hads1 == 300 : had1='lbd'
        elif hads1 == 310 : had1='lbd_b'

        if hads2 == 100 : had2='proton'
        elif hads2 == 105 : had2='neutron' # 'proton','neutron','helium','deuterium','lead'

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

                num_tmp[j] = fnc.ratio_y(had1,had2,zz,xbs,yy,f_prm,pt_pp,0.)
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

                num_tmp[j] = fnc.ratio_y(had1,had2,zz,xbs,yy,f_prm,pt_pp,0.)
                j+=1


        elif su2=='yes' and charm=='yes':

            j=0
            for nup, ndo, nst, nsea, ast, bup,bdo, pp in zip(df_prm.NUP,df_prm.NDO,df_prm.NST,df_prm.NSEA,df_prm.AST,df_prm.BUP,df_prm.BDO,df_prm.PP2):

                f_prm=arange(0.,18.,1.)
                f_prm[0] = nup
                f_prm[1] = nup#par[0]  #do
                f_prm[2] = nst  #par[0]      #st
                f_prm[3] = nsea  #par[0]      #upb
                f_prm[4] = nsea #par[0]      #dob
                f_prm[5] = ndo # par[0]      #stb
                #¯¯¯¯¯¯
                f_prm[6] = 0 #par[0]      #aup
                f_prm[7] = 0 #par[0]      #ado
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

                num_tmp[j] = fnc.ratio_y(had1,had2,zz,xbs,yy,f_prm,pt_pp,0.)
                j+=1

        mins_lp[i]=min(num_tmp)# fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,0)
        maxx_lp[i]=max(num_tmp)# fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,0)

        i+=1

    dati_lp['mins'] = mins_lp
    dati_lp['maxx'] = maxx_lp

    dati_lp.to_csv(r'fit_parameters/bands_/dati_lprot_su_charm_'+str(su2)+'_'+str(charm)+'_scale_'+str(sep)+'_PDF_'+str(pdf_name)+'_IC_'+str(IC)+'_def.csv',index=False)

    #return dati_lp


# In[ ]:





# In[20]:


def grids_ln_bands(df_prm,df,su2,charm,sep,IC,pdf_name):
    df_prm=df_prm.loc[(df_prm.index<cut_dfs)]
    
    yy=0.4
    dati_lp = pd.read_csv("fit_parameters/lneutr_point.csv")
    if sep == 28.6:dati_lp=dati_lp.loc[(dati_lp['xb']>0.05)]
    #sep=sep**2

    fnc = polarization(0.27,sep)
    fnc.mass = 1.115 
    fnc.frag2 = 'dss'
    fnc.g_k = 'PV17'
    fnc.mdl_den = 'pwr_lw_star'
    fnc.mdl_num = 'gauss'
    fnc.charm =charm
    if charm=='yes':nf=4
    elif charm=='no':nf=3
    fnc.nf= nf
    fnc.IC_num=IC
    fnc.pdf_name=pdf_name

    dati_lp=grids_ln(df,su2,charm,sep,IC,pdf_name)
    num_tmp=zeros(len(df_prm.NUP))
    mins_lp = zeros(len(dati_lp.hads1))
    maxx_lp = zeros(len(dati_lp.hads1))
    i=0 

    for hads1,hads2,zz,xbs in zip(dati_lp.hads1,dati_lp.hads2,dati_lp.z1,dati_lp.xb):

        if hads1 == 300 : had1='lbd'
        elif hads1 == 310 : had1='lbd_b'

        if hads2 == 100 : had2='proton'
        elif hads2 == 105 : had2='deuterium' # 'proton','neutron','helium','deuterium','lead'

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

                num_tmp[j] = fnc.ratio_y(had1,had2,zz,xbs,yy,f_prm,pt_pp,0.)
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

                num_tmp[j] = fnc.ratio_y(had1,had2,zz,xbs,yy,f_prm,pt_pp,0.)
                j+=1


        elif su2=='yes' and charm=='yes':

            j=0
            for nup, ndo, nst, nsea, ast, bup,bdo, pp in zip(df_prm.NUP,df_prm.NDO,df_prm.NST,df_prm.NSEA,df_prm.AST,df_prm.BUP,df_prm.BDO,df_prm.PP2):

                f_prm=arange(0.,18.,1.)
                f_prm[0] = nup
                f_prm[1] = nup#par[0]  #do
                f_prm[2] = nst  #par[0]      #st
                f_prm[3] = nsea  #par[0]      #upb
                f_prm[4] = nsea #par[0]      #dob
                f_prm[5] = ndo # par[0]      #stb
                #¯¯¯¯¯¯
                f_prm[6] = 0 #par[0]      #aup
                f_prm[7] = 0 #par[0]      #ado
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

                num_tmp[j] = fnc.ratio_y(had1,had2,zz,xbs,yy,f_prm,pt_pp,0.)
                j+=1

        mins_lp[i]=min(num_tmp)# fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,0)
        maxx_lp[i]=max(num_tmp)# fnc.ratio(had1,had2,zs1,zs2,f_prm,pt_pp,0)

        i+=1

    dati_lp['mins'] = mins_lp
    dati_lp['maxx'] = maxx_lp

    dati_lp.to_csv(r'fit_parameters/bands_/dati_deuterium_su_charm_'+str(su2)+'_'+str(charm)+'_scale_'+str(sep)+'_PDF_'+str(pdf_name)+'_IC_'+str(IC)+'_def.csv',index=False)
    
    #return dati_lp


# In[19]:

sep = 28.6 # 28.6, 44.7, 63.2, 104.9, 140.7

sep2 = 63.2 # 28.6, 44.7, 63.2, 104.9, 140.7
pdf_name1='CT14IC'# 'NNPDF40'#, 'CT14IC', 'CT10'
pdf_name2='NNPDF40'# 'NNPDF40'#, 'CT14IC', 'CT10'

IC=0 # 0 noIC, 2 BHPS



su2_1= 'no'
charm1='no'
#
su2_2= 'no'
charm2='yes'
#
su2_3= 'yes'
charm3='yes'

#grids_lp_bands(dfs1,df1,su2_1,charm1,sep,IC,pdf_name)

start = time.time()



p1 = multiprocessing.Process(target=grids_ln_bands,args=(dfs1,df1,su2_1,charm1,sep,IC,pdf_name1))
p2 = multiprocessing.Process(target=grids_ln_bands,args=(dfs2,df2,su2_2,charm2,sep,IC,pdf_name1))
p3 = multiprocessing.Process(target=grids_ln_bands,args=(dfs3,df3,su2_3,charm3,sep,IC,pdf_name1))

p4 = multiprocessing.Process(target=grids_ln_bands,args=(dfs1,df1,su2_1,charm1,sep2,IC,pdf_name1))
p5 = multiprocessing.Process(target=grids_ln_bands,args=(dfs2,df2,su2_2,charm2,sep2,IC,pdf_name1))
p6 = multiprocessing.Process(target=grids_ln_bands,args=(dfs3,df3,su2_3,charm3,sep2,IC,pdf_name1))

p7 = multiprocessing.Process(target=grids_ln_bands,args=(dfs2,df2,su2_2,charm2,sep,2,pdf_name1))
p8 = multiprocessing.Process(target=grids_ln_bands,args=(dfs2,df2,su2_2,charm2,sep,IC,pdf_name2))

p9 = multiprocessing.Process(target=grids_ln_bands,args=(dfs3,df3,su2_3,charm3,sep,2,pdf_name1))
p10 = multiprocessing.Process(target=grids_ln_bands,args=(dfs3,df3,su2_3,charm3,sep,IC,pdf_name2))

#p9 = multiprocessing.Process(target=grids_lp_bands,args=(dfs3,df3,su2_3,charm3,sep,IC,pdf_name2))



p1.start()
p2.start()
p3.start()
p4.start()
p5.start()
p6.start()
p7.start()
p8.start()
p9.start()
p10.start()


p1.join()
p2.join()
p3.join()
p4.join()
p5.join()
p6.join()
p7.join()
p8.join()
p9.join()
p10.join()


end = time.time()

mins=(end -start)/60
print('¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯¯')
print('time passed:')

print(str(end - start) +'   ' + 'sec')
print(str((end - start)/60) +'   ' + 'min')




