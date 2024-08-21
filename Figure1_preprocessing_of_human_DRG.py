#!/usr/bin/env python
# coding: utf-8

# # import function and sample information

# In[1]:


import scanpy as sc
import anndata as ad
import scvelo as scv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import sparse
from anndata import AnnData
import os

#pd.set_option('display.max_rows', 100000)
pd.set_option('display.max_columns', 100000)
#pd.set_option('display.width', 100000)


# In[2]:


rootdir = '/WorkSpace/DRG'
datadir = f'{rootdir}/01Data'
workdir = f'{rootdir}/02Results'
os.chdir(rootdir)


# In[3]:


import sys
import importlib
sys.path.append(rootdir)
import SCFunc
importlib.reload(SCFunc)
from SCFunc import *


# In[4]:


sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=100, dpi_save=300, fontsize=10, facecolor='white') 
sc.settings.figdir='./'


# In[5]:


import collections
MA={'neuron_crest_cell':'SOX10,SOX9,PAX3',
    'neuron_progenitor':'NEUROG1,NEUROG2,SIX1,EYA1,EYA2,IRX1',
    'primary_sensory_neuron':'ISL1,RBFOX3,POU4F1,UCHL1,TUBB3,ELAVL3,MAP2,POU4F2,POU4F3',
    'schwanncellprogenitor':'SOX2,FOXD3,FABP7,FABP4',
    'satellite_glial_cells':'APOE,KCNJ10,DBI,PLP1,GLUL,S100B,VIM,SLC1A3,GJA1',
    'Schwann_cell':'MPZ,PLLP,PRX,EDNRB,MBP,PMP22,MAG,GJB1',
    'remak_schwann_cell':'NCAM1',
    'GLIA_lineage':'ERBB2,ERBB3,NRG1,EGR2,SLC22A16,POU3F2,NFKB1,HDAC1,HDAC2,JAG1,NFATC4',
    'vascular_endothelial_cell':'CLDN5,?LY6C1,PLVAP,PECAM1,EMCN,ECSCR,CDH5,IGFBP7',
    'vascular_smooth_muscle_cells':'TPM2,ACTA2,DES',
    'vascular_endothelial_capillary':'S100A8',
    'fibroblast':'DCN,APOD,MGP',
    'immune_cells':'CD74',
    'red_blood_cell':'HBA1',
    'monocytic_cell':'LYZL2,LYZL2?LYZ2',
    'pericyte':'KCNJ8',
    'macrophage':'AIF1,CD68,MRC1,CD40,ICAM1,CIITA',
    'T-cells':'CD3G',
    'mesenchymal':'CDH2,FN1,ENG,ENG,NT5E,?CD73,THY1,?CD90,MME,?CD10,ANPEP,?CD13,NGFR,?CD271,SPARC,ATL1,ACTA2,COL1A1,COL3A1',
    'epithelial':'CDH1,TJP1,?ZO1,CLDN1,OCLN,COL4A1',
    'boundary_cap_cell':'EGR2,PRSS56,NTN5,SEMA6A,CDH7,SEMA3B,SEMA3G,CXCR4,HEY2,WIF1,HEYL,NPR3',
    'skin':'KRT8'}

MB={'neuron_crest_cell':'SOX10,SOX9,PAX3',
    'neuron_progenitor':'NEUROG1,NEUROG2,SIX1,EYA1,EYA2,IRX1',
    'primary_sensory_neuron':'NEUROD1,ISL1,RBFOX3,POU4F1,UCHL1,TUBB3,ELAVL3,MAP2,POU4F2,POU4F3',
    'boundary_cap_cell':'EGR2,PRSS56,NTN5,SEMA6A,CDH7,SEMA3B,SEMA3G,CXCR4,HEY2,WIF1,HEYL,NPR3',
    'peptidergic_neuron':'TAC1,?CGRP,CALCA,FAM19A1,TAFA1,FAM19A1/TAFA1,NEFH,NTRK1,MET',
    'non-peptidergic_neuron':'CD55,PLXNC1,?MRGPRD,P2RX3,TRPV1,SCN10A',
    'ALTMR_NF':'NEFH,LDHB,?NEACB2,CALB1,PVALB,CNTNAP2,SENP6,NTRK2,NTRK3,SHOX2',
    'proprioceptor':'NTRK3,?PARVB,RUNX3',
    'mechanoreceptor':'NTRK2,RET,MAFA',
    'cLTMR':'PIEZO2,TH,SLC17A8,TAFA4,GFRA2',
    'thermoceptor':'TRPV1,TRPM8,TRPA1',
    'pruriceptors':'PIEZO2,TH,P2RX3,SST',
    'nociceptor':'SCN10A,SCN9A,NTRK1,RUNX1,TRPV1,TRPA1'}

MB1={'neuron_subtypes':'ISL1,RBFOX3,POU4F1,UCHL1,TUBB3,STMN2,POU4F2,SIX1,MAP2,CACNA1H',
    'peptidergic_neuron':'TAC1,CALCA,?CGRP,CALCA,FAM19A1,TAFA1,FAM19A1/TAFA1,NEFH,NTRK1,MET',
    'non-peptidergic_neuron':'CD55,PLXNC1,MRGPRD,P2RX3,TRPV1,SCN10A',
    'ALTMR_NF':'NEFH,LDHB,?NEACB2,CALB1,PVALB,CNTNAP2,SENP6,NTRK2,NTRK3,SHOX2',
    'proprioceptor':'NTRK3,PARVB,RUNX3',
    'mechanoreceptor':'NTRK2,RET,MAFA',
    'cLTMR':'PIEZO2,TH,SLC17A8,TAFA4,GFRA2',
    'thermoceptor':'TRPV1,TRPM8,TRPA1',
    'pruriceptors':'PIEZO2,TH,P2RX3,SST',
    'nociceptor':'SCN10A,SCN9A,NTRK1,RUNX1,TRPV1,TRPA1'}


MC ={  'neuron_crest_cell': 'SOX10',
        'sensory_neuron_progenitor': 'NEUROG1,NEUROG2,SIX1,EYA1,EYA2,IRX1,NTRK1,NTRK2,NTRK3',
        'primary_sensory_neuron': 'ISL1,RBFOX3,POU4F1,UCHL1,TUBB3,ELAVL3,',
        'satellite_glial_cells': 'APOE,FABP7,KCNJ10,KIR4,DBI,PLP1,ERBB3',
        'Schwann_cell': 'MPZ,PLLP,PRX,SOX2,EDNRB',
        'GLIA': 'MBP',
        'vascular_endothelial_cell': 'CLDN5,LY6C1,PLVAP,PECAM1,CD31,EMCN,ECSCR,CDH5,LGFBP7',
        'vascular_smooth_muscle_cells': 'TPM2,ACTA2,DES',
        'vascular_endothelial_cells_capillary': 'S100A8',
        'glia_progenitor': 'NIFA',
        'fibroblast': 'DCN,APOD,MGP',
        'macrophage': 'ALF1,IBA1,CD68,MRC1',
        'immune_cells': 'CD74',
        'red_blood_cell': 'HBA1',
        'monocytic_cell': 'LYZ2',
        'pericyte': 'KCNJ8',
        'mesenchymal': 'CD34',
        'T-cells': 'CD3G',
        'connective_tissue': 'COL1A1'}

MD = {  'neuron_crest_cell': 'SOX10',
        'neuron_progenitor': 'NEUROG1,NEUROG2,SIX1,EYA1,EYA2,IRX1',
        'unspecific_sensory_neuron': 'ELAVL3',
        'peptidergic_neuron': 'TAC1,CGRP,CALCA,CGRP/CALCA,FAM19A1,PHETA1,FAM19A1/PHETA1,NEFH,NTRK1,PLXNC1,MET,SCN10A,TRPV1,RUNX1,TNFRSF1A,POU4F1',
        'non-peptidergic_neuron': 'CD55,PLXNC1,MRGPRD,P2RX3,CALCA,NTRK1',
        'ALTMR(NF)': 'NEFH,LDHB,NEACB2,CALB1,PV,PVALB,PV/PVALB,FAM19A1,PHETA1,FAM19A1/PHETA1,CNTNAP2,SSP1,NTRK2,NTRK3,SHOX2,SPP1',
        'NF': 'LDHB,NEACB2,CALB1,PVALB,PHETA1,CNTNAP2,SSP1,NTRK2,NTRK3',
        'TH': 'TH',
        'proprioceptor': 'NTRK3,PARVB,RUNX3,PVALB',
        'mechanoreceptor': 'NTRK2,RET,MAFA',
        'cLTMR(TH)': 'PIEZO2,TH,VGLUT3,FAM19A4,TAFA4,GFRA2,SLC17A8,FAM19A4,ZNF521,TLX3',
        'thermoceptor': 'TRPV1,TRPM8,TRPA1',
        'pruriceptors': 'PIEZO2,TH,P2X3,SST,NPPB',
        'nociceptor': 'SCN10A,SCN9A,NTRK1,RUNX1,TRPV1,TRPA1,MRGPRD'}

ME = {  'Neurons': 'SNAP25,MAP2,SYT1,SYP,RBFOX3',
        'Excitatory_Neuron': 'SLC17A7,CAMK2A,NRGN',
        'Inhibitory_Neuron/Interneuron': 'GAD1,GAD2,SLC32A1,SLC6A1',
        'EN_NRGN':  'NRGN,THY1',
        'Astrocyte': 'AQP4,GFAP,FGFR3,SLC1A2,SLC1A3,ALDH1A1,ALDH1L1',
        'Oligodendrocyte': 'MOG,MOBP,MBP,PLP1,OPALIN,HAPLN2',
        'Oligodendrocyte_Precursor_cell': 'PDGFRA,CSPG4,BCAN,PCDH15,APOD',
        'newly_born_Oligodendrocyte': 'GPR17,NKX2-2,BMP4,VCAN',
        'Microglia': 'PTPRC,P2RY12,C3,HLA-DRA,CD74,CD68,CX3CR1,C1QB',
        'Endothelial_cell': 'PECAM1,CLDN5,VWF,NOSTRIN,FLT1',
        'Pericyte': 'DCN,COL1A2,PDGFRB',
        'T_cell': 'IL7R,CD3E',
        'Gender': 'XIST,RPS4Y1,DDX3Y'}
            
markdict = collections.OrderedDict({"neuron crest cell": ["SOX10", "PAX3"],
                                    "sensory neuron progenitor": ["NEUROG1", "NEUROG2", "SIX1", "EYA1", "EYA2", "ISL1", "RBFOX3", "POU4F1", "ELAVL3"],
                                    "precursor1": ["NTRK3", "RUNX3"],
                                    "precursor2": ["NTRK1", "RUNX1"],
                                    "proprioceptor": ["NTRK3", "RUNX3"],
                                    "mechanoreceptor": ["NTRK2", "RET"],
                                    "nociceptor": ["NTRK1", "RUNX1", "TRPA1", "SCN10A", "TRPM8"],
                                    "TH": ["TH"],
                                    "schwann cell progenitor": ["SOX2", "FOXD3"],
                                    "schwann cell": ["MBP", "MAG", "MPZ", "PLLP"],
                                    "satellite glia cell": ["FABP7", "APOE", "GJA1", "SLC1A3", "S100B"],
                                    "endothelial cell": ["PECAM1", "EMCN"],
                                    "vascular smooth muscle cell": ["TPM2", "ACTA2"],
                                    "muscle cell": ["PAX7", "MYH3"],
                                    "macrophage": ["MRC1"],
                                    "red blood cell": ["HBA1"]})

MG = ['PVALB', 'KCNS1','ASIC1','RUNX3','NTRK3','NTRK2','TRPM8','SCN10A',
          'CPNE6', 'PENK','TRPA1', 'CHRNA3', 'SST','IL31RA','NPPB','GFRA2']

MABCDE = GetMk(MA, MB, MC, MD, ME)
MAB= GetMk(MA,MB)
MA = GetMk(MA)
MB = GetMk(MB)
ME = GetMk(ME)


MS = {
   'nociceptor': ["SCN10A", "SCN9A", "NTRK1", "RUNX1", "RET", "MET", "TAC1", "CALCA", 
                   "TAFA1", "NEFH", "NTRK1", "PLXNC1", "CD55", "PLXNC1", "MRGPRD", "P2RX3", 
                   "CALCA", "NTRK1", "PIEZO2", "TH", "VGLUT3", "TAFA4", "SLC17A8", "TRPV1", 
                   "TRPV2", "TRPV3", "TRPV4", "KCNK2", "KCNK4", "TRPM8", "TRPA1", "KCNK4", "KCNK2",
                   "KCNK4", "TRPV4", "TRPA1", "TRPV1", "KCNC4", "KCND3", "TRPA1", "MRGA1", "MRGA3", 
                   "MRGA4", "TRPV1", "TRPV2", "TRPM8", "TRPC3", "SCN11A", "MRGPRD", "P2X3", "CALCA",
                   "OPRM1", "ASIC3"],
    "proprioceptor": ["WHRN", "PVALB", "RUNX3", "NTRK3", "ETV1", "LMCD1", "CHAD", "FXYD7", "PIEZO1",
                      "PIEZO2", "CDH13", "SEMA5A", "CRTAC1", "HEG1", "NXPH1", "PCDH8", "WNT7A", "HEG1",
                      "POU4F3", "CALB2", "NXPH1"],
    "mechanoreceptor": ["MAFA", "RET", "GFRA2", "MAFA", "RET", "GFRA2", "NTRK2", "MAFA", "RET", "GFRA2",
                        "NTRK3", "MAFB", "C-MAF", "ASIC1A", "ASIC1B", "ASIC2A", "ASIC2B", "ASIC3", "ASIC4",
                        "PIEZO1", "PIEZO2", "KCNK2", "KCNK4", "KCNK10", "CACNA1H", "KCNQ4", "SCN11A", "TRPV1",
                        "TRPV2", "TRPV4", "TRPC1", "TRPC3", "TRPC5", "TRPC6", "TRPP1", "TRPP2", "TRPP3",
                        "TRPA1", "TRPM3", "TRPM4", "TRPM7"],
    'neuronglia': MA.V,
    'sensory_neuron_progenitor': ["NEUROG1", "NEUROG2", "SIX1", "EYA1", "EYA2", "IRX1"]
}

Recepters = {
    'nociceptor':{
        'putative C-LTMR': ['SLC17A8','GFRA2'], #'DFNA25', 'VGLUT3', 'GDNFRB', 'NRTNR-ALPHA', 'NTNRA', 'RETL2', 'TRNR2'],
        'pruriotgen-receptor' : ['SST','RET'],
        'Cold-nociceptor' : ['TRPMB'],
        'peptidergic' : ['TAFA1','TAC1'],
        'non-peptidergic' : ['CD55','RUNX1'],
    },

    'mechanoreceptor' :{
        'Immature LTMRs': ['GFRA2','NTRK2','SHOX2','NTRK3'],
        'ADelta-LTMR'   : ['NTRK2','MAFA','SHOX2','NTRK3','RUNX3','RUNX1','POU4F2'],
        'ABeta-SA-LTMR' : ['RET','GFRA2','NTRK2','MAFA'],
        'ABeta-RA-LTMR' : ['RET','GFRA2','NTRK2','MAFA'],
    },

    'proprioceptor' : {
        'Immature PNs': ['NTRK3','RUNX3'],
        'IA-PN' : ['VSTM2B','HTRA1'],
        'IB-PN' : ['SYNPR','DOC2B'],
        'II-PN' : ['GABRA2','ZEB2','TAC1','EPHA5'],
    },
}
MR = sum([sum(v.values(),[]) for v in Recepters.values()], [])
   
MS['receptors'] = MS['proprioceptor'] + MS['proprioceptor'] + MS['proprioceptor']


# In[6]:


sampleinfor = pd.DataFrame({ 'sampleid': ["DRG1", "DRG3", "DRG4", "DRG5", "DRG8", "DRG9", "DRG10", "S1", "S2",
                                           "S3", "S4", "S5", "S6", "S7",  "S8", 'HE10_DRG', 'DRG6', 'DRG7'],
                             "species" : ['human']*16 + ['macaca']*2,
                             'AGE'   :   ['GW8+3', 'GW10', 'GW7+2', 'GW9+5', 'GW14', 'GW17_1', 'GW21', 'GW7+4','GW8', 'GW9',
                                          'GW12_1', 'GW12_2', 'GW15_1', 'GW15_2', 'GW17_2', 'GW10_2', '13day', '6year'],
                             'GW'  :    ['GW8', 'GW10', 'GW7', 'GW10', 'GW14', 'GW17', 'GW21', 'GW7','GW8', 'GW9',
                                          'GW12', 'GW12', 'GW15', 'GW15', 'GW15', 'GW10', '13day', '6year'],
                             'SID'  :    ['GW8_2', 'GW10_2', 'GW7_1', 'GW10_1', 'GW14', 'GW17', 'GW21', 
                                          'GW7_2', 'GW8_1', 'GW9', 'GW12_1', 'GW12_2', 'GW15_1', 'GW15_2',
                                          'GW15_3', 'GW10_0', 'mm_13day', 'mm_6year'],
                             'GW_Rep':   ['8.0-8.3', '9.5-10', '7.2-7.4', '9.5-10', '14', '17', '21', 
                                          '7.2-7.4', '8.0-8.3', '9', '12', '12', '15', '15',
                                          '15', '10.2', 'mm13d', 'mm6y'],
                             'Time' :    [8.3, 10, 7.2, 9.5, 14, 17, 21, 7.4,8, 9, 12, 12, 15, 15, 15, 10.2, 2.1, 312.6]})

sampleinfor = sampleinfor.sort_values(by=['species', 'Time', 'SID']).reset_index(drop=True).astype('category')
for i in sampleinfor.columns:
    sampleinfor[i] = pd.Categorical(sampleinfor[i], categories = sampleinfor[i].drop_duplicates().tolist())
sampleinfor


# # merge adata

# In[ ]:


# AGGR
# %%bash
# OU=~/00DataBase/02AGGR
# cd $OU
# ~/cellranger-6.0.2/bin/cellranger aggr \
#     --id=V6_02_aggr \
#     --csv=~/00DataBase/02AGGR/V6.0.2_aggr.info.csv  \
#     --normalize=mapped \
#     --maxjobs 60 \
#     --localcores 60


# In[7]:


def humanAGGR(sinfo, IN, label='sampleid', mapid = {'1': 'AHC13', '2': 'AHC24', '3': 'AHC26', '4': 'AHC27'}):
    AGGR = IOs().getMX(IN)
    AGGRidx = AGGR.obs.index.to_series().str.split('-')
    AGGR.obs.index = AGGRidx.str[1].replace(mapid).str.cat(AGGRidx.str[0], sep=':' )
    AGGR.var.genome = 'GRCh38'

    #adata = ad.concat([AGGR, HE10_DRG])
    adata = AGGR
    adata.var.index.name = 'Gene'
    adata.var['mt'] = adata.var_names.str.startswith('MT-') 
    adata.var['ribo'] = adata.var_names.str.contains('^RP[SL]', regex=True)
    adata.var['hb'] = adata.var_names.str.contains('^HB[^(P)]', regex=True)
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)

    adata.obs[label] = adata.obs.index.str.split(':').str[0]
    adata.obs.index.name = 'CellID'
    adata.obs['barcordID'] = adata.obs.index
    adata.obs = adata.obs.reset_index().merge(sinfo, on=label, how="left").set_index('CellID')
    #adata.obs[label] = pd.Categorical(adata.obs[label], categories = sinfo[label].cat.categories)
    return(adata)

sinfo = sampleinfor[((sampleinfor.species =='human') & (sampleinfor.sampleid !='HE10_DRG'))].copy()
sinfo.SID = sinfo.SID.cat.remove_unused_categories()

groupby  = 'SID'
mapid = ["DRG1", "DRG3", "DRG4", "DRG5", "DRG8", "DRG9", "DRG10", "S1", "S2", "S3", "S4", "S5", "S6", "S7", "S8"]
mapid = { str(i+1): mapid[i] for i in range(len(mapid)) }
AGGRh5 = f'{datadir}/filtered_feature_bc_matrix/'
adatahs = humanAGGR(sinfo, AGGRh5, mapid=mapid)
# adatahs = sc.read_h5ad(f'{datadir}/hs.adata.raw.h5ad')
adatahs


# # set scruble parameters

# ##  parallele in parametering of scruble

# In[ ]:


import multiprocessing
def getscore(adatahs, expected_doublet_rate, groupby='SID',  ):
    iadata = Generalists(adatahs).GetScore(groupby, method=['scruble'], expected_doublet_rate=expected_doublet_rate)  
    return iadata.obs.groupby(groupby)['doubletP'].sum()

with multiprocessing.Pool(5) as pool:
    results = [pool.apply_async(getscore, args=(adatahs, _c, )).get()
               for _c in np.arange(0.06, 0.15, 0.01)]


# In[ ]:


from joblib import Parallel, delayed
def getscore(adatahs, groupby='SID', expected_doublet_rate=0.06):
    iadata = Generalists(adatahs).GetScore(groupby, method=['scruble'], 
                                           expected_doublet_rate=expected_doublet_rate)  
    dP = iadata.obs['doubletP']
    dP.name = '%.3f'%expected_doublet_rate
    return dP
results = Parallel( n_jobs=11, verbose=20 )( delayed(getscore)(adatahs.copy(), expected_doublet_rate=_c) 
                            for _c in np.arange(0.06, 0.161, 0.005))
results
results1 = pd.concat( [adatahs.obs]+results, axis=1)
results1.to_csv('scrublet.doubletP.Yen.detail.csv')


# In[ ]:


results1 = pd.concat( [adatahs.obs]+results, axis=1)
results1 = results1[ [groupby] + list(map(lambda x:'%.3f'%x, np.arange(0.06, 0.161, 0.005)))]
results2 = pd.concat([results1.groupby(groupby).size(), results1.groupby(groupby).sum(1)], axis=1)
results2.to_csv('scrublet.doubletP.Yen.counts.csv')
results2


# In[ ]:


import matplotlib.pyplot as plt
ax=results2.iloc[:,1:].T.plot(kind='line')
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))
leg = plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5), numpoints=1,  title='SID')
plt.xlabel( 'expected_doublet_rate')
plt.xticks( range(0, results2.shape[1]-1), results2.iloc[:,1:].T.index, rotation=90)
plt.ylabel( 'doublet_counts')
plt.savefig('scrublet.doubletP.Yen.counts.pdf',bbox_inches='tight')
plt.show()
plt.close()


# In[ ]:


adatahs1 = Generalists(adatahs).GetScore(groupby, expected_doublet_rate=0.1, threshold_method='Mean')


# In[ ]:


for i in np.arange(0.1,0.61,0.025):
    Generalists(adatahs).GetScore(groupby, expected_doublet_rate=0.1, 
                                  method=['scruble'], 
                                  scthre={'hs_GW21' : i}, specified=['hs_GW21'])


# # Human

# In[8]:


###########QC
def Filter1(adata, groupby='SID'):
    adata = adata.copy()
    sc.pp.filter_cells(adata, min_counts=800)
    sc.pp.filter_cells(adata, max_counts=40000)
    sc.pp.filter_cells(adata, min_genes=800)
    sc.pp.filter_cells(adata, max_genes=8000)
    adata = adata[((adata.obs.pct_counts_mt <= 5) & (adata.obs.doubletP ==False)) ,]
    sc.pp.filter_genes(adata, min_cells=11)
    # Plot(adata).QCplot(groupby=groupby, header='after', outdir='./scafterQC')
    Plot(adata).QCplot(groupby=groupby)
    return adata

###########QC
def Filter2(adata, groupby='SID'):
    adata = adata.copy()
    sc.pp.filter_cells(adata, min_counts=600)
    #sc.pp.filter_cells(adata, max_counts=40000)
    sc.pp.filter_cells(adata, min_genes=600)
    #sc.pp.filter_cells(adata, max_genes=8000)
    adata = adata[((adata.obs.pct_counts_mt <= 5) & (adata.obs.doubletP ==False)) ,]
    sc.pp.filter_genes(adata, min_cells=11)
    Plot(adata).QCplot(groupby=groupby, header='after', outdir='./scafterQC')
    return adata

def Filter9(adata, groupby='SID'):
    adata = adata.copy()
    sc.pp.filter_cells(adata, min_counts=1000)
    sc.pp.filter_cells(adata, max_counts=40000)
    sc.pp.filter_cells(adata, min_genes=800)
    sc.pp.filter_cells(adata, max_genes=8000)
    adata = adata[((adata.obs.pct_counts_mt <= 5) & (adata.obs.doubletP ==False)) ,]
    sc.pp.filter_genes(adata, min_cells=11)
    Plot(adata).QCplot(groupby=groupby, header='after', outdir='./scafterQC')
    return adata

def CellCounts(adata, adataraw, groupby='SID', out='QC.cells.counts.mean.gene.xls'):
    QC = pd.concat([ adataraw.obs.groupby(groupby).size(), adataraw.obs.groupby(groupby)['n_genes_by_counts'].mean(),
                     adataraw.obs.groupby('SID')['doubletP'].sum(),
                     adata.obs.groupby(groupby).size(), adata.obs.groupby(groupby)['n_genes_by_counts'].mean(),], axis=1)
    QC.columns = ['raw_cells', 'raw_meangene', 'doubletP', 'keep_cells', 'keep_meangene']
    QC['filter_cells']=QC.raw_cells-QC.keep_cells
    QC.to_csv(out, sep='\t')
    print(QC)
    return(QC)


# # add infor

# In[9]:


groupby='SID'
'''
adatahs = Generalists(adatahs).GetScore(groupby, expected_doublet_rate=0.14,
                                        threshold_method='Minimum',
                                        scthre={'GW7_1':0.3796, 'GW7_2':0.3474, 'GW8_1':0.3668, 'GW8_2':0.3278,
                                                'GW9':0.4, 'GW10_1':0.2727, 'GW10_2':0.3555, 'GW12_1':0.6096, 
                                                'GW12_2': 0.3763, 'GW14':0.3545, 'GW15_1':0.3097, 'GW15_2':0.4177,
                                                'GW15_3': 0.2713, 'GW17':0.5, 'GW21':0.2555})
'''
adatahs = Generalists(adatahs).GetScore(groupby, expected_doublet_rate=0.15, #test 0.06 new 0.15 old 0.14
                                        threshold_method='Minimum',
                                        ccgenecsv=f'{datadir}/human_mouse_cell_cycle_genes.csv',)

# adatahs.write('hs.adata.scanpy.raw.h5ad')

adata = adatahs.copy()
groupby = 'SID'
Plot(adata).QCplot(groupby=groupby)
# Plot(adata).QCplot(groupby=groupby, header='rawQC', outdir='./rawQC')


# In[10]:


# Unsupervise(adatahs).SplitCheck(split_key=groupby, outdir='./rawCheck')
Unsupervise(adatahs).SplitCheck(split_key=groupby)


# # raw data overview

# - **RawCheck**

# In[ ]:


groupby='SID'
adataR = Preprocession(adatahs.copy()).NormScale(batch_key=groupby, isnormal=True, isscale=True,
                                               dropMT=True, dropRibo=False,
                                               minnbat=4, min_mean=0.04, min_disp=0.4, max_mean=6)

adataR = Unsupervise(adataR).RawCheck(batch_key=groupby)


# In[9]:


# adataI = sc.read(f'{OUTR}/Scanoramav1/60/annotation/adata.Scanorama.60.h5ad')
# adataR.obs.loc[ (adataR.obs_names.intersection(adataI.obs_names)), 'CellType'] = \
#     adataI.obs.loc[ (adataR.obs_names.intersection(adataI.obs_names)), 'CellType']

# my_cmap=PlotUti().continuous_cmap(["lightgrey", 'yellow', 'red', 'darkred'])
# sc.pl.umap( adataR, color='CellType', color_map=my_cmap, na_color='black')
# sc.pl.umap( adataR, color='SID', color_map=my_cmap)

# sc.pl.umap( adataR, color='CellType', color_map=my_cmap, na_color='black')
# sc.pl.umap(adataR, color=['PIEZO1', 'PIEZO2'], color_map=my_cmap)


# In[ ]:


Plot(adatahs).QCplot(groupby='SID', header='before')


# In[ ]:


# os.chdir(f'{workdir}/RawView')
Plot(adataR).QCmap( clims= {'n_genes_by_counts':8000, 'total_counts':10000, 'pct_counts_mt':0.4, 'pct_counts_ribo':30, 'pct_counts_hb':0.3},
                    # save='.%s.%s.%s.mindist%s.QCFeature.pdf'%('quickcheck', 50, 2, 0.3)
                  )


# In[ ]:


# adataR.write('./RawView/adata.hs.RawView.h5ad')


# # clean Data

# In[10]:


import sys
import importlib
sys.path.append('/share/home/zhonw/JupyterCode')
import SCFunc
importlib.reload(SCFunc)
from SCFunc import *

adataR = Filter1(adatahs.copy())
adataR = adataR[(adataR.obs['SID'] !='GW15_3'),:].copy()

for i in sinfo.columns:
    adataR.obs[i] = adataR.obs[i].cat.remove_unused_categories(inplace=True)

# IOs().saveadata(adataR, '14samples')
CellCounts(adataR, adatahs.copy())
# Plot(adata).QCplot(groupby=groupby, header='14sampleQC', outdir='./14sampleQC')
Plot(adata).QCplot(groupby=groupby)


# In[11]:


print(adata)
Plot(adatahs).QCplot(groupby='SID', header='before', outdir='./scbeforeQC')


# In[ ]:


adataR=sc.read('adata.14samples.raw.h5ad')
adataR


# In[13]:


Plot(adataR).QCplot(groupby='SID', header='after', outdir='./scafterQC')


# # get adata

# In[7]:


OUTR = '/gpfs2/wulab17/WorkSpace/11Project/02DRG/01Analysis/Human/20220421/Scanoramav1/60'
os.makedirs(OUTR, exist_ok=True)
os.chdir(OUTR)

adataraw = sc.read(f'{OUTR}/annotation/adata.Scanorama.60.h5ad')
#adataraw = sc.read('/gpfs2/wulab17/WorkSpace/11Project/02DRG/01Analysis/Human/20220421/Scanoramav1/60/SubTypes/progenitor2receptors/Scanorama/60/adata.Scanorama.60.h5ad')
adataraw = adataraw.raw.to_adata().copy()

hsadata = '/gpfs2/wulab17/WorkSpace/11Project/02DRG/01Analysis/Human/20220421/hs.adata.raw.h5ad'
hsadata = sc.read_h5ad(hsadata)

adataraw.raw = hsadata[adataraw.obs_names, adataraw.var_names]
adataraw.X = adataraw.raw.X.copy()

adataraw.obs.loc[adataraw.obs['Clusters'].isin(['56']), 'CellType'] = 'satellite glia cell'
#adataraw.obs.loc[adataraw.obs['leiden_2.00'].isin(['12']), 'CellType'] = 'precursor1'
sc.pl.umap(adataraw, color=['CellType','Clusters'] )


# In[6]:


OUTR = '/gpfs2/wulab17/WorkSpace/11Project/02DRG/01Analysis/Human/20220421/RawView'
os.makedirs(OUTR, exist_ok=True)
os.chdir(OUTR)
#adatahs.write('hs.adata.raw.h5ad')
adataR = sc.read_h5ad('adata.hs.RawView.h5ad')
print(adataR.raw.X.max, adataR.X.max)
adataR


# In[9]:


adataR = adataR.raw.to_adata().copy()
del adataR.uns['log1p']


# In[ ]:


BATCH='Scanorama'
Res=5
groupby='SID'
CellType='CellType'
MARKERS = MA.V

louvain=f'louvain_{Res :.2f}'
leiden=f'leiden_{Res :.2f}'
COLS = [groupby, louvain, leiden]
    
min_dist=0.3

OUT = f'{OUTR}/{BATCH}'

adata = Preprocession(adataR.copy()).NormScale(batch_key=groupby, isnormal=False, isscale=False,
                                               dropMT=True, dropRibo=False,
                                               minnbat=4, min_mean=0.05, min_disp=0.5, max_mean=6)

for NPCS in range(50, 61,10):
    iOUT = '%s/%s'%(OUT,NPCS)
    os.makedirs(iOUT, exist_ok=True)
    os.chdir(iOUT)
    
    adataI = Integretion(adata).SCANO(NPCS=NPCS, batch_key=groupby, batch_size=8000, 
                                      n_neighbors=40, knn=40, alpha=0.1, sigma=100)
                
    adataI = Unsupervise(adataI).MAP(batch_key=groupby, NPCS=NPCS, method=BATCH,Res=Res, min_dist=min_dist)
    
    sc.pl.umap(adataI, color=['GW', 'leiden', 'CellType'],
               save=f'.{BATCH}.{NPCS}.GW.leiden.umap.png')

    my_cmap=PlotUti().continuous_cmap(["lightgrey", 'yellow', 'red','darkred'])
    sc.pl.umap(adataI, color=[i for i in MARKERS if i in adataI.raw.var.index], ncols=9, color_map=my_cmap, 
               show=False, save=f'.{BATCH}.{NPCS}.gene.classic.umap.min_dist{min_dist}.png')
    
    adataI = Unsupervise(adataI).Clust(clust=['leiden', 'louvain'])
    adataI = Unsupervise(adataI).Add23d(method='umap', npcs=3, min_dist=min_dist)

    AA = Ply(adataI).SData(obms='X_umap_3d', method='umap', groups=[groupby, leiden, louvain])
    Ply(adataI).Scatter3ds(AA, header = 'umap.%s.%s.%s'%(BATCH, NPCS, Res))

    os.chdir(iOUT)
    IOs().saveadata(adataI, f'{BATCH}.{NPCS}', todense=True)


# - **RawCheck**

# In[ ]:


groupby='SID'

adata = Preprocession(adataR.copy()).NormScale(batch_key=groupby, isnormal=True, isscale=True, dropMT=True, 
                                               minnbat=3, min_mean=0.035, min_disp=0.4, max_mean=6)
Unsupervise(adata).RawCheck(batch_key=groupby, outdir='./Nointegrate')


# - **Scanorama**

# In[ ]:


BATCH='Scanorama'
Res=5
groupby='SID'
CellType='CellType'
MARKERS = MA.V()

louvain=f'louvain_{Res :.2f}'
leiden=f'leiden_{Res :.2f}'
COLS = [groupby, louvain, leiden]
    
min_dist=0.3

OUT = f'{OUTR}/{BATCH}'

adata = Preprocession(adataR.copy()).NormScale(batch_key=groupby, isnormal=True, isscale=False,
                                               dropMT=True, dropRibo=False,
                                               minnbat=4, min_mean=0.04, min_disp=0.4, max_mean=6)

for NPCS in range(50, 61,10):
    iOUT = '%s/%s'%(OUT,NPCS)
    os.makedirs(iOUT, exist_ok=True)
    os.chdir(iOUT)
    
    adataI = Integretion(adata).SCANO(NPCS=NPCS, batch_key=groupby, batch_size=8000, 
                                      n_neighbors=40, knn=40, alpha=0.1, sigma=100)
                
    adataI = Unsupervise(adataI).MAP(batch_key=groupby, NPCS=NPCS, method=BATCH,Res=Res, min_dist=min_dist)
    
    sc.pl.umap(adataI, color=['GW','leiden'], save=f'.{BATCH}.{NPCS}.GW.leiden.umap.png')

    my_cmap=PlotUti().continuous_cmap(["lightgrey", 'yellow', 'red','darkred'])
    sc.pl.umap(adataI, color=[i for i in MARKERS if i in adataI.raw.var.index], ncols=9, color_map=my_cmap, 
               show=False, save=f'.{BATCH}.{NPCS}.gene.classic.umap.min_dist{min_dist}.png')
    
    adataI = Unsupervise(adataI).Clust(clust=['leiden', 'louvain'])
    adataI = Unsupervise(adataI).Add23d(method='umap', npcs=3, min_dist=min_dist)

    AA = Ply(adataI).SData(obms='X_umap_3d', method='umap', groups=[groupby, leiden, louvain])
    Ply(adataI).Scatter3ds(AA, header = 'umap.%s.%s.%s'%(BATCH, NPCS, Res))

    muscle=['DES','TTN','MYL2','MYH3','PAX7','MYOG','PITX1','PITX3','PAX3','PITX2','MYF5','MYF6','MYOD',
            'CXCR4','SDC3','ITGA7','FGFR4','ACTA2','RGS5','TAGLN','MCAM','ANPEP','CNN1','S100A4','TNC']
    mtgene = adataI.raw.var_names[adataI.raw.var_names.str.startswith('MT-')]
    ribogene = adataI.raw.var_names[adataI.raw.var_names.str.contains('^RP[SL]', regex=True)]
    QCgene= list(mtgene) + list(ribogene)

    geneList = MA.V() + MB.V() + QCgene + muscle + sum(markdict.values(),[])

    os.makedirs(f'{iOUT}/gene_3d', exist_ok=True)
    os.chdir(f'{iOUT}/gene_3d')
    px_cmap=[(0,"lightgrey"), (0.33,'yellow'), (0.67,'red'), (1,'darkred')]

    AA = Ply(adataI).SData(obms='X_umap_3d', method='umap', groups=geneList)

    Ply(adataI).Scatter3ds(AA, colormap=px_cmap, show=False, header = 'umap.%s.%s.%s'%(BATCH, NPCS, Res))
    os.makedirs(f'{iOUT}/cluters', exist_ok=True)
    os.chdir(f'{iOUT}/cluters')
    for i in range(2,9):
        for c in ['leiden', 'louvain']:
            iclust = f'{c}_{i}.00'
            sc.pl.umap(adataI, color=iclust, legend_fontsize =10, show=False, save=f'.{iclust}.2d.pdf')
            AA = Ply(adataI).SData(obms='X_umap', method='umap', groups=[iclust])
            Ply(adataI).Scatter2ds(AA, show=False, header= 'umap.%s.%s.%s'%(BATCH, NPCS, iclust))
            AA = Ply(adataI).SData(obms='X_umap_3d', method='umap', groups=[iclust])
            Ply(adataI).Scatter3ds(AA, show=False, header= 'umap.%s.%s.%s'%(BATCH, NPCS, iclust))
    
    os.chdir(iOUT)
    IOs().saveadata(adataI, f'{BATCH}.{NPCS}', todense=True)


#


from scipy.sparse import csr_matrix
import dbmap as dm
data = csr_matrix(adataI.X) 
diff = dm.diffusion.Diffusor(ann_dist='euclidean', #'cosine', 'jaccard', 'hamming' 
                             knn_dist='euclidean',
                             kernel_use= 'decay_adaptive', #'decay_adaptive',
                             n_jobs=15,
                             n_neighbors=40, 
                             n_components=120,
                             transitions=False, norm=False).fit(data)

min_dist = 0.3
db = diff.transform(data)
db = np.array(db)
res = diff.return_dict()
adataI.obsm['X_db'] = db

import matplotlib.pyplot as plt
plt.plot(range(0, len(res['EigenValues'])), res['EigenValues'], marker='o')


# In[ ]:


# I Diffusion graph layout with UMAP 
import umap
db_umap_emb = umap.UMAP(n_components=2, min_dist = 0.5, n_neighbors=30).fit_transform(db)
adataI.obsm['X_dbmap'] = db_umap_emb
sc.pl.embedding(adataI, basis ='X_dbmap', color='CellType') 

# II Diffusion graph layout with PaCMAP on the diffusion basis
db_pac_emb = dm.pacmapper.PaCMAP(n_dims=2, n_neighbors=50, MN_ratio=3, FP_ratio=.5) 
db_pac_emb_fit = db_pac_emb.fit_transform(db, init='random')
adataI.obsm['X_db_pacmap'] = db_pac_emb_fit
sc.pl.embedding(adataI, basis ='X_db_pacmap', color='CellType')

# III
sc.pp.neighbors(adataI, n_neighbors=15, use_rep='X_db', metric='euclidean')
sc.tl.umap(adataI, min_dist=0.3)
sc.pl.umap(adataI, color='CellType')

# IV
import pymde #pytorch 3.10 error
#https://topxometry.readthedocs.io/en/latest/welcome.html
mde = pymde.preserve_neighbors(db, embedding_dim=2, verbose=True, constraint=pymde.Standardized())
adataI.obsm['X_db_pymde'] = mde.embed(verbose=True).numpy()
#mde.embed(snapshot_every=1)
#mde.play(savepath='test.pymde.gif')
sc.pl.embedding(adataI, basis ='X_db_pymde', color='CellType') 


# In[ ]:


mde = pymde.preserve_neighbors(db, embedding_dim=3, verbose=True, constraint=pymde.Standardized())
adataI.obsm['X_db_pymde_3d'] = mde.embed(verbose=True).numpy()
AA = Ply(adataI).SData(obms='X_db_pymde_3d', method='umap', groups=['CellType'])
Ply(adataI).Scatter3ds(AA, show=True)


# In[ ]:


db_umap_emb_3d = umap.UMAP(n_components=3, min_dist = min_dist).fit_transform(db)
adataI.obsm['X_dbmap_3d'] = db_umap_emb_3d
AA = Ply(adataI).SData(obms='X_dbmap_3d', method='umap', groups=['CellType'])
Ply(adataI).Scatter3ds(AA, show=True)


# In[ ]:


db_pac_emb_3d = dm.pacmapper.PaCMAP(n_dims=3, n_neighbors=50, MN_ratio=5, FP_ratio=.5) 
adataI.obsm['X_db_pacmap_3d'] = db_pac_emb_3d.fit_transform(db, init='random')
AA = Ply(adataI).SData(obms='X_db_pacmap_3d', method='umap', groups=['CellType'])
Ply(adataI).Scatter3ds(AA, show=True)


# In[ ]:




