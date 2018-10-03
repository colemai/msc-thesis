
# coding: utf-8

# # Create NT files from CTD csvs
# <b>Author</b>: Ian Coleman <br/>
# <b>Function</b>: Takes local CTD csvs and turns each into a .nt file

# In[1]:


import pandas as pd
import numpy as np
import scipy as sp
import subprocess
import math


# ## Functions

# In[2]:


def convert_df_nt (df, output_file, subj_url, subj_col, obj_url, obj_col, pred_col, odd_url=0):
    """
    Input:
        DF: some rows and columns of a dataframe
        STR: name for the output file, include filetype .nt
        STR: subj_url is the url to be used for all subjects
        STR: subj_col is the column from which to get the id for the subj
        STR: obj_url is the url to be used for all objects
        STR: obj_col is the column from which to get the id for the obj
        STR: OPTIONAL odd_url is for when a subset obj/subj of require a different url 
    Output:
        NT file
    """
    f = open(output_file,'w')
    for index, row in df.iterrows():
        subj = '<' + subj_url + str(row[subj_col]) + '> '
        pred = '<' + 'http://ian.ie/' + str(row[pred_col]) + '> '
        if 'OMIM' in str(row[obj_col]):
            row[obj_col] = str(row[obj_col]).replace('OMIM:', '')
            obj = '<' + odd_url + str(row[obj_col]) + '> '
        else:
            obj = '<' + obj_url + str(row[obj_col]) + '> '
        f.write(subj + pred + obj + '.' + '\n')
    f.close()


# ## Download Databases

# In[3]:


# subprocess.call('pip3 install wget', shell=True)
# subprocess.call('wget http://ctdbase.org/reports/CTD_chem_gene_ixns.csv.gz', shell=True)
# subprocess.call('wget http://ctdbase.org/reports/CTD_chemicals_diseases.csv.gz', shell = True)
# subprocess.call('wget http://ctdbase.org/reports/CTD_chem_pathways_enriched.csv.gz', shell = True)
# subprocess.call('wget http://ctdbase.org/reports/CTD_genes_diseases.csv.gz', shell = True)
# subprocess.call('wget http://ctdbase.org/reports/CTD_genes_pathways.csv.gz', shell = True)
# subprocess.call('wget http://ctdbase.org/reports/CTD_diseases_pathways.csv.gz', shell = True)
# subprocess.call('wget http://ctdbase.org/reports/CTD_pheno_term_ixns.csv.gz', shell = True)


# In[5]:


# Move all the csvs to a subfolder and unzip them
# subprocess.call('mkdir csvs', shell=True)
# subprocess.call('mv *.gz csvs/', shell=True)
# subprocess.call('gunzip csvs/*.gz', shell=True)


# In[6]:


# too ambitious?? 
# TODO attempt to make one func to import and process all ctd databases
# def ctd_to_rdf(csv, output_file, subj_url, subj_col, obj_url, obj_col, pred_col):
#     """
#     """
#     df = pd.read_csv(csv, skiprows=27)
#     df = df.drop(0)
#     convert_df_nt(df, output_file, subj_url, subj_col, obj_url, obj_col, pred_col)


# ## CHEM-GENE 

# In[38]:


# Read in CTD sample, skipping the intro rows8888
df_cg = pd.read_csv('csvs/CTD_chem_gene_ixns.csv', skiprows=27)
df_cg = df_cg.drop(0)
df_cg = df_cg.rename(columns={'# ChemicalName': 'ChemicalName'}) # rename of a column
print('Progress: 1')

# In[ ]:


# Split the interactionActions into separate predicates RUN THIS ONLY ONCE
s = df_cg['InteractionActions'].str.split('|').apply(pd.Series, 1).stack()
s.index = s.index.droplevel(-1)
s.name = 'InteractionActions'
df_cg = df_cg.join(s.apply(lambda x: pd.Series(x.split('|'))))
print('Progress: 2')

# In[ ]:


# Make the new column prettier
df_cg = df_cg.rename(columns={0: 'Predicate'})
df_cg['Predicate'] = df_cg.Predicate.str.replace('^', '_')
df_cg['Predicate'] = df_cg.Predicate.str.replace(' ', '_')
print('Progress: 3')

# In[ ]:


# Need to change float to int for the url to work
df_cg['GeneID'] = df_cg.GeneID.astype(int)


# In[ ]:


subj_url = 'http://identifiers.org/ctd.chemical/' 
subj_col = 'ChemicalID'
obj_url = 'http://identifiers.org/ctd.gene/' 
obj_col = 'GeneID'
pred_col = 'Predicate'

convert_df_nt(df_cg, '/scratch/dragon/amd/colemai/output_cg.nt', subj_url, subj_col, obj_url, obj_col, pred_col)

print('Progress: 4')
# ## Chem-Disease

# In[14]:


# Read in CTD sample, skipping the intro rows
df_cd = pd.read_csv('csvs/CTD_chemicals_diseases.csv', skiprows=27)
df_cd = df_cd.drop(0)


# In[15]:


'OMIM' in (df_cd.loc[65,'DiseaseID'])


# In[16]:


# Process DiseaseID so as to be usable in url
df_cd['DiseaseID'] = df_cd['DiseaseID'].str.replace('MESH:', '')


# In[17]:

print('Progress: 5')
# Create Predicate Column
def cd_predicate(r):
    """
    Create predicate from directevidence if available
    """
    if r['DirectEvidence'] == "nan":
        return 'associated_by_inference_via_' + str(r.InferenceGeneSymbol)
    else:
        return 'associated_directly_with'

print('starting pred')    
df_cd['DirectEvidence'] = df_cd.DirectEvidence.astype(str) 
df_cd['Predicate'] = df_cd.apply(cd_predicate, axis=1)
print('ending pred')

# In[18]:


subj_url = 'http://identifiers.org/ctd.chemical/' 
subj_col = 'ChemicalID'
obj_url = 'http://identifiers.org/mesh/'
obj_url_2 = 'http://identifiers.org/omim/' # to use when CTD gives an omim disease id
obj_col = 'DiseaseID'
pred_col = 'Predicate'

convert_df_nt(df_cd, '/scratch/dragon/amd/colemai/output_cd.nt', subj_url, subj_col, obj_url, obj_col, pred_col, obj_url_2)


# ## Gene Disease

# In[19]:

print('Progress: 6')
# Read in CTD sample, skipping the intro rows
df_gd = pd.read_csv('csvs/CTD_genes_diseases.csv', skiprows=27)
df_gd = df_gd.drop(0)


# In[21]:


# Must make some quick refinements to ensure resulting URLs work
df_gd['GeneID'] = df_gd['GeneID'].fillna(0).astype(int)
df_gd['GeneID'] = df_gd.GeneID.astype(int) 
df_gd['DirectEvidence'] = df_gd.DirectEvidence.astype(str) 
df_gd['DiseaseID'] = df_gd['DiseaseID'].str.replace('MESH:', '')

# Create Predicate Column
def gd_predicate(r):
    """
    Create predicate
    """
    if r['DirectEvidence'] == "nan":
        return 'associated_by_inference_via_' + str(r.InferenceChemicalName).replace(' ', '_')
    else:
        return 'associated_directly_with'
    
df_gd['Predicate'] = df_gd.apply(gd_predicate, axis=1)


# In[24]:
print('Progress: 7')

subj_url = 'http://identifiers.org/ctd.gene/' 
subj_col = 'GeneID'
obj_url = 'http://identifiers.org/mesh/' 
obj_col = 'DiseaseID'
pred_col = 'Predicate'

convert_df_nt(df_gd, '/scratch/dragon/amd/colemai/output_gd.nt', subj_url, subj_col, obj_url, obj_col, pred_col)


# ## Gene Pathway

# In[25]:


# Read in CTD sample, skipping the intro rows
df_gp = pd.read_csv('csvs/CTD_genes_pathways.csv', skiprows=27)
df_gp = df_gp.drop(0)

print('Progress: 8')
# In[27]:


# Must make some quick refinements to ensure resulting URLs work
df_gp['GeneID'] = df_gp['GeneID'].fillna(0).astype(int)
df_gp['GeneID'] = df_gp.GeneID.astype(int) 
df_gp['PathwayID'] = df_gp['PathwayID'].str.replace('REACT:', '')

# Create Predicate Column
def gp_predicate(r):
    return 'associated_directly_with'
    
df_gp['Predicate'] = df_gp.apply(gp_predicate, axis=1)


# In[28]:


subj_url = 'http://identifiers.org/ctd.gene/' 
subj_col = 'GeneID'
obj_url = 'http://identifiers.org/reactome/' 
obj_col = 'PathwayID'
pred_col = 'Predicate'

convert_df_nt(df_gp, '/scratch/dragon/amd/colemai/output_gp.nt', subj_url, subj_col, obj_url, obj_col, pred_col)


# ## Disease Pathway

# In[29]:


# Read in CTD sample, skipping the intro rows
df_dp = pd.read_csv('csvs/CTD_diseases_pathways.csv', skiprows=27)
df_dp = df_dp.drop(0)


# In[31]:


# Must make some quick refinements to ensure resulting URLs work
df_dp['PathwayID'] = df_dp['PathwayID'].str.replace('REACT:', '')
df_dp['DiseaseID'] = df_dp['DiseaseID'].str.replace('MESH:', '')


# Create Predicate Column
def dp_predicate(r):
    return 'associated_by_inference_via_' + str(r.InferenceGeneSymbol)
    
df_dp['Predicate'] = df_dp.apply(dp_predicate, axis=1)


# In[32]:


subj_url = 'http://identifiers.org/mesh/' 
subj_col = 'DiseaseID'
obj_url = 'http://identifiers.org/reactome/' 
obj_col = 'PathwayID'
pred_col = 'Predicate'

convert_df_nt(df_dp, '/scratch/dragon/amd/colemai/output_dp.nt', subj_url, subj_col, obj_url, obj_col, pred_col)


# ## Phenotype Chemical
# I'm going to comment this section out as this is the Y data and shouldn't be in KG ( I think )

# In[33]:


# Read in CTD sample, skipping the intro rows
# df_pc = pd.read_csv('csvs/CTD_pheno_term_ixns.csv', skiprows=27, nrows = 100)
# df_pc = df_pc.drop(0)
# df_pc[10:20]


# In[34]:


# Split the interactionActions into separate predicates RUN THIS ONLY ONCE
# s = df_pc['interactionactions'].str.split('|').apply(pd.Series, 1).stack()
# s.index = s.index.droplevel(-1)
# s.name = 'interactionactions'
# df_pc = df_pc.join(s.apply(lambda x: pd.Series(x.split('|'))))
# df_pc = df_pc.rename(columns={0: 'Predicate'})
# df_pc['Predicate'] = df_pc.Predicate.str.replace('^', '_')
# df_pc['Predicate'] = df_pc.Predicate.str.replace(' ', '_')


# In[35]:


# subj_url = 'http://identifiers.org/ctd.chemical/'  
# subj_col = 'chemicalid'
# obj_url = 'http://identifiers.org/go/' 
# obj_col = 'phenotypeid'
# pred_col = 'Predicate'

# convert_df_nt(df_pc, '/scratch/dragon/amd/colemai/output_pc.nt', subj_url, subj_col, obj_url, obj_col, pred_col)


# ## Merge NT files

# In[36]:


subprocess.call('cat *.nt > /scratch/dragon/amd/colemai/master.nt', shell=True)

