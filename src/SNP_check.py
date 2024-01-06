# using bgen_reader

import os
from bgen_reader import read_bgen
import pandas as pd
import numpy as np
import argparse

case = '40'
SNP_file = '../datasets/6p/ukb_poi' + case +'_6p.regenie'
case_file = '../datasets/final/cleaned_case.csv'
control_file = '../datasets/final/cleaned_control.csv'
snp = pd.read_csv(SNP_file, sep='\t')
df_case = pd.read_csv(case_file)
df_control = pd.read_csv(control_file)
print(snp)
print(df_case)
print(df_control)

parser = argparse.ArgumentParser()
parser.add_argument('--index','-i', default=0)
args = parser.parse_args()

i = int(args.index)
chr_ = str(snp['CHROM'][i])
bgen = read_bgen(filepath='/dssg/home/acct-bmelgn/share/ukb22828_c' + chr_ + '_b0_v3.bgen',samples_filepath='/dssg/home/acct-bmelgn/share/ukb22828_c' + chr_ + '_b0_v3.sample',verbose=True)
if(chr_ == '19'):
    bgen = read_bgen(filepath='/dssg/home/acct-bmelgn/bmelgn-4/chr19_ukb/ukb22828_c19_b0_v3.bgen',samples_filepath='/dssg/home/acct-bmelgn/share/ukb22828_c' + chr_ + '_b0_v3.sample',verbose=True)
variants = bgen['variants']
samples = bgen['samples']
genotype = bgen["genotype"]
snp_id = snp['ID'][i]
print(snp_id)
variant = variants[variants["rsid"] == snp_id].compute()
variant_idx = variant.index.values[0]

# add new column: SNP annotation
df_case[snp_id + 'genotype1'] = 4
df_control[snp_id + 'genotype1'] = 4
df_case[snp_id + 'genotype2'] = 4
df_control[snp_id + 'genotype2'] = 4
df_case[snp_id + 'genotype3'] = 4
df_control[snp_id + 'genotype3'] = 4
df_case[snp_id] = 3
df_control[snp_id] = 3
case_sample_error = []
for j in range(df_case.shape[0]):
    sample_id = df_case['eid'][j]
    try:
        sample_idx = samples[samples == str(sample_id)].index.values[0]
        p = genotype[variant_idx].compute()["probs"][sample_idx]
        p_label = np.argmax(p)
        df_case.loc[df_case['eid'] == sample_id,snp_id] = p_label
        df_case.loc[df_case['eid'] == sample_id,snp_id + 'genotype1'] = p[0]
        df_case.loc[df_case['eid'] == sample_id,snp_id + 'genotype2'] = p[1]
        df_case.loc[df_case['eid'] == sample_id,snp_id + 'genotype3'] = p[2]
    except IndexError:
        #print(sample_id)
        case_sample_error.append(sample_id)
control_sample_error = []
for j in range(df_control.shape[0]):
    sample_id = df_control['eid'][j]
    try:
        sample_idx = samples[samples == str(sample_id)].index.values[0]
        p = genotype[variant_idx].compute()["probs"][sample_idx]
        p_label = np.argmax(p)
        df_control.loc[df_control['eid'] == sample_id,snp_id] = p_label
        df_control.loc[df_control['eid'] == sample_id,snp_id + 'genotype1'] = p[0]
        df_control.loc[df_control['eid'] == sample_id,snp_id + 'genotype2'] = p[1]
        df_control.loc[df_control['eid'] == sample_id,snp_id + 'genotype3'] = p[2]
    except IndexError:
        #print(sample_id)
        control_sample_error.append(sample_id)

# rearrange df_case and df_control
for idx in case_sample_error:
    df_case = df_case.drop(df_case[df_case['eid'] == idx].index)
df_case = df_case.reset_index(drop=True)
print(df_case)

for idx in control_sample_error:
    df_control = df_control.drop(df_control[df_control['eid'] == idx].index)
df_control = df_control.reset_index(drop=True)
print(df_control)

df_case.to_csv('../datasets/final_6p/snp_case_' + str(i+1) + '.csv',index=None)
df_control.to_csv('../datasets/final_6p/snp_control_' + str(i+1) + '.csv',index=None)