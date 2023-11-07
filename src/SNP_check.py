import os
from bgen_reader import read_bgen
import pandas as pd
import numpy as np

case = '40'
SNP_file = '../datasets/ukb_poi' + case +'_8p.regenie'
case_file = '../datasets/final/cleaned_case.csv'
control_file = '../datasets/final/cleaned_control.csv'
snp = pd.read_csv(SNP_file, sep='\t')
df_case = pd.read_csv(case_file)
df_control = pd.read_csv(control_file)
print(snp)
print(df_case)
print(df_control)

chr_ = ''
for i in range(snp.shape[0]):
    chr = str(snp['CHROM'][i])
    if(chr != chr_):
        chr_ = chr
        bgen = read_bgen(filepath='/dssg/home/acct-bmelgn/share/ukb22828_c' + chr_ + '_b0_v3.bgen',samples_filepath='/dssg/home/acct-bmelgn/share/ukb22828_c' + chr_ + '_b0_v3.sample',verbose=True)
        variants = bgen['variants']
        samples = bgen['samples']
        genotype = bgen["genotype"]
    snp_id = snp['ID'][i]
    variant = variants[variants["rsid"] == snp_id].compute()
    variant_idx = variant.index.values[0]

    # add new column: SNP annotation
    df_case[snp_id] = 0
    df_control[snp_id] = 0
    case_sample_error = []
    for j in range(df_case.shape[0]):
        sample_id = df_case['eid'][j]
        try:
            sample_idx = samples[samples == str(sample_id)].index.values[0]
            p = genotype[variant_idx].compute()["probs"][sample_idx]
            p_label = np.argmax(p)
            df_case.loc[df_case['eid'] == sample_id,snp_id] = p_label
        except IndexError:
            print(sample_id)
            case_sample_error.append(sample_id)
    control_sample_error = []
    for j in range(df_control.shape[0]):
        sample_id = df_control['eid'][j]
        try:
            sample_idx = samples[samples == str(sample_id)].index.values[0]
            p = genotype[variant_idx].compute()["probs"][sample_idx]
            p_label = np.argmax(p)
            df_control.loc[df_control['eid'] == sample_id,snp_id] = p_label
        except IndexError:
            print(sample_id)
            control_sample_error.append(sample_id)

    # rearrange df_case and df_control
    for id in case_sample_error:
        df_case = df_case.drop(df_case[df_case['eid'] == id].index)
    df_case = df_case.reset_index(drop=True)
    print(df_case)

    for id in control_sample_error:
        df_control = df_control.drop(df_control[df_control['eid'] == id].index)
    df_control = df_control.reset_index(drop=True)
    print(df_control)

df_case.to_csv('../datasets/final/snp_case.csv',index=None)
df_control.to_csv('../datasets/final/snp_control.csv',index=None)