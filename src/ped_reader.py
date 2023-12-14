import pandas as pd
import numpy as np
import linecache

'''
# check .ped file columns
file_path = '../ukb_22828_c19_b0_v3.ped'
with open(file_path) as r:
    line = r.readline()
    print(len(line.split())) # 4174040
'''
'''
Family ID
Individual ID
Paternal ID
Maternal ID
Sex (1=male; 2=female; other=unknown)
Phenotype
SNPs
'''
'''
# check .map file lines
file_path2 = '../ukb_22828_c19_b0_v3.map'
with open(file_path2) as r2:
    line = r2.readline()
    print(line)
    # all lines count (wc -l filename): 2087017 (*2=4174034)
'''

lead_snp = 'rs2384686'
chr_n = '19'
file_path2 = '../ukb_22828_c' + chr_n + '_b0_v3.map'
with open(file_path2) as r2:
    line = r2.readline()
    count = 1
    #print(line)
    while line:
        if(line.split()[1] == lead_snp):
            print(count) # 1964131 => 
            final_count = count
        line = r2.readline()
        count += 1

case_file = '../datasets/final/cleaned_case.csv'
control_file = '../datasets/final/cleaned_control.csv'

df_case = pd.read_csv(case_file)
df_control = pd.read_csv(control_file)
sample_file = '../../../share/ukb22828_c' + chr_n + '_b0_v3.sample'

# generate index of case and control samples
df_case['sample_index'] = -1
df_control['sample_index'] = -1

df_sample = pd.read_csv(sample_file,sep=' ')
print(df_sample)

for i in range(df_case.shape[0]):
    try:
        df_case['sample_index'][i] = df_sample.loc[df_sample['ID_1'] == df_case['eid'][i]].index[0] - 1
    except IndexError:
        pass
print(df_case)
for i in range(df_case.shape[0]):
    try:
        df_control['sample_index'][i] = df_sample.loc[df_sample['ID_1'] == df_control['eid'][i]].index[0] - 1
    except IndexError:
        pass


df_case['A1'] = '-'
df_control['A1'] = '-'
df_case['A2'] = '-'
df_control['A2'] = '-'

# generate phenotype
file_path = '../ukb_22828_c' + chr_n + '_b0_v3.ped'
with open(file_path) as r:
    count = 0
    line = r.readline()
    while line:
        if(count in df_case['sample_index'].tolist()):
            #print(count)
            text = line
            A1 = text.split()[2*final_count - 2 + 6]
            A2 = text.split()[2*final_count - 2 + 6 + 1]
            df_case.loc[df_case['sample_index'] == count,'A1'] = A1
            df_case.loc[df_case['sample_index'] == count,'A2'] = A2
            #print(df_case)
            #df_case.to_csv('../datasets/final/ped_snp_case_' + lead_snp + '.csv',index=None)
        if(count in df_control['sample_index'].tolist()):
            text = line
            A1 = text.split()[2*final_count - 2 + 6]
            A2 = text.split()[2*final_count - 2 + 6 + 1]
            df_control.loc[df_control['sample_index'] == count,'A1'] = A1
            df_control.loc[df_control['sample_index'] == count,'A2'] = A2
        line = r.readline()
        count += 1

df_case.to_csv('../datasets/final/ped_snp_case_' + lead_snp + '.csv',index=None)
df_control.to_csv('../datasets/final/ped_snp_control_' + lead_snp + '.csv',index=None)
