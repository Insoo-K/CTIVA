import numpy as np
import pandas as pd
import os
from statsmodels.formula.api import ols
from statsmodels.stats.anova import anova_lm


raw_data = []
with open('GSE9782-GPL96_series_matrix.txt','r') as f:
    while True :
        line = f.readline()
        raw_data.append(line)
        if not line : break

header_array = []
data_array = []
for sample in raw_data:
    if sample.strip().startswith('!'):
        p_sample = sample.replace('\n','')
        p_sample = p_sample.replace('!','')
        p_samples = p_sample.split('\t')
        header_array.append(p_samples[0])
        data_array.append(p_samples[1:])

sample_head_array = []
sample_array = []

for i in range(len(header_array)):
    head = header_array[i]
    check = head.split('_')[0]
    if check == 'Sample':
        sample_head_array.append(head)
        sample_array.append(data_array[i])

sample_array = np.array(sample_array)
sample_preprocess_array = sample_array.copy()
for i in range(len(sample_head_array)):
    head = sample_head_array[i]
    if head == 'Sample_characteristics_ch1' :
        for j in range(len(sample_array[i])):
            if sample_array[i,j].replace('"', '') != '':
                sample_head_array[i] = 'Sample_' + sample_array[i,j].replace('"','')
                break

sample_head_dict = {}
for i in range(len(sample_head_array)):
    sample_head_dict[sample_head_array[i]] = i

for i in range(len(sample_array[16])):
    preprocess_data = sample_preprocess_array[16,i].split('=')[1]
    preprocess_data = preprocess_data.replace('"', '')
    preprocess_data = preprocess_data.replace(' ', '')
    sample_preprocess_array[16,i] = preprocess_data

for i in range(9,31):
    sample_head = 'Sample_'+ sample_array[i,100].replace('"','')
    sample_head_array[i] = sample_head.split(' = ')[0]

sample_array = np.transpose(sample_array)
sample_preprocess_array = sample_array.copy()

sample_dataframe = pd.DataFrame(sample_array, columns=sample_head_array)
sample_dataframe.to_csv('SE9782-GPL96_preprocessed.csv',sep=',',index=False)


for i in range(len(sample_array)):
    characterisitcs = sample_array[i]

test_array = []
for sample in list(np.unique(sample_array[:,18])):
    head = sample_head_array[18].replace('Sample_','')
    if head in sample :
        test_array.append(True)
    else : test_array.append(False)

f_surv= open('GSE-surv.txt','w',encoding='utf-8')
event1_censor_check_true = sample_head_array[18].replace('Sample_', '')
event1_censor_check_false = sample_head_array[19].replace('Sample_', '')
event2_censor_check_true = sample_head_array[22].replace('Sample_', '')
event2_censor_check_false = sample_head_array[23].replace('Sample_', '')
for i in range(len(sample_array)):
    if event1_censor_check_true in sample_array[i,18] :
        event1_censor = '1'
        event1_time = sample_array[i,19].split(' = ')[1].replace('"', '')
        for j in range(20,24):
            if event2_censor_check_true in sample_array[i,j] :
                event2_censor = '1'
                event2_time = sample_array[i,j+1].split(' = ')[1].replace('"', '')
                break
            elif event2_censor_check_false in sample_array[i,j] :
                event2_censor = '0'
                event2_time = sample_array[i,j].split(' = ')[1].replace('"', '')
                break
    elif event1_censor_check_false in sample_array[i,18] :
        event1_censor = '0'
        event1_time = sample_array[i,18].split(' = ')[1].replace('"', '')
        for j in range(19,24):
            if event2_censor_check_true in sample_array[i,j] :
                event2_censor = '1'
                event2_time = sample_array[i,j+1].split(' = ')[1].replace('"', '')
                break
            elif event2_censor_check_false in sample_array[i,j] :
                event2_censor = '0'
                event2_time = sample_array[i,j].split(' = ')[1].replace('"', '')
                break
    elif event2_censor_check_true in sample_array[i,18]:
        event1_censor = '0'
        event1_time = '0'
        event2_censor = '1'
        event2_time = sample_array[i, 19].split(' = ')[1].replace('"', '')
    elif event2_censor_check_false in sample_array[i,18]:
        event1_censor = '0'
        event1_time = '0'
        event2_censor = '0'
        event2_time = sample_array[i,18].split(' = ')[1].replace('"', '')
    else :
        event1_censor = '0'
        event1_time = '0'
        for j in range(19,24):
            if event2_censor_check_true in sample_array[i,j] :
                event2_censor = '1'
                event2_time = sample_array[i,j+1].split(' = ')[1].replace('"', '')
                break
            elif event2_censor_check_false in sample_array[i,j] :
                event2_censor = '0'
                event2_time = sample_array[i,j].split(' = ')[1].replace('"', '')
                break
    event1_time = str(np.round(float(event1_time) / 365 , 4))
    event2_time = str(np.round(float(event2_time) / 365 , 4))
    line = event1_time + '\t' + event1_censor + '\t' + event2_time + '\t' + event2_censor + '\n'
    f_surv.write(line)
f_surv.close()

pGX_dict = {}
sex_dict = {}
age_dict = {}
pGX_raw_array = np.unique(sample_array[:,16])
sex_raw_array = np.unique(sample_array[:,11])

for i,pGX_raw in enumerate(pGX_raw_array):
    pGX_response = pGX_raw.split(' = ')[1].replace('"', '')
    pGX_dict[pGX_response] = str(i+1)

for i,sex_raw in enumerate(sex_raw_array):
    sex = sex_raw.split(' = ')[1].replace('"', '')
    sex_dict[sex] = str(i+1)

age_dict[64] = str(1)

f_covariate = open('GSE-covariate.txt','w',encoding='utf-8')
for i in range(len(sample_array)):
    pGX_raw = sample_array[i,16].split(' = ')[1].replace('"', '')
    age_raw = sample_array[i,13].split(' = ')[1].replace('"', '')
    sex_raw = sample_array[i,11].split(' = ')[1].replace('"', '')

    pGX = pGX_dict[pGX_raw]
    sex = sex_dict[sex_raw]
    if int(age_raw) > 64 : age = '2'
    else : age = '1'
    for j in range(18,24):
        if 'Number_of_Prior_Lines' in sample_array[i,j]:
            PL = sample_array[i,j].split(' = ')[1].replace('"','')
    line = '"%s"\t"%s"\t"%s"\t%s"\n' % (sex, age, pGX,PL)
    f_covariate.write(line)
f_covariate.close()

os.system('./CTIVA/GAIT/GAIT GSE-surv.txt GSE-covariate.txt')

GSE_pred = pd.read_csv('pred_T.txt',sep='\t',header=None).values
GSE_interval_time = np.abs(GSE_pred[:,1]-GSE_pred[:,0])

covariate_array = []
for i in range(len(sample_array)):
    pGX_raw = sample_array[i,16].split(' = ')[1].replace('"', '')
    age_raw = sample_array[i,13].split(' = ')[1].replace('"', '')
    sex_raw = sample_array[i,11].split(' = ')[1].replace('"', '')
    pGX = pGX_dict[pGX_raw]
    sex = sex_dict[sex_raw]
    if int(age_raw) > 64 : age = '2'
    else : age = '1'
    for j in range(18,24):
        if 'Number_of_Prior_Lines' in sample_array[i,j]:
            PL = sample_array[i,j].split(' = ')[1].replace('"','')
    covariate_array.append([sex,age,pGX,PL])

covariate_array = np.array(covariate_array)
data = np.concatenate((covariate_array, GSE_interval_time.reshape(-1,1)),axis=1)
columns = ['sex', 'age', 'pGX_Response','Number_of_Prior_Lines','Time']
data = pd.DataFrame(data,columns=columns)
data = data.astype(float)

model = ols('Time ~ C(sex)',data = data).fit()
print(anova_lm(model))
model = ols('Time ~ C(age)', data = data).fit()
print(anova_lm(model))
model = ols('Time ~ C(pGX_Response)', data = data).fit()
print(anova_lm(model))
model = ols('Time ~ C(Number_of_Prior_Lines)', data = data).fit()
print(anova_lm(model))