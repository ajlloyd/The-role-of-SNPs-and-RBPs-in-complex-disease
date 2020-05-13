import numpy as np
import os
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

#load in data chunks:
mylist = []
for idx, chunk in enumerate(pd.read_csv("./all_diseases.tsv", error_bad_lines=False, sep=r'\s*,\s*', delimiter="\t", low_memory=False, chunksize=100000)):
    print(f"Chunk {idx} loaded")
    mylist.append(chunk)
df = pd.concat(mylist, axis= 0)

# SNP in/out counts from the whole database:
full_inout_count = df["Inside"].value_counts().values
full_inout_count = np.flip(full_inout_count,0)
print("FULL COUNTS = ", full_inout_count)

# SNP in out counts from the CVD only database:
#is_cvd = df["DISEASE/TRAIT"].str.contains("Coronary artery disease|Cardiovascular disease|Ischemic stroke|Hypertension|Stroke|Coronary heart disease|myocardial infarction|Coronary artery disease|Sudden cardiac arrest")
is_cvd = df["DISEASE/TRAIT"].str.contains("Coronary artery disease|Cardiovascular disease|Ischemic stroke|Hypertension|Stroke|Resting heart rate|Coronary heart disease|Coronary artery disease (myocardial infarction, percutaneous transluminal coronary angioplasty, coronary artery bypass grafting, angina or chromic ischemic heart disease)|Myocardial infarction|Heart rate in heart failure with reduced ejection fraction|Cardiovascular disease risk factors|Heart rate|Sudden cardiac arrest|Aortic root size|Heart rate variability traits (pvRSA/HF)|Lipoprotein phospholipase A2 activity in cardiovascular disease|Large artery stroke|Coronary artery calcification|Coronary artery disease or large artery stroke|Idiopathic dilated cardiomyopathy|Cardiac troponin-I levels|Echocardiographic traits|Myocardial infarction (early onset)|Cardiac Troponin-T levels|Ischemic stroke (large artery atherosclerosis)|Cardiac hypertrophy|Coronary artery disease or ischemic stroke|Cardiovascular risk factors|White matter hyperintensities in ischemic stroke|Ischemic stroke (cardioembolic)|Coronary heart disease (SNP X SNP interaction)|Electrocardiographic traits|Ischemic stroke (small-vessel)|Cardiometabolic and hematological traits|Left ventricular internal dimension in diastole|Hypertension (SNP x SNP interaction)|Congenital left-sided heart lesions|Conotruncal heart defects (inherited effects)|Left ventricular internal dimension in systole|Nonobstructive coronary artery disease|Coronary artery disease and LDL cholesterol levels (multivariate analysis)|Left ventricle wall thickness|Mitral valve prolapse|Electrocardiographic conduction measures|Congenital heart disease (maternal effect)|Conotruncal heart defects (maternal effects)|Coronary artery disease and total cholesterol levels (multivariate analysis)|Carotid intima media thickness in rheumatoid arthritis|Small vessel stroke|Conotruncal heart defects|Lateral ventricular volume in normal aging|Mortality in heart failure|Coronary atherosclerosis (increased number of diseased vessels) (traffic exposure interaction)|Stroke (ischemic)|Coronary artery disease and triglyceride levels (multivariate analysis)|Left ventricle diastolic internal dimension|Heart failure|Coronary artery aneurysm in Kawasaki disease|Left ventricular obstructive tract defect (inherited effect)|Chronic obstructive pulmonary disease or coronary artery disease (pleiotropy)|Coronary artery disease and HDL cholesterol levels (multivariate analysis)|Aortic valve stenosis|Ischemic stroke (undetermined subtype)|Cardiovascular risk factors (age interaction)|Ventricular ectopy|Total ventricular volume|Left  ventricle systolic dysfunction|Cardiac structure and function|Postoperative myocardial infarction after cardiac surgery|Left ventricular QRS voltage|Thrombin-antithrombin complex levels in ischemic stroke|Intracranial, abdominal aortic or thoracic aortic aneurysm (pleiotropy)|Carotid plaques in rheumatoid arthritis|Congenital heart malformation|Ischemic heart disease in rheumatoid arthritis|Resistant hypertension|Pulmonary arterial hypertension|Congenital heart disease (inherited effect)|Vascular brain injury|Left ventricular mass|Postoperative stroke after cardiac surgery|Creatinine levels in ischemic stroke|Peak velocity of the mitral A-wave|Atrioventricular conduction|Left ventricular obstructive tract defect (maternal effect)|Thoracic aortic aneurysms and dissections|Supraventricular ectopy|Left ventricular fractional shortening|Congenital heart disease|Ischemic stroke (small artery occlusion)|Incident myocardial infarction|Cardiovascular event in rheumatoid arthritis|Cardiac repolarization|Stroke (pediatric)")
df_cvd = df[is_cvd]

#cvd_inout_count = df_cvd["Inside"].value_counts().values
cvd_inout_count = df_cvd["DELTA"] <= 50
cvd_inout_count = cvd_inout_count.value_counts().values
cvd_inout_count = np.flip(cvd_inout_count,0)
print("CVD COUNTS = ", cvd_inout_count)

# Difference
difference = full_inout_count - cvd_inout_count
#difference = np.flip(difference, 0)
print("DIFFERENCE:", difference)

# Contingency Table and Fisher:
full = np.c_[cvd_inout_count,difference]
print(full.T)

oddsratio, pvalue = fisher_exact(full)
print(oddsratio, pvalue)

