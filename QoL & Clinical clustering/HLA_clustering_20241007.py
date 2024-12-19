import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import pandas as pd
import seaborn as sns
mpl.rcParams['pdf.fonttype'] = 42
from scipy.stats import chi2_contingency

file_path2 = 'clinical_significant_data_CeD_clean9.csv'
CeD_cluster = pd.read_csv(file_path2)
columns_to_drop = ['PC4', 'PC5', 'PC6', 'PC7', 'PC8', 'PC9', 'PC10', 'PC11', 'PC12', 'PC13', 'PC14', 'PC15', 'PC16',
                    'tag.genotype', 'hibag.A.1', 'hibag.A.2', 'hibag.B.1', 'hibag.B.2', 'hibag.C.1', 'hibag.C.2', 
                    'hibag.DPB1.1', 'hibag.DPB1.2', 'hibag.DQA1.1', 'hibag.DQA1.2', 'hibag.DQB1.1', 'hibag.DQB1.2', 
                    'hibag.DRB1.1', 'hibag.DRB1.2', 'hibag.genotype', 'hla.la.A.1', 'hla.la.A.2', 'hla.la.B.1', 
                    'hla.la.B.2', 'hla.la.C.1', 'hla.la.C.2', 'hla.la.DPB1.1', 'hla.la.DPB1.2', 'hla.la.DQA1.1', 
                    'hla.la.DQA1.2', 'hla.la.DQB1.1', 'hla.la.DQB1.2', 'hla.la.DRB1.1', 'hla.la.DRB1.2', 'hla.la.genotype']
CeD_cluster = CeD_cluster.drop(columns=columns_to_drop)
#EHR data
CeD_cluster.loc[~CeD_cluster[['244.0:Hypothyroidism', '244.4:Hypothyroidism NOS', '245.0:Thyroiditis', '245.21:Chronic lymphocytic thyroiditis', '245.2:Chronic thyroiditis', '261.0:Vitamin deficiency', '261.4:Vitamin D deficiency',  '280.1:Iron deficiency anemias, unspecified or not due to blood loss', '340.0:Migraine', '535.0:Gastritis and duodenitis', '537.0:Other disorders of stomach and duodenum', '558.0:Noninfectious gastroenteritis', '561.0:Symptoms involving digestive system', '561.1:Diarrhea', '564.0:Functional digestive disorders', '564.1:Irritable Bowel Syndrome']].isna().all(axis=1),['244.0:Hypothyroidism', '244.4:Hypothyroidism NOS', '245.0:Thyroiditis', '245.21:Chronic lymphocytic thyroiditis', '245.2:Chronic thyroiditis', '261.0:Vitamin deficiency', '261.4:Vitamin D deficiency', '280.1:Iron deficiency anemias, unspecified or not due to blood loss', '340.0:Migraine', '535.0:Gastritis and duodenitis', '537.0:Other disorders of stomach and duodenum', '558.0:Noninfectious gastroenteritis', '561.0:Symptoms involving digestive system', '561.1:Diarrhea', '564.0:Functional digestive disorders', '564.1:Irritable Bowel Syndrome']]=CeD_cluster.loc[~CeD_cluster[['244.0:Hypothyroidism', '244.4:Hypothyroidism NOS', '245.0:Thyroiditis', '245.21:Chronic lymphocytic thyroiditis', '245.2:Chronic thyroiditis', '261.0:Vitamin deficiency', '261.4:Vitamin D deficiency', '280.1:Iron deficiency anemias, unspecified or not due to blood loss', '340.0:Migraine', '535.0:Gastritis and duodenitis', '537.0:Other disorders of stomach and duodenum', '558.0:Noninfectious gastroenteritis', '561.0:Symptoms involving digestive system', '561.1:Diarrhea', '564.0:Functional digestive disorders', '564.1:Irritable Bowel Syndrome']].isna().all(axis=1),['244.0:Hypothyroidism', '244.4:Hypothyroidism NOS', '245.0:Thyroiditis', '245.21:Chronic lymphocytic thyroiditis', '245.2:Chronic thyroiditis', '261.0:Vitamin deficiency', '261.4:Vitamin D deficiency', '280.1:Iron deficiency anemias, unspecified or not due to blood loss', '340.0:Migraine', '535.0:Gastritis and duodenitis', '537.0:Other disorders of stomach and duodenum', '558.0:Noninfectious gastroenteritis', '561.0:Symptoms involving digestive system', '561.1:Diarrhea', '564.0:Functional digestive disorders', '564.1:Irritable Bowel Syndrome']].fillna(0)

file_path2 = '20241002_HLA.csv'
HLA = pd.read_csv(file_path2)
print(HLA.columns.to_list())

filtered_HLA = HLA[['person_id', 'Impact label', 'genetic_sex', 'diagnosis domain', 'ancestry_pred']]
#filtered_HLA['status']=filtered_HLA['Impact label'].replace(['high','moderate','low','none'],['true CeD','true CeD','false CeD','false CeD'])
CeD_cluster = pd.merge(filtered_HLA, CeD_cluster, on='person_id', how='left')
file_path3 = 'NA_counts.csv'
# Assuming CeD_cluster is your DataFrame and you want to filter where 'CeD' is equal to 2
filtered_data = CeD_cluster[CeD_cluster['CeD'] == 2]
answers = pd.read_csv('survey_answers.csv')
filtered_data=pd.merge(filtered_data,answers[['person_id','Overall Health: Average Pain 7 Days']],on='person_id').drop('Average.Pain.7.Days',axis=1)
# Counting non-NA values for each column in the filtered DataFrame
filtered_data.replace(['PMI: Skip','No matching concept'],[np.nan,np.nan],inplace=True)
for i in ['6.or.More.Drinks Last year','Average.Daily.Drink.Count','Education.Level','General.Doctor.Visits','Respected.By.Provider','Spoken.To.Medical.Specialist','How.often.listening.to.doctor','Annual.Income','Health.Insurance','Smoke.Frequency']:
	filtered_data[i].replace(0,np.nan,inplace=True)

filtered_data['Health.Insurance'].replace([0,1,2],[np.nan,0,1],inplace=True)
filtered_data['Spoken.To.Medical.Specialist'].replace([0,1,2],[np.nan,0,1],inplace=True)
non_na_counts = filtered_data.count()
# Filter out columns where the non-NA count is less than 100
columns_to_keep = non_na_counts[non_na_counts >= 965].index
print(columns_to_keep)
# Update filtered_data to only include columns with non-NA count >= 100
filtered_data = filtered_data[columns_to_keep]
df = pd.read_csv(file_path3)
# Create a filter to exclude rows where 'To_remove' equals 'drop'
filtered_df = df[df['drop'] == 'drop']
 # Extract the values from the 'columns' column and store them in a list
columns_to_drop3 = filtered_df['Column'].tolist()

#consolidate survey/ehr results
def consolidate(columns=[]):
  result=filtered_data[columns].sum(axis=1)
  result[result>1]=1
  filtered_data.drop(columns,axis=1,inplace=True)
  return result


filtered_data['IBS']=consolidate(['564.1:Irritable Bowel Syndrome','irritable bowel syndrome (IBS)'])
filtered_data['Hypothyroidism']=consolidate([ '244.0:Hypothyroidism', 'hypothyroidism','244.4:Hypothyroidism NOS',])
filtered_data['Vitamin D deficiency']=consolidate([ '261.4:Vitamin D deficiency', 'vitamin D deficiency'])
filtered_data['Anemia']=consolidate([  '280.1:Iron deficiency anemias, unspecified or not due to blood loss', 'anemia'])
filtered_data['Migraine']=consolidate(['340.0:Migraine',  'migraine headaches'])
filtered_data['Thyroiditis']=consolidate(['245.0:Thyroiditis','245.2:Chronic thyroiditis','245.21:Chronic lymphocytic thyroiditis'])
filtered_data['Familial CeD']=consolidate(['celiac disease- Daughter', 'celiac disease- Father', 'celiac disease- Grandparent', 'celiac disease- Mother', 'celiac disease- Sibling', 'celiac disease- Son'])


filtered_data.drop(['PH_raw', 'Global_raw','MH_raw','PROMIS..Phyisical.health', 'PROMIS..Mental.health'],axis=1,inplace=True)
filtered_data.drop(filtered_data.columns[filtered_data.columns.str.contains('-diagnosis')],axis=1,inplace=True)
filtered_data.drop(filtered_data.columns[filtered_data.columns.str.contains('Daughter|Father|Grandparent|Mother|Son|Sibling')&~filtered_data.columns.str.contains('celiac')],axis=1,inplace=True)
file_path4 = '20240610_survey_cats_fdr.csv'
surveys = pd.read_csv(file_path4)



print(filtered_data.columns)
# Ensure filtered_data is correctly defined earlier in your code or loaded from your data source
count_cols=filtered_data[['261.0:Vitamin deficiency','535.0:Gastritis and duodenitis', '537.0:Other disorders of stomach and duodenum','558.0:Noninfectious gastroenteritis',        '561.0:Symptoms involving digestive system', '561.1:Diarrhea',        '564.0:Functional digestive disorders',"Crohn's disease",'acid reflux', 'allergies', 'an eating disorder','anxiety reaction/panic disorder', 'asthma', 'autism spectrum disorder','chronic fatigue', 'depression', 'endometriosis', 'fibroids', 'fibromyalgia', 'hyperthyroidism', 'memory loss or impairment','neuropathy', 'osteoporosis', 'peptic (stomach) ulcers', 'polycystic ovarian syndrome','reactions to anesthesia (such as hyperthermia)','restless leg syndrome', 'skin condition(s) (e.g., eczema, psoriasis)','spine, muscle, or bone disorders (non-cancer)', 'systemic lupus','type 1 diabetes', 'ulcerative colitis', 'vitamin B deficiency','Chronic sinus infections', 'Lyme disease','Recurrent urinary tract infections (UTI)/bladder infections','Recurrent yeast infections', 'West Nile virus','IBS', 'Hypothyroidism','Vitamin D deficiency', 'Anemia', 'Migraine', 'Thyroiditis','Familial CeD']].sum(axis=0)
drop_cols=count_cols[count_cols<=20]
filtered_data.drop(drop_cols.index,inplace=True,axis=1)
filtered_data['Impact label_code'] = filtered_data['Impact label'].replace(['none','low','moderate','high'],np.array(range(4)))
#filtered_data['ancestry_pred_code'] = filtered_data['ancestry_pred'].replace(['eur','afr','eas','sas','mid','amr'],np.array(range(6))+1)
filtered_data['diagnosis domain_code'] = filtered_data['diagnosis domain'].replace(['survey','ehr','survey+ehr'],np.array(range(3))+1)
filtered_data['is_gen_female_code'] = filtered_data['genetic_sex'].replace(['male','female'],np.array(range(2)))
filtered_data['is_hispanic'] = filtered_data['ethnicity'].replace(3,1)
columns_to_drop = ['Impact label', 'genetic_sex', 'diagnosis domain', 'ancestry_pred','ethnicity']




filtered_data.drop(['CeD','is_hispanic','ALP(U/L)', 'BMI(kg/m2)', 'T(°C)', 'weight(kg)', 'Ca(md/dl)',
        'Cholesterol(mg/dl)', 'SCrea(mg/dl)', 'RBC(x10^6/ul)', 'Glu(mg/dl)',
        'HR(bpm)', 'Hb(g/l)', 'WBC(10^3/ul)', 'MCHC(g/dl)', 'PLT(10^3/ul)',
        'SBP(mmHg)', 'race', 'sex_at_birth','6.or.More.Drinks Last year', 'Average.Daily.Drink.Count', 'Education.Level', 'General.Doctor.Visits', 'Respected.By.Provider', 'Spoken.To.Medical.Specialist', 'Annual.Income', 'Health.Insurance', 	'Average.Fatigue.7.Days', 'Emotional.Problem.7.Days', 'Everyday.Activities', 'General.Mental.Health','General.Physical.Health',
       'General.Quality', 'Social.Satisfaction' ],axis=1,inplace=True)

print(filtered_data.drop(['person_id'],axis=1).columns)
#filtered_data['ancestry_pred'].replace(['sas','eas','mid'],['other','other','other'],inplace=True)
total_ct=filtered_data.value_counts('Impact label')
cond_list = ['261.0:Vitamin deficiency', 
        '537.0:Other disorders of stomach and duodenum',
        '558.0:Noninfectious gastroenteritis',
        '561.0:Symptoms involving digestive system', '561.1:Diarrhea',
        '564.0:Functional digestive disorders', "Crohn's disease",
        'acid reflux', 'allergies', 'an eating disorder',
        'anxiety reaction/panic disorder', 'asthma', 'autism spectrum disorder',
        'chronic fatigue', 'depression','hyperthyroidism', 'memory loss or impairment',
        'neuropathy', 'osteoporosis', 'peptic (stomach) ulcers',
        'reactions to anesthesia (such as hyperthermia)',
        'restless leg syndrome', 'skin condition(s) (e.g., eczema, psoriasis)',
        'spine, muscle, or bone disorders (non-cancer)', 'systemic lupus',
        'type 1 diabetes', 'ulcerative colitis', 'vitamin B deficiency',
        'Chronic sinus infections', 'Lyme disease',
        'Recurrent urinary tract infections (UTI)/bladder infections',
        'Recurrent yeast infections',
         'IBS', 'Hypothyroidism',
        'Vitamin D deficiency', 'Anemia', 'Migraine', 'Thyroiditis',
        'Familial CeD','is_gen_female_code']
fem_spec_cond_list=['endometriosis', 'fibroids', 'fibromyalgia','polycystic ovarian syndrome']
percent_table=filtered_data[cond_list+['Impact label']].fillna(0).groupby('Impact label').agg('mean').T*100
filtered_fem=filtered_data[filtered_data['is_gen_female_code']==1]
fem_total = filtered_fem.value_counts('Impact label')
fem_spec_percents=filtered_fem[fem_spec_cond_list+['Impact label']].fillna(0).groupby('Impact label').agg('mean').T*100
fem_spec_percents['note']='fem specific percent'
total_percent=pd.concat([percent_table,fem_spec_percents])
total_percent.loc['total']=total_ct
total_percent.loc['fem total']=fem_total

counts=filtered_data[cond_list+fem_spec_cond_list+['Impact label']].groupby('Impact label').agg('sum').T
#counts.drop('other',axis=1,inplace=True)


for i in cond_list:
  df=counts.loc[i,['none','low','moderate','high']]
  df_alt=total_ct-counts.loc[i,['none','low','moderate','high']]
  test_df=pd.concat([df,df_alt],axis=1)
  res=chi2_contingency(test_df)
  total_percent.loc[i,'pvalue']=res.pvalue
  total_percent.loc[i,'X2']=res.statistic
  counts.loc[i,'pvalue']=res.pvalue
  counts.loc[i,'X2']=res.statistic


for i in fem_spec_cond_list:
  df=counts.loc[i,['none','low','moderate','high']]
  df_alt=fem_total-counts.loc[i,['none','low','moderate','high']]
  test_df=pd.concat([df,df_alt],axis=1)
  res=chi2_contingency(test_df)
  total_percent.loc[i,'pvalue']=res.pvalue
  total_percent.loc[i,'X2']=res.statistic
  counts.loc[i,'pvalue']=res.pvalue
  counts.loc[i,'X2']=res.statistic

#ancestry specific values
for i in cond_list:
 for j in ['none','low','moderate']:
  df=counts.loc[i,[j,'high']]
  df_alt=total_ct[[j,'high']]-counts.loc[i,[j,'high']]
  test_df=pd.concat([df,df_alt],axis=1)
  res=chi2_contingency(test_df)
  OR=df[j]*df_alt['high']/(df_alt[j]*df['high'])
  total_percent.loc[i,j+' pvalue']=res.pvalue
  total_percent.loc[i,j+' OR']=OR
  counts.loc[i,j+' pvalue']=res.pvalue
  counts.loc[i,j+' OR']=OR

for i in fem_spec_cond_list:
 for j in ['none','low','moderate']:
  df=counts.loc[i,[j,'high']]
  df_alt=fem_total[[j,'high']]-counts.loc[i,[j,'high']]
  test_df=pd.concat([df,df_alt],axis=1)
  res=chi2_contingency(test_df)
  OR=df[j]*df_alt['high']/(df_alt[j]*df['high'])
  total_percent.loc[i,j+' pvalue']=res.pvalue
  total_percent.loc[i,j+' OR']=OR
  counts.loc[i,j+' pvalue']=res.pvalue
  counts.loc[i,j+' OR']=OR


total_percent.to_csv('HLA_percent_v2.csv')
counts.to_csv('HLA_count_v2.csv')
filtered_data.to_csv('HLA_clust_intermediate_20241007.csv')

import pandas as pd
import numpy as np
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import matplotlib.pyplot as plt

# Load data
data = pd.read_csv('HLA_clust_intermediate.csv')

# Drop the 'Unnamed' column if it exists and other specified columns
columns_to_drop = ['Unnamed: 0', 'person_id', 'PC1', 'PC2', 'PC3']
data = data.drop(columns=[col for col in columns_to_drop if col in data.columns])


# Define numerical and categorical variables

numerical_columns = ['age','PROMIS.PH', 'PROMIS.MH', 'PROMIS.total','Overall Health: Average Pain 7 Days']


# Ensure 'Cluster' is treated as a categorical column
if 'ancestry_pred' not in categorical_columns:
    categorical_columns.append('ancestry_pred')

# Exclude 'Cluster' from numerical columns if mistakenly included
if 'ancestry_pred' in numerical_columns:
    numerical_columns.remove('ancestry_pred')

# Standardize numerical data
scaler = StandardScaler()
data_standardized = data[numerical_columns].copy()  # Create a copy to avoid modifying original data
data_standardized[numerical_columns] = scaler.fit_transform(data_standardized[numerical_columns])

# Add 'Cluster' back to the standardized data for grouping
data_standardized['ancestry_pred'] = data['ancestry_pred']

# Compute medians for numerical variables grouped by 'Cluster'
cluster_medians = data_standardized.groupby('ancestry_pred')[numerical_columns].median()
cluster_medians.rename({'age':'Age','PROMIS.PH':'PROMIS-Physical Health', 'PROMIS.MH':'PROMIS-Mental Health', 'PROMIS.total':'PROMIS-Overall','Overall Health: Average Pain 7 Days':'Pain'},axis=1,inplace=True)
# Transpose the result for better visualization (clusters on x-axis, variables on y-axis)
cluster_medians_transposed = cluster_medians.T

cluster_medians2 = data.groupby('ancestry_pred')[numerical_columns].median()
cluster_medians2.rename({'age':'Age','PROMIS.PH':'PROMIS-Physical Health', 'PROMIS.MH':'PROMIS-Mental Health', 'PROMIS.total':'PROMIS-Overall','Overall Health: Average Pain 7 Days':'Pain'},axis=1,inplace=True)
cluster_medians2_transposed = cluster_medians2.T
cluster_medians2_transposed.to_csv('numeric_medians.csv')
# Plotting the heatmap
plt.figure(figsize=(12, 10))
sns.heatmap(cluster_medians_transposed, cmap='coolwarm', cbar_kws={'label': 'Median Standardized Value'})
plt.title('Heatmap of Median Values for Standardized Numerical Variables by ancestry')
plt.xlabel('ancestry')
plt.ylabel('Numerical Variables')
plt.xticks(rotation=0)  # No need to rotate x-ticks as they are cluster numbers
plt.savefig('HLA figures 20240926/HLA_clust_intermediate.pdf',bbox_inches='tight',format='pdf')


##I have only run upto here


