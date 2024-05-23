### Preprocessing of data on ICU length of stay

# This code preprocesses from the MIMIC-IV data for the application to ICU
# length of stay in the paper "The Multivariate Bernoulli detector: Change point
# estimation in discrete survival analysis" by Willem van den Boom, Maria De
# Iorio, Fang Qian and Alessandra Guglielmi.

# The script is a modified version of the code by Tomer Mei available at
# https://github.com/tomer1812/DiscreteTimeSurvivalPenalization/blob/89693ebcb0f7db098b22c33b203769bc99e656e7/notebooks/mimiciv.ipynb

import pandas as pd
import sys
import os
sys.path.append('../')
import numpy as np
from sklearn.preprocessing import StandardScaler



# The following constants have been copied from
# https://github.com/tomer1812/DiscreteTimeSurvivalPenalization/blob/89693ebcb0f7db098b22c33b203769bc99e656e7/src/constants.py
ADMISSION_TIME_COL = 'admittime'
DISCHARGE_TIME_COL = 'dischtime'
DEATH_TIME_COL = 'deathtime'
ED_REG_TIME = 'edregtime'
ED_OUT_TIME = 'edouttime'
AGE_COL = 'anchor_age' 
GENDER_COL = 'gender'

YEAR_GROUP_COL = 'anchor_year_group'
SUBSET_YEAR_GROUP = '2017 - 2019'
SUBJECT_ID_COL = 'subject_id'
ADMISSION_ID_COL = 'hadm_id'
ADMISSION_TYPE_COL = 'admission_type'
CHART_TIME_COL = 'charttime'
STORE_TIME_COL = 'storetime'
LOS_EXACT_COL = 'LOS exact'
LOS_DAYS_COL = 'LOS days'
ADMISSION_LOCATION_COL = 'admission_location'
DISCHARGE_LOCATION_COL = 'discharge_location'
RACE_COL = 'race'
INSURANCE_COL = 'insurance'
ADMISSION_TO_RESULT_COL = 'admission_to_result_time'
ADMISSION_AGE_COL = 'admission_age'
ADMISSION_YEAR_COL = 'admission_year'
ADMISSION_COUNT_COL = 'admissions_count'
ITEM_ID_COL = 'itemid'
NIGHT_ADMISSION_FLAG = 'night_admission' 
MARITAL_STATUS_COL = 'marital_status'
STANDARDIZED_AGE_COL = 'standardized_age'
COEF_COL = '   coef   '
STDERR_COL = ' std err '
DIRECT_IND_COL = 'direct_emrgency_flag'
PREV_ADMISSION_IND_COL = 'last_less_than_diff'
ADMISSION_COUNT_GROUP_COL = ADMISSION_COUNT_COL + '_group'


DISCHARGE_REGROUPING_DICT = {
    'HOME': 'HOME',
    'HOME HEALTH CARE': 'HOME',
    'SKILLED NURSING FACILITY': 'FURTHER TREATMENT',
    'DIED': 'DIED',
    'REHAB': 'HOME',
    'CHRONIC/LONG TERM ACUTE CARE': 'FURTHER TREATMENT',
    'HOSPICE': 'FURTHER TREATMENT',
    'AGAINST ADVICE': 'CENSORED',
    'ACUTE HOSPITAL': 'FURTHER TREATMENT',
    'PSYCH FACILITY': 'FURTHER TREATMENT',
    'OTHER FACILITY': 'FURTHER TREATMENT',
    'ASSISTED LIVING': 'HOME',
    'HEALTHCARE FACILITY': 'FURTHER TREATMENT',
}

RACE_REGROUPING_DICT = {
    'WHITE': 'WHITE',
    'UNKNOWN': 'OTHER',
    'BLACK/AFRICAN AMERICAN': 'BLACK',
    'OTHER': 'OTHER',
    'ASIAN': 'ASIAN',
    'WHITE - OTHER EUROPEAN': 'WHITE',
    'HISPANIC/LATINO - PUERTO RICAN': 'HISPANIC',
    'HISPANIC/LATINO - DOMINICAN': 'HISPANIC',
    'ASIAN - CHINESE': 'ASIAN',
    'BLACK/CARIBBEAN ISLAND': 'BLACK',
    'BLACK/AFRICAN': 'BLACK',
    'BLACK/CAPE VERDEAN': 'BLACK',
    'PATIENT DECLINED TO ANSWER': 'OTHER',
    'WHITE - BRAZILIAN': 'WHITE',
    'PORTUGUESE': 'HISPANIC', 
    'ASIAN - SOUTH EAST ASIAN': 'ASIAN',
    'WHITE - RUSSIAN': 'WHITE',
    'ASIAN - ASIAN INDIAN': 'ASIAN',
    'WHITE - EASTERN EUROPEAN': 'WHITE',
    'AMERICAN INDIAN/ALASKA NATIVE': 'OTHER',
    'HISPANIC/LATINO - GUATEMALAN': 'HISPANIC',
    'HISPANIC/LATINO - MEXICAN': 'HISPANIC',
    'HISPANIC/LATINO - SALVADORAN': 'HISPANIC',
    'SOUTH AMERICAN': 'HISPANIC',
    'NATIVE HAWAIIAN OR OTHER PACIFIC ISLANDER': 'OTHER',
    'HISPANIC/LATINO - COLUMBIAN': 'HISPANIC',
    'HISPANIC/LATINO - CUBAN': 'HISPANIC',
    'ASIAN - KOREAN': 'ASIAN',
    'HISPANIC/LATINO - HONDURAN': 'HISPANIC',
    'HISPANIC/LATINO - CENTRAL AMERICAN': 'HISPANIC',
    'UNABLE TO OBTAIN': 'OTHER',
    'HISPANIC OR LATINO': 'HISPANIC'
}

# 'MCH': 'Mean Cell Hemoglobin', 
# 'MCHC': 'Mean Cell Hemoglobin Concentration',
table1_rename_columns = {
    'AnionGap': 'Anion Gap', 
    'Bicarbonate': 'Bicarbonate', 
    'CalciumTotal': 'Calcium Total', 
    'Chloride': 'Chloride', 
    'Creatinine': 'Creatinine',
    'Glucose': 'Glucose', 
    'Magnesium': 'Magnesium', 
    'Phosphate': 'Phosphate', 
    'Potassium': 'Potassium', 
    'Sodium': 'Sodium',
    'UreaNitrogen': 'Urea Nitrogen', 
    'Hematocrit': 'Hematocrit', 
    'Hemoglobin': 'Hemoglobin', 
    'MCH': 'MCH', 
    'MCHC': 'MCHC', 
    'MCV': 'MCV',
    'PlateletCount': 'Platelet Count', 
    'RDW': 'RDW', 
    'RedBloodCells': 'Red Blood Cells', 
    'WhiteBloodCells': 'White Blood Cells',
    NIGHT_ADMISSION_FLAG: 'Night Admission', 
    GENDER_COL: 'Sex', 
    DIRECT_IND_COL: 'Direct Emergency',
    PREV_ADMISSION_IND_COL: 'Previous Admission This Month', 
    ADMISSION_AGE_COL: 'Admission Age', 
    INSURANCE_COL: 'Insurance', 
    MARITAL_STATUS_COL: 'Marital Status',
    RACE_COL: 'Race', 
    ADMISSION_COUNT_GROUP_COL: 'Admissions Number', 
    LOS_DAYS_COL: 'LOS (days)', 
    DISCHARGE_LOCATION_COL: 'Discharge Location'
}

table1_rename_sex = {0: 'Male', 1: 'Female'}
table1_rename_race = {'ASIAN': 'Asian', 'BLACK': 'Black', 'HISPANIC': 'Hispanic', 'OTHER': 'Other', 'WHITE': 'White'}
table1_rename_marital = {'SINGLE': 'Single', 'MARRIED': 'Married', 'DIVORCED': 'Divorced', 'WIDOWED': 'Widowed'}
table1_rename_yes_no = {0: 'No', 1: 'Yes'}
table1_rename_normal_abnormal = {0: 'Normal', 1: 'Abnormal'}
table1_rename_discharge = {1: 'Home', 2: 'Further Treatment', 3: 'Died', 0: 'Censored'} 


rename_beta_index = {
 'AdmsCount 2': 'Admissions Number 2',
 'AdmsCount 3up': 'Admissions Number 3+',
 'AnionGap': 'Anion Gap',
 'Bicarbonate': 'Bicarbonate',
 'CalciumTotal': 'Calcium Total',
 'Chloride': 'Chloride',
 'Creatinine': 'Creatinine',
 'Ethnicity BLACK': 'Ethnicity Black',
 'Ethnicity HISPANIC': 'Ethnicity Hispanic',
 'Ethnicity OTHER': 'Ethnicity Other',
 'Ethnicity WHITE': 'Ethnicity White',
 'Glucose': 'Glucose',
 'Hematocrit': 'Hematocrit',
 'Hemoglobin': 'Hemoglobin',
 'Insurance Medicare': 'Insurance Medicare',
 'Insurance Other': 'Insurance Other',
 'MCH': 'MCH',
 'MCHC': 'MCHC',
 'MCV': 'MCV',
 'Magnesium': 'Magnesium',
 'Marital MARRIED': 'Marital Married',
 'Marital SINGLE': 'Marital Single',
 'Marital WIDOWED': 'Marital Widowed',
 'Phosphate': 'Phosphate',
 'PlateletCount': 'Platelet Count',
 'Potassium': 'Potassium',
 'RDW': 'RDW',
 'RedBloodCells': 'Red Blood Cells',
 'Sodium': 'Sodium',
 'UreaNitrogen': 'Urea Nitrogen',
 'WhiteBloodCells': 'White Blood Cells',
 'direct emrgency flag': 'Direct Emergency',
 'gender': 'Sex',
 'last less than diff': 'Recent Admission',
 'night admission': 'Night Admission',
 'standardized age': 'Standardized Age',
}

beta_units = {
    'Admissions Number 2': '2',
    'Admissions Number 3+': '3+',
    'Anion Gap': 'Abnormal',
    'Bicarbonate': 'Abnormal',
    'Calcium Total': 'Abnormal',
    'Chloride': 'Abnormal',
    'Creatinine': 'Abnormal',
    'Ethnicity Black': 'Black',
    'Ethnicity Hispanic': 'Hispanic',
    'Ethnicity Other': 'Other',
    'Ethnicity White': 'White',
    'Glucose': 'Abnormal',
    'Hematocrit': 'Abnormal',
    'Hemoglobin': 'Abnormal',
    'Insurance Medicare': 'Medicare',
    'Insurance Other': 'Other',
    'MCH': 'Abnormal',
    'MCHC': 'Abnormal',
    'MCV': 'Abnormal',
    'Magnesium': 'Abnormal',
    'Marital Married': 'Married',
    'Marital Single': 'Single',
    'Marital Widowed': 'Widowed',
    'Phosphate': 'Abnormal',
    'Platelet Count': 'Abnormal',
    'Potassium': 'Abnormal',
    'RDW': 'Abnormal',
    'Red Blood Cells': 'Abnormal',
    'Sodium': 'Abnormal',
    'Urea Nitrogen': 'Abnormal',
    'White Blood Cells': 'Abnormal',
    'Direct Emergency': 'Yes',
    'Sex': 'Female',
    'Recent Admission': 'Yes',
    'Night Admission': 'Yes',
    'Standardized Age': '',
}



# The required .csv.gz files with MIMIC-IV data are available from
# https://physionet.org/content/mimiciv/2.0/hosp/#files-panel
# (login required for access).
DATA_DIR = 'mimiciv/2.0/'


# Load Data
patients_file = os.path.join(DATA_DIR, 'hosp', 'patients.csv.gz')
admissions_file = os.path.join(DATA_DIR, 'hosp', 'admissions.csv.gz')
lab_file = os.path.join(DATA_DIR, 'hosp', 'labevents.csv.gz')
lab_meta_file = os.path.join(DATA_DIR, 'hosp', 'd_labitems.csv.gz')

patients_df = pd.read_csv(patients_file, compression='gzip')

COLUMNS_TO_DROP = ['dod']
patients_df.drop(COLUMNS_TO_DROP, axis=1, inplace=True)

admissions_df = pd.read_csv(admissions_file, compression='gzip', parse_dates=[ADMISSION_TIME_COL,
                            DISCHARGE_TIME_COL, DEATH_TIME_COL, ED_REG_TIME, ED_OUT_TIME])

COLUMNS_TO_DROP = ['hospital_expire_flag', 'edouttime', 'edregtime', 'deathtime', 'language']
admissions_df.drop(COLUMNS_TO_DROP, axis=1, inplace=True)

admissions_df = admissions_df.merge(patients_df, on=[SUBJECT_ID_COL])


# Calculate Age at Admission and Group of Admission Year

# Based on mimic IV example https://mimic.mit.edu/docs/iv/modules/hosp/patients/


# Diff column first
admissions_df[ADMISSION_YEAR_COL] = (admissions_df[ADMISSION_TIME_COL].dt.year - admissions_df['anchor_year'])

# Age at admission calculation
admissions_df[ADMISSION_AGE_COL] = (admissions_df[AGE_COL] + admissions_df[ADMISSION_YEAR_COL])

# Admission year group lower bound calculation
admissions_df[ADMISSION_YEAR_COL] = admissions_df[ADMISSION_YEAR_COL] + admissions_df[YEAR_GROUP_COL].apply(lambda x: int(x.split(' ')[0]))


# Calculating LOS (exact, days resolution) and night admission indicator
admissions_df[LOS_EXACT_COL] = (admissions_df[DISCHARGE_TIME_COL] - admissions_df[ADMISSION_TIME_COL])
admissions_df[NIGHT_ADMISSION_FLAG] = ((admissions_df[ADMISSION_TIME_COL].dt.hour >= 20) |                                        (admissions_df[ADMISSION_TIME_COL].dt.hour < 8) ).values
admissions_df[LOS_DAYS_COL] = admissions_df[LOS_EXACT_COL].dt.ceil('1d')


# Taking only SPECIFIC_ADMISSION_TYPE admissions from now on
SPECIFIC_ADMISSION_TYPE = ['DIRECT EMER.', 'EW EMER.']

admissions_df = admissions_df[admissions_df[ADMISSION_TYPE_COL].isin(SPECIFIC_ADMISSION_TYPE)]


# add direct emergency if needed
if 'DIRECT EMER.' in SPECIFIC_ADMISSION_TYPE:
    admissions_df[DIRECT_IND_COL] = (admissions_df[ADMISSION_TYPE_COL] == 'DIRECT EMER.').astype(int)


# Counting SPECIFIC_ADMISSION_TYPE admissions to each patient 
number_of_admissions = admissions_df.groupby(SUBJECT_ID_COL)[ADMISSION_ID_COL].nunique()
number_of_admissions.name = ADMISSION_COUNT_COL

admissions_df = admissions_df.merge(number_of_admissions, on=SUBJECT_ID_COL)


# Add recurrent admissions group per patient according to last admission
ADMISSION_COUNT_BINS = [1, 1.5, 2.5, 5000]
ADMISSION_COUNT_LABELS = ['1', '2', '3up']

admissions_df[ADMISSION_COUNT_GROUP_COL] = pd.cut(admissions_df[ADMISSION_COUNT_COL], 
                                                  bins=ADMISSION_COUNT_BINS, 
                                                  labels=ADMISSION_COUNT_LABELS, 
                                                  include_lowest=True)


# Adds last admission with previous admission in past month indicator
indicator_diff = pd.to_timedelta('30d')

tmp_admissions = admissions_df[admissions_df[ADMISSION_COUNT_COL] > 1]

ind_ser = tmp_admissions.sort_values(by=[SUBJECT_ID_COL, ADMISSION_TIME_COL]).groupby(
    SUBJECT_ID_COL).apply(
    lambda tmp_df: (tmp_df[ADMISSION_TIME_COL] - tmp_df[DISCHARGE_TIME_COL].shift(1)) <= indicator_diff)

ind_ser.index = ind_ser.index.droplevel(1)
ind_ser.name = PREV_ADMISSION_IND_COL
ind_ser = ind_ser.iloc[ind_ser.reset_index().drop_duplicates(subset=[SUBJECT_ID_COL], keep='last').index]
ind_ser

admissions_df = admissions_df.merge(ind_ser.astype(int), left_on=SUBJECT_ID_COL, right_index=True, how='outer')
admissions_df[PREV_ADMISSION_IND_COL].fillna(0, inplace=True)


# Keep only last admission per patient
only_last_admission = admissions_df.sort_values(by=[ADMISSION_TIME_COL]).drop_duplicates(subset=[SUBJECT_ID_COL], keep='last')


# Take only patients with last admission after MINIMUM YEAR
MINIMUM_YEAR = 2014
only_last_admission = only_last_admission[only_last_admission[ADMISSION_YEAR_COL] >= MINIMUM_YEAR]

pids = only_last_admission[SUBJECT_ID_COL].drop_duplicates()
adm_ids = only_last_admission[ADMISSION_ID_COL].drop_duplicates()


# Load relevant lab tests
LOAD_SPECIFIC_COLUMNS = [SUBJECT_ID_COL, ADMISSION_ID_COL, ITEM_ID_COL, 'storetime', 'flag']

print("Reading in labevents.csv.gz...")
chunksize = 10 ** 6
full_df = pd.DataFrame()
with pd.read_csv(lab_file, chunksize=chunksize, compression='gzip', parse_dates=[STORE_TIME_COL], usecols=LOAD_SPECIFIC_COLUMNS) as reader:
    for chunk in reader:
        tmp_chunk = chunk[chunk[SUBJECT_ID_COL].isin(pids) & chunk[ADMISSION_ID_COL].isin(adm_ids)]
        tmp_adms = only_last_admission[only_last_admission[SUBJECT_ID_COL].isin(pids) & only_last_admission[ADMISSION_ID_COL].isin(adm_ids)]
        tmp_chunk = tmp_chunk.merge(tmp_adms, on=[SUBJECT_ID_COL, ADMISSION_ID_COL])
        full_df = pd.concat([full_df, tmp_chunk])
        print(len(full_df))


# Continue only with included patients_df and admissions_df and full_df
pids = full_df[SUBJECT_ID_COL].drop_duplicates().values
adms_ids = full_df[ADMISSION_ID_COL].drop_duplicates().values
patients_df = patients_df[patients_df[SUBJECT_ID_COL].isin(pids)]
admissions_df = admissions_df[admissions_df[ADMISSION_ID_COL].isin(adms_ids)]


# Regrouping discharge location
discharge_regrouping_df = pd.Series(DISCHARGE_REGROUPING_DICT).to_frame()
discharge_regrouping_df.index.name = 'Original Group'
discharge_regrouping_df.columns = ['Regrouped']

admissions_df[DISCHARGE_LOCATION_COL].replace(DISCHARGE_REGROUPING_DICT, inplace=True)
full_df[DISCHARGE_LOCATION_COL].replace(DISCHARGE_REGROUPING_DICT, inplace=True)


# Regroup Race
race_regrouping_df = pd.Series(RACE_REGROUPING_DICT).to_frame()
race_regrouping_df.index.name = 'Original Group'
race_regrouping_df.columns = ['Regrouped']
race_regrouping_df

admissions_df[RACE_COL].replace(RACE_REGROUPING_DICT, inplace=True)
full_df[RACE_COL].replace(RACE_REGROUPING_DICT, inplace=True)


# Taking only results 24 hours from admission
full_df[ADMISSION_TO_RESULT_COL] = (full_df[STORE_TIME_COL] - full_df[ADMISSION_TIME_COL])

full_df = full_df[full_df[ADMISSION_TO_RESULT_COL] <= pd.to_timedelta('1d')]

full_df.sort_values(by=[ADMISSION_TIME_COL, STORE_TIME_COL]).drop_duplicates(subset=[SUBJECT_ID_COL, ADMISSION_ID_COL, ITEM_ID_COL], 
    inplace=True, keep='last')


# Most common lab tests upon arrival
lab_meta_df = pd.read_csv(lab_meta_file, compression='gzip')

threshold = 25000

common_tests = full_df.groupby(ITEM_ID_COL)[ADMISSION_ID_COL].nunique().sort_values(ascending=False)
included_in_threshold = common_tests[common_tests > threshold].to_frame().merge(lab_meta_df, on=ITEM_ID_COL)

full_df = full_df[full_df[ITEM_ID_COL].isin(included_in_threshold[ITEM_ID_COL].values)]

minimal_item_id = included_in_threshold.iloc[-1][ITEM_ID_COL]

pids = full_df[full_df[ITEM_ID_COL] == minimal_item_id][SUBJECT_ID_COL].drop_duplicates().values
adms_ids = full_df[full_df[ITEM_ID_COL] == minimal_item_id][ADMISSION_ID_COL].drop_duplicates().values
patients_df = patients_df[patients_df[SUBJECT_ID_COL].isin(pids)]
admissions_df = admissions_df[admissions_df[ADMISSION_ID_COL].isin(adms_ids)]
full_df = full_df[full_df[SUBJECT_ID_COL].isin(pids)]
full_df = full_df[full_df[ADMISSION_ID_COL].isin(adms_ids)]

full_df['flag'].fillna('normal', inplace=True)
full_df['flag'].replace({'normal': 0, 'abnormal':1}, inplace=True)

full_df = full_df.sort_values(by=[ADMISSION_TIME_COL, STORE_TIME_COL]).drop_duplicates(
    subset=[SUBJECT_ID_COL, ADMISSION_ID_COL, ITEM_ID_COL], 
    keep='last')

tmp = full_df[[SUBJECT_ID_COL, ADMISSION_ID_COL, ITEM_ID_COL, 'flag']]
fitters_table = pd.pivot_table(tmp, values=['flag'], index=[SUBJECT_ID_COL, ADMISSION_ID_COL], 
                               columns=[ITEM_ID_COL], aggfunc=np.sum)

fitters_table = fitters_table.droplevel(1, axis=0).droplevel(0, axis=1)

dummies_df = full_df.drop_duplicates(subset=[SUBJECT_ID_COL]).set_index(SUBJECT_ID_COL)

del full_df
del admissions_df
del patients_df


# Standardize age
scaler = StandardScaler()
dummies_df[STANDARDIZED_AGE_COL] = scaler.fit_transform(dummies_df[[AGE_COL]])

J_DICT = {'HOME': 1, 'FURTHER TREATMENT': 2, 'DIED': 3, 'CENSORED': 0} 
GENDER_DICT = {'F': 1, 'M': 0}

dummies_df[GENDER_COL] = dummies_df[GENDER_COL].replace(GENDER_DICT)


# Table 1
included_in_threshold['label'] = included_in_threshold['label'].apply(lambda x: x.replace(' ', '')).apply(lambda x: x.replace(',', ''))
RENAME_ITEMS_DICT = included_in_threshold[[ITEM_ID_COL, 'label']].set_index(ITEM_ID_COL).to_dict()['label']

table1 = pd.concat([
    fitters_table.copy(),
    dummies_df[[NIGHT_ADMISSION_FLAG,
                GENDER_COL, 
                DIRECT_IND_COL,
                PREV_ADMISSION_IND_COL,
                ADMISSION_AGE_COL]].astype(int),
    dummies_df[[INSURANCE_COL,
                MARITAL_STATUS_COL,
                RACE_COL,
                ADMISSION_COUNT_GROUP_COL]],
    dummies_df[LOS_DAYS_COL].dt.days,
    dummies_df[DISCHARGE_LOCATION_COL].dropna().replace(J_DICT).astype(int)
], axis=1)
    
table1.rename(RENAME_ITEMS_DICT, inplace=True, axis=1)  
table1.dropna(inplace=True)
table1

ADMINISTRATIVE_CENSORING = 28
censoring_index = table1[table1[LOS_DAYS_COL] > ADMINISTRATIVE_CENSORING].index
table1.loc[censoring_index, DISCHARGE_LOCATION_COL] = 0
table1.loc[censoring_index, LOS_DAYS_COL] = ADMINISTRATIVE_CENSORING + 1

table1[GENDER_COL].replace(table1_rename_sex, inplace=True)
table1[RACE_COL].replace(table1_rename_race, inplace=True)
table1[MARITAL_STATUS_COL].replace(table1_rename_marital, inplace=True)
table1[DIRECT_IND_COL].replace(table1_rename_yes_no, inplace=True)
table1[NIGHT_ADMISSION_FLAG].replace(table1_rename_yes_no, inplace=True)
table1[PREV_ADMISSION_IND_COL].replace(table1_rename_yes_no, inplace=True)
table1[DISCHARGE_LOCATION_COL].replace(table1_rename_discharge, inplace=True)
table1[ADMISSION_COUNT_GROUP_COL].replace({'3up': '3+'}, inplace=True)
table1.rename(table1_rename_columns, inplace=True, axis=1)

table1.dropna(inplace=True)

fitters_table = pd.concat([
    fitters_table.copy(),
    pd.get_dummies(dummies_df[INSURANCE_COL], prefix='Insurance', drop_first=True),
    pd.get_dummies(dummies_df[MARITAL_STATUS_COL], prefix='Marital', drop_first=True),
    pd.get_dummies(dummies_df[RACE_COL], prefix='Ethnicity', drop_first=True),
    pd.get_dummies(dummies_df[ADMISSION_COUNT_GROUP_COL], prefix='AdmsCount', drop_first=True),
    dummies_df[[NIGHT_ADMISSION_FLAG, 
                GENDER_COL, 
                DIRECT_IND_COL,
                PREV_ADMISSION_IND_COL]].astype(int),
    dummies_df[STANDARDIZED_AGE_COL],
    dummies_df[LOS_DAYS_COL].dt.days,
    dummies_df[DISCHARGE_LOCATION_COL].dropna().replace(J_DICT).astype(int)
], axis=1)

fitters_table.dropna(inplace=True)
fitters_table = fitters_table[fitters_table.index.isin(table1.index)]

fitters_table.reset_index(inplace=True)
fitters_table.rename({DISCHARGE_LOCATION_COL: 'J', LOS_DAYS_COL: 'X', SUBJECT_ID_COL: 'pid'}, inplace=True, axis=1)
fitters_table.rename(RENAME_ITEMS_DICT, inplace=True, axis=1)

fitters_table = fitters_table[fitters_table['X'] > 0]
fitters_table.loc[fitters_table.X > ADMINISTRATIVE_CENSORING, 'J'] = 0
fitters_table.loc[fitters_table.X > ADMINISTRATIVE_CENSORING, 'X'] = ADMINISTRATIVE_CENSORING + 1
fitters_table['J'] = fitters_table['J'].astype(int)

fitters_table.to_csv(path_or_buf="mimic.csv", index=False)
print(len(fitters_table))
