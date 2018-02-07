#!/usr/bin/python
import shutil

import pandas as pd
import numpy as np
import os
import glob
import low_high_conversion
import kosmoliaptsis
import sqlalchemy
from sqlalchemy import create_engine


def convert_to_sql():
    csv_database = create_engine('sqlite:///csv_database.db')

    for filename in glob.glob("/CsvFiles/*.csv"):
        df = pd.read_csv(filename)
        df = df.rename(columns={c: c.replace(' ', '') for c in df.columns})
        df.to_sql('table', csv_database, if_exists='append')


# Convert SAS to Csv file
def convert_to_csv(folder_path):

    print("Begin: convert_to_csv on ", folder_path)

    if folder_path != '':
        folder_path = folder_path + "/"
    for filename in glob.glob(folder_path + '*.sas7bdat'):
        df = pd.read_sas(filename, encoding='latin1')
        base = os.path.splitext(os.path.basename(filename))[0]
        if not os.path.exists("CsvFiles"):
            os.makedirs("CsvFiles")

        df.to_csv("CsvFiles/" + base + ".csv", index=False, encoding="latin1")

    print("End: convert_to_csv")


# Merge two tables with a merge key
def merge(first_table, second_table, merge_key):

    print("Begin: Merging ", first_table, " with ", second_table, " on ", merge_key)

    ftable = pd.read_csv(first_table, encoding='latin1')
    stable = pd.read_csv(second_table, encoding='latin1')
    merged = ftable.merge(stable, on=merge_key, suffixes=('', '__y'))

    # Drop Overlapping columns
    colsToDrop = [col for col in merged.columns if col[-3:] == "__y"]
    print("Dropping columns: ",  end='')
    for col in colsToDrop:
        print(col,  end='')

    print("")
    merged = merged.drop([col for col in merged.columns if col[-3:] == "__y"], axis=1)


    first_base = os.path.splitext(os.path.basename(first_table))[0]
    second_base = os.path.splitext(os.path.basename(second_table))[0]
    if not os.path.exists("DataFiles"):
        os.makedirs("DataFiles")

    oFile = "DataFiles/" + first_base + "_" + second_base + ".csv"
    merged.to_csv(oFile, index=False, encoding="latin1")

    print("Merged file at ", oFile)
    print("End: merge")


# Add two variables, days to graft failure date and total follow-up dates for prediction
def add_followup_dates(folder_path):

    print("Begin: add_followup_dates on ", folder_path)

    path = folder_path.split(".")[0]
    df = pd.read_csv(folder_path, encoding='latin1')
    df['DY_GRFAIL'] = (pd.to_datetime(df['TFL_GRAFT_DT']) - pd.to_datetime(df['REC_TX_DT'])).astype('timedelta64[D]')
    df['DY_TXFL'] = (pd.to_datetime(df['TFL_LAFUDATE']) - pd.to_datetime(df['REC_TX_DT'])).astype('timedelta64[D]')
    df.to_csv(path + "_with_variables.csv", index=False, encoding="latin1")

    # Drop patient graft failure with null value
    df2 = df[pd.notnull(df['TFL_GRAFT_DT'])]
    df2.to_csv(path+"dataset_with_valid_graft.csv", index=False, encoding="latin1")

    print("End: add_followup_dates")


# Create a dataframe of donor and receipent HLA information from main dataset
def get_hla_types(file_path):

    print("Begin: get_hla_types")

    donor_type = file_path.split("_")[-3].split(".")[0]

    df = pd.read_csv(file_path, encoding='latin1')
    df = df[['PERS_ID', 'DON_RACE', 'DON_A1', 'DON_A2', 'DON_B1', 'DON_B2', 'DON_DR1', 'DON_DR2', 'PERS_ID', 'CAN_RACE', 'REC_A1', 'REC_A2', 'REC_B1',
             'REC_B2', 'REC_DR1', 'REC_DR2']]

    # Drops records with all HLAs null
    df = df.dropna(subset=df.columns[2:8], how='all')

    # Decode race numbers to strings for further steps
    df['CAN_RACE'] = df['CAN_RACE'].map(
        {8: 'White', 16: 'Black or African American', 32: 'American Indian or Alaska Native', 64: 'Asian',
         128: 'Native Hawaiian or Other Pacific Islander', 256: 'Arab or Middle Eastern', 512: 'Indian Sub-continent',
         1024: 'Unknown (for Donor Referral only)', 2000: 'Hispanic/Latino'})
    df['DON_RACE'] = df['DON_RACE'].map(
        {8: 'White', 16: 'Black or African American', 32: 'American Indian or Alaska Native', 64: 'Asian',
         128: 'Native Hawaiian or Other Pacific Islander', 256: 'Arab or Middle Eastern', 512: 'Indian Sub-continent',
         1024: 'Unknown (for Donor Referral only)', 2000: 'Hipanic/Latino'})
    if not os.path.exists("DataFiles/donor_" + donor_type):
        os.makedirs("DataFiles/donor_" + donor_type)
    df.to_csv("DataFiles/donor_" + donor_type + "/low_hla_types_dataset.csv", index=False, encoding="latin1")


# Make DR file compatible with Kosmoliaptsis dictionary by removing ':' from HLAs
def format_file(file_name):
    lines = []
    with open(file_name, 'r') as input:
        lines = input.readlines()
    conversion = ':'
    new_text = ''
    output_lines = []
    for line in lines:
        temp = line[:]
        for c in conversion:
            temp = temp.replace(c, new_text)
        output_lines.append(temp)

    with open(file_name, 'w') as output:
        for line in output_lines:
            output.write(line + "\n")

    print("End: get_hla_types")


# Generate required input files for immunogenicity calculation, A-B for class I and DR for class II
def generate_kosmoliaptsis_input(input_filepath):

    print("Begin:generate_kosmoliaptsis_input")

    donor_type = input_filepath.split('/')[1]
    path = "DataFiles/" + donor_type + "/"
    if os.path.isfile(path + "A-B.csv"):
        os.remove(path + "A-B.csv")
    if os.path.isfile(path + "DR.csv"):
        os.remove(path + "DR.csv")
    # AB
    csv1 = pd.read_csv(input_filepath)
    csv2 = pd.read_csv(input_filepath)
    df1 = pd.DataFrame(csv1, columns=['PERS_ID', 'DON_A1', 'DON_A2', 'DON_B1', 'DON_B2'])
    df2 = pd.DataFrame(csv2, columns=['PERS_ID', 'REC_A1', 'REC_A2', 'REC_B1', 'REC_B2', 'REC_Cw1','REC_Cw2'])
    for i, c in enumerate(df1):
        if (i > 0) & (i < 5):
            df1 = pd.DataFrame(csv1, columns=[c, 'PERS_ID'])
            df1.merge(df2, on="PERS_ID").to_csv(path+"A-B.csv", mode='a', index=False, header=False, encoding="latin1")

    df_ab = pd.read_csv(path+"A-B.csv", header=None)
    df_ab.columns = ['Donor_Allele', 'ID', 'Recipient_HLA-A-1', 'Recipient_HLA-A-2', 'Recipient_HLA-B-1','Recipient_HLA-B-2','Recipient_HLA-C-1','Recipient_HLA-C-2']
    df_ab['Recipient_HLA-C-1'] = 'xxx'
    df_ab['Recipient_HLA-C-2'] = 'xxx'

    # Drop records with null value
    df_ab = df_ab.dropna(subset=['Donor_Allele'])
    df_ab = df_ab.dropna(subset=df_ab.columns[2:5], how='all')
    df_ab.to_csv(path + "A-B.csv", index=False, encoding="latin1")

    # DR
    csv3 = pd.read_csv(input_filepath)
    csv4 = pd.read_csv(input_filepath)
    df3 = pd.DataFrame(csv3, columns=['PERS_ID', 'DON_DR1', 'DON_DR2'])
    df4 = pd.DataFrame(csv4, columns=['PERS_ID', 'REC_DR1', 'REC_DR2'])
    for i, c in enumerate(df3):
        if (i > 0) & (i < 3):
            df3 = pd.DataFrame(csv3, columns=[c, 'PERS_ID'])
            df3.merge(df4, on="PERS_ID").to_csv(path+"DR.csv", mode='a', index=False, header=False, encoding="latin1")

    df_dr = pd.read_csv(path+"DR.csv", header=None)
    cols = ['Donor-A-Chain','Donor-B-Chain', 'ID', 'RecDRA1-1', 'RecDRA1-2', 'RecDRB1-1', 'RecDRB1-2','RecDRB345-1', 'RecDRB345-2', 'RecDQA1-1', 'RecDQA1-2', 'RecDQB1-1', 'RecDQB1-2','RecDPA1-1',
                     'RecDPA1-2', 'RecDPB1-1', 'RecDPB1-2']

    # Make the dataframe well-formed for kosmoliaptsis algorithm
    new_df = pd.DataFrame(columns=cols)
    new_df['Donor-B-Chain'] = df_dr[0]
    new_df['ID'] = df_dr[1]
    new_df['RecDRB1-1'] = df_dr[2]
    new_df['RecDRB1-2'] = df_dr[3]
    new_df['Donor-A-Chain'] = 'xxx'
    new_df['RecDRA1-1'] = 'xxx'
    new_df['RecDRA1-2'] = 'xxx'
    new_df['RecDRB345-1'] = 'xxx'
    new_df['RecDRB345-2'] = 'xxx'
    new_df['RecDQA1-1'] = 'xxx'
    new_df['RecDQA1-2'] = 'xxx'
    new_df['RecDQB1-1'] = 'xxx'
    new_df['RecDQB1-2'] = 'xxx'
    new_df['RecDPA1-1'] = 'xxx'
    new_df['RecDPA1-2'] = 'xxx'
    new_df['RecDPB1-1'] = 'xxx'
    new_df['RecDPB1-2'] = 'xxx'

    # Drop records with null value
    new_df = new_df.dropna(subset=['Donor-B-Chain'])
    new_df = new_df.dropna(subset=new_df.columns[5:6], how='all')

    new_df.to_csv(path + "DR.csv", index=False, encoding="latin1")
    format_file(path + "DR.csv")


# Calculate total immunogenicity from class I and class II results
def calculate_total_immunogenicity(input_file_classi, input_file_classii):
    donor_type = input_file_classi.split('/')[1]
    path = "DataFiles/"+donor_type+"/"
    # Class I
    df_class_i = pd.read_csv(input_file_classi)
    df_class_i = df_class_i.sort_values(['Donor_ID'], ascending=[True])
    df_class_i = df_class_i.groupby(['Donor_ID']).mean()
    # ClassII
    df_class_ii = pd.read_csv(input_file_classii)
    df_class_ii = df_class_ii.sort_values(['Donor_ID'], ascending=[True])
    df_class_ii = df_class_ii.groupby(['Donor_ID']).mean()

    # Rename class I and II variables
    df_class_i = df_class_i.rename(columns={'EMS': 'AB_EMS', 'HMS': 'AB_HMS', 'Number_AA_MM': 'AB_AMS'})
    df_class_ii = df_class_ii[['EMS_beta_chain', 'HMS_beta_chain', 'Number_AA_MM_beta']]
    df_class_ii = df_class_ii.rename(columns={'EMS_beta_chain': 'DR_EMS', 'HMS_beta_chain': 'DR_HMS', 'Number_AA_MM_beta': 'DR_AMS'})

    result_df = pd.concat([df_class_i, df_class_ii], axis=1)
    result_df['Avg_EMS'] = result_df[['AB_EMS', 'DR_EMS']].mean(axis=1)
    result_df['Avg_HMS'] = result_df[['AB_HMS', 'DR_HMS']].mean(axis=1)
    result_df['Avg_AMS'] = result_df[['AB_AMS', 'DR_AMS']].mean(axis=1)
    result_df.index.name = 'PERS_ID'
    result_df.to_csv(path + "Immunogenicity_dataset.csv", encoding="latin1")

    print("End: get_total_immunogenicity")


# Decode categorical variables and dates
def get_transformed_data(file_path):

    print("Begin: data transformation on ", file_path)

    donor_type = file_path.split('/')[1]
    path = "DataFiles/" + donor_type + "/"
    df = pd.read_csv(file_path, encoding='latin1')

    # Rename or Drop special variables
    df = df.rename(
        columns={'PERS_NEXTTX': 'PERS_NEXTTX_DT', 'PERS_RELIST': 'PERS_RELIST_DT', 'PERS_RETX': 'PERS_RETX_DT',
                 'TFL_ENDTXFU': 'TFL_ENDTXFU_DT', 'TFL_LAFUDATE': 'TFL_LAFUDATE_DT'})
    df = df.drop(['CAN_DGN_OSTXT', 'REC_DGN_OSTXT', 'DON_CANCER_OTHER_OSTXT'], axis=1)

    # Transform Date variables
    for col in df.columns:
        if col[-3:] == "_DT":
            df[col] = df[col].fillna(0)
            df[col] = pd.to_datetime(df[col])
            df[col] = df[col].apply(lambda x: x.toordinal())

    # One-Hot-Encoding for two-value categorical variables
    char_cols = df.dtypes.pipe(lambda x: x[x == 'object']).index
    label_mapping = {}
    for c in char_cols:
        if len(df[c].value_counts()) < 3:
            df[c], label_mapping[c] = pd.factorize(df[c])

    # Dummy coding for more than two-value categorical variables
    df = pd.get_dummies(df)
    df = df.fillna(-1)
    df.to_csv(path+"FINAL_DONOR_DECEASED_DATASET_ENCODING.csv")

    print("End: get_transformed_data")


def generate_donor_deceased_data():

    merge("DataFiles/cand_kipa_tx_ki.csv", "CSVFiles/donor_deceased.csv", "DONOR_ID")

    # Add prediction variables
    add_followup_dates('DataFiles/cand_kipa_tx_ki_donor_deceased.csv')

    # Extract only HLA type variables of donor and recipient
    get_hla_types("DataFiles/cand_kipa_tx_ki_donor_deceased_with_variables.csv")
    low_high_conversion.main("DataFiles/donor_deceased/low_hla_types_dataset.csv")

    # Prepare input files for Kosmoliaptsis algorithm
    generate_kosmoliaptsis_input("DataFiles/donor_deceased/high_hla_types_dataset.csv")

    # Calculate Immunogenicity for class I, II and finally their total
    kosmoliaptsis.select_kosmo_class_type(1, "DataFiles/donor_deceased/A-B.csv")
    kosmoliaptsis.select_kosmo_class_type(2, "DataFiles/donor_deceased/DR.csv")
    calculate_total_immunogenicity("DataFiles/donor_deceased/Immunogenicity_I.csv","DataFiles/donor_deceased/Immunogenicity_II.csv")

    # Concatenate Immunogenicity data with main dataset
    merge("DataFiles/cand_kipa_tx_ki_donor_deceased_with_variables.csv", "DataFiles/donor_deceased/Immunogenicity_dataset.csv", "PERS_ID")
    shutil.move("DataFiles/cand_kipa_tx_ki_donor_deceased_with_variables_Immunogenicity_dataset.csv", "DataFiles/donor_deceased/FINAL_DONOR_DECEASED_DATASET.csv")


if __name__ == "__main__":

    # Convert SAS data to CSV
    convert_to_csv("SASFiles")

    # Make main dataset for donor deceased
    merge("CsvFiles/cand_kipa.csv", "CSVFiles/tx_ki.csv", "PX_ID")
    generate_donor_deceased_data()

    # Data Transformation
    get_transformed_data("DataFiles/donor_deceased/FINAL_DONOR_DECEASED_DATASET.csv")

