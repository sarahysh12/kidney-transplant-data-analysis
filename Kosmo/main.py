import pandas as pd
import os
import glob
import kosmoliaptsis

def generate_kosmoliaptsis_input(input_filepath):
    if os.path.isfile("DataFiles/A-B.csv"):
        os.remove("DataFiles/A-B.csv")
    if os.path.isfile("DataFiles/DR.csv"):
        os.remove("DataFiles/DR.csv")
    # AB
    csv1 = pd.read_csv(input_filepath)
    csv2 = pd.read_csv(input_filepath)
    df1 = pd.DataFrame(csv1, columns=['PERS_ID', 'DON_A1', 'DON_A2', 'DON_B1', 'DON_B2'])
    df2 = pd.DataFrame(csv2, columns=['PERS_ID', 'REC_A1', 'REC_A2', 'REC_B1', 'REC_B2'])
    for i, c in enumerate(df1):
        if (i > 0) & (i < 5):
            df1 = pd.DataFrame(csv1, columns=[c, 'PERS_ID'])
            df1.merge(df2, on="PERS_ID").to_csv("DataFiles/A-B.csv", mode='a', index=False, header=False)

    df_ab = pd.read_csv("DataFiles/A-B.csv", header=None)
    ab_cols = ['Donor_Allele', 'ID', 'Recipient_HLA-A-1', 'Recipient_HLA-A-2', 'Recipient_HLA-B-1','Recipient_HLA-B-2', 'Recipient_HLA-C-1', 'Recipient_HLA-C-2']
    new_ab_df = pd.DataFrame(columns=ab_cols)
    new_ab_df['Donor_Allele'] = df_ab[0]
    new_ab_df['ID'] = df_ab[1]
    new_ab_df['Recipient_HLA-A-1'] = df_ab[2]
    new_ab_df['Recipient_HLA-A-2'] = df_ab[3]
    new_ab_df['Recipient_HLA-B-1'] = df_ab[4]
    new_ab_df['Recipient_HLA-B-2'] = df_ab[5]
    new_ab_df['Recipient_HLA-C-1'] = 'xxx'
    new_ab_df['Recipient_HLA-C-2'] = 'xxx'
    new_ab_df = new_ab_df.dropna(subset=['Donor_Allele'])
    new_ab_df = new_ab_df.dropna(subset=new_ab_df.columns[2:6], how='all')
    new_ab_df.to_csv("DataFiles/A-B.csv", index=False)
    print("finish AB")
    # DRDQ
    csv3 = pd.read_csv(input_filepath)
    csv4 = pd.read_csv(input_filepath)
    df3 = pd.DataFrame(csv3, columns=['PERS_ID', 'DON_DR1', 'DON_DR2'])
    df4 = pd.DataFrame(csv4, columns=['PERS_ID', 'REC_DR1', 'REC_DR2'])
    for i, c in enumerate(df3):
        if (i > 0) & (i < 3):
            df3 = pd.DataFrame(csv3, columns=[c, 'PERS_ID'])
            df3.merge(df4, on="PERS_ID").to_csv("DataFiles/DR.csv", mode='a', index=False, header=False)

    df_dr = pd.read_csv("DataFiles/DR.csv", header=None)
    cols = ['Donor-A-Chain','Donor-B-Chain', 'ID', 'RecDRA1-1', 'RecDRA1-2', 'RecDRB1-1', 'RecDRB1-2','RecDRB345-1', 'RecDRB345-2', 'RecDQA1-1', 'RecDQA1-2', 'RecDQB1-1', 'RecDQB1-2','RecDPA1-1',
                     'RecDPA1-2', 'RecDPB1-1', 'RecDPB1-2']

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
    new_df = new_df.dropna(subset=['Donor-B-Chain'])
    new_df = new_df.dropna(subset=new_df.columns[5:6], how='all')
    print("finish Dr")
    new_df.to_csv("DataFiles/DR.csv", index=False)


def get_total_immunogenicity(input_file_classi, input_file_classii):
    #Class I
    df_class_i = pd.read_csv(input_file_classi)
    df_class_i = df_class_i.sort_values(['Donor_ID'], ascending=[True])
    df_class_i = df_class_i.groupby(['Donor_ID']).mean()
    #ClassII
    df_class_ii = pd.read_csv(input_file_classii)
    df_class_ii = df_class_ii.sort_values(['Donor_ID'], ascending=[True])
    df_class_ii = df_class_ii.groupby(['Donor_ID']).mean()

    df_class_i = df_class_i.rename(columns={'EMS': 'AB_EMS', 'HMS': 'AB_HMS', 'Number_AA_MM': 'AB_AMS'})
    df_class_ii = df_class_ii[['EMS_beta_chain', 'HMS_beta_chain', 'Number_AA_MM_beta']]
    df_class_ii = df_class_ii.rename(columns={'EMS_beta_chain': 'DR_EMS', 'HMS_beta_chain': 'DR_HMS', 'Number_AA_MM_beta': 'DR_AMS'})

    result_df = pd.concat([df_class_i, df_class_ii], axis=1)
    result_df['Avg_EMS'] = result_df[['AB_EMS', 'DR_EMS']].mean(axis=1)
    result_df['Avg_HMS'] = result_df[['AB_HMS', 'DR_HMS']].mean(axis=1)
    result_df['Avg_AMS'] = result_df[['AB_AMS', 'DR_AMS']].mean(axis=1)
    result_df.index.name = 'PERS_ID'
    result_df.to_csv("DataFiles/Immunogenicity_dataset.csv")

def format_file(file_name):
    lines = []
    with open('DataFiles/DR.csv', 'r') as input:
        lines = input.readlines()
    conversion = ':'
    newtext = ''
    outputLines = []
    for line in lines:
        temp = line[:]
        for c in conversion:
            temp = temp.replace(c, newtext)
        outputLines.append(temp)

    with open('DataFiles/DR.csv', 'w') as output:
        for line in outputLines:
            output.write(line + "\n")



if __name__ == "__main__":
    generate_kosmoliaptsis_input("DataFiles/high_hla_types_dataset.csv")
    format_file("DataFiles/DR.csv")
    kosmoliaptsis.select_kosmo_class_type(1, "DataFiles/A-B.csv")
    kosmoliaptsis.select_kosmo_class_type(2, "DataFiles/DR.csv")
    get_total_immunogenicity("DataFiles/Immunogenicity_I.csv","DataFiles/Immunogenicity_II.csv")
    