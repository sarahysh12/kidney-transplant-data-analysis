#!/usr/bin/python
import pandas as pd


def select_kosmo_class_type(class_type, input_filepath):

    donor_type = input_filepath.split('/')[1]
    path = "DataFiles/"+donor_type+"/"

    # Make Amino Acids dictionary using by both classes
    amino_acids = pd.read_csv('Dictionaries/AA_Code.csv')
    amino_acids_dic = make_amino_acid_dict(amino_acids)

    # Decide calculating immunogenicity of class I or II
    if class_type == 1:
        # Make specified dictionary and calculate class one's immunogenicity
        aligned_seq = pd.read_csv('Dictionaries/Aligned_Seq.csv')
        aligned_sequence_dic = make_aligned_seq_dict(aligned_seq)
        input_allele_c1 = pd.read_csv(input_filepath)
        dataframe_result_c1 = calculate_class_one(amino_acids_dic, aligned_sequence_dic, input_allele_c1)
        dataframe_result_c1.to_csv(path+"Immunogenicity_I.csv", index=False)
    else:
        # Make specified dictionary and calculate class two's immunogenicity
        input_allele_c2 = pd.read_csv(input_filepath)
        all_sequence_dict = make_allele_sequence_dict(['DPA1', 'DPB1', 'DQA1', 'DRB1345', 'DQB1'])
        dataframe_result_c2 = calculate_class_two(amino_acids_dic, all_sequence_dict, input_allele_c2)
        dataframe_result_c2.to_csv(path+"Immunogenicity_II.csv", index=False)


def calculate_class_two(amino_acids_dic, all_sequence_dict, input_allele):

    print ("start calculating class II")
    df = pd.DataFrame()

    # Calculate HMS, EMS and Number_AA_MM (alpha and beta) for each donor allele by comparing AAs of donor and AAs of two alleles in recipient
    for index, allele in input_allele.iterrows():
        if index % 100 == 0:
            print (index)

        result = []

        # Compare each AA of donor with other two (DR-1, DR-2) AAs of recipient in the same index
        for don in allele[0:2]:

            hms, ems = 0, 0
            mismatch_result = [] * 3

            if don != 'xxx':

                # Specify type of chain (We have only DRB1 in our dataset, so we won't calculate alpha chain at all)
                if "A" in don:
                    donor_type = 1  # alpha_chain
                else:
                    donor_type = 2  # beta_chain

                # Get the sequence dictionary
                allele_dic = get_allele_dict(all_sequence_dict, don.split("*")[0])

                # Convert donor allele to sequence
                donor_allele_seq = get_amino_acid_details(don, amino_acids_dic, allele_dic[0])

                if donor_allele_seq is None:
                    continue
                else:
                    for i, donor_aa in enumerate(donor_allele_seq):

                        not_similar = []
                        # Limitations for alpha and beta (number of alpha/beta allele types)
                        alpha_allele = 6
                        beta_allele = 7

                        # Compare each AA of donor with other two AAs of recipient in the same index
                        for column_index, a_rec in enumerate(allele[3:16]):
                            if a_rec == 'xxx' or a_rec == '':
                                if column_index in [0, 1, 6, 7, 10, 11]:  # alpha alleles
                                    alpha_allele = alpha_allele - 1
                                else:
                                    beta_allele = beta_allele - 1
                                continue
                            else:
                                # Get sequence dictionary and convert recipient allele to AA
                                rec_allele_dic = get_allele_dict(all_sequence_dict, a_rec.split("*")[0])
                                receiver_allele = get_amino_acid_details(a_rec, amino_acids_dic, rec_allele_dic[0])

                                if receiver_allele is None:
                                    continue
                                else:
                                    # Keep the differences for each comparison between one donor's AA and recipient's AAs
                                    if (i <= len(receiver_allele) - 1) and (donor_aa != receiver_allele[i]):
                                        subtraction = tuple(
                                            map(lambda x, y: round(abs(x - y), 2), donor_aa, receiver_allele[i]))
                                        if ("B" not in a_rec) and (donor_type == 1):
                                            not_similar.append(subtraction)
                                        elif ("A" not in a_rec) and (donor_type == 2):
                                            not_similar.append(subtraction)

                        # Store the lowest difference, if all two AAs of recipient are different from donor's AA
                        if donor_type == 1:
                            if len(not_similar) == alpha_allele:
                                mismatch_result.append((min(not_similar, key=lambda minHB: minHB[0])[0],
                                                        min(not_similar, key=lambda minPI: minPI[1])[1]))
                        else:
                            if len(not_similar) == beta_allele:
                                mismatch_result.append((min(not_similar, key=lambda minHB: minHB[0])[0],
                                                        min(not_similar, key=lambda minPI: minPI[1])[1]))

                    # Calculate HMS and EMS, if we even have only one difference in donor and recipient AAs
                    if len(mismatch_result) > 0:
                        for mm in mismatch_result:
                            hms += mm[0]
                            ems += mm[1]
                    else:
                        hms, ems = 0, 0

            result.append((hms, ems, len(mismatch_result)))

        # Format the output file
        if len(result) == 2:
            df = df.append({'Donor_alpha': allele[0], 'Donor_beta': allele[1], 'Donor_ID': allele[2],
                            'EMS_alpha_chain': round(result[0][1], 2), 'EMS_beta_chain': round(result[1][1], 2),
                            'EMS_total': round((result[0][1] + result[1][1]), 2),
                            'HMS_alpha_chain': round(result[0][0], 2), 'HMS_beta_chain': round(result[1][0], 2),
                            'HMS_total': round((result[0][0] + result[1][0]), 2), 'Number_AA_MM_alpha': result[0][2],
                            'Number_AA_MM_beta': result[1][2], 'Number_AA_MM_total': (result[0][2] + result[1][2])},
                           ignore_index=True)

    return df


def calculate_class_one(amino_acids_dic, aligned_sequence_dic, input_allele):

    print ("start calculating class I")
    df = pd.DataFrame()

    # Calculate HMS, EMS and Number_AA_MM for each donor allele by comparing 276 AAs of donor and 276 AAs of four alleles in recipient
    for index, allele in input_allele.iterrows():
        if index % 100 == 0:
            print (index)

        mismatch_result = [] * 3

        # Convert donor allele to sequence of 276 AAs
        donor_allele = get_amino_acid_details(allele[0], amino_acids_dic, aligned_sequence_dic)
        if donor_allele is None:
            continue

        hms, ems = 0, 0

        for i, donor_aa in enumerate(donor_allele):

            not_similar = []
            all_allele_diff = 6

            # Compare each AA of donor (of 276) with other four AAs of recipient in the same index (out of 276)
            for a_rec in allele[2:8]:

                # Change the limitation according to the number of HLAs we have in dataset for recipient
                if a_rec == 'xxx' or a_rec == '':
                    all_allele_diff = all_allele_diff - 1
                    continue
                else:
                    # Convert recipient allele to sequence of 276 AAs
                    receiver_allele = get_amino_acid_details(a_rec, amino_acids_dic, aligned_sequence_dic)

                    if receiver_allele is None:
                        continue
                    elif donor_aa != receiver_allele[i]:
                        # Keep the differences for each comparison between one donor's AA and recipient's AAs (4 AA)
                        not_similar.append(tuple(map(lambda x, y: round(abs(x - y), 2), donor_aa, receiver_allele[i])))

            # Store the lowest difference in mismatch_result, if all four AAs of recipient are different from donor's AA
            if len(not_similar) == all_allele_diff:
                mismatch_result.append(
                    (min(not_similar, key=lambda minHB: minHB[0])[0], min(not_similar, key=lambda minPI: minPI[1])[1]))

        # Calculate HMS and EMS, if we even have only one difference (in mismatch_result) in donor and recipient AAs
        if len(mismatch_result) > 0:
            for mm in mismatch_result:
                hms += mm[0]
                ems += mm[1]
        df = df.append({'Donor': allele[0], 'Donor_ID': allele[1], 'HMS': round(hms, 2), 'EMS': round(ems, 2),
                        'Number_AA_MM': len(mismatch_result)}, ignore_index=True)
    return df


# Select the proper dictionary according to the type of allele
def get_allele_dict(all_sequence_dict, allele_type):
    if allele_type == "DRB1" or allele_type == "DRB3" or allele_type == "DRB4" or allele_type == "DRB5":
        return all_sequence_dict["DRB1345"]
    else:
        return all_sequence_dict[allele_type]


# Return a sequence of 276 Amino Acids(AA) for each HLA allele where each AA replaced with a pair of HB and PI
def get_amino_acid_details(allele, amino_acids, seq_dict):
    person_allele = []
    if allele not in seq_dict:
        return None
    else:
        for seq in seq_dict[allele]:
            for aa in seq:
                person_allele.append(((amino_acids[aa])[0], (amino_acids[aa])[1]))
        return person_allele


# Create different dictionaries for Class II which maps HLA alleles(e.g. DR1) to a sequence of Amino Acids
def make_allele_sequence_dict(files_name):
    all_dictionaries = {}
    for file_n in files_name:
        seq_dataframe = pd.read_csv("Dictionaries/Class_Two_Seq/" + file_n + '_seq.csv')
        seq_dic = {}
        for aseq in seq_dataframe.columns:
            if aseq not in seq_dic:
                seq_dic[aseq] = []
            seq_dic[aseq].append(seq_dataframe[aseq])
        if file_n not in all_dictionaries:
            all_dictionaries[file_n] = []
        all_dictionaries[file_n].append(seq_dic)
    return all_dictionaries


# Create a dictionary for Class I which maps one HLA allele (A/B/C) to a sequence of 276 Amino Acids
def make_aligned_seq_dict(aseq_dataframe):
    aligned_seq_dic = {}
    for aseq in aseq_dataframe.columns:
        if aseq not in aligned_seq_dic:
            aligned_seq_dic[aseq] = []
        aligned_seq_dic[aseq].append(aseq_dataframe[aseq])
    return aligned_seq_dic


# Create a dictionary for 21 amino acids which maps each amino acid to its HB and PI values
def make_amino_acid_dict(aa_dataframe):
    aa_dict = {}
    for index, aa in aa_dataframe.iterrows():
        if aa[1] not in aa_dict:
            aa_dict[aa[1]] = []
        aa_dict[aa[1]].append(aa[2])  # HB
        aa_dict[aa[1]].append(aa[3])  # PI
    return aa_dict