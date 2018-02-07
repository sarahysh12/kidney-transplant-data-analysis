import pandas as pd


def select_kosmo_class_type(class_type, input_filepath):
    amino_acids = pd.read_csv('Dictionaries/AA_Code.csv')
    amino_acids_dic = make_amino_acid_dict(amino_acids)

    if class_type == 1:
        aligned_seq = pd.read_csv('Dictionaries/Aligned_Seq.csv')
        aligned_sequence_dic = make_aligned_seq_dict(aligned_seq)
        input_allele_c1 = pd.read_csv(input_filepath)
        dataframe_result_c1 = calculate_class_one(amino_acids_dic, aligned_sequence_dic, input_allele_c1)
        dataframe_result_c1.to_csv("DataFiles/Immunogenicity_I.csv", index=False)
    else:
        input_allele_c2 = pd.read_csv(input_filepath)
        all_sequence_dict = make_allele_sequence_dict(['DPA1', 'DPB1', 'DQA1', 'DRB1345', 'DQB1'])
        dataframe_result_c2 = calculate_class_two(amino_acids_dic, all_sequence_dict, input_allele_c2)
        dataframe_result_c2.to_csv("DataFiles/Immunogenicity_II.csv", index=False)


def calculate_class_two(amino_acids_dic, all_sequence_dict, input_allele):
    df = pd.DataFrame()
    for index, allele in input_allele.iterrows():
        if index%1000==0:
            print(index)
        result = []
        for don in allele[0:2]:
            hms, ems = 0, 0
            mismatch_result = [] * 3
            if don != 'xxx':
                if "A" in don:
                    donor_type = 1  # alpha_chain
                else:
                    donor_type = 2  # beta_chain
                allele_dic = get_allele_dict(all_sequence_dict, don.split("*")[0])
                donor_allele_seq = get_amino_acid_details(don, amino_acids_dic, allele_dic[0])
                if donor_allele_seq is None:
                    continue
                else:
                    for i, donor_aa in enumerate(donor_allele_seq):
                        not_similar = []
                        alpha_allele = 6
                        beta_allele = 7
                        for column_index, a_rec in enumerate(allele[3:16]):
                            if a_rec == 'xxx' or a_rec == '':
                                if column_index in [0, 1, 6, 7, 10, 11]:  # alpha alleles
                                    alpha_allele = alpha_allele - 1
                                else:
                                    beta_allele = beta_allele - 1
                                continue
                            else:
                                rec_allele_dic = get_allele_dict(all_sequence_dict, a_rec.split("*")[0])
                                receiver_allele = get_amino_acid_details(a_rec, amino_acids_dic, rec_allele_dic[0])
                                if receiver_allele is None:
                                    continue
                                else:
                                    if (i <= len(receiver_allele) - 1) and (donor_aa != receiver_allele[i]):
                                        subtraction = tuple(
                                            map(lambda x, y: round(abs(x - y), 2), donor_aa, receiver_allele[i]))
                                        if ("B" not in a_rec) and (donor_type == 1):  # calculating alpha esm,hms
                                            not_similar.append(subtraction)
                                        elif ("A" not in a_rec) and (donor_type == 2):
                                            not_similar.append(subtraction)

                        if donor_type == 1:  # calculating alpha esm,hms
                            if len(not_similar) == alpha_allele:
                                mismatch_result.append((min(not_similar, key=lambda minHB: minHB[0])[0],
                                                        min(not_similar, key=lambda minPI: minPI[1])[1]))
                        else:
                            if len(not_similar) == beta_allele:
                                mismatch_result.append((min(not_similar, key=lambda minHB: minHB[0])[0],
                                                        min(not_similar, key=lambda minPI: minPI[1])[1]))

                    if len(mismatch_result) > 0:
                        for mm in mismatch_result:
                            hms += mm[0]
                            ems += mm[1]
                    else:
                        hms, ems = 0, 0

            result.append((hms, ems, len(mismatch_result)))
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
    df = pd.DataFrame()
    for index, allele in input_allele.iterrows():
        if index%1000==0:
            print(index)
        mismatch_result = [] * 3
        donor_allele = get_amino_acid_details(allele[0], amino_acids_dic, aligned_sequence_dic)
        if donor_allele is None:
            continue
        hms, ems = 0, 0
        for i, donor_aa in enumerate(donor_allele):

            not_similar = []
            all_allele_diff = 6
            for a_rec in allele[2:8]:
                if a_rec == 'xxx' or a_rec == '':
                    all_allele_diff = all_allele_diff - 1
                    continue
                else:
                    receiver_allele = get_amino_acid_details(a_rec, amino_acids_dic, aligned_sequence_dic)
                    if receiver_allele is None:
                        continue
                    elif donor_aa != receiver_allele[i]:
                        not_similar.append(tuple(map(lambda x, y: round(abs(x - y), 2), donor_aa, receiver_allele[i])))

            if len(not_similar) == all_allele_diff:
                mismatch_result.append(
                    (min(not_similar, key=lambda minHB: minHB[0])[0], min(not_similar, key=lambda minPI: minPI[1])[1]))
        if len(mismatch_result) > 0:
            for mm in mismatch_result:
                hms += mm[0]
                ems += mm[1]
        df = df.append({'Donor': allele[0], 'Donor_ID': allele[1], 'HMS': round(hms, 2), 'EMS': round(ems, 2),
                        'Number_AA_MM': len(mismatch_result)}, ignore_index=True)
    return df


def get_allele_dict(all_sequence_dict, allele_type):
    if allele_type == "DRB1" or allele_type == "DRB3" or allele_type == "DRB4" or allele_type == "DRB5":
        return all_sequence_dict["DRB1345"]
    else:
        return all_sequence_dict[allele_type]


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


def get_amino_acid_details(allele, amino_acids, seq_dict):
    person_allele = []
    if allele not in seq_dict:
        return None
    else:
        for seq in seq_dict[allele]:
            for aa in seq:
                person_allele.append(((amino_acids[aa])[0], (amino_acids[aa])[1]))
        return person_allele


def make_aligned_seq_dict(aseq_dataframe):
    aligned_seq_dic = {}
    for aseq in aseq_dataframe.columns:
        if aseq not in aligned_seq_dic:
            aligned_seq_dic[aseq] = []
        aligned_seq_dic[aseq].append(aseq_dataframe[aseq])
    return aligned_seq_dic


def make_amino_acid_dict(aa_dataframe):
    aa_dict = {}
    for index, aa in aa_dataframe.iterrows():
        if aa[1] not in aa_dict:
            aa_dict[aa[1]] = []
        aa_dict[aa[1]].append(aa[2])  # HB
        aa_dict[aa[1]].append(aa[3])  # PI
    return aa_dict