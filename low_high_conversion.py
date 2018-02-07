#!/usr/bin/python
import pandas as pd
import os
import re
import math


def is_low_allele_valid(allele):
    return True


def is_high_allele_valid(allele):
    return True


# Extract two High HLAs from WHO assigned type value (e.g. A23(9) converts to A23 and A9 )
def get_parenthesis_tuple(allele_str):
    first_allele = allele_str[:allele_str.find("(")]
    inside_paren = allele_str[allele_str.find("(") + 1:allele_str.find(")")]
    hla_type = ''.join(re.findall(r'[A-Za-z]+', first_allele))
    second_allele = hla_type + inside_paren
    return first_allele, second_allele


# Mapping HLA allele containing parentheses(broad antigens) to WHO assigned type
def make_allele_mapping_dict(dict_file):
    dictionary_reference = pd.read_csv(dict_file)
    mapping_dict = {}

    for i, w in dictionary_reference.iterrows():
        allele_sets = [allele.strip() for allele in w[2].split('/')]
        for allele_set in allele_sets:
            if "(" in allele_set:
                first_allele, second_allele = get_parenthesis_tuple(allele_set)
                if second_allele not in mapping_dict:
                    mapping_dict[second_allele] = set()
                mapping_dict[second_allele].add(first_allele)
    return mapping_dict


# Creating main HLA Dictionary according to Expert assigned type
def make_allele_conversion_dict(dict_file):
    ret = {}
    dictionary_reference = pd.read_csv(dict_file)

    for index, row in dictionary_reference.iterrows():
        # Create an item for all HLA types with 'blank' prefix in dictionary for mapping of low HLAs with no value
        # (e.g blankA, blankB, blankDR)
        blank_modified_allele = 'blank' + (row[0].split('*'))[0]
        if blank_modified_allele not in ret:
            ret[blank_modified_allele] = []
        ret[blank_modified_allele].append(row[0])

        # If two low resolution HLAs have the same high resolution equivalent which is shown with '/' in HLA Dictionary
        if '/' in row[1]:
            alleles = [stripped_allele.strip() for stripped_allele in row[1].split('/')]
        else:
            alleles = [row[1].strip()]

        for allele in alleles:
            if is_low_allele_valid(allele) and is_high_allele_valid(allele):
                modified_allele = allele
                if modified_allele.lower() == 'blank':
                    modified_allele = modified_allele.lower() + row[0][0]
                if modified_allele not in ret:
                    ret[modified_allele] = []
                ret[modified_allele].append(row[0])
    return ret


# Creating a dictionary which maps A~B~DR haplotype to frequency
def make_freq_dic_of_race(dict_file):
    ret = {}
    dictionary_reference = pd.read_csv(dict_file)

    for index, row in dictionary_reference.iterrows():
        if row[0] not in ret:
            ret[row[0]] = {}
        if row[2] not in ret[row[0]]:
            ret[row[0]][row[2]] = {}
        if row[4] not in ret[row[0]][row[2]]:
            ret[row[0]][row[2]][row[4]]= row[6]
        else:
            # Summing across the various C and DRB3/4/5 types because the focus is on A~B~DR
            sum_freq = ret[row[0]][row[2]][row[4]] + row[6]
            ret[row[0]][row[2]][row[4]] = sum_freq

    return ret


# Making frequency dictionary(NMDP) for different races
def make_haplo_frequency_dict(dict_folder):
    freq_races_dic = {}
    for filename in os.listdir(dict_folder):
        if filename.endswith(".csv"):
            if filename not in freq_races_dic:
                freq_races_dic[filename] = make_freq_dic_of_race(dict_folder + "/" + filename)
    return freq_races_dic


# Translating races from our dataset to names using in NMDP dataset
def get_race(raw_race):
    race_dict = {"Black or African American": "AFA", "Arab or Middle Eastern": "CAU", "Hispanic/Latino": "HIS",
                 "American Indian or Alaska Native": "NAM", "Asian": "API",
                 "Native Hawaiian or Other Pacific Islander": "API", "Indian Sub-continent": "API", "White": "NAM"}
    if raw_race in race_dict:
        return race_dict[raw_race]
    else:
        return None


# Replacing Null value HLAs with blank+type
def replace_blank_in_haplo(haplo):
    new_haplo = []
    blank_addition = ['A', 'B', 'DRB1']

    for i in range(len(haplo)):
        if not haplo[i].strip():
            new_haplo.append('blank' + blank_addition[i])
        else:
            new_haplo.append(haplo[i])
    return new_haplo


# Getting high resolution equivalents of one HLA from hla dictionary OR mapping dict for special cases
def get_high_hla_list(low_hla, mapping_dict, hla_dict):
    high_list = []
    if low_hla in mapping_dict:
        mapped_list = []
        for mapped_hla in list(mapping_dict[low_hla]) + [low_hla]:
            if mapped_hla in hla_dict:
                mapped_list.append(mapped_hla)
    else:
        mapped_list = [low_hla] if low_hla in hla_dict else []
    for low_hla in mapped_list:
        high_list.extend(hla_dict[low_hla])
    return list(set(high_list))


# Getting maximum frequency among all possible equivalent of one haplotype
def max_haplo_freq(one_haplo, mapping_dict, hla_dic, freq_dic, race):
    possib = []
    one_haplo = replace_blank_in_haplo(one_haplo)

    # Getting all records containing frequency for this haplotype from nested NMDP dictionary
    for a in get_high_hla_list(one_haplo[0], mapping_dict, hla_dic):
        if a not in freq_dic[race + ".csv"]:
            continue
        for b in get_high_hla_list(one_haplo[1], mapping_dict, hla_dic):
            if b not in freq_dic[race + ".csv"][a]:
                continue
            for dr in get_high_hla_list(one_haplo[2], mapping_dict, hla_dic):
                if dr not in freq_dic[race + ".csv"][a][b]:
                    continue
                possib.append(([a,b,dr], freq_dic[race + ".csv"][a][b][dr]))
    if not possib:
        return None

    # Returning the record with the highest frequency
    return max(possib, key=lambda arr: arr[1])


# Replacing the Null value HLA with another one too, if one of the two HLAs of same type does not have value
def compare_haplo_blanks(two_haplo):
    for i in range(0, 5, 2):
        if two_haplo[i] is ' ' and two_haplo[i + 1] is not ' ':
            two_haplo[i] = two_haplo[i + 1]
        elif two_haplo[i + 1] is ' ' and two_haplo[i] is not ' ':
            two_haplo[i + 1] = two_haplo[i]
    return two_haplo


# Converting a low resolution dataframe to high resolution format
def get_high_resolution(my_dataset, result_columns, allele_mapping_dict, allele_conversion_dict, haplo_freq_dict):
    return_list = []
    column_names = ['A', 'A', 'B', 'B', 'DR', 'DR']

    for index, row in my_dataset.iterrows():
        if index%100==0:
            print (index)

        # Making HLAs well-formed
        haplos = []
        for i, allele in enumerate(row[2:8]):
            if not math.isnan(allele):
                haplos.append(str(column_names[i]) + str(int(allele)))
            else:
                haplos.append("")

        # Returning empty record, if no mapping found in dictionary for none of the HLAs in haplotype
        if any(hapl and hapl not in allele_conversion_dict and hapl not in allele_mapping_dict for hapl in haplos):
            return_list.append([row[0], row[1]] + [None] * 6)
            continue

        # Separating two haplotypes
        two_haplos = []
        two_haplos.append(compare_haplo_blanks(haplos)[::2])
        two_haplos.append(compare_haplo_blanks(haplos)[1::2])
        race = get_race(row[1])

        # Returning empty record, if race does not exist in or categories
        if race is None:
            return_list.append([row[0], row[1]] + [None] * 6)
            continue

        candidates = []
        for i1 in range(2):
            for i2 in range(2):
                # Generating all combinations of A, B and DR
                # A
                first_haplo = [two_haplos[0][0]]
                second_haplo = [two_haplos[1][0]]
                # B
                first_haplo.append(two_haplos[i1][1])
                second_haplo.append(two_haplos[(i1 + 1) % 2][1])
                # DR
                first_haplo.append(two_haplos[i2][2])
                second_haplo.append(two_haplos[(i2 + 1) % 2][2])

                # Finding maximum frequency for each haplotype
                max_haplo_tuple_first = max_haplo_freq(first_haplo, allele_mapping_dict, allele_conversion_dict, haplo_freq_dict, race)
                max_haplo_tuple_second = max_haplo_freq(second_haplo, allele_mapping_dict, allele_conversion_dict, haplo_freq_dict, race)

                # Getting total frequency of a haplotype pair
                if max_haplo_tuple_first and max_haplo_tuple_second:
                    if max_haplo_tuple_first[0] == max_haplo_tuple_second[0]:
                        candidates.append(
                            (
                                        max_haplo_tuple_first[0],
                                        max_haplo_tuple_second[0],
                                        max_haplo_tuple_first[1] * max_haplo_tuple_second[1]
                            )
                        )
                    else:
                        candidates.append(
                            (
                                        max_haplo_tuple_first[0],
                                        max_haplo_tuple_second[0],
                                        2 * max_haplo_tuple_first[1] * max_haplo_tuple_second[1]
                            )
                        )

        # Returning the found high resolution haplotype with highest frequency
        if len(candidates) == 0:
            return_list.append([row[0], row[1]] + [None] * 6)
        else:
            best_candidate = list(max(candidates, key=lambda arr: arr[2]))
            return_list.append([row[0], row[1], best_candidate[0][0], best_candidate[1][0], best_candidate[0][1]
                                   , best_candidate[1][1], best_candidate[0][2], best_candidate[1][2]])
    return pd.DataFrame(return_list, columns=result_columns)


def main(input_filepath):

    donor_type = input_filepath.split('/')[1]
    # Making Dicitonaries
    allele_mapping = make_allele_mapping_dict('Dictionaries/HLAConversionDataDictionary.csv')
    allele_dict = make_allele_conversion_dict('Dictionaries/HLAConversionDataDictionary.csv')
    freq_dict = make_haplo_frequency_dict('Dictionaries/freq_dic')
    print("finished with dictionaries")

    input_data_frame = pd.read_csv(input_filepath, encoding='latin1')

    # Donor Conversion
    donor_columns = ['PERS_ID', 'DON_RACE', 'DON_A1', 'DON_A2', 'DON_B1', 'DON_B2', 'DON_DR1', 'DON_DR2']
    donor_df = get_high_resolution(input_data_frame.iloc[:, 0:8], donor_columns, allele_mapping, allele_dict,
                                   freq_dict)
    # Recipient Conversion
    recipient_columns = ['PERS_ID', 'REC_RACE', 'REC_A1', 'REC_A2', 'REC_B1', 'REC_B2', 'REC_DR1', 'REC_DR2']
    recipient_df = get_high_resolution(input_data_frame.iloc[:, 8:16], recipient_columns, allele_mapping, allele_dict,
                                       freq_dict)

    result_df = pd.concat([donor_df, recipient_df], axis=1)
    result_df.to_csv("DataFiles/"+ donor_type + "/high_hla_types_dataset.csv", index=False, encoding='latin1')
