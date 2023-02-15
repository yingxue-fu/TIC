#!/usr/bin/python

# JUMP_isotope_ditsribution calculation
# Created by Surendhar Reddy Chepyala, Modified by Ji-Hoon Cho and Yingxue Fu

import os
import pickle, numpy as np, pandas as pd
from collections import defaultdict

# load preComputedIsotopes file
preComputedIsotopes = os.path.join(os.path.dirname(os.path.realpath(__file__)), "isotopeMassIntensity.pkl")

# Open the default elementary dictionary
with open(preComputedIsotopes, 'rb') as f:
    iso_mass_inten_dict = pickle.load(f)


# calculate isotope distribution for an individual element of number 'M'
def element_isoDistr_1toM(iso_mass_inten_dict, element, M, large_num_to_store, inten_threshold_trim):
    
    # for range 2~M
    for i in range(2, M):
        element_intensity_temp = np.array(iso_mass_inten_dict[element]['Intensity'][1]).reshape(-1, 1) @ np.array(
            iso_mass_inten_dict[element]['Intensity'][i - 1]).reshape(1, -1)
        element_mass_temp = np.asmatrix(iso_mass_inten_dict[element]['Mass'][1]).T + np.asmatrix(
            iso_mass_inten_dict[element]['Mass'][i - 1])
        element_intensity_temp = np.array([element_intensity_temp[::-1, :].diagonal(i).sum() for i in
                                           range(-element_intensity_temp.shape[0] + 1,
                                                 element_intensity_temp.shape[1])])
        element_mass_temp = np.array([np.asarray(element_mass_temp[::-1, :].diagonal(i))[0][0] for i in
                                      range(-element_mass_temp.shape[0] + 1, element_mass_temp.shape[1])])
        element_mass_temp_1 = np.array(element_mass_temp);
        element_mass_temp_1 = element_mass_temp_1[np.array(element_intensity_temp) > 1e-5]
        element_intensity_temp_1 = np.array(element_intensity_temp);
        element_intensity_temp_1 = element_intensity_temp_1[np.array(element_intensity_temp) > 1e-5]
        # element_intensity_temp_1
        iso_mass_inten_dict[element]['Intensity'][i] = element_intensity_temp[element_intensity_temp > inten_threshold_trim]  
        # element_mass_temp_1
        iso_mass_inten_dict[element]['Mass'][i] = element_mass_temp[element_intensity_temp > inten_threshold_trim]  
    
    # generating element istopic peaks for 1000, 10000, 100000... so on upto 10E8
    large_num_array = []
    for n in (large_num_to_store):
        large_num_array.append(np.repeat(n, 9).tolist())
    large_num_array = [j for i in large_num_array for j in i]
    if len(large_num_array) > 0:
        iso_inten_temp = iso_mass_inten_dict[element]['Intensity'][large_num_array[0]]
        iso_mass_temp = iso_mass_inten_dict[element]['Mass'][large_num_array[0]]
        current_num = large_num_array[0]
        inten_threshold_trim2 = 1e-10
        for n in large_num_array:
            current_num += n
            element_intensity_temp = np.array(iso_mass_inten_dict[element]['Intensity'][n]).reshape(-1, 1) * np.array(
                iso_inten_temp).reshape(1, -1)
            element_mass_temp = np.array(iso_mass_inten_dict[element]['Mass'][n]).reshape(-1, 1) + np.array(
                iso_mass_temp).reshape(1, -1)
            iso_inten_temp = np.array([element_intensity_temp[::-1, :].diagonal(i).sum() for i in
                                       range(-element_intensity_temp.shape[0] + 1, element_intensity_temp.shape[1])])
            iso_mass_temp = np.array([np.asarray(element_mass_temp[::-1, :].diagonal(i))[0] for i in
                                      range(-element_mass_temp.shape[0] + 1, element_mass_temp.shape[1])])
            iso_mass_temp = iso_mass_temp[iso_inten_temp > inten_threshold_trim2]
            iso_inten_temp = iso_inten_temp[iso_inten_temp > inten_threshold_trim2]
            if current_num in large_num_to_store:
                iso_mass_inten_dict[element]['Intensity'][current_num] = iso_inten_temp[
                    iso_inten_temp > inten_threshold_trim]
                iso_mass_inten_dict[element]['Mass'][current_num] = iso_mass_temp[iso_inten_temp > inten_threshold_trim]
    return iso_mass_inten_dict


# Creating a dictionary with isotopic peak intensity and mass for the mono elemnts ( with cutoff )  
def isotope_distribution_indElement(elemInfo_dict, iso_mass_inten_dict, inten_threshold_trim):
    large_num_to_store= np.power(10, range(1,5)).tolist()
    for element in list(elemInfo_dict.keys()):
        iso_mass_inten_dict[element] = {'Mass': {1: list(elemInfo_dict[element].keys())},
                                        'Intensity': {1: list(elemInfo_dict[element].values())}}
        iso_mass_inten_dict = element_isoDistr_1toM(iso_mass_inten_dict, element, 11, large_num_to_store, inten_threshold_trim)
    return iso_mass_inten_dict


"""std_aa_comp{}
A dictionary with elemental compositions of the twenty standard amino acid residues, 
amino acid modifications, selenocysteine, pyrrolysine, and standard H- and -OH terminal groups.
"""
std_aa_comp = {
    'A': {'H': 5, 'C': 3, 'O': 1, 'N': 1},
    'C': {'H': 5, 'C': 3, 'O': 1, 'N': 1, 'S': 1, },  # Static modification of cysteine alkylation (addition of C2H3NO)
    'D': {'H': 5, 'C': 4, 'O': 3, 'N': 1},
    'E': {'H': 7, 'C': 5, 'O': 3, 'N': 1},
    'F': {'H': 9, 'C': 9, 'O': 1, 'N': 1},
    'G': {'H': 3, 'C': 2, 'O': 1, 'N': 1},
    'H': {'H': 7, 'C': 6, 'O': 1, 'N': 3, },
    'I': {'H': 11, 'C': 6, 'O': 1, 'N': 1},
    'K': {'H': 12, 'C': 6, 'O': 1, 'N': 2},
    'L': {'H': 11, 'C': 6, 'O': 1, 'N': 1},
    'M': {'H': 9, 'C': 5, 'O': 1, 'N': 1, 'S': 1, },
    'N': {'H': 6, 'C': 4, 'O': 2, 'N': 2},
    'P': {'H': 7, 'C': 5, 'O': 1, 'N': 1},
    'Q': {'H': 8, 'C': 5, 'O': 2, 'N': 2},
    'R': {'H': 12, 'C': 6, 'O': 1, 'N': 4},
    'S': {'H': 5, 'C': 3, 'O': 2, 'N': 1},
    'T': {'H': 7, 'C': 4, 'O': 2, 'N': 1},
    'V': {'H': 9, 'C': 5, 'O': 1, 'N': 1},
    'W': {'H': 10, 'C': 11, 'O': 1, 'N': 2, },
    'Y': {'H': 9, 'C': 9, 'O': 2, 'N': 1},
    'U': {'H': 5, 'C': 3, 'O': 1, 'N': 1, 'Se': 1},
    'O': {'H': 19, 'C': 12, 'O': 2, 'N': 3},
    'H-': {'H': 1},
    '-OH': {'O': 1, 'H': 1},
}


# Module to convert the amino acid seq to chemical element dictionary    
def pepSeq_to_chemComp(pep_seq, Charge, aa_comp, TMT_ver):
    chem_comp = defaultdict(int)
    for aa in pep_seq:
        if aa in aa_comp:
            for elem, cnt in aa_comp[aa].items():
                chem_comp[elem] += cnt
        else:
            print('No information for %s in `aa_comp`' % aa)
            
    chem_comp['H'] += 2 + (1 * Charge)
    chem_comp['O'] += 1

    # Add TMT balancer
    count_K = pep_seq.count('K')
    
    if TMT_ver == 'None':
        chem_comp = chem_comp
    elif TMT_ver == 'TMTpro':
        chem_comp['C'] += 8 * (count_K + 1)
        chem_comp['N'] += 1 * (count_K + 1)
        chem_comp['H'] += 25 * (count_K + 1)
        chem_comp['O'] += 3 * (count_K + 1)

        chem_comp['x'] += 7 * (count_K + 1)  # x = C13
        chem_comp['y'] += 2 * (count_K + 1)  # y = N15
    elif TMT_ver == 'TMT':
        chem_comp['C'] += 8 * (count_K + 1)
        chem_comp['N'] += 1 * (count_K + 1)
        chem_comp['H'] += 20 * (count_K + 1)
        chem_comp['O'] += 2 * (count_K + 1)

        chem_comp['x'] += 4 * (count_K + 1)  # x = C13
        chem_comp['y'] += 1 * (count_K + 1)  # y = N15
    return chem_comp


# 
def gen_array_combi(n, element):
    com_array = []
    c = 0
    if element in ['C', 'N', 'H', 'O']:
        while (n > 0):
            r = n % 10
            if ((c > 1) & (r > 1)):
                com_array.append(np.repeat(int((r * (10 ** c)) / r), r).tolist())
            elif (r > 0):
                com_array.append([r * (10 ** c)])
            c += 1
            n //= 10
    elif element == 'S':
        while (n > 0):
            r = n % 10
            if ((c >= 1) & (r > 1)):
                com_array.append(np.repeat(int((r * (10 ** c)) / r), r).tolist())
            elif (r > 0):
                com_array.append([r * (10 ** c)])
            c += 1
            n //= 10
    else:
        while (n > 0):
            r = n % 10
            if ((c >= 0) & (r > 0)):
                com_array.append(np.repeat(int((r * (10 ** c)) / r), r).tolist())
            elif (r > 0):
                com_array.append([r * (10 ** c)])
            c += 1
            n //= 10
    com_array = [j for i in com_array for j in i]
    return com_array


# calculate element isotope distribution if its size greater than default one
def iso_distri_largeNum(element, count, iso_mass_inten_dict):
    
    # 
    gen_iso_combi_array = gen_array_combi(count, element)
    iso_inten_temp = iso_mass_inten_dict[element]['Intensity'][gen_iso_combi_array[0]]
    iso_mass_temp = iso_mass_inten_dict[element]['Mass'][gen_iso_combi_array[0]]
    # 
    for n in gen_iso_combi_array[1:]:
        element_intensity_temp = np.array(iso_mass_inten_dict[element]['Intensity'][n]).reshape(-1, 1) * \
                                 np.array(iso_inten_temp).reshape(1, -1)
        element_mass_temp = np.array(iso_mass_inten_dict[element]['Mass'][n]).reshape(-1, 1) + \
                            np.array(iso_mass_temp).reshape(1, -1)
        iso_inten_temp = np.array([element_intensity_temp[::-1, :].diagonal(i).sum() for i in
                                   range(-element_intensity_temp.shape[0] + 1, element_intensity_temp.shape[1])])
        iso_mass_temp = np.array([np.asarray(element_mass_temp[::-1, :].diagonal(i))[0] for i in
                                  range(-element_mass_temp.shape[0] + 1, element_mass_temp.shape[1])])
        
    pep_iso_distr_df = pd.DataFrame(iso_inten_temp, columns=['isotope_inten'])
    pep_iso_distr_df['isotope_mass'] = iso_mass_temp
    
    return pep_iso_distr_df


pd.set_option('mode.chained_assignment', None)


#
def iso_distri_combine_eleme(pep_iso_distr_df, next_elem_iso_inesity_distr, next_elem_iso_mass_distr):
    
    #
    peptide_intensity = np.array(pep_iso_distr_df.isotope_inten.values).reshape(-1, 1) * \
                        np.array(next_elem_iso_inesity_distr).reshape(1, -1)
    peptide_mass = np.array(np.asmatrix(pep_iso_distr_df.isotope_mass.values).T + \
                            np.asmatrix(next_elem_iso_mass_distr))
    
    # 
    pep_iso_distr_df = pd.DataFrame((np.concatenate((peptide_mass.reshape(-1, 1), peptide_intensity.reshape(-1, 1)), axis=1)),
                                    columns=["isotope_mass", "isotope_inten"])
    #
    pep_iso_distr_df = pep_iso_distr_df[pep_iso_distr_df.isotope_inten > 1e-10].sort_values(['isotope_inten'], ascending=[False])
    #
    if pep_iso_distr_df.shape[0] > 50000:
        pep_iso_distr_df = pep_iso_distr_df.head(50000)
        
    return pep_iso_distr_df


#
def iso_distri_neutral(iso_mass_inten_dict, chemical_com, isotope_cutoff, mass_tolerance, mass_calculation_method):
    # mass_calculation_method = 1 ('weighted') or 2 ('strongest')
    
    # default No. of elements that have pre-calculated isotope distribution distribution 
    default_elementDict_size = {'C': 200, 'N': 100, 'H': 300, 'O': 100, 'S': 10, 'P': 10, 'F': 10, 'Na': 10, 'K': 10, 'Si': 10, 
                                'Cl': 10, 'Mg': 10, 'Fe': 10, 'Ca': 10, 'Zn': 10, 'Br': 10, 'Pb': 10, 'Cu': 10, 'Al': 10, 'Cd': 10, 
                                'I': 10, 'Ti': 10, 'B': 10, 'Se': 10, 'Ni': 10, 'Mn': 10, 'As': 10, 'Li': 10, 'Mo': 10, 'Co': 10, 
                                'x': 10, 'y': 10}
    
    # if any(k > 200 for k in list(chemical_com.values())):
    if any(chemical_com[k] > default_elementDict_size[k] for k in set(chemical_com).intersection(default_elementDict_size)):  
        
        # calculate an initial distribution
        for element, count in list(chemical_com.items())[:1]:
            if count > default_elementDict_size[element]:
                pep_iso_distr_df = iso_distri_largeNum(element, count, iso_mass_inten_dict)
            else:  # if elemnt count is not grater than deafaults size
                pep_iso_distr_df = pd.DataFrame((iso_mass_inten_dict[element]['Intensity'][count]),
                                                columns=['isotope_inten'])
                pep_iso_distr_df['isotope_mass'] = (iso_mass_inten_dict[element]['Mass'][count])
        
        # combine other element to the initial distribution
        for element, count in list(chemical_com.items())[1:]:
            if count > default_elementDict_size[element]:
                pep_iso_distr_df_ = iso_distri_largeNum(element, count, iso_mass_inten_dict)
                pep_iso_distr_df = iso_distri_combine_eleme(pep_iso_distr_df, pep_iso_distr_df_.isotope_inten.values,
                                                            pep_iso_distr_df_.isotope_mass.values)
            else:
                pep_iso_distr_df = iso_distri_combine_eleme(pep_iso_distr_df,
                                                            iso_mass_inten_dict[element]['Intensity'][count],
                                                            iso_mass_inten_dict[element]['Mass'][count])
    else: # if any elemnt size is not greater than default elemnts size
        # calculate an initial distribution
        for element, count in list(chemical_com.items())[:1]:
            pep_iso_distr_df = pd.DataFrame((iso_mass_inten_dict[element]['Intensity'][count]),
                                            columns=['isotope_inten'])
            pep_iso_distr_df['isotope_mass'] = (iso_mass_inten_dict[element]['Mass'][count])
            
        # combine other element to the initial distribution
        for element, count in list(chemical_com.items())[1:]:  # chemical_com.items():
            pep_iso_distr_df = iso_distri_combine_eleme(pep_iso_distr_df,
                                                        iso_mass_inten_dict[element]['Intensity'][count],
                                                        iso_mass_inten_dict[element]['Mass'][count])

    pep_iso_distr_df.columns = ["pep_mass", "pep_intensity"]
    pep_iso_distr_df = pep_iso_distr_df.sort_values(['pep_intensity'], ascending=[False])
    pep_iso_distr_df_temp = pep_iso_distr_df[pep_iso_distr_df.pep_intensity > 1e-10].sort_values(['pep_intensity'],ascending=[False])
    pep_iso_distr_df_temp['groups'] = 0
    
    # define groups to aggregate close peaks
    i = 0
    while len(pep_iso_distr_df_temp.loc[pep_iso_distr_df_temp['groups'] == 0, 'pep_intensity']) > 0:
        i += 1
        maxIntensity_isopeak = (pep_iso_distr_df_temp.loc[pep_iso_distr_df_temp.groups == 0])['pep_mass'].values[0]
        lb = maxIntensity_isopeak - ((mass_tolerance / 1e6) * maxIntensity_isopeak)  # [0]
        ub = maxIntensity_isopeak + ((mass_tolerance / 1e6) * maxIntensity_isopeak)  # [0]
        pep_iso_distr_df_temp.loc[pep_iso_distr_df_temp['pep_mass'].between(lb, ub, inclusive="neither"), 'groups'] = i
    
    pep_iso_distr_df = pep_iso_distr_df_temp.copy()
    
    # aggregate peaks of the same group
    ## mass_calculation_method = 1 ('weighted') or 2 ('strongest')
    if mass_calculation_method == 1:
        pep_iso_distr_df['Rel_inten'] = pep_iso_distr_df['pep_intensity'] / pep_iso_distr_df.groupby('groups')['pep_intensity'].transform('sum')
        pep_iso_distr_df['weighted_mass'] = pep_iso_distr_df['pep_mass'] * pep_iso_distr_df['Rel_inten']
        pep_iso_distr_df = pep_iso_distr_df.groupby(['groups']).agg(isotope_mass=('weighted_mass', 'sum'),
                                                                    isotope_inten=('pep_intensity', 'sum'))
    elif mass_calculation_method == 2:
        pep_iso_distr_df = pep_iso_distr_df.groupby(['groups']).agg(isotope_mass=('pep_mass', 'first'),
                                                                    isotope_inten=('pep_intensity', 'sum'))
    
    # re-order based on isotope mass
    pep_iso_distr_df = pep_iso_distr_df.sort_values(['isotope_mass'],ascending=[True])
    
    # define mass positions
    pep_mono_mass = pep_iso_distr_df.isotope_mass.values[0]
    pep_iso_distr_df['Mass_position'] = round(pep_iso_distr_df['isotope_mass']-pep_mono_mass)
    
    return pep_iso_distr_df


def iso_distri(iso_mass_inten_dict, chemical_com, Charge, isotope_cutoff, mass_tolerance,
               mass_calculation_method, is_pep):
    #
    pep_iso_distr_df = iso_distri_neutral(iso_mass_inten_dict, chemical_com, isotope_cutoff, mass_tolerance, mass_calculation_method)
    #pep_iso_distr_df['isotope_mass'] = (pep_iso_distr_df['isotope_mass'].values + 1.007276466812000 * Charge) / abs(Charge)  # electron mass
    pep_iso_distr_df['isotope_mass'] = (pep_iso_distr_df['isotope_mass'].values - 0.000548579909460 * Charge*is_pep)/abs(Charge)
    #proton_mass = 1.007276466812
    #pep_iso_distr_df['isotope_mass'] = (pep_iso_distr_df['isotope_mass'].values + proton_mass * Charge)/abs(Charge)
    
    pep_iso_distr_df['isotope_inten'] = pep_iso_distr_df.isotope_inten / pep_iso_distr_df.isotope_inten.max()
    pep_iso_distr_df = pep_iso_distr_df[pep_iso_distr_df['isotope_inten'] > isotope_cutoff]
    pep_iso_distr_df['isotope_inten'] = (pep_iso_distr_df.isotope_inten / pep_iso_distr_df.isotope_inten.sum()) * 100

    return pep_iso_distr_df
