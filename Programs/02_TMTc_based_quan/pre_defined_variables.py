# -*- coding: utf-8 -*-
"""
Created on 08/09/2022

@author: yfu
"""

import numpy as np


############################# Pre-defined variables ###########################
TMTtag_modi_mass={'TMT': 229.162932,
                  'TMTpro': 304.207145}
electron_mass=0.000548579909065
proton_mass=1.007276466621
hydrogen_mass=1.00782503207

# Predifined sample mixing method
N_mix = ['sig126','sig128N','sig129N','sig130N','sig131N','sig131C','sig132N','sig134N','sig135N']
C_mix = ['sig127C','sig128C','sig129C','sig130C','sig131N','sig131C','sig132C','sig134C','sig135N']
mix_18 = ["sig126","sig127N","sig127C","sig128N","sig128C","sig129N","sig129C","sig130N","sig130C",
          "sig131N","sig131C","sig132N","sig132C","sig133N","sig133C","sig134N","sig134C","sig135N"]


# TMT(pro) tag balancer part impurities
TMT_balancer_impurity={
    "sig126": np.array([0.0321, 0.9602, 0.0077]),
    "sig127N": np.array([0.0301, 0.9567, 0.0132]),
    "sig127C": np.array([0.0321, 0.9602, 0.0077]),
    "sig128N": np.array([0.0301, 0.9567, 0.0132]),
    "sig128C": np.array([0.0301, 0.9567, 0.0132]),
    "sig129N": np.array([0.0193, 0.9610, 0.0197]),
    "sig129C": np.array([0.0193, 0.9610, 0.0197]),
    "sig130N": np.array([0.0295, 0.9404, 0.0301]),
    "sig130C": np.array([0.0295, 0.9404, 0.0301]),
    "sig131N": np.array([0.0127, 0.9481, 0.0392]),
    "sig131C": np.array([0.0102, 0.9366, 0.0532])}

TMTpro_balancer_impurity={
    "sig126": np.array([0.0321, 0.9602, 0.0077]),
    "sig127N": np.array([0.0301, 0.9567, 0.0132]),
    "sig127C": np.array([0.0321, 0.9602, 0.0077]),
    "sig128N": np.array([0.0301, 0.9567, 0.0132]),
    "sig128C": np.array([0.0301, 0.9567, 0.0132]),
    "sig129N": np.array([0.0193, 0.9610, 0.0197]),
    "sig129C": np.array([0.0193, 0.9610, 0.0197]),
    "sig130N": np.array([0.0295, 0.9404, 0.0301]),
    "sig130C": np.array([0.0295, 0.9404, 0.0301]),
    "sig131N": np.array([0.0127, 0.9481, 0.0392]),
    "sig131C": np.array([0.0102, 0.9366, 0.0532]),
    "sig132N": np.array([0.0074, 0.9356, 0.057]),
    "sig132C": np.array([0.0074, 0.9356, 0.057]),
    "sig133N": np.array([0., 0.9311, 0.0689]),
    "sig133C": np.array([0.0074, 0.9356, 0.057]),
    "sig134N": np.array([0., 0.9311, 0.0689]),
    "sig134C": np.array([0., 0.9311, 0.0689]),
    "sig135N": np.array([0., 0.9327, 0.0673])}

TMTtag_balancer_impurity_dict = {'TMT': TMT_balancer_impurity,
                                 'TMTpro': TMTpro_balancer_impurity}

# TMT(pro) whole tag impurities
TMT_whole_impurity={
    "sig126": np.array([0.0635, 0.8446, 0.0919]),
    "sig127N": np.array([0.0679, 0.8383, 0.0938]),
    "sig127C": np.array([0.0635, 0.8446, 0.0919]),
    "sig128N": np.array([0.0679, 0.8383, 0.0938]),
    "sig128C": np.array([0.0679, 0.8383, 0.0938]),
    "sig129N": np.array([0.056, 0.8575, 0.0865]),
    "sig129C": np.array([0.056, 0.8575, 0.0865]),
    "sig130N": np.array([0.0712, 0.8433, 0.0855]),
    "sig130C": np.array([0.0712, 0.8433, 0.0855]),
    "sig131N": np.array([0.0669, 0.8458, 0.0873]),
    "sig131C": np.array([0.0622, 0.8477, 0.0901])}

TMTpro_whole_impurity={
    "sig126": np.array([0.0635, 0.8446, 0.0919]),
    "sig127N": np.array([0.0679, 0.8383, 0.0938]),
    "sig127C": np.array([0.0635, 0.8446, 0.0919]),
    "sig128N": np.array([0.0679, 0.8383, 0.0938]),
    "sig128C": np.array([0.0679, 0.8383, 0.0938]),
    "sig129N": np.array([0.056, 0.8575, 0.0865]),
    "sig129C": np.array([0.056, 0.8575, 0.0865]),
    "sig130N": np.array([0.0712, 0.8433, 0.0855]),
    "sig130C": np.array([0.0712, 0.8433, 0.0855]),
    "sig131N": np.array([0.0669, 0.8458, 0.0873]),
    "sig131C": np.array([0.0622, 0.8477, 0.0901]),
    "sig132N": np.array([0.0621, 0.845, 0.0929]),
    "sig132C": np.array([0.0621, 0.845, 0.0929]),
    "sig133N": np.array([0.0715, 0.8471, 0.0814]),
    "sig133C": np.array([0.0621, 0.845, 0.0929]),
    "sig134N": np.array([0.0715, 0.8471, 0.0814]),
    "sig134C": np.array([0.0715, 0.8471, 0.0814]),
    "sig135N": np.array([0.07, 0.85, 0.08])}

TMTtag_whole_impurity_dict = {'TMT': TMT_whole_impurity,
                              'TMTpro': TMTpro_whole_impurity}

## difine average TMT whole tag impurity
TMTtag_whole_impurity_mean = {'TMT': np.array([0.066, 0.846, 0.088]),
                              'TMTpro': np.array([0.066, 0.846, 0.088])}


# TMTpro HCD Monoisotopic Reporter Mass (from thermo TMTpro user guide)
TMTtag_reportor_mass={
    "sig126": 126.127726,
    "sig127N": 127.124761,
    "sig127C": 127.131081,
    "sig128N": 128.128116,
    "sig128C": 128.134436,
    "sig129N": 129.131471,
    "sig129C": 129.13779,
    "sig130N": 130.134825,
    "sig130C": 130.141145,
    "sig131N": 131.13818,
    "sig131C": 131.1445,
    "sig132N": 132.141535,
    "sig132C": 132.147855,
    "sig133N": 133.14489,
    "sig133C": 133.15121,
    "sig134N": 134.148245,
    "sig134C": 134.154565,
    "sig135N": 135.1516}


CO_l_mass=27.99491461956 # light: C[12]O 
CO_h_mass=28.99826945736 # heavy: C[13]O

TMTpro_CO_mass={
    "sig126" : CO_h_mass,
    "sig127N": CO_h_mass,
    "sig127C": CO_l_mass,
    "sig128N": CO_l_mass,
    "sig128C": CO_l_mass,
    "sig129N": CO_l_mass,
    "sig129C": CO_l_mass,
    "sig130N": CO_l_mass,
    "sig130C": CO_l_mass,
    "sig131N": CO_l_mass,
    "sig131C": CO_h_mass,
    "sig132N": CO_h_mass,
    "sig132C": CO_h_mass,
    "sig133N": CO_h_mass,
    "sig133C": CO_l_mass,
    "sig134N": CO_l_mass,
    "sig134C": CO_l_mass,
    "sig135N": CO_l_mass}

TMT11_CO_mass={
    "sig126" : CO_h_mass,
    "sig127N": CO_h_mass,
    "sig127C": CO_h_mass,
    "sig128N": CO_l_mass,
    "sig128C": CO_h_mass,
    "sig129N": CO_h_mass,
    "sig129C": CO_h_mass,
    "sig130N": CO_h_mass,
    "sig130C": CO_l_mass,
    "sig131N": CO_l_mass,
    "sig131C": CO_l_mass}

TMTtag_CO_mass={'TMT': TMT11_CO_mass,
                'TMTpro': TMTpro_CO_mass}


# NO. of heavy isotopes in complementary ion
TMTproC_n_heavy={
    "sig126" : 'h8',
    "sig127N": 'h7',
    "sig127C": 'h8',
    "sig128N": 'h7',
    "sig128C": 'h7',
    "sig129N": 'h6',
    "sig129C": 'h6',
    "sig130N": 'h5',
    "sig130C": 'h5',
    "sig131N": 'h4',
    "sig131C": 'h3',
    "sig132N": 'h2',
    "sig132C": 'h2',
    "sig133N": 'h1',
    "sig133C": 'h2',
    "sig134N": 'h1',
    "sig134C": 'h1',
    "sig135N": 'h0'}

TMT11C_n_heavy={
    "sig126" : 'h4',
    "sig127N": 'h3',
    "sig127C": 'h3',
    "sig128N": 'h3',
    "sig128C": 'h2',
    "sig129N": 'h1',
    "sig129C": 'h1',
    "sig130N": 'h0',
    "sig130C": 'h1',
    "sig131N": 'h0',
    "sig131C": 'h0'}

TMTc_n_heavy_dict={'TMT': TMT11C_n_heavy,
                   'TMTpro': TMTproC_n_heavy}

min_n_peaks_fit={'TMT': 3,
                 'TMTpro': 5}