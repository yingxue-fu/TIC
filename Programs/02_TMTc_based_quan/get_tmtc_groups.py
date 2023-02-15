# -*- coding: utf-8 -*-

import pandas as pd
from pre_defined_variables import TMTc_n_heavy_dict


# 
def get_TMTc_groups(used_channels, TMT_ver):
    
    # TMT_ver = 'TMT' / 'TMTpro'
    if TMT_ver == 'TMT':
        res = pd.DataFrame({'TMTc_grp': [TMTc_n_heavy_dict[TMT_ver][channel] for channel in used_channels],
                            'TMT_channels': used_channels})
    elif TMT_ver == 'TMTpro':
        res = pd.DataFrame({'TMTc_grp': [TMTc_n_heavy_dict[TMT_ver][channel] for channel in used_channels],
                            'TMT_channels': used_channels})
        
    return res


# get TMTc groups of used channels
def get_combined_channels(used_channels, TMT_ver):
    
    res = get_TMTc_groups(used_channels, TMT_ver)
        
    #
    res = res.groupby('TMTc_grp', as_index = False).agg({'TMT_channels': '.'.join})
    res = res['TMT_channels'].tolist()
    
    return res



# merge sample ratios belonging to the same TMTc group
def get_combined_ratios(sample_ratios, used_channels, TMT_ver):
    
    res = get_TMTc_groups(used_channels, TMT_ver)
    res['Ratio'] = list(sample_ratios)
        
    #
    res = res.groupby('TMTc_grp', as_index = False).agg({'TMT_channels': '.'.join, 'Ratio': 'sum'})
    
    return res
