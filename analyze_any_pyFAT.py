#!/usr/local/bin/ python3
# -*- coding: utf-8 -*-
import numpy as np
import support_functions as sf
import matplotlib.font_manager as mpl_fm
from astropy.io import fits
import os
import pickle
def main():
    version = 'v0.1.8'
    directory = '/home/peter/FAT_Main/Test_Sets/Proper_Test_For_Beta/'
    missing_links = False
    read_all_input=True
    adddelt= False
    if read_all_input:
        adddelt =True
    if not adddelt:
        read_all_input=False
    mpl_fm.fontManager.addfont("/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf")
 
    #mpl_fm.fontManager.addfont("/home/pkamphuis/Fonts/Times_New_Roman.ttf")
    Input_File = {'pyFAT':{'version':version,
                           'program': 'pyFAT',
                            'Database': {'dir': f'{directory}'\
                                        ,'config':'FAT_defaults.yml'}
                          }
                  }
    #Input_File = {'pyFAT':{'version':version,
    #                       'program': 'pyFAT',
    #                        #'Database': {'dir': f'{directory}/Database_pyFAT_{version}/'\
    #                        #            ,'config':'FAT_defaults.yml'},
    #                        'LVHIS': {'dir': f'{directory}/LVHIS-26_pyFAT_{version}/'\
    #                                ,'config':'FAT_defaults.yml'}}}

    if read_all_input:
        Input_Parameters = sf.obtain_parameters(Input_File,\
            missing_links=missing_links)
        with open(f"{directory}/Input_Parameters_pyFAT{Input_File['pyFAT']['version']}.pkl",'wb') as tmp:
            pickle.dump(Input_Parameters,tmp)
    else:
        with open(f"{directory}/Input_Parameters_pyFAT{Input_File['pyFAT']['version']}.pkl", 'rb') as f:
            Input_Parameters = pickle.load(f)

    if adddelt:
        sf.add_deltas(Input_Parameters)
        with open(f"{directory}/Input_Delt_pyFAT{Input_File['pyFAT']['version']}.pkl",'wb') as tmp:
            pickle.dump(Input_Parameters,tmp)     
    else:
        with open(f"{directory}/Input_Delt_pyFAT{Input_File['pyFAT']['version']}.pkl", 'rb') as f:
            Input_Parameters = pickle.load(f)
    
     
    for database in Input_File['pyFAT']:
        if database in ['version','program']:
            continue
        program=Input_File['pyFAT']['program']
      
        sf.analyze(Input_File['pyFAT']['program'],database,Input_Parameters,basename=\
                f'{program}_{Input_File[program]["version"]}_Results',\
                main_directory=Input_File[program][database]['dir'])
            
        if os.path.isfile(f'{Input_File[program][database]["dir"]}/Timing_Result.txt'):
            sf.analyze_timing(Input_File[program][database]["dir"],Input_Parameters,\
                        basename=f'pyFAT_{Input_File[program]["version"]}_Results',\
                            database = database, program=program)


  

 

if __name__ == '__main__':
    main()
