#!/usr/local/bin/ python3
# -*- coding: utf-8 -*-
import numpy as np
import support_functions as sf
import matplotlib.font_manager as mpl_fm
from astropy.io import fits
import os
import pickle
def main():
    '''Compare the pyFAT and GDL runs.'''
    Pversion = 'v0.1.4'
    Gversion = 'v2.0.2'
    read_all_input = True
    adddelt = True
    analyze_individual = True
    directory = '/home/peter/FAT_Main/Analysis/Compare_GDL_pyFAT/'
    create_individual = ['pyFAT','GDL']
    mpl_fm.fontManager.addfont("/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf")
    #mpl_fm.fontManager.addfont("/home/pkamphuis/Fonts/Times_New_Roman.ttf")
    missing_links = True
    Input_File = {'pyFAT':{'version':Pversion,
                           'program': 'pyFAT',
                            'Database': {'dir': f'{directory}pyFAT_{Pversion}/Database/'\
                                        ,'config':'FAT_defaults.yml'},
                            'LVHIS': {'dir': f'{directory}pyFAT_{Pversion}/LVHIS-26/'\
                                    ,'config':'FAT_defaults.yml'}},
                  'GDL':{'version':Gversion,
                         'program': 'GDL',
                            'Database': {'dir': f'{directory}GDL_{Gversion}/Database/'\
                                        ,'config':'FAT_INPUT.config'},
                            'LVHIS': {'dir': f'{directory}GDL_{Gversion}/LVHIS-26/'\
                                    ,'config':'FAT_INPUT.config'}}}


    if read_all_input:
        Input_Parameters = sf.obtain_parameters(Input_File,\
            missing_links=missing_links)
        with open(f'Input_Parameters_GDL{Gversion}_pyFAT{Pversion}.pkl','wb') as tmp:
            pickle.dump(Input_Parameters,tmp)
    else:
        with open(f'Input_Parameters_GDL{Gversion}_pyFAT{Pversion}.pkl', 'rb') as f:
            Input_Parameters = pickle.load(f)

    if adddelt:
        sf.add_deltas(Input_Parameters)
        with open(f'Input_Parameters_delt_GDL{Gversion}_pyFAT{Pversion}.pkl','wb') as tmp:
            pickle.dump(Input_Parameters,tmp)
        print(Input_Parameters['ESO_223_G009_8.8Beams_3.0SNR'])
    else:
        with open(f'Input_Parameters_delt_GDL{Gversion}_pyFAT{Pversion}.pkl', 'rb') as f:
            Input_Parameters = pickle.load(f)
    if analyze_individual:
        for program in create_individual:
            for database in Input_File[program]:
                if database == 'version':
                    continue
                sf.analyze(program,database,Input_Parameters,basename=\
                    f'{program}_{Input_File[program]["version"]}_Results',\
                    main_directory=Input_File[program][database]['dir'])
                if program == 'pyFAT':
                    if os.path.isfile(f'{Input_File[program][database]["dir"]}/Timing_Result.txt'):
                        sf.analyze_timing(Input_File[program][database]["dir"],Input_Parameters,\
                            basename=f'pyFAT_{Input_File[program]["version"]}_Results',\
                             database = database, program=program)


    programs = [x for x in Input_File]
    if len(programs) > 1:
        #Compare GDL and pyFAT results
        sf.compare(Input_Parameters,programs=programs)


    #individual_results = sf.obtain_individual_results(Input_Parameters)

    exit()

    #First analyze the database.
    deltas = sf.analyze(Database_Directory,Database_config_file, basename=f'pyFAT_{version}_Results',missing_links=missing_links)
    deltas_LVHIS = sf.analyze(LVHIS_Directory,LVHIS_config_file, basename=f'pyFAT_{version}_Results',missing_links=missing_links,LVHIS=True)

    if os.path.isfile(f'{Database_Directory}/Timing_Result.txt'):
        sf.analyze_timing(Database_Directory,basename=f'pyFAT_{version}_Results', deltas = deltas)
    if os.path.isfile(f'{LVHIS_Directory}/Timing_Result.txt'):
        sf.analyze_timing(LVHIS_Directory,basename=f'pyFAT_{version}_Results', deltas = deltas_LVHIS)




if __name__ == '__main__':
    main()
