#!/usr/local/bin/ python3
# -*- coding: utf-8 -*-
import numpy as np
import support_functions as sf
import matplotlib.font_manager as mpl_fm
from astropy.io import fits
import os

def main():
    '''Analyze both the Database and LVHIS fits of the GDL runs.'''
    version = 'v0.1.1'

    mpl_fm.fontManager.addfont("/usr/share/fonts/truetype/msttcorefonts/Times_New_Roman.ttf")
    #mpl_fm.fontManager.addfont("/home/pkamphuis/Fonts/Times_New_Roman.ttf")

    LVHIS_Directory = f'/home/peter/FAT/LVHIS-26_pyFAT_v0.1.1/'
    missing_links = True
    LVHIS_config_file = 'FAT_defaults.yml'
    Database_config_file = 'FAT_defaults.yml'
    Database_Directory = '/home/peter/FAT/Database_pyFAT_v0.1.1/'
    #Database_config_file = 'Noconvolve.yml'
    #Database_Directory = '/home/peter/FAT_Main/Bad_Fits'

    #First analyze the database.
    deltas = sf.analyze(Database_Directory,Database_config_file, basename=f'pyFAT_{version}_Results',missing_links=missing_links)
    deltas_LVHIS = sf.analyze(LVHIS_Directory,LVHIS_config_file, basename=f'pyFAT_{version}_Results',missing_links=missing_links,LVHIS=True)

    if os.path.isfile(f'{Database_Directory}/Timing_Result.txt'):
        sf.analyze_timing(Database_Directory,basename=f'pyFAT_{version}_Results', deltas = deltas)
    if os.path.isfile(f'{LVHIS_Directory}/Timing_Result.txt'):
        sf.analyze_timing(LVHIS_Directory,basename=f'pyFAT_{version}_Results', deltas = deltas_LVHIS)




if __name__ == '__main__':
    main()
