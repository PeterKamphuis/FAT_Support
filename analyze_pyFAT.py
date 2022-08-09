#!/usr/local/bin/ python3
# -*- coding: utf-8 -*-
import numpy as np
import support_functions as sf

from astropy.io import fits
def main():
    '''Analyze both the Database and LVHIS fits of the GDL runs.'''
    version = 'v0.0.9'
    LVHIS_Directory = f'/home/peter/FAT/LVHIS-26_pyFAT_{version}'
    LVHIS_config_file = 'FAT_defaults.yml'
    Database_Directory = f'/home/peter/FAT/Database_pyFAT_{version}'
    Database_config_file = 'pyFAT_fit.yml'

    #First analyze the database.
    sf.analyze(Database_Directory,Database_config_file, basename=f'pyFAT_{version}_Results')
    #sf.analyze(LVHIS_Directory,LVHIS_config_file, basename=f'pyFAT_{version}_Results',LVHIS=True)




if __name__ == '__main__':
    main()
