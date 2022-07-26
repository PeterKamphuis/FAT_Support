#!/usr/local/bin/ python3
# -*- coding: utf-8 -*-
import numpy as np
import support_functions as sf

from astropy.io import fits
def main():
    '''Analyze both the Database and LVHIS fits of the GDL runs.'''
    version = 'v2.0.1'
    LVHIS_Directory = f'/home/peter/FAT/LVHIS-26_GDL_{version}'
    LVHIS_config_file = 'FAT_INPUT.config'
    Database_Directory = f'/home/peter/FAT/Database_GDL_{version}'
    Database_config_file = 'FAT_INPUT.config'

    #First analyze the database.
    sf.analyze(Database_Directory,Database_config_file, basename=f'GDL_{version}_Results',GDL=True)
    sf.analyze(LVHIS_Directory,LVHIS_config_file, basename=f'GDL_{version}_Results',LVHIS=True,GDL=True)




if __name__ == '__main__':
    main()
