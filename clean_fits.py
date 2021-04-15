#!/usr/local/bin/ python3

import os
import shutil
import glob
import numpy as np
#This script goes through a FAT fitted main directory and cleans all the files from previous fits.



maindir= '/home/peter/FAT/LVHIS-26_pyFAT_v0.0.2'
cat = 'Data_LVHIS.txt'
#clean all
clean = ['One_Step','Old_pyFAT', 'GDL','Finalmodel','Sofia_Output' ]
#clean GDL
clean = ['GDL']
#clean pyFAT
#clean = ['One_Step','Old_pyFAT']

#clean = ['GDL']

def catalogue(filename, debug = False):
    Catalogue = {}
    with open(filename,'r') as tmpfile:
        #Define the exsiting catalogue input()
        input_columns = [x.strip().upper() for x in tmpfile.readline().split('|')]
        Catalogue['ENTRIES'] = ['ENTRIES']
        Catalogue['ENTRIES'].extend(input_columns)
        for key in input_columns:
            Catalogue[key] = []

        for line in tmpfile.readlines():
            input = [x.strip() for x  in line.split('|')]
            for i,key in enumerate(input_columns):
                if key == 'DISTANCE':
                    Catalogue[key].append(float(input[i]))
                else:
                    Catalogue[key].append(input[i])
    if 'NUMBER' in Catalogue['ENTRIES']:
        Catalogue['NUMBER'] = np.array(Catalogue['NUMBER'],dtype=int)

    return Catalogue

def remove_these(maindir,list,dirs,files):
    for subs in list:
        for dir in dirs:
            try:
                print(f"removing {subs}/{dir}")
                shutil.rmtree(f"{subs}/{dir}")
            except:
                pass
        for file in files:
            if '*' in file:
                list_files = glob.glob(f"{subs}/{file}")
                for wc_files in list_files:
                    try:
                        print(f"removing {wc_files}")
                        os.remove(f"{wc_files}")
                    except:
                        pass
            else:
                try:
                    print(f"removing {subs}/{file}")
                    os.remove(f"{subs}/{file}")
                except:
                    pass

input_cat = catalogue(f"{maindir}/{cat}")
try:
    list_subfolders_with_paths = [f"{maindir}/{x}" for x in input_cat['DIRNAME']]
except KeyError:
    list_subfolders_with_paths = [f"{maindir}/{x}" for x in input_cat['DIRECTORYNAME']]

if 'One_Step' in clean:
    dirs = ['Logs','One_Step_Convergence','tmp_incl_check']
    files = ['Overview_Prev.png','*_FAT-Basic_Info.txt','*_FAT.fits','tmp_incl_check_In.def']
    remove_these(maindir,list_subfolders_with_paths,dirs,files)



if 'Old_pyFAT' in clean:
    dirs = ['Logs','Centre_Convergence','Extend_Convergence','Sofia_Output']
    files = ['Overview_Prev.png','*_FAT-Basic_Info.txt','*_FAT.fits']
    remove_these(maindir,list_subfolders_with_paths,dirs,files)

if 'GDL' in clean:
    dirs = ['Intermediate','No_Warp','Warp_Info','PV-Diagrams','Residuals','Moments']
    files = ['Prev_Log.txt','Log.txt','the_last_input.def', 'Overview_Prev.png' ,'BasicInfo*.txt','*_preprocessed*.fits']
    remove_these(maindir,list_subfolders_with_paths,dirs,files)

if 'Finalmodel' in clean:
    dirs = ['Finalmodel']
    files = ['Overview.png']
    remove_these(maindir,list_subfolders_with_paths,dirs,files)

if 'Sofia_Output' in clean:
    dirs = ['Sofia_Output']
    files = []
    remove_these(maindir,list_subfolders_with_paths,dirs,files)
