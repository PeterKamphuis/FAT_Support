#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


from support_functions import load_input_catalogue
from pyFAT_astro.Support.read_functions import tirific_template
from pyFAT_astro.Support.write_functions import tirific as write_template
from astropy.io import fits

import os
import subprocess
import numpy as np
import copy
#The LVHIS conversions are weird but it is not clear why there is an error in the radi
class InputError(Exception):
    pass


def main():

    #Fat in put list for galaxies to convert
    list_to_convert = 'Data_LVHIS.txt'
    GDL_list = True
    main_directory = '/home/peter/FAT_Main/GDL_release_check/LVHIS-26_GDL_v2.0.1'


    #first read in the list

    database_inp_catalogue = load_input_catalogue(f'{main_directory}/{list_to_convert}',GDL=GDL_list)


    for i,directory in enumerate(database_inp_catalogue['DIRECTORYNAME']):
        # first check the directory exists
        if not os.path.exists(f'{main_directory}/{directory}/'):
            print(f'We are skipping the directory {directory} as it is non-existent in {main_directory}.' )
            continue
        else:
            print(f'We are processing the directory {directory}.' )


        #Then load a tirific template
        Tirific_Template=tirific_template()

        if database_inp_catalogue['DISTANCE'][i] != -1:
            Tirific_Template['DISTANCE'] = f"{database_inp_catalogue['DISTANCE'][i]}"

        #read the RC file, this is a straightforward table
        header_originalVF = fits.getheader(f'{main_directory}/{directory}/OriginalVF.fits')

        rc_file =f'{main_directory}/{directory}/Parameters.RC'
        if os.path.exists(rc_file):
            load_RC_to_Template(Tirific_Template,filename=rc_file, header=header_originalVF)
        df_file =f'{main_directory}/{directory}/Parameters.DF'
        if os.path.exists(df_file):
            load_DF_to_Template(Tirific_Template,filename=df_file, header=header_originalVF,
                                    disk_to_load_to= 2)

        for key in Tirific_Template:
            try:
                if float(Tirific_Template[key]) == 10.:
                    Tirific_Template[key] = '0.'
            except:
                pass
        if 'BMAJ' in header_originalVF:
            Tirific_Template['BMAJ'] = str(header_originalVF['BMAJ']*3600.)
        if 'BMIN' in header_originalVF:
            Tirific_Template['BMIN'] = str(header_originalVF['BMIN']*3600.)
        if 'BPA' in header_originalVF:
            Tirific_Template['BPA'] = header_originalVF['BPA']


        write_template({'OUTPUTLOG': None,'FITTING_DIR':''},Tirific_Template,name =f'{main_directory}/{directory}/ModelInput.def')
        #os.remove(f'{main_directory}/{directory}/Model_Input.def')
        #os.remove(f'{main_directory}/{directory}/Model_Input.Compare.def')


def load_DF_to_Template(Tirific_Template,filename='DiskFit.df',disk_to_load_to = 1,\
                        variables_to_transfer=['RADI','VROT','VROT_ERR',\
                                        'INCL','INCL_ERR','PA','PA_ERR',\
                                        'XPOS','XPOS_ERR','YPOS','YPOS_ERR',\
                                        'VSYS','VSYS_ERR'], header =None):
    column_translation_dict = {'RADI': 0,
                            'VROT': 2,
                            '# VROT_ERR': 3,
                            'VRAD': 4,
                            '# VRAD_ERR': 5,
                            'VMT': 6,
                            '# VMT_ERR':7,
                            'VMR': 8,
                            '# VMR_ERR':9,
                            'NPTS':1}

    if disk_to_load_to == 1:
        ext=''
    else:
        ext = f'_{int(disk_to_load_to)}'
        if f'RADI' in variables_to_transfer and f'RADI{ext}' not in Tirific_Template:
            prev_key = 'Empty'
            for key in Tirific_Template:
                if key == f'VROT{ext}':
                    Tirific_Template.insert(prev_key,f'RADI{ext}','10.')
                    break
                else:
                    prev_key = copy.deepcopy(key)
    if 'INCL' in variables_to_transfer:
        try:
            Original_INCL =  Tirific_Template[f'INCL{ext}']
        except KeyError:
            Original_INCL = 'EMPTY'

    try:
        with open(filename,'r') as DF_File:
            DF_lines = DF_File.readlines()

        for i,line in enumerate(DF_lines):
            elements = line.split()
            start= " ".join([x.lower() for x in elements[:3]])
            if start == 'fitted velocity components':
                column_start = i+3
        df_upload=np.genfromtxt(filename, skip_header=column_start,dtype=str).transpose()
        nur = len(df_upload[0])
        rmax = df_upload[column_translation_dict['RADI']][-1]
        for i,variable in enumerate(variables_to_transfer):
            if variable[-4:] == '_ERR':
                if variable[:2] != '# ':
                    variable = f'# {variable}'
                    variables_to_transfer[i] = variable
                if f'{variable}{ext}' not in Tirific_Template:
                    stripped = variable[2:-4]
                    Tirific_Template.insert(f'{stripped}{ext}',f'{variable}{ext}','10.')
            if variable in column_translation_dict:
                if variable == 'RADI':
                    selection = [str(abs(float(e))*((abs(header['CDELT1'])\
                                +abs(header['CDELT2']))/2.)*3600.) for e in\
                                df_upload[column_translation_dict[variable]] ]

                else:
                    selection = df_upload[column_translation_dict[variable]]

                Tirific_Template[f'{variable}{ext}'] = combine_string(selection)
        #for the other parameters we need to read the header

        Trigger = False
        phi_m = 0.
        el_m = 0.
        for line in DF_lines:
            ind_val=line.split(' ')
            ind_val=[e.strip().lower() for e in ind_val if e.strip() != '']
            if len(ind_val) > 0:
                if ind_val[0] == 'best' and not Trigger:
                    Trigger = True
                if len(ind_val) > 1:
                        if ind_val[1] == 'pa,':
                            if 'PA' in variables_to_transfer:
                                Tirific_Template[f'PA{ext}'] = \
                                    combine_string(ind_val[4],length=nur)
                            if '# PA_ERR' in variables_to_transfer and Trigger:
                                Tirific_Template[f'# PA_ERR{ext}'] = \
                                    combine_string(ind_val[6],length=nur)
                            else:
                                Tirific_Template[f'# PA_ERR{ext}'] = \
                                   combine_string(0.,length=nur)
                        if ind_val[1] == 'incl':
                            if Trigger:
                                index= 3
                            else:
                                index = 4
                            if 'INCL' in variables_to_transfer:
                                Tirific_Template[f'INCL{ext}'] = \
                                    combine_string(ind_val[index],length=nur)
                            if '# INCL_ERR' in variables_to_transfer and Trigger:
                                Tirific_Template[f'# INCL_ERR{ext}'] = \
                                    combine_string(ind_val[5],length=nur)
                            else:
                                Tirific_Template[f'# INCL_ERR{ext}'] =  \
                                    combine_string(0.,length=nur)
                        if ind_val[0] == 'x,y':
                            if Trigger:
                                y_index= 7
                            else:
                                y_index= 5
                            ypos_df=header['CRVAL2']+header['CDELT2']*\
                                    (float(ind_val[y_index])-header['CRPIX2'])
                            if 'XPOS' in variables_to_transfer:
                                Tirific_Template[f'XPOS{ext}'] = \
                                    combine_string(header['CRVAL1']+\
                                        header['CDELT1']*(float(ind_val[4])-\
                                        header['CRPIX1'])/np.cos(np.radians(ypos_df)),length=nur)
                            if '# XPOS_ERR' in variables_to_transfer and Trigger:
                                Tirific_Template[f'# XPOS_ERR{ext}'] = \
                                    combine_string((float(ind_val[6].strip(','))\
                                        -header['CRPIX1'])*header['CDELT1'],length=nur)
                            else:
                                Tirific_Template[f'# XPOS_ERR{ext}'] = \
                                    combine_string(0.,length=nur)
                            if 'YPOS' in variables_to_transfer:
                                Tirific_Template[f'YPOS{ext}'] =\
                                    combine_string(ypos_df,length=nur)
                            if '# YPOS_ERR' in variables_to_transfer and Trigger:
                                Tirific_Template[f'# YPOS_ERR{ext}'] = \
                                    combine_string((float(ind_val[9].strip(','))\
                                    -header['CRPIX2'])*header['CDELT2'],length=nur)

                            else:
                                Tirific_Template[f'# YPOS_ERR{ext}'] =  \
                                   combine_string(0.,length=nur)
                        if ind_val[0] == 'vsys':
                            if 'VSYS' in variables_to_transfer:
                                Tirific_Template[f'VSYS{ext}'] = \
                                    combine_string(ind_val[2],length=nur)
                            if '# VSYS_ERR' in variables_to_transfer and Trigger:
                                Tirific_Template[f'# VSYS_ERR{ext}'] = \
                                    combine_string(ind_val[4],length=nur)
                            else:
                                Tirific_Template[f'# VSYS_ERR{ext}'] =  \
                                   combine_string(0.,length=nur)
                        if ind_val[0] == 'r_w':
                            r_warp = float(ind_val[3])
                            if Trigger:
                                r_warp_err= float(ind_val[5])
                            else:
                                r_warp_err= 0.
                        if len(ind_val) > 3:
                            if ind_val[2] == 'wphim:':
                                phi_m = float(ind_val[3])
                                if Trigger:
                                    phi_m_err = float(ind_val[5])
                                else:
                                    phi_m_err = 0.

                            if ind_val[2] == 'welm:':
                                el_m=float(ind_val[3])
                                if Trigger:
                                    el_m_err = float(ind_val[5])
                                else:
                                    el_m_err = 0.
                        if ind_val[1] == 'eps:':
                            calc_incl_df = float(np.degrees(np.arccos(np.sqrt(((1-float(ind_val[2]))**2-0.2**2)/0.96)))+2.)

        if phi_m != 0. and 'PA' in variables_to_transfer:
            phi_0= phi_m/(float(rmax)-r_warp)**2.
            if f'RADI{ext}' in Tirific_Template:
                radius =  np.array([float(x) for x in Tirific_Template[f'RADI{ext}'].split()],dtype=float)
            else:
                radius = np.array([float(x) for x in Tirific_Template[f'RADI'].split()],dtype=float)
            pa = np.array([float(x) for x in Tirific_Template[f'PA{ext}'].split()],dtype=float)
            for i,ring in enumerate(radius):
                if float(ring) > float(r_w):
                    pa[i] = pa[i]+phi_0*(ring-r_w)**2.
            Tirific_Template[f'PA{ext}'] = combine_string(pa)
        if 'INCL' in variables_to_transfer:
            if el_m != 0 :
                incl_out = []
                elm_0=el_m/(float(rmax)-r_w)**2.
                if f'RADI{ext}' in Tirific_Template:
                    radius =  np.array([float(x) for x in Tirific_Template[f'RADI{ext}'].split()],dtype=float)
                else:
                    radius = np.array([float(x) for x in Tirific_Template[f'RADI'].split()],dtype=float)
                incl = np.array([float(x) for x in Tirific_Template[f'INCL{ext}'].split()],dtype=float)
                for i,ring in enumerate(radius):
                    if float(ring) > float(r_w):
                        incl[i] = incl[i]+elm_0*(ring-r_w)**2
                Tirific_Template[f'INCL{ext}'] = combine_string(incl)
            else:
                if Tirific_Template[f'INCL{ext}'] == Original_INCL:
                    Tirific_Template[f'INCL{ext}'] = combine_string(calc_incl_df,length=nur)
                    if '# INCL_ERR' in variables_to_transfer:
                        Tirific_Template[f'# INCL_ERR{ext}'] =  \
                            combine_string(0.,length=nur)

    except FileNotFoundError:
        print("No DF file present")


def combine_string(value,length = 0):
    if length == 0:
        try:
            for x in value:
                pass
        except TypeError:
            value=[value]
        string = " ".join([str(x) for x in value])
    else:
        string = " ".join([str(x) for x in np.full(int(length), value)])
    return string
def load_RC_to_Template(Tirific_Template,filename='Rotcur.rc',disk_to_load_to = 1,\
                        variables_to_transfer=['RADI','VROT','VROT_ERR',\
                                        'INCL','INCL_ERR','PA','PA_ERR',\
                                        'XPOS','XPOS_ERR','YPOS','YPOS_ERR',\
                                        'VSYS','VSYS_ERR'], header =None):
    # first read the rotcur file
    column_translation_dict = {'RADI': 0,
                            'RING_WIDTH':1 ,
                            'VSYS': 2,
                            '# VSYS_ERR': 3,
                            'VROT': 4,
                            '# VROT_ERR': 5,
                            'VRAD': 6,
                            '# VRAD_ERR': 7,
                            'PA': 8,
                            '# PA_ERR': 9,
                            'INCL': 10,
                            '# INCL_ERR': 11,
                            'XPOS': 12,
                            '# XPOS_ERR': 13,
                            'YPOS': 14,
                            '# YPOS_ERR': 15,
                            'NPTS': 16,
                            'SIG':17}
    #rc_upload=np.loadtxt(filename, skiprows=11)
    rc_upload=np.genfromtxt(filename, skip_header=11).transpose()

    if disk_to_load_to == 1:
        ext=''
        Tirific_Template[f'NUR'] = len(rc_upload[0])
    else:
        ext = f'_{int(disk_to_load_to)}'
        if f'RADI' in variables_to_transfer and f'RADI{ext}' not in Tirific_Template:
            prev_key = 'Empty'
            for i,key in Tirific_Template:
                if key == f'VROT{ext}':
                    Tirific_Template.insert(f'RADI{ext}',prev_key,'10.')
                    break
                else:
                    prev_key = copy.deepcopy(key)

    for variable in variables_to_transfer:
        if variable[-4:] == '_ERR':
            if variable[:2] != '# ':
                variable = f'# {variable}'
            if f'{variable}{ext}' not in Tirific_Template:
                stripped = variable[2:-4]
                Tirific_Template.insert(f'{stripped}{ext}',f'{variable}{ext}','10.')
        if variable not in column_translation_dict:
            print(f'We do not know how to translate the variable {variable} to a RotCur column. Please use TiRiFiC input variables available in RotCur')
            continue
        if variable in ['XPOS','# XPOS_ERR','YPOS','# YPOS_ERR'] and not header:
            raise InputError('RotCur prints its central position in pixel values, we can not translate this without the header of the velocity field')
        if variable == 'XPOS':
            extra_var = header['CRVAL2']+header['CDELT2']*np.array(rc_upload[column_translation_dict['YPOS']])
            selection = header['CRVAL1']+header['CDELT1']*np.array(rc_upload[column_translation_dict['XPOS']])/np.cos(np.radians(extra_var))
        elif variable == '# XPOS_ERR':
            selection = abs(header['CDELT1']*np.array(rc_upload[column_translation_dict[variable]]))
        elif variable == '# YPOS_ERR':
            selection = abs(header['CDELT2']*np.array(rc_upload[column_translation_dict[variable]]))
        elif variable == 'YPOS':
            selection = header['CRVAL2']+header['CDELT2']*np.array(rc_upload[column_translation_dict[variable]])
        else:
            selection = rc_upload[column_translation_dict[variable]]

        Tirific_Template[f'{variable}{ext}'] = " ".join(str(e) for e in selection)











if __name__ == '__main__':
    main()
