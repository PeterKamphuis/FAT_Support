#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


from collections import OrderedDict #used in Proper_Dictionary
from datetime import datetime
from inspect import getframeinfo,stack
from scipy.optimize import curve_fit
from scipy import ndimage
from scipy import interpolate
from astropy.wcs import WCS
from astropy.io import fits
from omegaconf import OmegaConf
from pyFAT_astro.config.defaults import defaults
from pyFAT_astro.Support.support_functions import convertskyangle,columndensity,setup_configuration
from pyFAT_astro.Support.fat_errors import BadCatalogueError
import os
import signal
import traceback
import numpy as np
import copy
import warnings
import time
import sys

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.patches import Ellipse


import re
import subprocess
class SupportRunError(Exception):
    pass
class SmallSourceError(Exception):
    pass
class FileNotFoundError(Exception):
    pass
class ModelError(Exception):
    pass


# A class of ordered dictionary where keys can be inserted in at specified locations or at the end.
class Proper_Dictionary(OrderedDict):
    def insert(self, existing_key, new_key, key_value):
        done = False
        if new_key in self:
            self[new_key] = key_value
            done = True
        else:
            new_orderded_dict = self.__class__()
            for key, value in self.items():
                new_orderded_dict[key] = value
                if key == existing_key:
                    new_orderded_dict[new_key] = key_value
                    done = True
            if not done:
                new_orderded_dict[new_key] = key_value
                done = True
                print(
                    "----!!!!!!!! YOUR new key was appended at the end as you provided a non-existing key to add it after!!!!!!---------")
            self.clear()
            self.update(new_orderded_dict)

        if not done:
            print("----!!!!!!!!We were unable to add your key!!!!!!---------")



def analyze(program,database,Input_Parameters,basename='EMPTY',\
        main_directory='EMPTY'):
    LVHIS = False
    if database == 'LVHIS':
        LVHIS = True
    if not os.path.exists(f'''{main_directory}/{basename}'''):
        os.system(f'''mkdir {main_directory}/{basename}''')

    plot_overview(Input_Parameters,\
        filename=f'''{main_directory}/{basename}/Release_All'''\
        ,program=program,database=database)
    print(database)
    RCs = get_all_RCs(Input_Parameters,program=program,database=database)


    print(f'''{main_directory}/{basename}/RCs''')
    plot_RCs(RCs,database=database,\
        filename=f'''{main_directory}/{basename}/RCs''')

def compare(Input_Parameters,programs= None):




    print(Input_Parameters)

def get_all_RCs(Input_Parameters,program=None,database=None):
    RCs = {}
    for galaxy in Input_Parameters:
        if Input_Parameters[galaxy]['Model']['Cat_Type'] != database:
            continue
        print(f'Adding {galaxy} to the RCs')
        RCShape = Input_Parameters[galaxy]['Model']['RCShape']
        if RCShape not in RCs:
            RCs[RCShape] = {'MODEL': {'RADIUS':Input_Parameters[galaxy]\
                ['Model']['RADI'],
                'RC':Input_Parameters[galaxy]['Model']['VROT'],
                'DISTANCE': [Input_Parameters[galaxy]['Model']['Distance']],
                'STATUS': -1}}
            if float(RCs[RCShape]['MODEL']['DISTANCE'][0]) == 0.:
                RCs[RCShape]['MODEL']['DISTANCE'] = [0.5]

        if 'RADI_2' in Input_Parameters[galaxy]['Model']:
            RCs[RCShape]['MODEL_2'] = {'RADIUS':Input_Parameters[galaxy]\
                ['Model']['RADI_2'],
                'RC':Input_Parameters[galaxy]['Model']['VROT_2'],
                'DISTANCE': [Input_Parameters[galaxy]['Model']['Distance']],
                'STATUS': -1}
            if float(RCs[RCShape]['MODEL_2']['DISTANCE'][0]) == 0.:
                RCs[RCShape]['MODEL_2']['DISTANCE'] = [0.5]

        RCs[RCShape][galaxy] = {'RADIUS':Input_Parameters[galaxy][program]['RADI']\
            , 'RC':Input_Parameters[galaxy][program]['VROT'],
            'DISTANCE':[Input_Parameters[galaxy][program]['Distance']],\
            'STATUS': int(Input_Parameters[galaxy][program]['Result'][0])}
        if float(RCs[RCShape][galaxy]['DISTANCE'][0]) == 0.:
                    RCs[RCShape][galaxy]['DISTANCE'] = [0.5]
    return RCs

def obtain_individual_parameters(dict,LVHIS = False, Model =False):
    if Model:
        tmp = load_tirific(f"{dict['Directory']}/ModelInput.def"\
                            , LVHIS = LVHIS)
    else:
        tmp = load_tirific(get_model_filename(dict), LVHIS = LVHIS)
    for key in tmp:
        dict[key] = tmp[key]
    if Model:
        dict['FLUX'] = get_flux_values(dict['Directory'],input=True)
    else:
        if dict['Result'][0] == 2:
            dict['FLUX'] = get_flux_values(dict['Directory'])
        else:
            dict['FLUX'] = [float('NaN'),float('NaN'),float('NaN')]

def get_model_filename(galaxy_dictionary):
    if galaxy_dictionary['Result'][0] == 2 or \
        (galaxy_dictionary['Result'][0] == 1 and os.path.isfile(f'{galaxy_dictionary["Directory"]}/Finalmodel/Finalmodel.def')):
        filename = f"{galaxy_dictionary['Directory']}/Finalmodel/Finalmodel.def"
    elif galaxy_dictionary['Result'][0] == 1 :
        filename = f"{galaxy_dictionary['Directory']}/No_Warp/No_Warp.def"
    else:
        filename = 'EMPTY'
    return filename


def obtain_model_settings(galaxy,dict,LVHIS=False,LVHIS_Names = None):
    if not LVHIS:
        diameter_in_beams, SNR, RCshape = get_name_info(galaxy)
    else:
        diameter_in_beams = [float(dict['RADI'][-1])/float(dict['BMAJ'][0])*2.]
        if 'RADI_2' in dict:
            diameter_in_beams.append(float(dict['RADI_2'][-1])\
                            /float(dict['BMAJ'][0])*2.)
        else:
            diameter_in_beams.append(float(dict['RADI'][-1])\
                /float(dict['BMAJ'][0])*2.)

        SNR = dict['FLUX'][2]/dict['NOISE']
        RCshape = LVHIS_Names[galaxy]
    if 'SBR' in dict:
        RHI =get_RHI(sbr_profile=[dict[f'SBR'],dict[f'SBR_2']],radi=dict[f'RADI'],\
            systemic=dict[f'VSYS'][0],distance=dict['Distance'])
    else:
        RHI = get_diff_rmax(dict[f'RADI'], dict['RADI'],dict['BMAJ'][0])
        if 'RADI_2' in dict:
            RHI.append(get_diff_rmax(dict[f'RADI_2'], dict['RADI'],dict['BMAJ'][0]))
    dict['R_HI'] = RHI
    dict['Size_in_Beams']= diameter_in_beams
    dict['SNR']= SNR
    dict['RCShape']= RCshape

def obtain_overview_parameters(Input_Parameters,galaxy,key,inp_index,tmp_inp_cat\
        ,tmp_config, LVHIS =False):
    if galaxy not in Input_Parameters:
        Input_Parameters[galaxy] = {'pyFAT': {}, 'GDL':{}, 'Model': None}
    Input_Parameters[galaxy][key]['Directory'] =\
        f"{tmp_config['MAIN_DIRECTORY']}{galaxy}"
    #Some general info
    Input_Parameters[galaxy][key]['Distance'] = tmp_inp_cat['DISTANCE'][inp_index]
    Input_Parameters[galaxy][key]['Cat_Type'] = 'Database'
    if LVHIS:
        Input_Parameters[galaxy][key]['Cat_Type'] = 'LVHIS'
    Input_Parameters[galaxy][key]['Cube'] =f"{tmp_inp_cat['CUBENAME'][inp_index]}.fits"
    if 'Gauss' in Input_Parameters[galaxy][key]['Cube']:
        Input_Parameters[galaxy][key]['Corruption'] = 'Gauss'
    elif 'UC' in Input_Parameters[galaxy][key]['Cube']:
        Input_Parameters[galaxy][key]['Corruption'] = 'Uncorrupted'
    elif 'CS' in Input_Parameters[galaxy][key]['Cube']:
        Input_Parameters[galaxy][key]['Corruption'] = 'CASA'
    else:
        if LVHIS:
            Input_Parameters[galaxy][key]['Corruption'] = ['ROTCUR','DISK_FIT']
        else:
            Input_Parameters[galaxy][key]['Corruption'] = 'Real'
    cube = fits.open(f"{Input_Parameters[galaxy][key]['Directory']}/{Input_Parameters[galaxy][key]['Cube']}")
    fact = 1.
    try:
        if cube[0].header['CUNIT3'].upper() == 'M/S':
            fact = 1000.
    except KeyError:
        if cube[0].header['CDELT3'] > 500:
            fact = 1000.
    Input_Parameters[galaxy][key]['CHANNEL_WIDTH'] = cube[0].header['CDELT3']/fact
    try:
        Input_Parameters[galaxy][key]['NOISE'] = cube[0].header['FAT_NOISE']
    except KeyError:
        Input_Parameters[galaxy][key]['NOISE'] = np.std(cube[0].data[1,:,:])
    cube.close()

def obtain_parameters(Input_Files,missing_links=False):

    Input_Parameters = {}

    for key in Input_Files:

        if key == 'GDL':
            GDL = True
        else:
            GDL = False
        keyinputs = [x for x in Input_Files[key] if x != 'version']
        for key2 in keyinputs:
            print(key2)

            LVHIS = False
            if key2 == 'LVHIS':
                LVHIS = True

            tmp_config = load_config_file(f"{Input_Files[key][key2]['dir']}/{Input_Files[key][key2]['config']}")
            if LVHIS:
                LVHIS_Names = load_LVHIS_Name(tmp_config["MAIN_DIRECTORY"])
            else:
                LVHIS_Names = None
            tmp_inp_cat = load_input_catalogue(tmp_config['INPUT_CATALOGUE'],\
                                                GDL=GDL)
            tmp_out_cat = load_output_catalogue(
                            tmp_config['OUTPUT_CATALOGUE'],binary=GDL)
            if missing_links and not GDL:
                fix_links(tmp_config,tmp_out_cat)
            for i,galaxy in enumerate(tmp_out_cat['DIRECTORY_NAME']):
                if galaxy not in tmp_inp_cat['DIRECTORYNAME']:
                    print(f'''The galaxy {galaxy} is not in the input catalogue.
This should not be possible, Please fix your input ''')
                else:
                    inp_index = tmp_inp_cat['DIRECTORYNAME'].index(galaxy)
                obtain_overview_parameters(Input_Parameters,galaxy,key,inp_index,\
                    tmp_inp_cat,tmp_config,LVHIS = LVHIS)

                Input_Parameters[galaxy][key]['Result'] = get_result(\
                    tmp_out_cat,tmp_config,i,galaxy,\
                    binary = GDL)

                obtain_individual_parameters(\
                    Input_Parameters[galaxy][key],LVHIS=LVHIS)
                #And get the model
                if Input_Parameters[galaxy]['Model'] is None:
                    Input_Parameters[galaxy]['Model'] = {'Directory': \
                        f"{tmp_config['MAIN_DIRECTORY']}{galaxy}"}
                    obtain_overview_parameters(Input_Parameters,galaxy,'Model',inp_index,\
                        tmp_inp_cat,tmp_config,LVHIS = LVHIS)
                    obtain_individual_parameters(\
                        Input_Parameters[galaxy]['Model'],LVHIS=LVHIS, Model=True)
                    obtain_model_settings(galaxy,Input_Parameters[galaxy]['Model']\
                            ,LVHIS=LVHIS,LVHIS_Names = LVHIS_Names)
# Load the output parameters
    #Obtain the model parameters

    return Input_Parameters
def warn_with_traceback(message, category, filename, lineno, file=None, line=None):
    log = file if hasattr(file,'write') else sys.stderr
    traceback.print_stack(file=log)
    log.write(warnings.formatwarning(message, category, filename, lineno, line))

def add_deltas(Input_Parameters):
    warnings.showwarning = warn_with_traceback
    for galaxy in Input_Parameters:
        if 'Model' not in Input_Parameters[galaxy]:
            print(f'''We have no model values for the galaxy {galaxy}
skipping this galaxy''')
        keys = [x for x in Input_Parameters[galaxy] if x != 'Model']
        for key in keys:
            try:
                deltas = galaxy_deltas(Input_Parameters[galaxy][key],\
                        Input_Parameters[galaxy]['Model'])
            except Exception as e:
                print(f"The deltas failed in {galaxy}")
                traceback.print_exception(type(e),e,e.__traceback__)
                exit()

def plot_RCs(RCs, filename='RCs',database= None):

    LVHIS = False
    if database == 'LVHIS':
        LVHIS=True

    labelfont= {'family':'Times New Roman',
                'weight':'normal',
                'size':24}
    plt.rc('font',**labelfont)

    plotsize = 8
    no_plots = len([x for x in RCs])
    length = int(np.ceil(np.sqrt(no_plots)))
    #length = no_plots

    plt.figure(89,figsize=(plotsize*length,plotsize*length),dpi=300,facecolor = 'w', edgecolor = 'k')

    #this runs counter to figsize. How can a coding language be this illogical?
    gs = gridspec.GridSpec(length,length )
    gs.update(wspace=0.25, hspace=0.25)
    if LVHIS:
        linew = 3
    else:
        linew = 1
    for i,key in enumerate(RCs):

        print(f'Plotting the galaxy RC shape {key}')
        ax = plt.subplot(gs[i])
        ax.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=True,         # ticks along the top edge are off
        labelbottom=True,
        right = True,
        left= True,
        labelleft = True)
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(4)

            # increase tick width
        ax.tick_params(width=4)
        failed = 0
        tot = 0
        plot = True
        if LVHIS:
            galaxy = [x for x in RCs[key] if x not in ['MODEL','MODEL_2']]
            if RCs[key][galaxy[0]]['STATUS'] == 0.:
                plot = False

        if plot:
            for indi in RCs[key]:

                print(f'Plotting the actual galaxy {indi}')


                if indi not in ['MODEL','MODEL_2']:
                    tot += 1

                kpcradii = np.array(convertskyangle({'OUTPUTLOG': None,\
                    'TIMING':None,'DEBUG': False, 'DISTANCE':float(RCs[key][indi]['DISTANCE'][0]) }\
                    ,RCs[key][indi]['RADIUS']))
                print(f''' Plotting the RC {RCs[key][indi]['RC']} with:
    radi (arcsec) = {RCs[key][indi]['RADIUS']}
    radi (kpc) = {kpcradii}
    distance = {float(RCs[key][indi]['DISTANCE'][0])}''')

                if float(RCs[key][indi]['DISTANCE'][0]) == 0. :
                    print(f'There is something wrong with the distance in {key}')
                    exit()
                if kpcradii[2] < 0. :
                    print(f'There is something wrong with the distance in {key}')
                    exit()
                if indi == 'MODEL':
                    ax.plot(kpcradii, RCs[key][indi]['RC'], 'b',zorder= 2)
                    ax.plot(kpcradii, RCs[key][indi]['RC'], 'bo',zorder=2)
                elif indi == 'MODEL_2':
                    #print(np.isnan(np.array(RCs[key][indi]['RC'],dtype=float)),all(np.isnan(np.array(RCs[key][indi]['RC'],dtype=float))))
                    if not np.sum(RCs[key][indi]['RC']) == 0.:
                        ax.plot(kpcradii, RCs[key][indi]['RC'], 'r',zorder= 2)
                        ax.plot(kpcradii, RCs[key][indi]['RC'], 'ro',zorder=2)
                else:
                    if len(kpcradii) > len(RCs[key][indi]['RC']):
                        kpcradii=kpcradii[:len(RCs[key][indi]['RC'])]
                    if RCs[key][indi]['STATUS'] == 0:
                        failed += 1
                    elif RCs[key][indi]['STATUS'] == 1:
                        #ymin, ymax = ax.get_ylim()
                        ax.plot(kpcradii, RCs[key][indi]['RC'], 'k--',zorder= 1,linewidth=linew, alpha =0.5)
                        if LVHIS:
                            ax.plot(kpcradii, RCs[key][indi]['RC'], 'ko',zorder= 1,linewidth=linew, alpha =0.5)
                        #ax.set_ylim(ymin,ymax)
                    else:
                        ax.plot(kpcradii, RCs[key][indi]['RC'], 'k',zorder= 1 ,linewidth=linew, alpha =0.75)
                        if LVHIS:
                            ax.plot(kpcradii, RCs[key][indi]['RC'], 'ko',zorder= 1 ,linewidth=linew, alpha =0.75)
            ax.set_xlabel('Radius (kpc)', **labelfont)
            ax.set_ylabel('V$_{rot}$ (km s$^{-1}$)', **labelfont)
            if RCs[key][indi]['STATUS'] == 1 and LVHIS:
                print(f"Trying a marker")
                ax.scatter(0.95,0.95,marker='*',color='k', s=237,transform=ax.transAxes)
        else:
            ax.text(0.5,0.5,f'The galaxy {key} failed to fit ', transform=ax.transAxes,horizontalalignment= 'center', verticalalignment='top')
            ax.axis('off')

        ax.set_title(key)

        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin,ymax+(ymax-ymin)/10.)
        if not LVHIS:
            ax.text(0.95,0.95,f'Out of {tot} galaxies, {failed} failed to fit. ', transform=ax.transAxes,horizontalalignment= 'right', verticalalignment='top')



    if LVHIS:
        version= 'LVHIS'
    else:
        version='Database'

    plt.savefig(f'{filename}_RC_{version}.png', bbox_inches='tight')
    print(f'{filename}_RC_{version}.png')
    plt.close()


def analyze_timing(Database_Directory, Input_Parameters,\
                 database = None, program=None,basename = 'Analysis_Output'):
    '''This function only works for the pyFAT runs not GDL runs'''
    Individual_Times = read_timing_results(f'{Database_Directory}/Timing_Result.txt')

    plot_time_duration(f'{Database_Directory}/{basename}/Time_Statistic.png', \
       Individual_Times, Input_Parameters, database = database, program=program)


def create_patch(patch,order,Input_Parameters,program=None):

    x= np.array(get_all_specific(patch[1],order,Input_Parameters,\
        program=program)[0])
    y= np.array(get_all_specific(patch[2],order,Input_Parameters,\
        program=program)[0])

    if patch[0] == 'Ellipse':
        patch_fig = Ellipse(xy=[np.nanmean(x),np.nanmean(y)],\
                        width=np.nanstd(x) ,height=np.nanstd(y), angle=0,\
                        edgecolor='none', alpha=0.6, lw=4, facecolor='k', \
                        hatch = '////', zorder=-1)
    return patch_fig

def plot_time_duration(filename, Times_Dictionary,Input_Parameters, \
                            database = None,program =None ):

    duration= []
    prep_duration = []
    fit_duration = []
    finish_duration = []
    beams = []
    snr = []
    inclination= []
    flux = []
    RHI = []
    status = []
    #order = [x for x in Input_Parameters if \
    #        Input_Parameters[x][program]['Cat_Type'] == database]
    order = []

    for galaxy in  Times_Dictionary:
        duration.append(Times_Dictionary[galaxy]['Summed_Time'][2])
        prep_duration.append(Times_Dictionary[galaxy]['Ini_Time'][2])
        fit_duration.append(Times_Dictionary[galaxy]['Fit_Time'][2])
        finish_duration.append(Times_Dictionary[galaxy]['Shaker_Time'][2])
        order.append(galaxy)
    beams = get_all_specific('Size_in_Beams',order,Input_Parameters,\
                                program='Model')[0]
    snr = get_all_specific('SNR',order,Input_Parameters,\
                                program='Model')[0]
    inclination = get_all_specific('INCL',order,Input_Parameters,\
                                program='Model')[0]
    flux =  get_all_specific('FLUX',order,Input_Parameters,\
                                program='Model')[0]
    RHI = get_all_specific('R_HI',order,Input_Parameters,\
                                program='Model')[0]
    status = get_all_specific('Result',order,Input_Parameters,\
                                program=program)[0]

    plot_assembly = {'PLOT_1':{
                        'WINDOW_0': {'LOCATION': 0,
                                     'Y': [duration, 'Duration (hrs)' ],
                                     'X': [beams,'Maj Ax Beams'],
                                     'NO_MEAN': False,
                                     #'PATCH':  Ellipse(xy=[np.nanmean(np.array([x[0] for x in deltas['XPOS']])),\
                                    #                    np.nanmean(np.array([x[0] for x in deltas['YPOS']]))],\
                                    #           width=np.nanstd(np.array([x[0] for x in deltas['XPOS']])) ,\
                                    #           height=np.nanstd(np.array([x[0] for x in deltas['YPOS']])), angle=0,\
                                    #        edgecolor='none', alpha=0.6, lw=4, facecolor='k', hatch = '////', zorder=-1)
                                    },
                        'WINDOW_1': {'LOCATION': 1,
                                     'Y': [duration, 'Duration (hrs)' ],
                                     'X': [snr,'SNR'] #error from https://www.wolframalpha.com/widgets/view.jsp?id=8ac60957610e1ee4894b2cd58e753
                                     },
                        'WINDOW_2': {'LOCATION': 3,
                                     'Y': [prep_duration, 'Preparation Duration (hrs)' ],
                                     'X': [beams ,'Maj Ax Beams'],
                                     'NO_MEAN': False,
                                     #'PATCH':  Ellipse(xy=[np.nanmean(np.array([x[0] for x in deltas['XPOS']])),\
                                    #                    np.nanmean(np.array([x[0] for x in deltas['YPOS']]))],\
                                    #           width=np.nanstd(np.array([x[0] for x in deltas['XPOS']])) ,\
                                    #           height=np.nanstd(np.array([x[0] for x in deltas['YPOS']])), angle=0,\
                                    #        edgecolor='none', alpha=0.6, lw=4, facecolor='k', hatch = '////', zorder=-1)
                                    },
                        'WINDOW_3': {'LOCATION': 4,
                                     'Y': [prep_duration, 'Preparation Duration (hrs)' ],
                                     'X': [snr,'SNR'] #error from https://www.wolframalpha.com/widgets/view.jsp?id=8ac60957610e1ee4894b2cd58e753
                                     },
                        'WINDOW_4': {'LOCATION': 6,
                                     'Y': [fit_duration, 'Fit Duration (hrs)' ],
                                     'X': [beams ,'Maj Ax Beams'],
                                     'NO_MEAN': False,
                                     #'PATCH':  Ellipse(xy=[np.nanmean(np.array([x[0] for x in deltas['XPOS']])),\
                                    #                    np.nanmean(np.array([x[0] for x in deltas['YPOS']]))],\
                                    #           width=np.nanstd(np.array([x[0] for x in deltas['XPOS']])) ,\
                                    #           height=np.nanstd(np.array([x[0] for x in deltas['YPOS']])), angle=0,\
                                    #        edgecolor='none', alpha=0.6, lw=4, facecolor='k', hatch = '////', zorder=-1)
                                    },
                        'WINDOW_5': {'LOCATION': 7,
                                     'Y': [fit_duration, 'Fit Duration (hrs)' ],
                                     'X': [snr,'SNR'] #error from https://www.wolframalpha.com/widgets/view.jsp?id=8ac60957610e1ee4894b2cd58e753
                                     },
                        'WINDOW_6': {'LOCATION': 9,
                                     'Y': [finish_duration, 'Finish Duration (hrs)' ],
                                     'X': [beams ,'Maj Ax Beams'],
                                     'NO_MEAN': False,
                                     #'PATCH':  Ellipse(xy=[np.nanmean(np.array([x[0] for x in deltas['XPOS']])),\
                                    #                    np.nanmean(np.array([x[0] for x in deltas['YPOS']]))],\
                                    #           width=np.nanstd(np.array([x[0] for x in deltas['XPOS']])) ,\
                                    #           height=np.nanstd(np.array([x[0] for x in deltas['YPOS']])), angle=0,\
                                    #        edgecolor='none', alpha=0.6, lw=4, facecolor='k', hatch = '////', zorder=-1)
                                    },
                        'WINDOW_7': {'LOCATION': 10,
                                     'Y': [finish_duration, 'Finish Duration (hrs)' ],
                                     'X': [snr,'SNR'] #error from https://www.wolframalpha.com/widgets/view.jsp?id=8ac60957610e1ee4894b2cd58e753
                                     },
                        'WINDOW_8': {'LOCATION': 2,
                                     'Y': [duration, 'Duration (hrs)' ],
                                     'X': [inclination,'Central Inclination'] #error from https://www.wolframalpha.com/widgets/view.jsp?id=8ac60957610e1ee4894b2cd58e753
                                     },
                        'WINDOW_9': {'LOCATION': 5,
                                     'Y': [duration, 'Duration (hrs)' ],
                                     'X': [status,'Fit status'] #error from https://www.wolframalpha.com/widgets/view.jsp?id=8ac60957610e1ee4894b2cd58e753
                                     },
                        'WINDOW_10': {'LOCATION': 8,
                                     'Y': [duration, 'Duration (hrs)' ],
                                     'X': [flux,'Total Flux'] #error from https://www.wolframalpha.com/widgets/view.jsp?id=8ac60957610e1ee4894b2cd58e753
                                     },
                        'WINDOW_11': {'LOCATION': 11,
                                     'Y': [duration, 'Duration (hrs)' ],
                                     'X': [RHI,'R$_{HI}$'] #error from https://www.wolframalpha.com/widgets/view.jsp?id=8ac60957610e1ee4894b2cd58e753
                                     },

                                     }}

    for plot in plot_assembly:
        labelfont= {'family':'Times New Roman',
                    'weight':'normal',
                    'size':26}
        plt.rc('font',**labelfont)
        plt.figure(89,figsize=(18,24),dpi=300,facecolor = 'w', edgecolor = 'k')

        #this runs counter to figsize. How can a coding language be this illogical?
        gs = gridspec.GridSpec(4,3)
        gs.update(wspace=0.3, hspace=0.3)


        for i,key in enumerate(plot_assembly[plot]):
            if 'NO_MEAN' in plot_assembly[plot][key]:
                nomean = plot_assembly[plot][key]['NO_MEAN']
            else:
                nomean = False

            print(plot_assembly[plot][key]['X'][1],plot_assembly[plot][key]['Y'][1])

            ax,legend_items = make_plot(plot_assembly[plot][key]['X'][0],\
                           plot_assembly[plot][key]['Y'][0],\
                           xlabel = plot_assembly[plot][key]['X'][1],\
                           ylabel = plot_assembly[plot][key]['Y'][1],\
                           location = gs[plot_assembly[plot][key]['LOCATION']],\
                           No_Mean = nomean,no_error=[True, True])
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(4)

                    # increase tick width
                ax.tick_params(width=4)
            if 'PATCH' in plot_assembly[plot][key]:
                patch_fig = create_patch(Input_Parameters, plot_assembly[plot][key]['PATCH'])
                ax.add_patch(patch_fig)
            if i == 1:

                # make a legend
                labelfont= {'family':'Times New Roman',
                            'weight':'normal',
                            'size':18}
                plt.rc('font',**labelfont)
                chartBox = ax.get_position()
                ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*1.0, chartBox.height])
                ax.legend(handles=legend_items,loc='upper left', bbox_to_anchor=(1.25, 1.0), shadow=True, ncol=1)
                # Beams vs delta inclination
                labelfont= {'family':'Times New Roman',
                            'weight':'normal',
                            'size':26}
                #AND below it a corruption size indicator
                #unique_corruption = np.unique(deltas['CORRUPTION'])
                #if unique_corruption.size > 1:

                #    sizes


                plt.rc('font',**labelfont)




        plt.savefig(filename, bbox_inches='tight')
        plt.close()


def read_timing_results(filename):
    with open(filename) as file:
        lines=file.readlines()
    if lines[0].split()[0] == 'Timing' and  lines[0].split()[1] == 'results':
        Timing_Results = read_section_times(lines)
    elif lines[0].split()[0] == 'This' and  lines[0].split()[1] == 'file':
        Timing_Results = read_individual_times(lines)
    else:
        raise FileNotFoundError(f'The file {filename} is not recognized as a timing reselt file.')

    return Timing_Results

def read_section_times(lines):
    individual_times = {}
    current_galaxy = 'None'

    for line in lines:
        individual_words = line.split()
        if individual_words[-1][-1] == '.':
            individual_words[-1] = individual_words[-1][:-1]

        if individual_words[0] == 'The' and individual_words[1] == 'galaxy':
            if current_galaxy != 'None':
                summed_time = [initialize_time[0],run_time[1],\
                    np.nansum([initialize_time[2],prep_time[2],sofia_time[2],ig_time[2],\
                        fit_time[2]])]
                individual_times[current_galaxy] = {'Name': current_galaxy,
                        'Run_Time': run_time, 'Summed_Time': summed_time,  \
                        'Ini_Time': initialize_time,'Prep_Time': prep_time ,\
                        'Sofia_Time': sofia_time, 'IG_Time': ig_time,\
                        'Fit_Time':fit_time,'Shaker_Time': ts_time}

            full_directory = individual_words[4]
            current_galaxy = full_directory.split('/')[-2]
        elif individual_words[0] == 'Full' and individual_words[2] == 'time':
            run_time = [f'{individual_words[4]} { individual_words[5]}',\
                        f'{individual_words[7]} { individual_words[8]}']
            run_time.append(get_duration(run_time[0],run_time[1]))
        elif (individual_words[0] == 'Inititalization' or \
            individual_words[0] == 'Initialization')    and individual_words[1] == 'from':
            initialize_time = [f'{individual_words[2]} { individual_words[3]}',\
                        f'{individual_words[5]} { individual_words[6]}']
            initialize_time.append(get_duration(initialize_time[0],initialize_time[1]))
        elif individual_words[0] == 'Preparation' and individual_words[1] == 'time':
            prep_time = [f'{individual_words[3]} { individual_words[4]}',\
                        f'{individual_words[6]} { individual_words[7]}']
            prep_time.append(get_duration(prep_time[0],prep_time[1]))
        elif individual_words[3] == 'Sofia':
            sofia_time = [f'{individual_words[5]} { individual_words[6]}',\
                        f'{individual_words[8]} { individual_words[9]}']
            sofia_time.append(get_duration(sofia_time[0],sofia_time[1]))
        elif individual_words[3] == 'Initial':
            ig_time = [f'{individual_words[5]} { individual_words[6]}',\
                        f'{individual_words[8]} { individual_words[9]}']
            ig_time.append(get_duration(ig_time[0],ig_time[1]))
        elif individual_words[2] == 'TirShaking':
            ts_time = [f'{individual_words[3]} { individual_words[4]}',\
                        f'{individual_words[6]} { individual_words[7]}']
            ts_time.append(get_duration(ts_time[0],ts_time[1]))
        elif individual_words[3] == 'fitting':
            fit_time = [f'{individual_words[4]} { individual_words[5]}',\
                        f'{individual_words[7]} { individual_words[8]}']
            fit_time.append(get_duration(fit_time[0],fit_time[1]))

    return individual_times



def read_individual_times(lines):
    '''Read the times from the timing file and put in dictionary with duration in hrs'''

    individual_times = {}
    current_galaxy = 'None'
    for line in lines:
        individual_words = line.split()
        if individual_words[-1][-1] == '.':
            individual_words[-1] = individual_words[-1][:-1]

        if individual_words[0] == 'The' and individual_words[1] == 'galaxy':
            if current_galaxy != 'None':
                summed_time = [run_time[0],run_time[1],\
                            np.sum([prep_time[2],fit_time[2],ts_time[2]])]

                individual_times[current_galaxy] = {'Name': current_galaxy,
                        'Run_Time': run_time, 'Summed_Time': summed_time,  \
                        'Ini_Time': ini_time,'Prep_Time': prep_time ,\
                        'Sofia_Time': sofia_time, 'IG_Time': ig_time,\
                        'Fit_Time':fit_time,'Shaker_Time': ts_time}
            full_directory = individual_words[4]

            current_galaxy = full_directory.split('/')[-2]
            start_time = f'{individual_words[-2]} { individual_words[-1]}'
        elif individual_words[0] == 'Finished' and individual_words[1] == 'preparations':
            prep_time = [start_time,f'{individual_words[-2]} { individual_words[-1]}']
            prep_time.append(get_duration(start_time,prep_time[1]))
        elif individual_words[0] == 'Converged' and individual_words[1] == 'to':
            fit_time = [prep_time[1],f'{individual_words[-2]} { individual_words[-1]}']
            fit_time.append(get_duration(prep_time[1],fit_time[1]))
        elif individual_words[0] == 'It' and individual_words[1] == 'finished':
            run_time = [start_time,f'{individual_words[-2]} { individual_words[-1]}']
            run_time.append(get_duration(run_time[0],run_time[1]))
            ts_time = [fit_time[1],f'{individual_words[-2]} { individual_words[-1]}']
            ts_time.append(get_duration(ts_time[0],ts_time[1]))
            ini_time = ['Not Started', 'Not Completed', float('NaN')]
            sofia_time = ['Not Started', 'Not Completed', float('NaN')]
            ig_time =['Not Started', 'Not Completed', float('NaN')]

    return individual_times


def get_duration(start,end):
    '''Obtain the difference between two string times in hrs'''
    try:
        start_object = datetime.strptime(start, '%Y-%m-%d %H:%M:%S.%f')
        end_object = datetime.strptime(end, '%Y-%m-%d %H:%M:%S.%f')
        difference  =end_object - start_object
        diff_hr = difference.days*24+difference.seconds/3600.
    except:
        diff_hr = float('NaN')
    return diff_hr

def get_all_specific(select_parameter,order,Input_Parameters,program=None):
    out = [[],[]]
    if select_parameter[0] == 'D':
        counter =0
        while Input_Parameters[order[counter]][program]['deltas'] == None:
            counter += 1

        models = [x for x in Input_Parameters[order[counter]][program]['deltas']]
        if models[0] == 'ROTCUR' and len(models) == 1:
            models.append('DISK_FIT')
        parameter = select_parameter[1:]
        for galaxy in order:
            for i,model in enumerate(models):
                if parameter == 'CENTRAL':

                    try:
                        x,xerr = Input_Parameters[galaxy][program]['deltas'][model]\
                            ['XPOS']

                        y,yerr = Input_Parameters[galaxy][program]['deltas'][model]\
                            ['YPOS']
                        x = np.array(x,dtype=float)
                        xerr = np.array(xerr,dtype=float)
                        y = np.array(y,dtype=float)
                        yerr = np.array(yerr,dtype=float)
                        out[0].append(np.sqrt(x**2+y**2))
                        out[1].append(np.sqrt(x**2*xerr**2+y**2*yerr)**2/\
                                (x**2+y**2))
                        #print(out)
                    except:
                        out[0].append(float('NaN'))
                        out[1].append(float('NaN'))

                else:
                    try:
                        out[0].append(float(Input_Parameters[galaxy][program]['deltas'][model]\
                            [parameter][0]))
                        out[1].append(float(Input_Parameters[galaxy][program]['deltas'][model]\
                            [parameter][1]))
                    except:
                        out[0].append(float('NaN'))
                        out[1].append(float('NaN'))
    else:
        parameter  = select_parameter
        for galaxy in order:
            if select_parameter in ['Corruption','SNR']:
                out[0].append(Input_Parameters[galaxy][program][select_parameter])
            else:
                print(f'IP = {Input_Parameters[galaxy][program][select_parameter]}')
                out[0].append(float(Input_Parameters[galaxy][program][select_parameter][0]))
                try:
                    out[1].append(float(Input_Parameters[galaxy][program][f'{select_parameter}_ERR'][0]))
                except:
                    try:
                        out[1].append(float(Input_Parameters[galaxy][program][select_parameter][1]))
                    except:
                        out[1].append(float('NaN'))

    return out

def plot_difference(Input_Compare):
    plot_assembly = {'PLOT_1':{
                        'WINDOW_0': {'LOCATION': 0,
                                     'X': ['DXPOS', '$\Delta$ RA (beams)',program ],
                                     'Y': ['DYPOS','$\Delta$ DEC (beams)',program ],
                                     'NO_MEAN': True,
                                     'PATCH': ['Ellipse','DXPOS','DYPOS']
                                     },
                        'WINDOW_1': {'LOCATION': 1,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DCENTRAL','$\Delta$ Central (beams)',program ], #error from https://www.wolframalpha.com/widgets/view.jsp?id=8ac60957610e1ee4894b2cd58e753
                                     'NO_ERROR': [True,False]
                                     },
                        'WINDOW_2': {'LOCATION': 3,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DVSYS', r'$\Delta$ ${\rm V_{sys}}$ (channels)',program],
                                     'NO_ERROR': [True,False],
                                     },
                        'WINDOW_3': {'LOCATION': 4,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DINCL', '$\Delta$ $i$ ($^{\circ}$)',program],
                                     'NO_ERROR': [True,False]},
                        'WINDOW_4': {'LOCATION': 6,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DPA', '$\Delta$ PA ($^{\circ}$)',program],
                                     'NO_ERROR': [True,False]},
                        'WINDOW_5': {'LOCATION': 7,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DVROT', r'$\Delta$ V$_{\rm rot}$  (channels)',program],
                                     'NO_ERROR': [True,False]},
                        'WINDOW_6': {'LOCATION': 8,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DFLUX', r'$\Delta$ Tot Flux  (Jy/beam km/s)',program],
                                     'NO_ERROR': [True,False]}},
                'PLOT_2':{
                        'WINDOW_0': {'LOCATION': 0,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DSDIS',r'$\Delta$ Dispersion  (channels)',program],
                                     'NO_ERROR': [True,False]},
                        'WINDOW_1': {'LOCATION': 1,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DZ0',r'$\Delta$ Scaleheight (beams)',program],
                                     'NO_ERROR': [True,False]
                                     },
                        'WINDOW_2': {'LOCATION': 3,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DR_HI',r'$\Delta$ R$_{\rm HI}$  (beams)',program],
                                     'NO_ERROR': [True,False]},
                        'WINDOW_3': {'LOCATION': 4,
                                     'X': ['DR_HI', r'$\Delta$ R$_{\rm HI}$  (beams)',program],
                                     'Y': ['DINCL',r'$\Delta$ $i$ ($^{\circ}$)',program],

                                     },
                        'WINDOW_4': {'LOCATION': 6,
                                     'X': ['SNR', 'SNR','Model'],
                                     'Y': ['DFLUX', r'$\Delta$ Tot Flux  (Jy/beam km/s)',program],
                                      'NO_ERROR': [True,False]},
                        'WINDOW_5': {'LOCATION': 7,
                                     'X': ['SNR', 'SNR','Model'],
                                     'Y': ['DVROT', r'$\Delta$ V$_{\rm rot}$  (channels)',program],
                                      'NO_ERROR': [True,False]},
                        'WINDOW_6': {'LOCATION': 8,
                                     'X': ['SNR', 'SNR','Model'],
                                     'Y': ['DR_HI',r'$\Delta$ R$_{\rm HI}$  (beams)',program],
                                      'NO_ERROR': [True,False]}}
                }

    for plot in plot_assembly:
        labelfont= {'family':'Times New Roman',
                        'weight':'normal',
                        'size':26}
        plt.rc('font',**labelfont)
        plt.figure(89,figsize=(24,24),dpi=300,facecolor = 'w', edgecolor = 'k')

            #this runs counter to figsize. How can a coding language be this illogical?
        gs = gridspec.GridSpec(3,3 )
        gs.update(wspace=0.3, hspace=0.3)
            #First the RA and DEC
            #RAval = np.array(RAval,dtype=float)
            #DECval = np.array(DECval,dtype=float)
            #if LVHIS:
            #    coloring = [0 if x == 'ROTCUR' else 90. for x in deltas['INPUT_MODEL']]
            #
            #else:
        order = [x for x in Input_Parameters]
        symbol= get_all_specific('Corruption',order,Input_Parameters,\
                                            program='Model')
        if database == 'LVHIS':
            symbol= get_all_specific('Corruption',order,Input_Parameters,\
                                        program='Model')[0]
            symbol_string= 'Compared to: '

        else:
            symbol= get_all_specific('Corruption',order,Input_Parameters,\
                                        program='Model')[0]

            symbol_string= 'Corrupted with: '
        for i,key in enumerate(plot_assembly[plot]):
            if 'NO_MEAN' in plot_assembly[plot][key]:
                nomean = plot_assembly[plot][key]['NO_MEAN']
            else:
                nomean = False
            if 'NO_ERROR' in plot_assembly[plot][key]:
                noerr= plot_assembly[plot][key]['NO_ERROR']
            else:
                noerr = [False,False]
            values =[]
            for par in ['X','Y']:
                values.append(get_all_specific(plot_assembly[plot][key][par][0]
                ,order,Input_Parameters,program=plot_assembly[plot][key][par][2]))
            if noerr[0]:
                values[0] = values[0][0]
            if noerr[1]:
                values[1] = values[1][0]


            ax,legend_items = make_plot(values[0],values[1],\
                           xlabel = plot_assembly[plot][key]['X'][1],\
                           ylabel = plot_assembly[plot][key]['Y'][1],\
                           location = gs[plot_assembly[plot][key]['LOCATION']],\
                           color=coloring,color_scale= coloring_scale,status = result_status, \
                           symbol=symbol,symbol_string= symbol_string\
                           ,No_Mean = nomean,size=result_status, size_string = 'Status =  ',\
                           no_error =  noerr )
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(4)

                    # increase tick width
                ax.tick_params(width=4)
            if 'PATCH' in plot_assembly[plot][key]:
                patch_fig = create_patch(plot_assembly[plot][key]['PATCH'],\
                    order,Input_Parameters,program=program)
                ax.add_patch(patch_fig)

            if i == 1:

                # make a legend
                labelfont= {'family':'Times New Roman',
                            'weight':'normal',
                            'size':18}
                plt.rc('font',**labelfont)
                chartBox = ax.get_position()
                ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*1.0, chartBox.height])
                ax.legend(handles=legend_items,loc='upper left', bbox_to_anchor=(1.25, 1.0), shadow=True, ncol=1)
                # Beams vs delta inclination
                labelfont= {'family':'Times New Roman',
                            'weight':'normal',
                            'size':26}
                #AND below it a corruption size indicator
                #unique_corruption = np.unique(deltas['CORRUPTION'])
                #if unique_corruption.size > 1:

                #    sizes


                plt.rc('font',**labelfont)

        #make a color bar
        ax = plt.subplot(gs[5])
        #Make a color bar for the inlination

        img = plt.imshow(np.array([coloring_scale]), cmap="rainbow")
        plt.gca().set_visible(False)
        cax = plt.axes([0.63, 0.4, 0.01, 0.45])
        barr = plt.colorbar(orientation="vertical", cax=cax)

        #barr.set_ticks([np.nanmin(coloring), np.max()])

        barr.set_label('Inclination', rotation=270, verticalalignment='bottom')


        labelfont= {'family':'Times New Roman',
                    'weight':'normal',
                    'size':37}
        plt.rc('font',**labelfont)
        galaxies = len(result_status)
        succes_galaxies=0
        for x in result_status:
            if x > 0.1:
                succes_galaxies+=1


        plt.figtext(0.5,0.91,f'''Out of {galaxies} galaxies, {succes_galaxies} were succesfully fitted''', horizontalalignment='center')

        #print([b[0] for item in RAval for b in item[0]])
        #RAval = np.array(RAval)

        if LVHIS:
            version= 'LVHIS'
        else:
            version='Database'

        plt.savefig(f'{filename}_{plot[-1]}_{version}.png', bbox_inches='tight')
        plt.close()



def plot_overview(Input_Parameters,filename='Overview_Difference',LVHIS=False\
                    ,program = None, database = None):

    plot_assembly = {'PLOT_1':{
                        'WINDOW_0': {'LOCATION': 0,
                                     'X': ['DXPOS', '$\Delta$ RA (beams)',program ],
                                     'Y': ['DYPOS','$\Delta$ DEC (beams)',program ],
                                     'NO_MEAN': True,
                                     'PATCH': ['Ellipse','DXPOS','DYPOS']
                                     },
                        'WINDOW_1': {'LOCATION': 1,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DCENTRAL','$\Delta$ Central (beams)',program ], #error from https://www.wolframalpha.com/widgets/view.jsp?id=8ac60957610e1ee4894b2cd58e753
                                     'NO_ERROR': [True,False]
                                     },
                        'WINDOW_2': {'LOCATION': 3,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DVSYS', r'$\Delta$ ${\rm V_{sys}}$ (channels)',program],
                                     'NO_ERROR': [True,False],
                                     },
                        'WINDOW_3': {'LOCATION': 4,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DINCL', '$\Delta$ $i$ ($^{\circ}$)',program],
                                     'NO_ERROR': [True,False]},
                        'WINDOW_4': {'LOCATION': 6,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DPA', '$\Delta$ PA ($^{\circ}$)',program],
                                     'NO_ERROR': [True,False]},
                        'WINDOW_5': {'LOCATION': 7,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DVROT', r'$\Delta$ V$_{\rm rot}$  (channels)',program],
                                     'NO_ERROR': [True,False]},
                        'WINDOW_6': {'LOCATION': 8,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DFLUX', r'$\Delta$ Tot Flux  (Jy/beam km/s)',program],
                                     'NO_ERROR': [True,False]}},
                'PLOT_2':{
                        'WINDOW_0': {'LOCATION': 0,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DSDIS',r'$\Delta$ Dispersion  (channels)',program],
                                     'NO_ERROR': [True,False]},
                        'WINDOW_1': {'LOCATION': 1,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DZ0',r'$\Delta$ Scaleheight (beams)',program],
                                     'NO_ERROR': [True,False]
                                     },
                        'WINDOW_2': {'LOCATION': 3,
                                     'X': ['Size_in_Beams', 'Diameter (beams)','Model'],
                                     'Y': ['DR_HI',r'$\Delta$ R$_{\rm HI}$  (beams)',program],
                                     'NO_ERROR': [True,False]},
                        'WINDOW_3': {'LOCATION': 4,
                                     'X': ['DR_HI', r'$\Delta$ R$_{\rm HI}$  (beams)',program],
                                     'Y': ['DINCL',r'$\Delta$ $i$ ($^{\circ}$)',program],

                                     },
                        'WINDOW_4': {'LOCATION': 6,
                                     'X': ['SNR', 'SNR','Model'],
                                     'Y': ['DFLUX', r'$\Delta$ Tot Flux  (Jy/beam km/s)',program],
                                      'NO_ERROR': [True,False]},
                        'WINDOW_5': {'LOCATION': 7,
                                     'X': ['SNR', 'SNR','Model'],
                                     'Y': ['DVROT', r'$\Delta$ V$_{\rm rot}$  (channels)',program],
                                      'NO_ERROR': [True,False]},
                        'WINDOW_6': {'LOCATION': 8,
                                     'X': ['SNR', 'SNR','Model'],
                                     'Y': ['DR_HI',r'$\Delta$ R$_{\rm HI}$  (beams)',program],
                                      'NO_ERROR': [True,False]}}
                }

    for plot in plot_assembly:
        labelfont= {'family':'Times New Roman',
                    'weight':'normal',
                    'size':26}
        plt.rc('font',**labelfont)
        plt.figure(89,figsize=(24,24),dpi=300,facecolor = 'w', edgecolor = 'k')

        #this runs counter to figsize. How can a coding language be this illogical?
        gs = gridspec.GridSpec(3,3 )
        gs.update(wspace=0.3, hspace=0.3)
        #First the RA and DEC
        #RAval = np.array(RAval,dtype=float)
        #DECval = np.array(DECval,dtype=float)
        #if LVHIS:
        #    coloring = [0 if x == 'ROTCUR' else 90. for x in deltas['INPUT_MODEL']]
        #
        #else:
        order = [x for x in Input_Parameters if \
            Input_Parameters[x][program]['Cat_Type'] == database]
        print(database)


        result_status = get_all_specific('Result',order,Input_Parameters,\
                                        program=program)[0]
        coloring = get_all_specific('INCL',order,Input_Parameters,\
                                        program='Model')[0]

        coloring_scale=[0.,90.]
        colorbarlabel = 'Inclination'
        if database == 'LVHIS':
            symbol= get_all_specific('Corruption',order,Input_Parameters,\
                                        program='Model')[0]
            symbol_string= 'Compared to: '

        else:
            symbol= get_all_specific('Corruption',order,Input_Parameters,\
                                        program='Model')[0]

            symbol_string= 'Corrupted with: '
        print(f'We are createing plot {plot}')
        for i,key in enumerate(plot_assembly[plot]):
            if 'NO_MEAN' in plot_assembly[plot][key]:
                nomean = plot_assembly[plot][key]['NO_MEAN']
            else:
                nomean = False
            if 'NO_ERROR' in plot_assembly[plot][key]:
                noerr= plot_assembly[plot][key]['NO_ERROR']
            else:
                noerr = [False,False]
            values =[]
            for par in ['X','Y']:

                values.append(get_all_specific(plot_assembly[plot][key][par][0]
                ,order,Input_Parameters,program=plot_assembly[plot][key][par][2]))
            if noerr[0]:
                values[0] = values[0][0]
            if noerr[1]:
                values[1] = values[1][0]

            print(f' plotting {plot_assembly[plot][key]["X"][0]} vs {plot_assembly[plot][key]["Y"][0]} with the values:' )
            print(values[0],values[1])
            ax,legend_items = make_plot(values[0],values[1],\
                           xlabel = plot_assembly[plot][key]['X'][1],\
                           ylabel = plot_assembly[plot][key]['Y'][1],\
                           location = gs[plot_assembly[plot][key]['LOCATION']],\
                           color=coloring,color_scale= coloring_scale,status = result_status, \
                           symbol=symbol,symbol_string= symbol_string\
                           ,No_Mean = nomean,size=result_status, size_string = 'Status =  ',\
                           no_error =  noerr )
            for axis in ['top','bottom','left','right']:
                ax.spines[axis].set_linewidth(4)

                    # increase tick width
                ax.tick_params(width=4)
            if 'PATCH' in plot_assembly[plot][key]:
                patch_fig = create_patch(plot_assembly[plot][key]['PATCH'],\
                    order,Input_Parameters,program=program)
                ax.add_patch(patch_fig)

            if i == 1:

                # make a legend
                labelfont= {'family':'Times New Roman',
                            'weight':'normal',
                            'size':18}
                plt.rc('font',**labelfont)
                chartBox = ax.get_position()
                ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*1.0, chartBox.height])
                ax.legend(handles=legend_items,loc='upper left', bbox_to_anchor=(1.25, 1.0), shadow=True, ncol=1)
                # Beams vs delta inclination
                labelfont= {'family':'Times New Roman',
                            'weight':'normal',
                            'size':26}
                #AND below it a corruption size indicator
                #unique_corruption = np.unique(deltas['CORRUPTION'])
                #if unique_corruption.size > 1:

                #    sizes


                plt.rc('font',**labelfont)

        #make a color bar
        ax = plt.subplot(gs[5])
        #Make a color bar for the inlination

        img = plt.imshow(np.array([coloring_scale]), cmap="rainbow")
        plt.gca().set_visible(False)
        cax = plt.axes([0.63, 0.4, 0.01, 0.45])
        barr = plt.colorbar(orientation="vertical", cax=cax)

        #barr.set_ticks([np.nanmin(coloring), np.max()])

        barr.set_label('Inclination', rotation=270, verticalalignment='bottom')


        labelfont= {'family':'Times New Roman',
                    'weight':'normal',
                    'size':37}
        plt.rc('font',**labelfont)
        galaxies = len(result_status)
        succes_galaxies=0
        for x in result_status:
            if x > 0.1:
                succes_galaxies+=1
        print(f'Succes Galaxies {succes_galaxies}')

        plt.figtext(0.5,0.91,f'''Out of {galaxies} galaxies, {succes_galaxies} were succesfully fitted''', horizontalalignment='center')

        #print([b[0] for item in RAval for b in item[0]])
        #RAval = np.array(RAval)

        if LVHIS:
            version= 'LVHIS'
        else:
            version='Database'
        print(filename,plot[-1],version)
        plt.savefig(f'{filename}_{plot[-1]}_{version}.png', bbox_inches='tight')
        plt.close()



def make_plot(x_in,y_in, status= None, location = [0,1], symbol= None,symbol_string= '',
                    xlabel = '',ylabel = '', No_Mean = False, size = None,color_scale= [0.,1.],
                    size_string= '',color=None,no_error = [False,False]):
        '''
        print(f' Make plot Input is:
x_in = {x_in}
y_in = {y_in}
status= {status}, location = {location}
symbol = {symbol}, symbol_string = {symbol_string}
xlabel = {xlabel}, ylabel = {ylabel}
size = {size}, size_string= {size_string}
color_scale= {color_scale}, color = {color}
No_Mean = {No_Mean},no_error = {no_error}')
        '''
        #x = np.array([v[0] for v in x_in], dtype=float)


        if no_error[0]:
            x = np.array(x_in)
            x_err = np.array([0. for v in x], dtype=float)
        else:
            x = np.array(x_in[0])
            x_err = np.array(x_in[1])
        #y = np.array([v[0] for v in y_in], dtype=float)

        if no_error[1]:
            y = np.array(y_in)
            y_err = np.array([0. for v in y], dtype=float)
        else:
            y = np.array(y_in[0])
            y_err = np.array(y_in[1])

        if not status is None:
            status = np.array(status)
            succes = np.where(status > 0.1)[0]
            mean = np.nanmean(y[succes])
            stdev = np.nanstd(y[succes]-mean)
            stat_elements = np.unique(status)
            norm_elements  = copy.deepcopy(stat_elements)
            min = np.min(norm_elements)
            max = np.max(norm_elements)
            if min/max < 0.1:
                norm_elements = norm_elements+0.1*(np.max(norm_elements)+0.01*np.max(norm_elements))
                max = np.max(norm_elements)
            norm_elements = norm_elements/max
            transparency = copy.deepcopy(status)
            for i,els in enumerate(stat_elements):
                transparency[status == els] = norm_elements[i]
            transparency = np.array(transparency,dtype = float)
        else:
            stat_elements = np.array([0.])
            norm_elements = [1]
            print(y)
            mean = np.nanmean(y[:])
            stdev = np.nanstd(y[:]-mean)
            transparency =np.ones(len(x[:]))
        ax = plt.subplot(location)
        ax.tick_params(
        axis='both',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom=True,      # ticks along the bottom edge are off
        top=True,         # ticks along the top edge are off
        labelbottom=True,
        right = True,
        left= True,
        labelleft = True)



        if not size is None:
            size_types = np.unique(size)
            numerical_size= [x for x in range(len(size_types))]
            size_size = copy.deepcopy(size)
            for i,c_type in enumerate(size_types):
                size_size = [numerical_size[i] if x == c_type else x for x in size_size]
            size_size = np.array(size_size,dtype=float)
            size = np.array(size)
            size_legend_items= []
            if len(numerical_size) > 1:
                sizes = (np.array(numerical_size,dtype=float)+4)**3.
                for i,siz in enumerate(sizes):

                    if siz != 4**3:
                        lab_string = f'{size_string}{size_types[i]}'
                        tmp_fig = plt.figure(1,figsize=(1,1),dpi=30,facecolor = 'w', edgecolor = 'k')
                        tmp_plot = plt.scatter([0,1],[0,1], c = 'k', s=siz, marker = 'o',label = lab_string)
                        size_legend_items.append(tmp_plot)
                        plt.close(tmp_fig)

        else:
            size_size = np.zeros(len(x[:]))
            size_legend_items= 0

        symlist = ["o", "v", "^", "<",">","s","P","*","X","D","1","3"]
        alphabet = [f"${x}$" for x in map(chr,range(97,123))]
        symlist = symlist+alphabet


        if not symbol is None:
            symbol= np.array(symbol)
            req_no_elements = np.unique(symbol)
            while len(req_no_elements) > len(symlist):
                symlist=symlist+symlist
            symbol_use = [symlist[i] for i,shape in enumerate(req_no_elements)]
        else:
            req_no_elements = ['Unspecified']
            symbol = np.array(['Unspecified' for gh in x])
            symbol_use = ['o']
        if not color is None:
            color = np.array(color,dtype=float)
            color = color/(np.nanmax(color_scale))-np.nanmin(color_scale)
            cmap = plt.cm.get_cmap('rainbow')
            rgba_cols = [cmap(color)]

        else:
            color = np.zeros(len(x))
            cmap = plt.cm.get_cmap('gist_gray')
            rgba_cols = [cmap(color)]
        shape_legend_items = []
        proc_lab_string = []
        for i,shaped in enumerate(req_no_elements):
                proc = np.where(symbol == shaped)[0]
                if len(proc) > 0:
                    lab_string = f'{symbol_string}{shaped}'
                    if lab_string not in proc_lab_string:
                        tmp_fig = plt.figure(1,figsize=(2,2),dpi=30,facecolor = 'w', edgecolor = 'k')
                        tmp_plot = plt.scatter(x[proc],y[proc],cmap= 'rainbow',\
                            c = 'k', s=4**3, marker = symbol_use[i],alpha = 1,\
                            label = lab_string)
                        shape_legend_items.append(tmp_plot)
                        proc_lab_string.append(lab_string)
                        plt.close(tmp_fig)
                    siz = (size_size[proc]+4)**3.
                    #try:
                    #ax.scatter(x[proc[add]],y[proc[add]],cmap= 'rainbow', c = rgba_cols[0][proc[add]][:], s=2**4, marker = symbol_use[i],alpha = 1,label = lab_string)

                    ax.scatter(x[proc],y[proc],cmap= 'rainbow', c = rgba_cols[0][proc][:], s=siz, marker = symbol_use[i])

                    plt.errorbar(x[proc],y[proc],xerr=x_err[proc],yerr=y_err[proc], linestyle="None", ecolor = rgba_cols[0][proc][:])
                    #except:
                    #    ax.scatter(x[proc[add]],y[proc[add]],cmap= 'rainbow', c = rgba_cols[0][proc[add]][:], s=siz, marker = symbol_use[i],alpha = norm_elements[j],label = lab_string)
                    #    plt.errorbar(x[proc[add]],y[proc[add]],xerr=np.zeros(len(add)),yerr=y[proc[add],1], linestyle="None", ecolor = rgba_cols[0][proc[add]][:],alpha = norm_elements[j])

        if not No_Mean:
            xmin,xmax = ax.get_xlim()
            ax.plot([xmin-1,xmax+2.],[mean,mean], c = 'k', alpha= 0.5)
            ax.plot([xmin-1,xmax+2.],[mean-stdev,mean-stdev], 'k--', alpha= 0.5)
            ax.plot([xmin-1,xmax+2.],[mean+stdev,mean+stdev], 'k--', alpha= 0.5)
            ax.set_xlim(xmin,xmax)

        ax.text(0.95,0.95,f'Mean = {mean:.1f} $\pm$ {stdev:.1f} ', transform=ax.transAxes,horizontalalignment= 'right', verticalalignment='top')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        if size_legend_items == 0.:
            legend_items = None
        else:
            legend_items = shape_legend_items+size_legend_items
        return ax,legend_items
# Function to calculate the difference between model and fit
def get_diff(
        val,model, radii = [], model_radii = [], single = False, errors = []\
        ,second = [] ,second_model = [], second_errors = [],norm=1.):
    val=np.array(val)
    if len(model) == 1:
        model = [model[0] for x in model_radii]
    model= np.array(model)
    model_radii = np.array(model_radii)
    errors = [x  if not np.isnan(x) else 0. for x in errors]
    mean_error = np.median(errors)

    errors = np.array(errors)
    errors[errors >  20.*mean_error] = 20.*mean_error
    second = np.array(second)
    if len(second_model) == 1:
        second_model = [second_model[0] for x in model_radii]
    second_model = np.array(second_model)
    second_errors = [x  if not np.isnan(x) else 0. for x in second_errors]
    mean_error = np.median(errors)
    second_errors = np.array(second_errors)
    second_errors[second_errors >  20.*mean_error] = 20.*mean_error
    to_use = np.where(abs(val) > 0.)[0]

    if to_use[0] != 0:
        to_use =np.hstack((0,to_use))
    model_to_use = np.where(abs(model) > 0.)[0]

    if model_to_use[0] != 0:
        model_to_use =np.hstack((0,model_to_use))

    if len(model_radii) > 0 and len(radii) > 0:
        model_int_function = interpolate.interpolate.interp1d(model_radii[model_to_use],model[model_to_use],fill_value="extrapolate")
        model = model_int_function(radii)

        if len(second) > 0:
            if np.sum(second_model) == 0.:
                second_model = model
            second_to_use = np.where(abs(second) > 0.)[0]
            if second_to_use[0] != 0:
                second_to_use =np.hstack((0,second_to_use))
            second_model_to_use = np.where(abs(second_model) > 0.)[0]
            if second_model_to_use[0] != 0:
                second_model_to_use =np.hstack((0,second_model_to_use))
            model_int_function = interpolate.interpolate.interp1d(model_radii[second_model_to_use],second_model[second_model_to_use],fill_value="extrapolate")
            second_model = model_int_function(radii)

    difference = abs(val[to_use]-model[to_use])
    if len(second) > 0:
        difference = np.hstack((difference,abs(second[second_to_use]-second_model[second_to_use])))
        if len(second_errors) > 0:
            errors = np.hstack((errors[to_use],second_errors[second_to_use]))
        else:
            errors = np.hstack((errors[to_use],errors[second_to_use]))
    else:
        errors = errors[to_use]
    difference =difference/norm
    if len(errors) > 0 and np.sum(errors) != 0.:
        errors = errors/norm
        difference = difference/errors
        value = np.sum(difference)/np.sum(1./errors)
        error  = np.mean(errors)
    else:
        value = np.mean(difference)
        if value == difference[0]:
            # If we really no idea we assume a 10% error on the normalisation.
            error = 0.1
        else:
            error = np.std(difference)
    if False:
        print(f'In values = {val}')
        print(f'In model = {model}')
        print(f'In Errors = {errors}')
        if len(second) > 0:
            print(f'In values2 = {second}')
            print(f'In model2 = {second_model}')
        print(f'In difference = {difference}')
        print(f'In errors = {errors}')
        print(f'norm = {norm} , value = {value}, error = {error}')
        exit()
    return [value,error]

# decompose the database name into quantities
def get_name_info(name):
    try:
        beam=float(name.split('ba')[1].split('SNR')[0])
        RCShape = name.split('-')[0]
    except:
        in_front = name.split('Beams')[0]
        individual = in_front.split('_')
        beam = float(individual[-1])
        RCShape = "_".join([x for x in individual[:-1]])
    try:
        SNR=float(name.split('SNR')[1].split('bm')[0])
    except:
        SNR =float(name.split('_')[-1].split('SNR')[0])
    return [beam],SNR,RCShape

def  get_LVHIS_average_flux(directory):
    for file in os.listdir(f'{directory}/Moments'):
        if file.endswith("preprocessed_mom0.fits") or \
            file.endswith("preprocessed_mom0_small.fits"):
            map = fits.open(f'{directory}/Moments/{file}')
    pixels = map[0].data[map[0].data > 0.].size
    flux = np.sum(map)/(map[0].header['BMAJ']/np.mean([abs(map[0].header['CDELT1']),abs(map[0].header['CDELT2'])]))
    return flux/pixels

def get_totflux(Configuration,map_name,debug=False):
    image = fits.open(f"{Configuration['FITTING_DIR']}{map_name}")
    #We are taking these from the moment map so we have to divide out the km/s
    flux = float(np.nansum(image[0].data)/Configuration['BEAM_IN_PIXELS'][2]/Configuration['CHANNEL_WIDTH'])
    #Should this not have an additional channel width parameter
    error = np.sqrt((np.where(image[0].data> 0.)[0].size)/Configuration['BEAM_IN_PIXELS'][2])*Configuration['NOISE']
    image.close()
    return [flux,error]
get_totflux.__doc__ =f'''
 NAME:
    get_totflux

 PURPOSE:
    Get the total flux from a intensity map

 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration
    map_name = name of the intensity fits file

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    total flux in the map in Jy*km/s

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

def get_flux_values(directory, input =False):
    #determine the mask and cube
    masked_cube = get_masked_cube(directory,input=input)
    Jy_conversion = get_beam_in_pixels(masked_cube[0].header)
    if input:
        try:
            noise = masked_cube[0].header['FAT_NOISE']
        except:
            noise = float('NaN')
    else:
        noise = float('NaN')
    total_flux = float(np.nansum(masked_cube[0].data)/Jy_conversion)

    error_flux = float(np.sqrt(masked_cube[0].data[masked_cube[0].data != 0.].size\
                    /Jy_conversion)*noise)
    average_flux = float(total_flux/masked_cube[0].data[masked_cube[0].data != 0.].size\
                        *Jy_conversion)#Jy/beam
    return [total_flux,error_flux,average_flux]

def get_beam_in_pixels(Template_Header, beam= None):
    if beam is None:
        beam = [Template_Header["BMAJ"],Template_Header["BMIN"]]
    #  https://science.nrao.edu/facilities/vla/proposing/TBconv
    beamarea=(np.pi*abs((beam[0]*beam[1])))/(4.*np.log(2.))
    return beamarea/(abs(Template_Header['CDELT1'])*abs(Template_Header['CDELT2']))

def get_masked_cube(directory, input =False):
    for file in os.listdir(f'{directory}/Sofia_Output/'):
        if '_binmask' in file:
            mask = fits.open(f'{directory}/Sofia_Output/{file}')
            cube_end = file.split('_')[-1]
            if cube_end == 'binmask.fits':
                cube_end = f'''{file.split('_')[-2]}.fits'''
        if '_mask.fits' in file:
            mask = fits.open(f'{directory}/Sofia_Output/{file}')
            cube_end = '_FAT.fits'
    if input:
        for file in os.listdir(f'{directory}'):
            if file.endswith(cube_end):
                print(f'We are opening the Cube {file}')
                cube=fits.open(f'{directory}/{file}')
    else:
        cube = fits.open(f'{directory}/Finalmodel/Finalmodel.fits')
    cube[0].data[mask[0].data < 0.1] = 0.
    return cube

def removefunction_get_sdsdsdflux_values(directory):
    for file in os.listdir(f'{directory}/Sofia_Output/'):
        if '_binmask' in file:
            mask = fits.open(f'{directory}/Sofia_Output/{file}')
            cube_end = file.split('_')[-1]
            if cube_end == 'binmask.fits':
                cube_end = f'''{file.split('_')[-2]}.fits'''
        if '_mask.fits' in file:
            mask = fits.open(f'{directory}/Sofia_Output/{file}')
            cube_end = '_FAT.fits'
    print(f'This is the end {cube_end}')
    for file in os.listdir(f'{directory}'):
        if file.endswith(cube_end):
            print(f'We are opening the Cube {file}')
            cube=fits.open(f'{directory}/{file}')
    noise = np.nanstd(cube[0].data[0,:,:])
    model = fits.open(f'{directory}/Finalmodel/Finalmodel.fits')
    cube[0].data[mask[0].data < 0.1] = 0.
    model[0].data[mask[0].data < 0.1] = 0.
    beamarea=(np.pi*abs(cube[0].header['BMAJ']*cube[0].header['BMIN']))/(4.*np.log(2.))
    beam_in_pixels = beamarea/(abs(cube[0].header['CDELT1'])*abs(cube[0].header['CDELT2']))
    print(f'This is the beam area in pixels {beam_in_pixels}  from beam {cube[0].header["BMAJ"]} x {cube[0].header["BMIN"]} and pixel {np.mean([abs(cube[0].header["CDELT1"]),abs(cube[0].header["CDELT2"])])}')
    total_flux = [np.nansum(cube[0].data)/beam_in_pixels,np.nansum(model[0].data)/beam_in_pixels] #Jy
    print(f'We found these fluxes in {total_flux} Jy')
    error_flux = [np.sqrt(cube[0].data[mask[0].data > 0.1].size/beam_in_pixels)*noise,\
                  np.sqrt(model[0].data[mask[0].data > 0.1].size/beam_in_pixels)*noise ]
    average_flux = [x/cube[0].data[mask[0].data > 0.1].size*beam_in_pixels for x  in total_flux]#Jy/beam
    return total_flux,error_flux,average_flux


def get_diff_rmax(model_rad,fat_rad,bmaj):
    '''Get the difference between the maximum extend in the model and the fit'''
    model_last_ring = (model_rad[-1]-model_rad[-2])/2.
    model_max= float(model_rad[-1]) + model_last_ring
    fat_last_ring =  (fat_rad[-1]-fat_rad[-2])/2.
    fat_max= float(fat_rad[-1]) + fat_last_ring
    delta_max= (model_max-fat_max)/(bmaj)
    return [delta_max, np.max([model_last_ring,fat_last_ring])/(float(bmaj))]

def load_input_catalogue(filename,split_char='|', debug = False,column_check=True,GDL=False):
    Catalogue = Proper_Dictionary({})
    with open(filename,'r') as tmpfile:
        firstline = tmpfile.readline()
        all_columns_check = False
        required_columns= ['ID','DISTANCE','DIRECTORYNAME','CUBENAME']
        while not all_columns_check:
            input_columns = [x.strip().upper() for x in firstline.split(split_char)]
            if GDL:
                if 'NUMBER' in input_columns:
                    input_columns[input_columns.index('NUMBER')] = 'ID'
                if 'DIRNAME' in input_columns:
                    input_columns[input_columns.index('DIRNAME')] = 'DIRECTORYNAME'



            Catalogue['ENTRIES'] = ['ENTRIES']
            Catalogue['ENTRIES'].extend(input_columns)

            if column_check:
                for key in required_columns:
                    print(f'Checking key {key}')
                    if key not in Catalogue['ENTRIES']:
                        if split_char == '|':
                            print(f'Key {key} not found')
                            split_char=' '
                            all_columns_check = False
                            break
                        else:
                            raise BadCatalogueError(f'We can not find the column for {key} in your input catalogue')
                    else:
                        all_columns_check = True
                        continue



        for key in input_columns:
            Catalogue[key] = []

        for line in tmpfile.readlines():
            input = [x.strip() for x  in line.split(split_char)]
            if len(input) == len(input_columns):
                for i,key in enumerate(input_columns):
                    if key == 'DISTANCE':
                        Catalogue[key].append(float(input[i]))
                    else:
                        Catalogue[key].append(input[i])
            else:
                print(f'READ_CATALOGUE: Your line "{line}" in the input catalogue does not have correct number of columns, skipping it')

    #if 'NUMBER' in Catalogue['ENTRIES']:
    #    Catalogue['NUMBER'] = np.array(Catalogue['NUMBER'],dtype=int)

    return Catalogue
load_input_catalogue.__doc__ =f'''
 NAME:
    catalogue

 PURPOSE:
    Read the FAT input catalogue and write into the a dictionary

 CATEGORY:
    read_functions

 INPUTS:
    filename = name of the catalogue to read

 OPTIONAL INPUTS:
    debug = False

 OUTPUTS:
    Catalogue = dictionary with the read file

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''



def load_output_catalogue(filename, debug = False,binary = False):
    Results = {}
    with open(filename) as file:
        ini_line=file.readline()
        tmp =[x.upper().strip() for x in ini_line.split()]
        if binary:
            columns = [f'{tmp[0]}_{tmp[1]}',tmp[2],tmp[3],f'{tmp[4]}_{tmp[5]}_{tmp[6]}_{tmp[7]}']
        else:
            columns = [f'{tmp[0]}_{tmp[1]}',tmp[2],f'{tmp[3]}_{tmp[4]}_{tmp[5]}_{tmp[6]}']
        for val in columns:
            Results[val] = []
        for line in file.readlines():
            vals = [x.strip() for x in line.split()]
            for i in range(len(columns)):
                if i !=  len(columns)-1:
                    Results[columns[i]].append(vals[i])
                else:
                    Results[columns[i]].append(f'{" ".join(vals[i:])}')
    return Results


def load_config_file(filename, debug = False, old_style=False):
    file_extension = os.path.splitext(filename)
    Translation_Dictionary = {'MAIN_DIRECTORY':'MAINDIR',
                              'CATALOGUE_START_ID':'STARTGALAXY',
                              'CATALOGUE_END_ID':'ENDGALAXY',
                              'INPUT_CATALOGUE':'CATALOGUE',
                              'OUTPUT_CATALOGUE':'OUTPUTCATALOGUE',
                              'FIXED_PARAMETERS': ['FIX_INCLINATION','FIX_PA','FIX_SDIS','FIX_Z0']}

    if file_extension[1] in ['.txt','.config']:
        print('We detected an old style config file')
        old_style = True
    elif file_extension[1] in ['.yml','.yaml']:
        print('We detected an yaml config file')

    else:
        print('We do not recognize this type of file assuming it is a yaml file')
    if old_style:
        required_configuration_keys = [Translation_Dictionary[key] for key in Translation_Dictionary]
        required_configuration_keys = [x for sub_list in required_configuration_keys for x in sub_list]
    else:
        required_configuration_keys = [key for key in Translation_Dictionary]

    tmpfile = open(filename, 'r')
    Configuration = Proper_Dictionary({})
# Separate the keyword names
    if old_style:
        key_list = list(Translation_Dictionary.keys())
        value_list =  list(Translation_Dictionary.values())
        for tmp in tmpfile.readlines():
            if tmp[0] != '#':
            # python is really annoying with needing endlines. Let's strip them here and add them when writing
                add_key = tmp.split('=', 1)[0].strip().upper()
                if add_key in Translation_Dictionary['FIXED_PARAMETERS']:
                    if 'FIXED_PARAMETERS' not in Configuration:
                        Configuration['FIXED_PARAMETERS'] = []
                    invalid_input = True
                    inp = tmp.split('=', 1)[1].strip()
                    while invalid_input:
                        if inp.lower() == "true" or inp.lower() == "t" or inp.lower() == "y" or inp.lower() == "yes" or inp[0] == '1':
                            value = True
                            invalid_input = False
                        elif inp.lower() == "false" or inp.lower() == "f" or inp.lower() == "n" or inp.lower() == "no" or inp[0] == '0':
                            value = False
                            invalid_input = False
                        else:
                            inp = input("The parameter {} in the configuration file  must be true/false or yes/no. Please give the correct value. \n".format(add_key))
                    if value:
                        parameter = add_key.split('_')[1]
                        Configuration['FIXED_PARAMETERS'].append(parameter)
                elif add_key in value_list:
                    Configuration[key_list[value_list.index(add_key)]] = tmp.split('=', 1)[1].strip()
                else:

                    Configuration[add_key] = tmp.split('=', 1)[1].strip()

        if not os.path.exists(Configuration['MAIN_DIRECTORY']):
            indir = '/'.join(filename.split('/')[:-1])
            if Configuration['MAIN_DIRECTORY'] != indir and\
                os.path.exists(indir):
                Configuration['MAIN_DIRECTORY']=indir
            else:
                Configuration['MAIN_DIRECTORY'] = input(f'''
                    Your main directory ({Configuration['MAIN_DIRECTORY']}) does not exist.
                    Please provide the correct file name.
                    ''')
    else:
        cfg = OmegaConf.structured(defaults)

        # read command line arguments anything list input should be set in '' e.g. pyROTMOD 'rotmass.MD=[1.4,True,True]'
        yaml_config = OmegaConf.load(filename)
#merge yml file with defaults
        cfg = OmegaConf.merge(cfg,yaml_config)
        # translate into our dictionary
        if not os.path.exists(cfg.input.main_directory):
            indir = '/'.join(filename.split('/')[:-1])

            if cfg.input.main_directory != indir and\
                os.path.exists(indir):
                cfg.input.main_directory=indir
            else:
                cfg.input.main_directory = input(f'''
                    Your main directory ({cfg.input.main_directory}) does not exist.
                    Please provide the correct file name.
                    ''')
        if not os.path.exists(cfg.input.catalogue):
            indir = '/'.join(filename.split('/')[:-1])
            catalogue = cfg.input.catalogue.split('/')[-1]

            if cfg.input.catalogue != f"{indir}{catalogue}" and\
                os.path.exists(f"{indir}{catalogue}"):
                cfg.input.catalogue = f"{indir}{catalogue}"
            else:
                cfg.input.catalogue = input(f'''
                    Your catalogue ({cfg.input.catalogue}) does not exist.
                    Please provide the correct file name.
                    ''')

        Configuration_ini = setup_configuration(cfg)

        Configuration_ini['INPUT_CATALOGUE'] = Configuration_ini['CATALOGUE']
        Configuration= Configuration_ini
        #print(Configuration['MAIN_DIRECTORY'],Configuration['INPUT_CATALOGUE'])
        #Configuration['INPUT_CATALOGUE']= f'''{Configuration['MAIN_DIRECTORY']}{Configuration['INPUT_CATALOGUE']}'''

    while not os.path.exists(Configuration['INPUT_CATALOGUE']):
        file = Configuration['INPUT_CATALOGUE'].split('/')[-1]

        if os.path.exists(f"{Configuration['MAIN_DIRECTORY']}{file}"):
            Configuration['INPUT_CATALOGUE'] = f"{Configuration['MAIN_DIRECTORY']}{file}"
        else:
            Configuration['INPUT_CATALOGUE'] = input('''
                    Your input catalogue ({}) does not exist.
                    Please provide the correct file name.
                    '''.format(Configuration['INPUT_CATALOGUE']))
    #The output catalogue only needs to be in a valid directory as we create it
    output_catalogue_dir = Configuration['OUTPUT_CATALOGUE'].split('/')
    if len(output_catalogue_dir) > 1:
        check_dir = '/'.join(output_catalogue_dir[:-1])
        while not os.path.isdir(check_dir):
            if os.path.exists(f"{Configuration['MAIN_DIRECTORY']}{output_catalogue_dir[-1]}"):
                Configuration['OUTPUT_CATALOGUE'] = f"{Configuration['MAIN_DIRECTORY']}{output_catalogue_dir[-1]}"
                check_dir=f"{Configuration['MAIN_DIRECTORY']}"
            else:
                check_dir= input('''
                    The directory for your output catalogue ({}) does not exist.
                    Please provide the correct directory name.
                    '''.format(Configuration['OUTPUT_CATALOGUE']))
                Configuration['OUTPUT_CATALOGUE'] = check_dir+'/'+output_catalogue_dir[-1]



    return Configuration
load_config_file.__doc__ ='''
;+
; NAME:
;       config_file(input_parameters, start_dir)
;
; PURPOSE:
;       Read the FAT config file and write into the a dictionary
;
; CATEGORY:
;       read
;
;
; INPUTS:
;       Configuration
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;          result = dictionary with the config file
;
; OPTIONAL OUTPUTS:
;       -
;
; PROCEDURES CALLED:
;      split, strip, open
;
; EXAMPLE:
;
;
'''

def load_LVHIS_Name(directory):


    with open(f'{directory}/Galaxy_Names.txt') as tmpfile:
        tmp = tmpfile.readlines()

    name_dictionary = {}

    for line in tmp:
        line_split = [x.strip() for x in line.split()]
        if len(line_split) == 5:
            pass
        elif len(line_split) == 4:
            name_dictionary[line_split[0].upper()] = f'{line_split[1].upper()} {line_split[2]}'
        else:
            print(f' The line {line} does not adhere to our expectations.' )
    return name_dictionary
#Function for loading the variables of a tirific def file into a set of variables to be used
def load_tirific(filename,Variables = None, LVHIS=False):
    if Variables is None:
        Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',\
                     'VROT_ERR','Z0', 'Z0_ERR', 'SBR','SBR_ERR', 'INCL',\
                     'INCL_ERR','PA','PA_ERR','XPOS','XPOS_ERR','YPOS',\
                     'YPOS_ERR','VSYS','VSYS_ERR','SDIS','SDIS_ERR','VROT_2',\
                     'VROT_2_ERR',  'Z0_2','Z0_2_ERR','SBR_2','SBR_2_ERR',
                     'INCL_2', 'INCL_2_ERR','PA_2','PA_2_ERR','XPOS_2',\
                     'XPOS_2_ERR','YPOS_2','YPOS_2_ERR','VSYS_2','VSYS_2_ERR',\
                     'SDIS_2','SDIS_2_ERR','CONDISP','CFLUX','CFLUX_2']

    if LVHIS:
        Variables.append('RADI_2')
    Variables = np.array([e.upper() for e in Variables],dtype=str)
    if filename == 'EMPTY':
        template = ['RADI = 0. 100.  200.\n']
    else:
        with open(filename) as file:
            template = file.readlines()
    output = {}
    # Separate the keyword names
    for line in template:
        if line.count('=') > 1:
            print(f'This taco is not correct. \n You have lines in  def file {filename} where = occurs multiple times.')
            print(f'This is the offending line {line}')
            exit()
        var_concerned = str(line.split('=')[0].strip().upper())
        if len(var_concerned) <2:
            continue
        if var_concerned[0] == '#':
            var_concerned = var_concerned[1:].strip()
        if len(var_concerned) > 1:
            if var_concerned in Variables:

                output[var_concerned] = [float(x) for x in line.split('=')[1].rsplit()]
    for input in Variables:
        if input not in output:
            if input == 'RADI_2':
                pass
            else:
                output[input] = [float('NaN') for x in output['RADI']]

    return output



def remove_too_faint(tirific_dictionary):
    '''Remove all elements in the def dictionary that are too faint'''
    for key in  tirific_dictionary:
        if len(tirific_dictionary[key]) > 1:
            if key == 'RADI':
                sbr = [np.max([x,y]) for x,y in zip(tirific_dictionary['SBR'], tirific_dictionary['SBR_2'])]
            elif key[-1] == 2:
                sbr = tirific_dictionary['SBR_2']
            else:
                sbr = tirific_dictionary['SBR']
            current = sbr[-1]
            while current < 1e-12:
                sbr = sbr[:-1]
                current = sbr[-1]
            tirific_dictionary[key] = tirific_dictionary[key][:len(sbr)]
    return tirific_dictionary

def fix_links(database_config, database_out_catalogue):
    for i, galaxy in enumerate(database_out_catalogue['DIRECTORY_NAME']):
        if str_to_bool(database_out_catalogue['OS'][i]):
            if not os.path.isfile(f'{database_config["MAIN_DIRECTORY"]}/{galaxy}/Finalmodel/Finalmodel.def') and\
                os.path.isdir(f'{database_config["MAIN_DIRECTORY"]}/{galaxy}/Finalmodel/'):
                linkname = f'{database_config["MAIN_DIRECTORY"]}/{galaxy}/Fit_Tirific_OSC/Fit_Tirific_OSC'
                os.symlink(f'{linkname}.fits',f'{database_config["MAIN_DIRECTORY"]}/{galaxy}/Finalmodel/Finalmodel.fits')
                os.symlink(f'{linkname}.def',f'{database_config["MAIN_DIRECTORY"]}/{galaxy}/Finalmodel/Finalmodel.def')

def get_result(database_out_catalogue,database_config,index,galaxy,binary = False):
    status =0
    if binary:
        if int(database_out_catalogue['AC1'][index]) ==1:
            status =1
        if int(database_out_catalogue['AC2'][index]) ==1:
            status =2
    else:
        print(f'''What on earth {database_out_catalogue['OS'][index]}''')
        if str_to_bool(database_out_catalogue['OS'][index]):
            status = 2
        else:
            if os.path.isfile(f'{database_config["MAIN_DIRECTORY"]}/{galaxy}/Finalmodel/Finalmodel.def'):
                status = 1
    return [status,database_out_catalogue['COMMENTS_ON_FIT_RESULT'][index]]

def galaxy_deltas(output_dict,model_dict,LVHIS=False):
    #Get the deltas for an individual galaxy
    ext = ''
    if 'RADI_2' in model_dict:
        models = ['ROTCUR','DISK_FIT']
    elif np.sum(model_dict[f"VROT_2"]) == 0:
        models = ['ROTCUR']
    else:
        models = ['TIRIFIC']
    if output_dict['Result'][0] > 0:
        deltas ={}
        for model in models:
            if model == 'DISK_FIT':
                ext = '_2'
            deltas[model] = {}
        #First some odd deltas
            deltas[model]['MAX_EXTEND'] = get_diff_rmax(\
                model_dict[f'RADI{ext}'],output_dict['RADI'],output_dict['BMAJ'][0])
        #If the model has a SBR profile we can determine a Tru RHI
        if model == 'TIRIFIC':
            model_RHI = get_RHI(sbr_profile=[model_dict[f'SBR'],\
                model_dict[f'SBR_2']],radi=model_dict[f'RADI'],\
                systemic=model_dict[f'VSYS'][0],distance=model_dict['Distance'])
            fitted_RHI = get_RHI(sbr_profile=[output_dict[f'SBR'],\
                output_dict[f'SBR_2']],radi=output_dict[f'RADI'],\
                systemic=output_dict[f'VSYS'][0],distance=model_dict['Distance'])

            deltas[model]['R_HI'] = [(model_RHI[0]-fitted_RHI[0])/\
                (float(output_dict['BMAJ'][0])),(model_RHI[1]+fitted_RHI[1])/\
                (float(output_dict['BMAJ'][0]))]
        else:
            #If there is not SBR we use a proxy
            deltas[model]['R_HI'] = deltas[model]['MAX_EXTEND']


        deltas[model]['FLUX'] = [output_dict['FLUX'][0]-\
            model_dict['FLUX'][0],output_dict['FLUX'][1]+model_dict['FLUX'][1]]

        #Then we go through the parameters in the output dictionary
        for parameter in output_dict:
            print(f'Starting to check {parameter}')
            #If not a standard parameter we skip it
            try:
                length =  len(output_dict[parameter])
                number = float(output_dict[parameter][0])
            except:
                print(f"{parameter} is not a proper parameter varying per ring. Skipping it")
                continue

            if parameter[-1] == '2' or parameter == 'RADI' or \
                len(output_dict[parameter]) != len(output_dict['RADI']) \
                 or parameter[-4:] == '_ERR':
                 print(f"{parameter} is not a proper parameter varying per ring. Skipping it")
                 continue
            if np.sum(output_dict[parameter]) == 0. or \
                np.isnan(output_dict[parameter][0]):
                    raise ModelError(f'There is an error in the output {parameter} in {galaxy}')

            # Set the normalisations values
            if parameter in ['VROT','SDIS','VSYS']:
                normalisation = float(output_dict['CHANNEL_WIDTH'])
            elif parameter in ['XPOS','YPOS']:

                normalisation = float(output_dict['BMAJ'][0])/3600.
            else:
                normalisation = 1.
            if np.sum(model_dict[f"{parameter}{ext}"]) == 0.:
                if parameter not in ['Z0','SBR','SDIS']:
                    raise ModelError(f'We did not find model parameters for {parameter} in {model}')
                else:
                    deltas[model][parameter] = \
                        [float('NaN'),float('NaN')]
                continue
            if model in ['ROTCUR','DISK_FIT']:
                deltas[model][parameter] = get_diff(output_dict[parameter],\
                    model_dict[f'{parameter}{ext}'],radii = output_dict['RADI'],\
                    model_radii=model_dict[f'RADI{ext}'],errors =  \
                    output_dict[f'{parameter}_ERR'],norm = normalisation,\
                    second = output_dict[f'{parameter}_2'], \
                    second_model = model_dict[f'{parameter}{ext}'],second_errors = \
                    output_dict[f'{parameter}{ext}_ERR'])
            elif model == 'TIRIFIC':
                deltas[model][parameter] = get_diff(output_dict[parameter],\
                    model_dict[f'{parameter}'],radii = output_dict['RADI'],\
                    model_radii=model_dict['RADI'],errors =  \
                    output_dict[f'{parameter}_ERR'],norm = normalisation, \
                    second = output_dict[f'{parameter}_2'], \
                    second_model = model_dict[f'{parameter}_2'],second_errors = \
                    output_dict[f'{parameter}_2_ERR'])
            else:
                raise ModelError(f"The model {model} is not ours")
            print(f'For {parameter} we find a difference between mod and fit {deltas[model][parameter]}')
        output_dict['deltas'] = deltas
    else:
        output_dict['deltas'] = None



def removefunction_retrieve_deltas_and_RCs(database_config, database_inp_catalogue, database_out_catalogue,binary=False,LVHIS= False):
    '''for every galaxy read the Finalmodel and retrieve RCs and deltas'''
    if LVHIS:
        LVHIS_Names = load_LVHIS_Name(database_config["MAIN_DIRECTORY"])
    RCs = {}
    #deltas = {'NAME': [], 'BEAMS_ACROSS': [], 'SNR':[], 'MAX_EXTEND': [], 'R_HI':[]}
    deltas= {}
    for i, galaxy in enumerate(database_out_catalogue['DIRECTORY_NAME']):
        print(f'Processing directory {galaxy}')
        status  = get_result(database_out_catalogue,database_config,i,galaxy,\
                    binary = binary)
        status = status[0]
        # First read the the input model Class

        #Then read the model

        model_parameters = load_tirific(
            f'{database_config["MAIN_DIRECTORY"]}/{galaxy}/ModelInput.def',LVHIS=LVHIS)
        print(f'We have the status of {status} so let go')
        if status == 0:
            output_parameters =  load_tirific('EMPTY')
            total_flux = [float('NaN'),float('NaN')]
            error_flux =[float('NaN'),float('NaN')]
            average_flux = [float('NaN'),float('NaN')]


        else:
            output_parameters =  load_tirific(
                f'{database_config["MAIN_DIRECTORY"]}/{galaxy}/Finalmodel/Finalmodel.def')
            total_flux = get_flux_values(f'{database_config["MAIN_DIRECTORY"]}/{galaxy}/')

        output_parameters = remove_too_faint(output_parameters)
        if 'DISTANCE' in model_parameters:
            distance=model_parameters['DISTANCE']
        else:
            distance=[model_parameters['VSYS'][0]/69.7] #km/s/Mpc
        if len(distance) == 0.:
            distance=[model_parameters['VSYS'][0]/69.7]
        if distance[0] == 0.:
            distance=[model_parameters['VSYS'][0]/69.7]
        #write our RCs to a variable
        # First read the the input model Class
        if LVHIS:

            diameter_in_beams = [float(model_parameters['RADI'][-1])/float(model_parameters['BMAJ'][0])*2.]
            if 'RADI_2' in model_parameters:
                print('What happened')
                diameter_in_beams.append(float(model_parameters['RADI_2'][-1])/float(model_parameters['BMAJ'][0])*2.)
            else:
                diameter_in_beams.append(float(model_parameters['RADI'][-1])/float(model_parameters['BMAJ'][0])*2.)
            if not np.isnan(average_flux[0]):
                SNR = average_flux[0]/output_parameters['RMS']
            else:
                SNR = -1.
            RCshape = LVHIS_Names[galaxy]
        else:
            diameter_in_beams, SNR, RCshape = get_name_info(galaxy)
        if RCshape not in RCs:
            RCs[RCshape] = {'MODEL': {'RADIUS':model_parameters['RADI'], 'RC':model_parameters['VROT'], 'DISTANCE': distance, 'STATUS': -1}}
        else:
            if 'MODEL' not in RCs[RCshape]:
                RCs[RCshape]['MODEL'] = {'RADIUS':model_parameters['RADI'], 'RC':model_parameters['VROT'], 'DISTANCE': distance, 'STATUS': -1}
        if LVHIS:
            if 'MODEL_2' not in RCs[RCshape]:
                RCs[RCshape]['MODEL_2'] = {'RADIUS':model_parameters['RADI_2'], 'RC':model_parameters['VROT_2'], 'DISTANCE': distance, 'STATUS': -1}
        if len(output_parameters['VROT']) < len(output_parameters['RADI']):
            vrot = output_parameters['VROT_2']
        else:
            vrot = output_parameters['VROT']
        RCs[RCshape][galaxy] = {'RADIUS':output_parameters['RADI'], 'RC':vrot, 'DISTANCE':distance, 'STATUS': status}
        if binary:
            cubename = f'{database_inp_catalogue["CUBENAME"][database_inp_catalogue["DIRECTORYNAME"].index(galaxy)]}_preprocessed.fits'
        else:
            cubename = f'{database_inp_catalogue["CUBENAME"][database_inp_catalogue["DIRECTORYNAME"].index(galaxy)]}_FAT.fits'

        tmp_corruption = database_inp_catalogue["CUBENAME"][database_inp_catalogue["DIRECTORYNAME"].index(galaxy)].split('_')
        if len(tmp_corruption) > 1:
            if tmp_corruption[-1] in ['CS', 'Gauss','UC']:
                if tmp_corruption[-1] == 'CS':
                    corruption = 'Casa Sim'
                elif tmp_corruption[-1] == 'UC':
                    corruption = 'Uncorrupted'
                else:
                    corruption =  tmp_corruption[-1]
            else:
                corruption = 'Unspecified'
        else:
            corruption = 'Unspecified'
        if status > 0:
            hdr = fits.getheader(f'{database_config["MAIN_DIRECTORY"]}{galaxy}/{cubename}')
            channel_width = hdr['CDELT3']/1000.
        #and then calculate the deltas
        if LVHIS:
            input_models = ['ROTCUR','DISKFIT']
        else:
            input_models = ['TIRIFIC']
        ext=['','_2']
        for j,model in enumerate(input_models):
            if i == 0:
                deltas = {'NAME': [], 'BEAMS_ACROSS': [], 'SNR':[], \
                                 'MAX_EXTEND': [], 'R_HI':[],'STATUS': [],\
                                 'CENTRAL_INPUT_INCLINATION': [],'INPUT_MODEL': [],\
                                  'RCSHAPE': [], 'TOTAL_FLUX': [], 'CORRUPTION': []}
            deltas['NAME'].append(galaxy)
            deltas['BEAMS_ACROSS'].append(diameter_in_beams[j])
            deltas['SNR'].append(SNR)
            deltas['CORRUPTION'].append(corruption)
            deltas['INPUT_MODEL'].append(model)
            deltas['RCSHAPE'].append(RCshape)
            deltas['STATUS'].append(status)
            deltas['CENTRAL_INPUT_INCLINATION'].append(model_parameters[f'INCL{ext[j]}'][0])

            # First the delta in extend as it is slightly different
            if status > 0.:
                deltas['MAX_EXTEND'].append(
                    get_diff_rmax(model_parameters[f'RADI{ext[j]}'],output_parameters['RADI'],hdr['BMAJ']))
                if model == 'TIRIFIC':
                    print(f'In the galaxy {galaxy}')
                    model_RHI = get_RHI(sbr_profile=[model_parameters[f'SBR'],model_parameters[f'SBR_2']],\
                                radi=model_parameters[f'RADI'],systemic=model_parameters[f'VSYS'][0],\
                                distance=distance )
                    fitted_RHI = get_RHI(sbr_profile=[output_parameters[f'SBR'],output_parameters[f'SBR_2']],\
                                radi=output_parameters[f'RADI'],systemic=output_parameters[f'VSYS'][0],distance=distance )

                    deltas['R_HI'].append([(model_RHI[0]-fitted_RHI[0])/(hdr['BMAJ']*3600.),\
                                            (model_RHI[1]+fitted_RHI[1])/(hdr['BMAJ']*3600.)])

                else:
                    deltas['R_HI'].append(get_diff_rmax(model_parameters[f'RADI{ext[j]}'],output_parameters['RADI'],hdr['BMAJ']))

            else:
                deltas['MAX_EXTEND'].append([float('NaN'),float('NaN')])
                deltas['R_HI'].append([float('NaN'),float('NaN')])
            #for now leave R_HI
            deltas['TOTAL_FLUX'].append([total_flux[1]-total_flux[0], error_flux[1]+error_flux[0]])
            for key in output_parameters:
                print(f'Starting to check {key}')
                if key[-1] == '2' or key == 'RADI' or \
                    len(output_parameters[key]) < 2 or key[-4:] == '_ERR':
                    continue
                if i == 0:
                    deltas[key] = []
                if np.sum(output_parameters[key]) == 0. or np.isnan(output_parameters[key][0]):
                    if key in deltas:
                        deltas[key].append([float('NaN'),float('NaN')])
                        print(f'For {key} in {galaxy} we find no values so we set the difference to NaN')
                else:
                    if key in ['VROT','SDIS','VSYS']:
                        normalisation = channel_width
                    elif key in ['XPOS','YPOS']:
                        normalisation = hdr['BMAJ']
                    else:
                        normalisation = 1.
                    if model =='ROTCUR':
                        if np.sum(model_parameters[f'{key}']) == 0.:
                            if key not in ['Z0','Z0_2','SBR','SBR_2','SDIS','SDIS_2']:
                                raise ModelError(f'We did not find model parameters for {key} in {model}')
                            else:
                                deltas[key].append([float('NaN'),float('NaN')])
                        else:
                            deltas[key].append(get_diff(output_parameters[key],model_parameters[key],\
                                radii = output_parameters['RADI'], model_radii= model_parameters['RADI'],
                                errors =  output_parameters[f'{key}_ERR'],norm = normalisation,
                                second = output_parameters[f'{key}_2'], second_model = model_parameters[f'{key}'],
                                second_errors = output_parameters[f'{key}_ERR']))
                    elif model =='DISKFIT':
                        if np.sum(model_parameters[f'{key}_2']) == 0.:
                            deltas[key].append([float('NaN'),float('NaN')])
                        else:
                            deltas[key].append(get_diff(output_parameters[key],model_parameters[f'{key}_2'],\
                                radii = output_parameters['RADI'], model_radii= model_parameters['RADI_2'],
                                errors =  output_parameters[f'{key}_ERR'],norm = normalisation,
                                second = output_parameters[f'{key}_2'], second_model = model_parameters[f'{key}_2'],
                                second_errors = output_parameters[f'{key}_2_ERR'] ))
                    else:
                        if np.sum(model_parameters[f'{key}']) == 0.:
                            if key not in ['Z0','Z0_2','SBR','SBR_2','SDIS','SDIS_2']:
                                raise ModelError(f'We did not find model parameters for {key} in {model}')
                            else:
                                deltas[key].append([float('NaN'),float('NaN')])
                        else:
                            deltas[key].append(get_diff(output_parameters[key],model_parameters[key],\
                                radii = output_parameters['RADI'], model_radii= model_parameters['RADI'],
                                errors =  output_parameters[f'{key}_ERR'],norm = normalisation,
                                second = output_parameters[f'{key}_2'], second_model = model_parameters[f'{key}_2'],
                                second_errors = output_parameters[f'{key}_2_ERR'] ))

                    print(f'For {key} in {galaxy} we find a difference between mod and fit {deltas[key][-1]}')
        # We also need the flux difference
        #deltas[key].append(sf.get_flux_diff(


    return deltas, RCs


def get_DHI(Configuration,Model='Finalmodel' ,debug=False):
    #Get the sbrs
    radi,sbr,sbr_2,systemic = load_tirific(Configuration,f"{Configuration['FITTING_DIR']}{Model}/{Model}.def",Variables = ['RADI','SBR','SBR_2','VSYS'],debug=debug)
    #convert to solar_mass/pc^2
    sbr_msolar = columndensity(Configuration,sbr*1000.,systemic=systemic[0],arcsquare=True,solar_mass_output=True)
    sbr_2_msolar = columndensity(Configuration,sbr_2*1000.,systemic=systemic[0],arcsquare=True,solar_mass_output=True)
    # interpolate these to ~1" steps
    new_radii = np.linspace(0,radi[-1],int(radi[-1]))
    new_sbr_msolar = np.interp(new_radii,radi,sbr_msolar)
    new_sbr_2_msolar = np.interp(new_radii,radi,sbr_2_msolar)

    index_1 = np.where(new_sbr_msolar > 1.)[0]
    index_2 = np.where(new_sbr_2_msolar > 1.)[0]
    if index_1.size > 0 and index_2.size > 0:
        DHI = float(new_radii[index_1[-1]]+new_radii[index_2[-1]])
    elif index_1.size > 0:
        DHI = float(new_radii[index_1[-1]])
    elif index_2.size > 0:
        DHI = float(new_radii[index_2[-1]])
    else:
        DHI = float('NaN')
    return DHI
get_DHI.__doc__ =f'''
 NAME:
    get_DHI

 PURPOSE:
    get the DHI as determined by the SBR profiles in the fit from the Tirific Template

 CATEGORY:
    read_functions

 INPUTS:
    Configuration = Standard FAT configuration

 OPTIONAL INPUTS:
    debug = False

    Model = 'Finalmodel'
    location of the def file to get DHI from. it should be in the fitting dir in the {{Model}}/{{Model}}.def

 OUTPUTS:

 OPTIONAL OUTPUTS:

 PROCEDURES CALLED:
    Unspecified

 NOTE:
'''

#Calculate the actual number of rings in the model from ring size and the size in beams:
def get_RHI(sbr_profile=[0.,0.],radi= [0.],systemic = 100.,distance=1.):

    sbr_msolar = columndensity({'OUTPUTLOG': None, 'DISTANCE': distance,'DEBUG':False, 'TIMING': True},np.array(sbr_profile[0],dtype=float)*1000.\
                        ,systemic=systemic,arcsquare=True,solar_mass_output=True)
    sbr_2_msolar = columndensity({'OUTPUTLOG': None, 'DISTANCE': distance,'DEBUG':False,'TIMING': True},np.array(sbr_profile[1],dtype=float)*1000.\
                        ,systemic=systemic,arcsquare=True,solar_mass_output=True)
    # interpolate these to ~1" steps
    new_radii = np.linspace(0,radi[-1],int(radi[-1]))
    while len(radi) > sbr_msolar.size:
        sbr_msolar = np.concatenate([sbr_msolar,[0.]])

    new_sbr_msolar = np.interp(new_radii,radi,sbr_msolar)

    while len(radi) > sbr_2_msolar.size:
        sbr_2_msolar = np.concatenate([sbr_2_msolar,[0.]])

    new_sbr_2_msolar = np.interp(new_radii,radi,sbr_2_msolar)

    index_1 = np.where(new_sbr_msolar > 1.)[0]
    index_2 = np.where(new_sbr_2_msolar > 1.)[0]
    if index_1.size > 0 and index_2.size > 0:
        DHI = float(new_radii[index_1[-1]]+new_radii[index_2[-1]])
    elif index_1.size > 0:
        DHI = float(new_radii[index_1[-1]])*2.
    elif index_2.size > 0:
        DHI = float(new_radii[index_2[-1]])*2.
    else:
        DHI = float('NaN')

    return [DHI/2.,(radi[-1]-radi[-2])/10.]

#Really python bool(str) is always true
def str_to_bool(inp):
    if inp == "":
        return False
    elif inp.lower() == "true" or inp.lower() == "t" or inp.lower() == "y" or inp.lower() == "yes":
        return True
    elif inp.lower() == "false" or inp.lower() == "f" or inp.lower() == "n" or inp.lower() == "no":
        return False
    else:
        print(f"Error: the input must be true/false or yes/no. {inp} is not a clear boolean")
