#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


from collections import OrderedDict #used in Proper_Dictionary
from inspect import getframeinfo,stack
from scipy.optimize import curve_fit
from scipy import ndimage
from scipy import interpolate
from astropy.wcs import WCS
from astropy.io import fits

import os
import signal
import traceback
import numpy as np
import copy
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt


import re
import subprocess
class SupportRunError(Exception):
    pass
class SmallSourceError(Exception):
    pass
class FileNotFoundError(Exception):
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



#Calculate the actual number of rings in the model from ring size and the size in beams:

def calc_rings(Configuration,size_in_beams = 0., ring_size  = 0.,debug=False):
    if ring_size == 0.:
        ring_size = Configuration['RING_SIZE']
    if size_in_beams == 0.:
        size_in_beams = Configuration['SIZE_IN_BEAMS']
    if debug:
        print_log(f'''CALC_RINGS: Calculating the number of rings in the model.
{'':8s} size in beams = {size_in_beams}
{'':8s} ring_size = {ring_size}
''',Configuration['OUTPUTLOG'],debug = debug)
    est_rings = round((size_in_beams)/(ring_size)+2)
    if est_rings > 20 and Configuration['MAX_RINGS'] > 25:
        Configuration['OUTER_RINGS_DOUBLED'] = True
        no_rings = 2+10+round((est_rings-10-2)/2.)
    else:
        Configuration['OUTER_RINGS_DOUBLED'] = False
        no_rings = est_rings
    if debug:
        print_log(f'''CALC_RINGS: The model will have {no_rings} Rings.
''',Configuration['OUTPUTLOG'],debug = False)
    return int(no_rings)

def load_output_catalogue(filename, debug = False):
    Results = {}
    with open(filename) as file:
        ini_line=file.readline()
        tmp =[x.upper().strip() for x in ini_line.split()]
        columns = [f'{tmp[0]}_{tmp[1]}',tmp[2],tmp[3],f'{tmp[4]}_{tmp[5]}_{tmp[6]}_{tmp[7]}']
        for val in columns:
            Results[val] = []
        for line in file.readlines():
            vals = [x.strip() for x in line.split()]
            Results[columns[0]].append(vals[0])
            Results[columns[1]].append(vals[1])
            Results[columns[2]].append(vals[2])
            Results[columns[3]].append(f'{" ".join(vals[3:])}')
    return Results

def load_catalogue(filename, debug = False):
    Catalogue = Proper_Dictionary({})
    tmpfile = open(filename,'r')
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
load_catalogue.__doc__ ='''
;+
; NAME:
;       catalogue(filename)
;
; PURPOSE:
;       Read the FAT input catalogue and write into the a dictionary
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
;          result = dictionary with the read file
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



#Function to read FAT configuration file into a dictionary
def load_config_file(filename, debug = False):
    tmpfile = open(filename, 'r')
    Configuration = Proper_Dictionary({})
    boolean_keys = ['NEW_OUTPUT', 'HANNING','FIX_INCLINATION','FIX_PA','FIX_SDIS','FIX_Z0','WARP_OUTPUT']
    string_keys = ['OUTPUTLOG', 'OUTPUTCATALOGUE','MAINDIR','CATALOGUE']
    integer_keys = ['STARTGALAXY','ENDGALAXY','MAPS_OUTPUT','OPT_PIXELBEAM','FINISHAFTER']
    # Separate the keyword names
    for tmp in tmpfile.readlines():
        if tmp[0] != '#':
        # python is really annoying with needing endlines. Let's strip them here and add them when writing
            add_key = tmp.split('=', 1)[0].strip().upper()
            if add_key in boolean_keys:
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
                Configuration[add_key] = value
            elif add_key in string_keys:
                Configuration[add_key] = tmp.split('=', 1)[1].strip()
            elif add_key in integer_keys:
                Configuration[add_key] = int(tmp.split('=', 1)[1].strip())
            else:
                Configuration[add_key] = float(tmp.split('=', 1)[1].strip())



    #Make the input idiot safe
    if Configuration['MAINDIR'][-1] != '/':
        Configuration['MAINDIR'] = Configuration['MAINDIR']+'/'

    while not os.path.isdir(Configuration['MAINDIR']):
        Configuration['MAINDIR'] = input('''
                    Your main fitting directory ({}) does not exist.
                    Please provide the correct directory.
                    '''.format(Configuration['MAINDIR']))
    while not os.path.exists(Configuration['CATALOGUE']):
        Configuration['CATALOGUE'] = input('''
                    Your input catalogue ({}) does not exist.
                    Please provide the correct file name.
                    '''.format(Configuration['CATALOGUE']))
    #The output catalogue only needs to be in a valid directory as we create it
    output_catalogue_dir = Configuration['OUTPUTCATALOGUE'].split('/')
    if len(output_catalogue_dir) > 1:
        check_dir = '/'.join(output_catalogue_dir[:-1])
        while not os.path.isdir(check_dir):
            check_dir= input('''
                    The directory for your output catalogue ({}) does not exist.
                    Please provide the correct directory name.
                    '''.format(Configuration['OUTPUTCATALOGUE']))
            Configuration['OUTPUTCATALOGUE'] = check_dir+'/'+output_catalogue_dir[-1]


    required_configuration_keys = ['FIX_INCLINATION','FIX_PA','FIX_SDIS','FIX_Z0','HANNING','STARTGALAXY', 'ENDGALAXY', 'TESTING', 'START_POINT','RING_SIZE', 'FINISHAFTER', 'CATALOGUE', 'MAINDIR', 'OUTPUTCATALOGUE', 'OUTPUTLOG', 'NEW_OUTPUT', 'OPT_PIXELBEAM', 'MAPS_OUTPUT','WARP_OUTPUT']

    for key in required_configuration_keys:
        if key not in Configuration:
            if key == 'STARTGALAXY':
                Configuration[key] = 0
            if key == 'FINISHAFTER':
                Configuration[key] = 2
            if key == 'TESTING':
                Configuration[key] = 0
            if key == 'START_POINT': #Previously calle allnew
                if 'ALLNEW' not in Configuration:
                    Configuration[key] = 1
            if key == 'ENDGALAXY':
                Configuration[key] = -1
            if key == 'NEW_OUTPUT':   # Called newresult in the gdl code
                Configuration[key] = True
            if key == 'HANNING':
                if 'VELOCITY_RESOLUTION' in Configuration:
                    Configuration[key] = Configuration['VELOCITY_RESOLUTION']
                else:
                    Configuration[key] = False
            if key == 'RING_SIZE': #Previosuly called RINGSPACING in
                Configuration[key] = 1.1
            if key == 'OPT_PIXELBEAM':
                Configuration[key] = 4
            if key == 'MAPS_OUTPUT': # Previously called bookkeeping
                Configuration[key] = 3
            if key == 'WARP_OUTPUT':
                Configuration[key] = False
            if key == 'OUTPUTLOG':
                Configuration[key] = None
    if Configuration['RING_SIZE'] < 0.5:
        Configuration['RING_SIZE'] = 0.5
    if Configuration['MAPS_OUTPUT'] == 5:
        Configuration['MAPS_OUTPUT'] = 4
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


#Function for loading the variables of a tirific def file into a set of variables to be used
def load_tirific(filename,Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT','VROT_ERR',
                 'Z0', 'Z0_ERR', 'SBR','SBR_ERR', 'INCL','INCL_ERR','PA','PA_ERR','XPOS','YPOS','VSYS','SDIS','SDIS_ERR'
                 ,'VROT_2','VROT_2_ERR',  'Z0_2','Z0_2_ERR','SBR_2','SBR_2_ERR',
                 'INCL_2', 'INCL_2_ERR','PA_2','PA_2_ERR','XPOS_2','YPOS_2','VSYS_2','SDIS_2','SDIS_2_ERR','CONDISP','CFLUX','CFLUX_2'],
                 unpack = True , debug = False ):
    if debug:
        print_log(f'''LOAD_TIRIFIC: Starting to extract the following paramaters:
{'':8s}{Variables}
''',None,screen=True, debug = True)
    Variables = np.array([e.upper() for e in Variables],dtype=str)

    tmp = open(filename, 'r')

    numrings = [int(e.split('=')[1].strip()) for e in tmp.readlines() if e.split('=')[0].strip().upper() == 'NUR']
    tmp.seek(0)
    outputarray=np.zeros((numrings[0],len(Variables)),dtype=float)
    unarranged = tmp.readlines()
    # Separate the keyword names
    for line in unarranged:
        if line.count('=') > 1:
            print(f'This taco is not correct. \n You have lines in  def file {filename} where = occurs multiple times.')
            print(f'This is the offending line {line}')
            exit()
        var_concerned = str(line.split('=')[0].strip().upper())

        if debug:
            print_log(f'''LOAD_TIRIFIC: extracting line
{'':8s}{var_concerned}.
''',None,screen=True, debug = True)
        if len(var_concerned) < 1:
            var_concerned = 'xxx'
        varpos = np.where(Variables == var_concerned)[0]

        if varpos.size > 0:

            try:
                tmp =  np.array(line.split('=')[1].rsplit(),dtype=float)
            except ValueError:
                if var_concerned == 'CONDISP':
                    tmp = line.split('=')[1].rsplit()
                    tmp = np.array([tmp[0][1:]],dtype=float)
            if len(outputarray[:,0]) < len(tmp):
                tmp_out=outputarray
                outputarray = np.zeros((len(tmp), len(Variables)), dtype=float)
                outputarray[0:len(tmp_out),:] = tmp_out
            outputarray[0:len(tmp),int(varpos)] = tmp[0:len(tmp)]
        else:
            if var_concerned[0] == '#':
                varpos = np.where(var_concerned[2:].strip() == Variables)[0]
                if debug:
                    print_log(f'''LOAD_TIRIFIC: comparing {var_concerned[2:].strip()} to the variables.
{'':8s}Found {varpos}.
''',None,screen=True, debug = True)
                if varpos.size > 0:
                    tmp = np.array(line.split('=')[1].rsplit(),dtype=float)
                    if len(outputarray[:, 0]) < len(tmp):
                        tmp_out = outputarray
                        outputarray = np.zeros((len(tmp), len(Variables)), dtype=float)
                        outputarray[0:len(tmp_out), :] = tmp_out
                    outputarray[0:len(tmp),int(varpos)] = tmp[:]

    if unpack:
        return (*outputarray.T,)
    else:
        return outputarray

#batch convert types
def convert_type(array, type = 'float',debug = False):
    if debug:
        print_log(f'''CONVERT_TYPE: Start.
''',None,debug =True)
    if type =='int':
        return  [int(x) for x in array]
    elif type =='str':
        return  [str(x) for x in array]
    else:
        return  [float(x) for x in array]

convert_type.__doc__ = '''
;+
; NAME:
;       def convert_type(array, type = 'float'):
;
; PURPOSE:
;       Convert a list of variables from one type to another and be able to have them unpack into single varaiable again.
;
; CATEGORY:
;       Support
;
; INPUTS:
;     Array to convert
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
;
; MODULES CALLED:
;
;
; EXAMPLE:
;

'''
def get_name_info(name):
    try:
        beam=float(name.split('ba')[1].split('SNR')[0])
        RCShape = name.split('-')[0]
    except:
        if name[0] == 'C':
            RCShape= name.split('_')[0]
            beam =float(name.split('_')[1].split('Beams')[0])
        else:
            RCShape= name.split('_')[0]+'_'+ name.split('_')[1]
            beam =float(name.split('_')[2].split('Beams')[0])
    try:
        SNR=float(name.split('SNR')[1].split('bm')[0])
    except:
        if name[0] == 'C':
            SNR =float(name.split('_')[2].split('SNR')[0])
        else:
            SNR =float(name.split('_')[3].split('SNR')[0])
    return beam,SNR,RCShape


def get_flux_from_info(name):
    with open(name+'-Info.txt') as file:
        for line in file.readlines():
            tmp = line.split()
            if tmp[0] == 'HI_Mass':
                Mass = float(tmp[1])
            if tmp[0] == 'At':
                Distance = float(tmp[4])
    return Mass/(2.36E5*Distance**2) #Flux in Jy*km/s


#Function to convert column densities
# levels should be n mJy/beam when flux is given
def columndensity(levels,systemic = 100.,beam=[1.,1.],channel_width=1.,column= False,arcsquare=False,solar_mass_input =False,solar_mass_output=False, debug = False):
    if debug:
        print_log(f'''COLUMNDENSITY: Starting conversion from the following input.
{'':8s}Levels = {levels}
{'':8s}Beam = {beam}
{'':8s}channel_width = {channel_width}
''',None,debug =True)
    beam=np.array(beam)
    f0 = 1.420405751786E9 #Hz rest freq
    c = 299792.458 # light speed in km / s
    pc = 3.086e+18 #parsec in cm
    solarmass = 1.98855e30 #Solar mass in kg
    mHI = 1.6737236e-27 #neutral hydrogen mass in kg
    if debug:
                print_log(f'''COLUMNDENSITY: We have the following input for calculating the columns.
{'':8s}COLUMNDENSITY: level = {levels}, channel_width = {channel_width}, beam = {beam}, systemic = {systemic})
''',None,debug=debug)
    if systemic > 10000:
        systemic = systemic/1000.
    f = f0 * (1 - (systemic / c)) #Systemic frequency
    if arcsquare:
        HIconv = 605.7383 * 1.823E18 * (2. *np.pi / (np.log(256.)))
        if column:
            # If the input is in solarmass we want to convert back to column densities
            if solar_mass_input:
                levels=levels*solarmass/(mHI*pc**2)
            #levels=levels/(HIconv*channel_width)
            levels = levels/(HIconv*channel_width)
        else:

            levels = HIconv*levels*channel_width
            if solar_mass_output:
                levels=levels*mHI/solarmass*pc*pc
    else:
        if beam.size <2:
            beam= [beam,beam]
        b=beam[0]*beam[1]
        if column:
            if solar_mass_input:
                levels=levels*solarmass/(mHI*pc**2)
            TK = levels/(1.823e18*channel_width)
            levels = TK/(((605.7383)/(b))*(f0/f)**2)
        else:
            TK=((605.7383)/(b))*(f0/f)**2*levels
            levels = TK*(1.823e18*channel_width)
    if ~column and solar_mass_input:
        levels = levels*mHI*pc**2/solarmass
    return levels
columndensity.__doc__ = '''
;+
; NAME:
;       columndensity(levels,systemic = 100.,beam=[1.,1.],channel_width=1.,column= False,arcsquare=False,solar_mass_input =False,solar_mass_output=False)
;
; PURPOSE:
;       Convert the various surface brightnesses to other values
;
; CATEGORY:
;       Support
;
; INPUTS:
;       levels = the values to convert
;       systemic = the systemic velocity of the source
;        beam  =the beam in arcse
;       channelwidth = width of a channel in km/s
;     column = if true input is columndensities
;
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
;
; MODULES CALLED:
;
;
; EXAMPLE:
;

'''
        # a Function to convert the RA and DEC into hour angle (invert = False) and vice versa (default)
def convertRADEC(RAin,DECin,invert=False, colon=False,debug = False):
    if debug:
            print_log(f'''CONVERTRADEC: Starting conversion from the following input.
    {'':8s}RA = {RAin}
    {'':8s}DEC = {DECin}
''',None,debug =True)
    RA = copy.deepcopy(RAin)
    DEC = copy.deepcopy(DECin)
    if not invert:
        try:
            _ = (e for e in RA)
        except TypeError:
            RA= [RA]
            DEC =[DEC]
        for i in range(len(RA)):
            xpos=RA
            ypos=DEC
            xposh=int(np.floor((xpos[i]/360.)*24.))
            xposm=int(np.floor((((xpos[i]/360.)*24.)-xposh)*60.))
            xposs=(((((xpos[i]/360.)*24.)-xposh)*60.)-xposm)*60
            yposh=int(np.floor(np.absolute(ypos[i]*1.)))
            yposm=int(np.floor((((np.absolute(ypos[i]*1.))-yposh)*60.)))
            yposs=(((((np.absolute(ypos[i]*1.))-yposh)*60.)-yposm)*60)
            sign=ypos[i]/np.absolute(ypos[i])
            if colon:
                RA[i]="{}:{}:{:2.2f}".format(xposh,xposm,xposs)
                DEC[i]="{}:{}:{:2.2f}".format(yposh,yposm,yposs)
            else:
                RA[i]="{}h{}m{:2.2f}".format(xposh,xposm,xposs)
                DEC[i]="{}d{}m{:2.2f}".format(yposh,yposm,yposs)
            if sign < 0.: DEC[i]='-'+DEC[i]
        if len(RA) == 1:
            RA = str(RA[0])
            DEC = str(DEC[0])
    else:
        if isinstance(RA,str):
            RA=[RA]
            DEC=[DEC]

        xpos=RA
        ypos=DEC

        for i in range(len(RA)):
            # first we split the numbers out
            tmp = re.split(r"[a-z,:]+",xpos[i])
            RA[i]=(float(tmp[0])+((float(tmp[1])+(float(tmp[2])/60.))/60.))*15.
            tmp = re.split(r"[a-z,:'\"]+",ypos[i])
            if float(tmp[0]) != 0.:
                DEC[i]=float(np.absolute(float(tmp[0]))+((float(tmp[1])+(float(tmp[2])/60.))/60.))*float(tmp[0])/np.absolute(float(tmp[0]))
            else:
                DEC[i] = float(np.absolute(float(tmp[0])) + ((float(tmp[1]) + (float(tmp[2]) / 60.)) / 60.))
                if tmp[0][0] == '-':
                    DEC[i] = float(DEC[i])*-1.
        if len(RA) == 1:
            RA= float(RA[0])
            DEC = float(DEC[0])
        else:
            RA =np.array(RA,dtype=float)
            DEC = np.array(DEC,dtype=float)
    return RA,DEC


# function for converting kpc to arcsec and vice versa

def convertskyangle(angle, distance=1., unit='arcsec', distance_unit='Mpc', physical=False,debug = False):
    if debug:
            print_log(f'''CONVERTSKYANGLE: Starting conversion from the following input.
    {'':8s}Angle = {angle}
    {'':8s}Distance = {distance}
''',None,debug =True)
    try:
        _ = (e for e in angle)
    except TypeError:
        angle = [angle]

        # if physical is true default unit is kpc
    angle = np.array(angle)
    if physical and unit == 'arcsec':
        unit = 'kpc'
    if distance_unit.lower() == 'mpc':
        distance = distance * 10 ** 3
    elif distance_unit.lower() == 'kpc':
        distance = distance
    elif distance_unit.lower() == 'pc':
        distance = distance / (10 ** 3)
    else:
        print('CONVERTSKYANGLE: ' + distance_unit + ' is an unknown unit to convertskyangle.\n')
        print('CONVERTSKYANGLE: please use Mpc, kpc or pc.\n')
        raise SupportRunError('CONVERTSKYANGLE: ' + distance_unit + ' is an unknown unit to convertskyangle.')
    if not physical:
        if unit.lower() == 'arcsec':
            radians = (angle / 3600.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'arcmin':
            radians = (angle / 60.) * ((2. * np.pi) / 360.)
        elif unit.lower() == 'degree':
            radians = angle * ((2. * np.pi) / 360.)
        else:
            print('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.\n')
            print('CONVERTSKYANGLE: please use arcsec, arcmin or degree.\n')
            raise SupportRunError('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.')


        kpc = 2. * (distance * np.tan(radians / 2.))
    else:
        if unit.lower() == 'kpc':
            kpc = angle
        elif unit.lower() == 'mpc':
            kpc = angle / (10 ** 3)
        elif unit.lower() == 'pc':
            kpc = angle * (10 ** 3)
        else:
            print('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.\n')
            print('CONVERTSKYANGLE: please use kpc, Mpc or pc.\n')
            raise SupportRunError('CONVERTSKYANGLE: ' + unit + ' is an unknown unit to convertskyangle.')

        radians = 2. * np.arctan(kpc / (2. * distance))
        kpc = (radians * (360. / (2. * np.pi))) * 3600.
    if len(kpc) == 1:
        kpc = float(kpc[0])
    return kpc

def finish_current_run(Configuration,current_run,debug=False):
    print_log(f"FINISH_CURRENT_RUN: Is Tirific Running? {Configuration['TIRIFIC_RUNNING']}. \n",Configuration['OUTPUTLOG'],debug=debug,screen=True)
    if Configuration['TIRIFIC_RUNNING']:
        try:
            os.kill(Configuration['TIRIFIC_PID'], signal.SIGKILL)
            print_log(f"FINISH_CURRENT_RUN: We killed PID = {Configuration['TIRIFIC_PID']}. \n",Configuration['OUTPUTLOG'],debug=debug,screen=True)
        except:
            try:
                current_run.kill()
                print_log(f"FINISH_CURRENT_RUN: We killed the current run although we failed on the PID = {Configuration['TIRIFIC_PID']}. \n",Configuration['OUTPUTLOG'],debug=debug,screen=True)
            except AttributeError:
                print_log(f"FINISH_CURRENT_RUN: We failed to kill the current run even though we have tirific running",Configuration['OUTPUTLOG'],debug=debug,screen=True)
                raise TirificRunError('FINISH_CURRENT_RUN: Despite having an initialized tirific we could not kill it.')
        Configuration['TIRIFIC_RUNNING'] = False
        Configuration['TIRIFIC_PID'] = 'Not Initialized'
    else:
        print_log(f"FINISH_CURRENT_RUN: No run is initialized.",Configuration['OUTPUTLOG'],debug=debug,screen=True)


def gaussian_function(axis,peak,center,sigma):
    return peak*np.exp(-(axis-center)**2/(2*sigma**2))

def fit_gaussian(x,y, covariance = False,debug = False):
    if debug:
            print_log(f'''FIT_GAUSSIAN: Starting to fit a Gaussian.
{'':8s}x = {x}
{'':8s}y = {y}
''',None,debug =True)
    # First get some initial estimates
    est_peak = np.nanmax(y)
    est_center = float(x[np.where(y == est_peak)])
    est_sigma = np.nansum(y*(x-est_center)**2)/np.nansum(y)
    gauss_parameters, gauss_covariance = curve_fit(gaussian_function, x, y,p0=[est_peak,est_center,est_sigma])
    if covariance:
        return gauss_parameters, gauss_covariance
    else:
        return gauss_parameters
#Put template values in a list !!!!!!!! This is very similar to load_template in read_funtcions maybe use one?
def get_from_template(Tirific_Template,Variables, debug = False):
    out = []
    if debug:
        print(f'''{'':8s}GET_FROM_TEMPLATE: Trying to get the following profiles {Variables}
''')
    for key in Variables:
        out.append([float(x) for x  in Tirific_Template[key].split()])
    #Because lists are stupid i.e. sbr[0][0] = SBR[0], sbr[1][0] = SBR_2[0] but  sbr[:][0] = SBR[:] not SBR[0],SBR_2[0] as logic would demand
    if debug:
        print(f'''{'':8s}GET_FROM_TEMPLATE: We extracted the following profiles from the Template.
{'':8s}GET_FROM_TEMPLATE: {out}
''' )
    #Beware that lists are stupid i.e. sbr[0][0] = SBR[0], sbr[1][0] = SBR_2[0] but  sbr[:][0] = SBR[:] not SBR[0],SBR_2[0] as logic would demand
    # However if you make a np. array from it make sure that you specify float  or have lists of the same length else you get an array of lists which behave just as dumb
    return out

# Function to calculate the difference between model and fit
def get_diff(val,model, radii = [], model_radii = [], single = False, errors = [],second = [] ,second_model = [], second_errors = [],norm=1.):
    to_use = np.where(val > 0.)[0]
    model_to_use = np.where(model > 0.)[0]
    if len(model_radii) > 0 and len(radii) > 0:
        model_int_function = interpolate.interpolate.interp1d(model_radii[model_to_use],model[model_to_use],fill_value="extrapolate")
        model = model_int_function(radii)
        if len(second) > 0:
            if np.sum(second_model) == 0.:
                second_model = model
            second_to_use = np.where(second > 0.)[0]
            second_model_to_use = np.where(second_model > 0.)[0]
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
        error = np.std(difference)
    if value > 1e8:
        print(f'In values = {val}')
        print(f'In model = {model}')
        if len(second) > 0:
            print(f'In values2 = {second}')
            print(f'In model2 = {second_model}')
        print(f'In difference = {difference}')
        print(f'In errors = {errors}')
        print(f'norm = {norm} , value = {value}, error = {error}')
        exit()
    return value,error

# Function to get the amount of inner rings to fix
def get_inner_fix(Configuration,Tirific_Template, debug =False):
    if debug:
        print_log(f'''GET_INNER_FIX: Attempting to get the inner rings to be fixed.
''',Configuration['OUTPUTLOG'], debug = debug, screen = True)
    sbr_av = np.array([(float(x)+float(y))/2. for x,y  in zip(Tirific_Template['SBR'].split(),Tirific_Template['SBR_2'].split())],dtype = float)
    column_levels = columndensity(sbr_av, arcsquare = True, debug = debug)
    column_levels[0]= 1e21
    tmp = np.where(column_levels > 1e20)[0]
    return set_limits(int(np.floor(tmp[-1]/1.5-1)), 4, int(Configuration['NO_RINGS']*0.9))

def get_usage_statistics(process_id, debug = False):
    result = subprocess.check_output(['top',f'-p {process_id}','-d 1','-n 1'])
    #result = subprocess.check_output(['ps','u'])
    lines = result.decode('utf8').split('\n')
    column_names = [x.upper() for x in lines[6].strip().split()]
    if debug:
        print(f'''{'':8s}GET_usage_statistics: We extracted the following column names {column_names}
''')
    CPU = float(0.)
    mem=float(0.)
    column_var = [x for x in lines[7].strip().split()]
    if debug:
        print(f'''{'':8s}GET_usage_statistics: We extracted the following variables {column_var}
''')
    CPU = float(column_var[column_names.index('%CPU')])
    mem = float(column_var[column_names.index('RES')])/1024**2
    try:
        if int(column_var[column_names.index('PID')]) == int(process_id):
            CPU = float(column_var[column_names.index('%CPU')])
            mem = float(column_var[column_names.index('RES')])/1024**2
    except:
        #if the PID is not numeric it got merged with the crap in shiftcentercounter
        try:
            if column_var[column_names.index('COMMAND')-1] == 'tirific':
                CPU = float(column_var[column_names.index('%CPU')-1])
                mem = float(column_var[column_names.index('RES')-1])/1024**2
        except:
            pass

    return CPU,mem
# A simple function to return the line numbers in the stack from where the functions are called
def linenumber(debug=False):
    line = []
    for key in stack():
        if key[1] == 'FAT.py':
            break
        if key[3] != 'linenumber' and key[3] != 'print_log':
            file = key[1].split('/')
            line.append(f"In the function {key[3]} at line {key[2]} in file {file[-1]}")
    if len(line) > 0:
        if debug:
            line = ', '.join(line)+f'\n{"":8s}'
        else:
            line = f'{"":8s}'
    else:
        for key in stack():
            if key[1] == 'FAT.py':
                line = f"{'('+str(key[2])+')':8s}"
                break
    return line
def make_plot(x,y, color= None, status= None, location = [0,1], symbol= None,
                    xlabel = '',ylabel = '', colorbarlabel = '', legend = False, No_Mean = False):
        x = np.array(x, dtype=float)
        y = np.array(y, dtype=float)
        ax = plt.subplot(location)

        if not status is None:
            status = np.array(status)

            succes = np.where(status > 0.1)[0]
            mean = np.nanmean(y[succes,0])
            stdev = np.nanstd(y[succes,0]-mean)
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
            mean = np.nanmean(y[:,0])
            stdev = np.nanstd(y[:,0]-mean)
            transparency =np.ones(len(x[:,0]))
        symlist = ["o", "v", "^", "<",">","s","P","*","X","D","1","3","$a$","$b$","$c$","$d$","$e$","$f$","$g$","$h$"]


        if not symbol is None:
            symbol= np.array(symbol)
            req_no_elements = np.unique(symbol)
            symbol_use = [symlist[i] for i,shape in enumerate(req_no_elements)]
        else:
            req_no_elements = ['Unspecified']
            symbol = np.array(['Unspecified' for gh in x[:,0]])
            symbol_use = ['o']
        if not color is None:
            color = np.array(color,dtype=float)
            color = color/90.
            cmap = plt.cm.get_cmap('rainbow')
            rgba_cols = [cmap(color)]

        else:
            color = np.zeros(len(x[:,0]))
        for i,shaped in enumerate(req_no_elements):
            proc = np.where(symbol == shaped)[0]
            for j,transparency_val in enumerate(norm_elements):
                add = np.where(transparency[proc] == transparency_val)[0]
                if len(add) > 0:
                    lab_string = f'RC Shape {shaped}, fit quality {stat_elements[j]}'
                    siz = (stat_elements[j]+2)**4.
                    try:
                        ax.scatter(x[proc[add],0],y[proc[add],0],cmap= 'rainbow', c = rgba_cols[0][proc[add]][:], s=siz, marker = symbol_use[i],alpha = norm_elements[j],label = lab_string)
                        plt.errorbar(x[proc[add],0],y[proc[add],0],xerr=x[proc[add],1],yerr=y[proc[add],1], linestyle="None", ecolor = rgba_cols[0][proc[add]][:],alpha = norm_elements[j])
                    except:
                        ax.scatter(x[proc[add]],y[proc[add],0],cmap= 'rainbow', c = rgba_cols[0][proc[add]][:], s=siz, marker = symbol_use[i],alpha = norm_elements[j],label = lab_string)
                        plt.errorbar(x[proc[add]],y[proc[add],0],xerr=np.zeros(len(add)),yerr=y[proc[add],1], linestyle="None", ecolor = rgba_cols[0][proc[add]][:],alpha = norm_elements[j])


        if not No_Mean:
            xmin,xmax = ax.get_xlim()
            ax.plot([xmin-1,xmax+2.],[mean,mean], c = 'k', alpha= 0.5)
            ax.plot([xmin-1,xmax+2.],[mean-stdev,mean-stdev], 'k--', alpha= 0.5)
            ax.plot([xmin-1,xmax+2.],[mean+stdev,mean+stdev], 'k--', alpha= 0.5)
            ax.set_xlim(xmin,xmax)

        ax.text(0.95,0.95,f'Mean = {mean:.1f} $\pm$ {stdev:.1f} ', transform=ax.transAxes,horizontalalignment= 'right', verticalalignment='top')
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        return ax

def print_log(log_statement,log, screen = False,debug = False):
    log_statement = f"{linenumber(debug=debug)}{log_statement}"
    if screen or not log:
        print(log_statement)
    if log:
        with open(log,'a') as log_file:
            log_file.write(log_statement)

print_log.__doc__ = '''
;+
; NAME:
;       print_log
;
; PURPOSE:
;       Print statements to log if existent and screen if Requested
;
; CATEGORY:
;       Support
;
; CALLING SEQUENCE:
;       cleanup(Configuration)
;
; INPUTS:
;     Configuration = Structure that has the current fitting directory and expected steps
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
;
; MODULES CALLED:
;       os
;
; EXAMPLE:
;

'''

def rename_fit_products(Configuration,stage = 'initial', fit_stage='Undefined_Stage',debug = False):
    extensions = ['def','log','ps','fits']
    for filetype in extensions:
        if filetype == 'log':
            if os.path.exists(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.{filetype}"):
                os.system(f"cp {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.{filetype} {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}_Prev.{filetype} ")
        else:
            if filetype == 'def':
                if fit_stage == 'Extent_Convergence' and stage == 'run_ec':
                    Loopnr = f"{Configuration['EC_LOOPS']}"
                elif fit_stage == 'Centre_Convergence' and stage == 'run_cc' :
                    Loopnr = f"{Configuration['CC_LOOPS']}"
                else:
                    Loopnr = 'before_'+stage
                if os.path.exists(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.{filetype}"):
                    os.system(f"mv {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.{filetype} {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}_{Loopnr}.{filetype}")

            elif os.path.exists(f"{Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.{filetype}"):
                os.system(f"mv {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}.{filetype} {Configuration['FITTING_DIR']}{fit_stage}/{fit_stage}_Prev.{filetype}")


rename_fit_products.__doc__ = '''
;+
; NAME:
;       rename_fit_products(Configuration,stage)
;
; PURPOSE:
;       rename the tirific product from the previous stage.
;
; CATEGORY:
;       Support
;
;
; INPUTS:
;     Configuration, and the cube header
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
;
; MODULES CALLED:
;
;
; EXAMPLE:
;

'''

def remove_inhomogeneities(Configuration,fits_map,inclination=30., pa = 90. , center = [0.,0.],WCS_center = True, debug=False):
    if debug:
        print_log(f'''REMOVE_INHOMOGENEITIES: These are the values we get as input
{'':8s}Inclination = {inclination}
{'':8s}pa = {pa}
{'':8s}center = {center}
''',Configuration['OUTPUTLOG'],debug = True)
    map = fits_map[0].data
    # first rotate the pa
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        map_wcs = WCS(fits_map[0].header)
    # convert the boundaries to real coordinates
    if WCS_center:
        x, y = map_wcs.wcs_world2pix(*center, 0.)
    else:
        x = center[0]
        y = center[1]
    rot_map = rotateImage(map,pa-90,[x,y],debug=debug)
    # deproject
    dep_map = deproject(copy.deepcopy(rot_map),inclination,center = y, debug=debug)

    if debug:
        fits.writeto(f"{Configuration['FITTING_DIR']}rot_map.fits",rot_map,fits_map[0].header,overwrite = True)
        fits.writeto(f"{Configuration['FITTING_DIR']}dep_map.fits",dep_map,fits_map[0].header,overwrite = True)

    angles = np.linspace(5.,360.,71)
    minimum_map = copy.deepcopy(dep_map)
    for angle in angles:
        rot_dep_map =  rotateImage(copy.deepcopy(dep_map),angle,[x,y],debug=debug)

        #tmp = np.where(rot_dep_map < minimum_map)[0]
        minimum_map[rot_dep_map < minimum_map] =rot_dep_map[rot_dep_map < minimum_map]
    clean_map = rotateImage(deproject(copy.deepcopy(minimum_map),inclination,center = y,invert= True,debug=debug),-1*(pa-90),[x,y],debug=debug)

    if debug:
        fits.writeto(f"{Configuration['FITTING_DIR']}minimum_map.fits",minimum_map,fits_map[0].header,overwrite = True)
        fits.writeto(f"{Configuration['FITTING_DIR']}clean_map.fits",clean_map,fits_map[0].header,overwrite = True)
    fits_map[0].data = clean_map
    return fits_map

def deproject(map,angle, center = 0., invert = False,debug=False):
    axis = range(len(map[:,0]))-center
    if invert:
        newaxis = axis/np.cos(np.radians(angle))
    else:
        newaxis = axis*np.cos(np.radians(angle))
    for x in range(len(map[0,:])):
        profile = copy.deepcopy(map[:,x])
        new_profile = np.interp(np.array(newaxis,dtype=float),np.array(axis,dtype=float),np.array(profile,dtype=float))
        map[:,x] = new_profile
    return map


#function to rotate a cube without losing info
def rotateImage(Cube, angle, pivot, debug = False):
    padX = [int(Cube.shape[1] - pivot[0]), int(pivot[0])]
    padY = [int(Cube.shape[0] - pivot[1]), int(pivot[1])]
    imgP = np.pad(Cube, [padY, padX], 'constant')
    imgR = ndimage.rotate(imgP, angle, axes=(1, 0), reshape=False)
    return imgR[padY[0]: -padY[1], padX[0]: -padX[1]]

def sbr_limits(Configuration,hdr, systemic= 100. , debug = False):
    radii = set_rings(Configuration,debug=debug)
    if debug:
        print_log(f'''SBR_LIMITS: Got {len(radii)} radii
''',Configuration['OUTPUTLOG'], debug=debug,screen =True)
    level = hdr['FATNOISE']*1000
    bm = [hdr['BMAJ']*3600.,hdr['BMIN']*3600.]
    noise_in_column = columndensity(level,beam = bm,systemic = systemic,channel_width=hdr['CDELT3']/1000.)
    J2007col=9.61097e+19
    ratio=(noise_in_column/J2007col)**0.5
    beamsolid=(np.pi*bm[0]*bm[1])/(4.*np.log(2.))
    ringarea= [0 if radii[0] == 0 else np.pi*((radii[0]+radii[1])/2.)**2]
    ringarea = np.hstack((ringarea,
                         [np.pi*(((y+z)/2.)**2-((y+x)/2.)**2) for x,y,z in zip(radii,radii[1:],radii[2:])],
                         [np.pi*((radii[-1]+0.5*(radii[-1]-radii[-2]))**2-((radii[-1]+radii[-2])/2.)**2)]
                         ))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sbr_ring_limits=9e-4*(ringarea/beamsolid)**(-0.82)*ratio
    if ringarea[0] == 0.:
         sbr_ring_limits[0]=0.
    if len(Configuration['LIMIT_MODIFIER']) == 1:

        sbr_ring_limits= sbr_ring_limits*Configuration['LIMIT_MODIFIER']
    else:
        mod_list = list(Configuration['LIMIT_MODIFIER'])
        while len(mod_list) < len(sbr_ring_limits):
            mod_list.append(Configuration['LIMIT_MODIFIER'][-1])
        Configuration['LIMIT_MODIFIER'] = np.array(mod_list,dtype=float)
        sbr_ring_limits=[x*y for x,y in zip(sbr_ring_limits,Configuration['LIMIT_MODIFIER'])]
    if debug:
        print_log(f'''SBR_LIMITS: Retrieved these radii and limits:
{'':8s}{radii}
{'':8s}{sbr_ring_limits}
''',Configuration['OUTPUTLOG'], debug=False,screen =True)
    return radii,sbr_ring_limits


sbr_limits.__doc__ = '''
;+
; NAME:
;       sbr_limits(Configuration,hdr)
;
; PURPOSE:
;       Create the radii at which to evaluate  the model and a corresponding array that has the sbr reliability sbr_limits
;
; CATEGORY:
;       Support
;
;
; INPUTS:
;     Configuration, and the cube header
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
;
; MODULES CALLED:
;
;
; EXAMPLE:
;

'''

def set_limits(value,minv,maxv,debug = False):
    if value < minv:
        return minv
    elif value > maxv:
        return maxv
    else:
        return value

set_limits.__doc__ = '''
;+
; NAME:
;       set_limits(value,min,max)
;
; PURPOSE:
;       Make sure Value is between min and max else set to min when smaller or max when larger.
;
; CATEGORY:
;       Support
;
;
; INPUTS:
;     value
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
;
; MODULES CALLED:
;
;
; EXAMPLE:
;

'''
#simple function keep track of how to modify the edge limits
def set_limit_modifier(Configuration,Inclination, debug= False):
    if debug:
        print_log(f'''SET_LIMIT_MODIFIER: Checking the limit modifier.
''', Configuration['OUTPUTLOG'], debug=debug)
    if not Inclination.shape:
        Inclination = [Inclination]
    modifier_list = []
    for inc in Inclination:
        if 40 < inc < 50:
            modifier_list.append(set_limits(1.+(50-(inc)*0.05),1,2.5))
        elif inc < 40:
            modifier_list.append(set_limits(np.sin(np.radians(75.))/np.sin(np.radians(inc)),1.,2.5))
        else:
            modifier_list.append(1.)
    if Configuration['OUTER_RINGS_DOUBLED']:
        if len(modifier_list) > 10:
            modifier_list[10:]= np.sqrt(modifier_list[10:])
    Configuration['LIMIT_MODIFIER'] = np.array(modifier_list,dtype=float)
    print_log(f'''SET_LIMIT_MODIFIER: We updated the LIMIT_MODIFIER to {Configuration['LIMIT_MODIFIER']}.
''', Configuration['OUTPUTLOG'], debug=debug)


def set_ring_size(Configuration, debug = False, size_in_beams = 0., check_set_rings = True):

    if size_in_beams == 0.:
        check_set_rings = False
        size_in_beams =  Configuration['SIZE_IN_BEAMS']
    ring_size = Configuration['RING_SIZE']
    no_rings = calc_rings(Configuration,ring_size=ring_size,size_in_beams=size_in_beams,debug=debug)
    if debug:
        print_log(f'''SET_RING_SIZE: Starting with the following parameters.
{'':8s}SIZE_IN_BEAMS = {size_in_beams}
{'':8s}RING_SIZE = {ring_size}
''', Configuration['OUTPUTLOG'],debug=debug)

    while ring_size > 0.5 and  no_rings < Configuration['MINIMUM_RINGS']:
        previous_ringsize = ring_size
        ring_size = set_limits(ring_size/1.5,0.5,float('NaN'),debug=debug)
        no_rings = calc_rings(Configuration,ring_size=ring_size,size_in_beams=size_in_beams,debug=debug)
        print_log(f'''SET_RING_SIZE: Because we had less than four rings we have reduced the ring size from {previous_ringsize} to {ring_size}
''',Configuration['OUTPUTLOG'],debug = debug)

    while no_rings < Configuration['MINIMUM_RINGS'] and size_in_beams !=  Configuration['MAX_SIZE_IN_BEAMS']:
        size_in_beams = set_limits(size_in_beams+1.*ring_size,1, Configuration['MAX_SIZE_IN_BEAMS'])
        no_rings = calc_rings(Configuration,ring_size=ring_size,size_in_beams=size_in_beams,debug=debug)
        print_log(f'''SET_RING_SIZE: The initial estimate is too small to fit adding a ring to it.
''',Configuration['OUTPUTLOG'],debug = False)

    if check_set_rings:
        return size_in_beams,ring_size
    else:
        if debug:
            print_log(f'''SET_RING_SIZE: Setting the following parameters.
{'':8s}SIZE_IN_BEAMS = {size_in_beams}
{'':8s}RING_SIZE = {ring_size}
{'':8s}NO_RINGS = {no_rings}
''', Configuration['OUTPUTLOG'],debug=False)
        Configuration['NO_RINGS'] = int(no_rings)
        Configuration['SIZE_IN_BEAMS'] = size_in_beams
        Configuration['RING_SIZE'] = ring_size
        if Configuration['NO_RINGS'] < Configuration['MINIMUM_RINGS']:
            print_log(f'''SET_RING_SIZE: With a ring size of {Configuration['RING_SIZE']} we still only find {Configuration['NO_RINGS']}.
    {"":8s}SET_RING_SIZE: This is not enough for a fit.
    ''',Configuration['OUTPUTLOG'],debug = False)
            raise SmallSourceError('This source is too small to reliably fit.')


def set_rings(Configuration,ring_size = 0.,size_in_beams = 0. , debug = False):
    if ring_size == 0.:
        ring_size = Configuration['RING_SIZE']
    if size_in_beams == 0.:
        size_in_beams = Configuration['SIZE_IN_BEAMS']
    if debug:
        print_log(f'''SET_RINGS: Starting with the following parameters.
{'':8s}SIZE_IN_BEAMS = {size_in_beams}
{'':8s}RING_SIZE = {ring_size}
''', Configuration['OUTPUTLOG'],debug=debug)
    no_rings = calc_rings(Configuration,debug=debug)
    if debug:
        print_log(f'''SET_RINGS: We find {no_rings} rings.
''', Configuration['OUTPUTLOG'],debug=False)
    #Configuration['NO_RINGS'] = Configuration['SIZE_IN_BEAMS']/Configuration['RING_SIZE']
    if Configuration['OUTER_RINGS_DOUBLED']:
        print_log(f'''SET_RINGS: This is a large galaxy (Size = {size_in_beams}) Therefore we use twice the ring size in the outer parts.
''',Configuration['OUTPUTLOG'],debug = Configuration['DEBUG'])
        radii = [0.,1./5.*Configuration['BMMAJ']]
        radii = np.hstack((radii,(np.linspace(Configuration['BMMAJ']*ring_size,Configuration['BMMAJ']*10.*ring_size, \
                                        10)+1./5*Configuration['BMMAJ'])))
        radii = np.hstack((radii,(np.linspace(Configuration['BMMAJ']*11.*ring_size, \
                                              Configuration['BMMAJ']*(size_in_beams), \
                                              no_rings-12) \
                                              +1./5*Configuration['BMMAJ'])))
    else:
        radii = [0.,1./5.*Configuration['BMMAJ']]
        radii = np.hstack((radii,(np.linspace(Configuration['BMMAJ']*ring_size,Configuration['BMMAJ']*size_in_beams, \
                                        no_rings-2)+1./5.*Configuration['BMMAJ'])))
    if debug:
        print_log(f'''SET_RINGS: Got the following radii.
{'':8s}{radii}
{'':8s}We should have {Configuration['NO_RINGS']} rings and have {len(radii)} rings.
''', Configuration['OUTPUTLOG'],debug=False)
    #Configuration['NO_RINGS'] = len(radii)
    #Configuration['SIZE_IN_BEAMS']= int((radii[-1]-1./5.*bmaj)/bmaj)
    return np.array(radii,dtype = float)
set_rings.__doc__ = '''
;+
; NAME:
;       set_rings(Configuration,hdr)
;
; PURPOSE:
;       make an array that chas all the ring radii in ARCSEC
;
; CATEGORY:
;       Support
;
;
; INPUTS:
;     Configuration, and the cube header
;
; OPTIONAL INPUTS:
;
;
; KEYWORD PARAMETERS:
;       -
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;       -
;
; MODULES CALLED:
;
;
; EXAMPLE:
;

'''

def sofia_output_exists(Configuration,Fits_Files, debug = False):
    if debug:
        print_log(f'''SOFIA_OUTPUT_EXISTS: Starting check
''', Configuration['OUTPUTLOG'],debug = debug)
    req_files= ['MOMENT1','MOMENT0','MOMENT2','MASK']
    for file in req_files:
        if os.path.exists(Configuration['FITTING_DIR']+'Sofia_Output/'+Fits_Files[file]):
            continue
        else:
            log_statement = f"CHECK_SOFIA_OUTPUT: The file {Configuration['FITTING_DIR']+'Sofia_Output/'+Fits_Files[file]} is not found."
            print_log(log_statement, Configuration['OUTPUTLOG'],debug = debug)
            raise FileNotFoundError(log_statement)

    if not os.path.exists(Configuration['FITTING_DIR']+'Sofia_Output/'+Configuration['BASE_NAME']+'_cat.txt'):
        log_statement = f"CHECK_SOFIA_OUTPUT: The file {Configuration['FITTING_DIR']+'Sofia_Output/'+Configuration['BASE_NAME']+'_cat.txt'} is not found."
        print_log(log_statement, Configuration['OUTPUTLOG'],debug = debug)
        raise FileNotFoundError(log_statement)


sofia_output_exists.__doc__ =f'''
Simple function to make sure all sofia output is present as expeceted
'''
