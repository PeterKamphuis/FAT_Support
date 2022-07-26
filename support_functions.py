#!/usr/local/bin/ python3
# This module contains a set of functions and classes that are used in several different Python scripts in the Database.


from collections import OrderedDict #used in Proper_Dictionary
from inspect import getframeinfo,stack
from scipy.optimize import curve_fit
from scipy import ndimage
from scipy import interpolate
from astropy.wcs import WCS
from astropy.io import fits
from omegaconf import OmegaConf
from pyFAT_astro.config.defaults import defaults
from pyFAT_astro.Support.support_functions import convertskyangle,columndensity
from pyFAT_astro.Support.fat_errors import BadCatalogueError
import os
import signal
import traceback
import numpy as np
import copy
import warnings
import time

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



def analyze(Database_Directory,config_file, basename = 'Analysis_Output', LVHIS= False,GDL=False):
        #load the configuration file
    database_config = load_config_file(
        f'{Database_Directory}/{config_file}')
    #load the input catalogue
    database_inp_catalogue = load_input_catalogue(database_config['INPUT_CATALOGUE'],GDL=GDL)
    #end the results
    database_out_catalogue = load_output_catalogue(
        database_config['OUTPUT_CATALOGUE'],binary=GDL)

    ####

    deltas, RCs = retrieve_deltas_and_RCs(database_config,
            database_inp_catalogue, database_out_catalogue,binary=GDL,LVHIS=LVHIS)
    if not os.path.exists(f'''{database_config['MAIN_DIRECTORY']}/{basename}'''):
        os.system(f'''mkdir {database_config['MAIN_DIRECTORY']}/{basename}''')
    plot_overview(database_config,deltas,filename=f'''{database_config['MAIN_DIRECTORY']}/{basename}/Release_All''',LVHIS=LVHIS)
    plot_RCs(RCs,LVHIS=LVHIS,filename=f'''{database_config['MAIN_DIRECTORY']}/{basename}/RCs''')


def plot_RCs(RCs, filename='RCs',LVHIS =False):

    labelfont= {'family':'Times New Roman',
                'weight':'normal',
                'size':24}
    plt.rc('font',**labelfont)

    plotsize = 8
    no_plots = len([x for x in RCs])

    length = int(np.ceil(np.sqrt(no_plots)))

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
        failed = 0
        tot = 0
        for indi in RCs[key]:
            print(f'Plotting the actual galaxy {indi}')


            if indi not in ['MODEL','MODEL_2']:
                tot += 1
            kpcradii = np.array(convertskyangle({'OUTPUTLOG': None},RCs[key][indi]['RADIUS'],distance=float(RCs[key][indi]['DISTANCE'][0])))
            print(f''' Plotting the RC {RCs[key][indi]['RC']} with:
radi (arcsec) = {RCs[key][indi]['RADIUS']}
radi (kpc) = {kpcradii}
distance = {float(RCs[key][indi]['DISTANCE'][0])}''')
            if indi == 'MODEL':
                ax.plot(kpcradii, RCs[key][indi]['RC'], 'b',zorder= 2)
                ax.plot(kpcradii, RCs[key][indi]['RC'], 'bo',zorder=2)
            elif indi == 'MODEL_2':
                #print(np.isnan(np.array(RCs[key][indi]['RC'],dtype=float)),all(np.isnan(np.array(RCs[key][indi]['RC'],dtype=float))))
                if not np.sum(RCs[key][indi]['RC']) == 0.:
                    ax.plot(kpcradii, RCs[key][indi]['RC'], 'r',zorder= 2)
                    ax.plot(kpcradii, RCs[key][indi]['RC'], 'ro',zorder=2)
            else:
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
        ax.set_title(key)

        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin,ymax+(ymax-ymin)/10.)
        if not LVHIS:
            ax.text(0.95,0.95,f'Out of {tot} galaxies, {failed} failed to fit. ', transform=ax.transAxes,horizontalalignment= 'right', verticalalignment='top')
        else:

            if RCs[key][indi]['STATUS'] == 1:
                print(f"Trying a marker")
                ax.scatter(0.95,0.95,marker='*',color='k', s=137,transform=ax.transAxes)
    if LVHIS:
        version= 'LVHIS'
    else:
        version='Database'

    plt.savefig(f'{filename}_RC_{version}.png', bbox_inches='tight')

    plt.close()


def plot_overview(config,deltas,filename='Overview_Difference',LVHIS=False):
    plot_assembly = {'PLOT_1':{
                        'WINDOW_0': {'LOCATION': 0,
                                     'X': [deltas['XPOS'], '$\Delta$ RA (beams)' ],
                                     'Y': [deltas['YPOS'],'$\Delta$ DEC (beams)'],
                                     'NO_MEAN': True,
                                     'PATCH':  Ellipse(xy=[np.nanmean(np.array([x[0] for x in deltas['XPOS']])),\
                                                        np.nanmean(np.array([x[0] for x in deltas['YPOS']]))],\
                                               width=np.nanstd(np.array([x[0] for x in deltas['XPOS']])) ,\
                                               height=np.nanstd(np.array([x[0] for x in deltas['YPOS']])), angle=0,\
                                            edgecolor='none', alpha=0.6, lw=4, facecolor='k', hatch = '////', zorder=-1) },
                        'WINDOW_1': {'LOCATION': 1,
                                     'X': [deltas['BEAMS_ACROSS'], 'Diameter (beams)'],
                                     'Y': [[[np.sqrt(float(x[0])**2+float(y[0])**2),\
                                                 np.sqrt((float(x[0])**2*float(x[1])**2+\
                                                 float(y[0])**2*float(y[1])**2)/\
                                                 (float(x[0])**2+float(y[0])**2))] \
                                                 for x,y in zip(deltas['XPOS'],\
                                                 deltas['YPOS'])],\
                                                 '$\Delta$ Central (beams)'] #error from https://www.wolframalpha.com/widgets/view.jsp?id=8ac60957610e1ee4894b2cd58e753
                                     },
                        'WINDOW_2': {'LOCATION': 3,
                                     'X': [deltas['BEAMS_ACROSS'], 'Diameter (beams)'],
                                     'Y': [deltas['VSYS'], r'$\Delta$ ${\rm V_{sys}}$ (channels)']},
                        'WINDOW_3': {'LOCATION': 4,
                                     'X': [deltas['BEAMS_ACROSS'], 'Diameter (beams)'],
                                     'Y': [deltas['INCL'], '$\Delta$ $i$ ($^{\circ}$)']},
                        'WINDOW_4': {'LOCATION': 6,
                                     'X': [deltas['BEAMS_ACROSS'], 'Diameter (beams)'],
                                     'Y': [deltas['PA'], '$\Delta$ PA ($^{\circ}$)']},
                        'WINDOW_5': {'LOCATION': 7,
                                     'X': [deltas['BEAMS_ACROSS'], 'Diameter (beams)'],
                                     'Y': [deltas['VROT'], r'$\Delta$ V$_{\rm rot}$  (channels)']},
                        'WINDOW_6': {'LOCATION': 8,
                                     'X': [deltas['BEAMS_ACROSS'], 'Diameter (beams)'],
                                     'Y': [deltas['TOTAL_FLUX'], r'$\Delta$ Tot Flux  (Jy/beam km/s)']}},
                'PLOT_2':{
                        'WINDOW_0': {'LOCATION': 0,
                                     'X': [deltas['BEAMS_ACROSS'], 'Diameter (beams)'],
                                     'Y': [deltas['SDIS'],r'$\Delta$ Dispersion  (channels)']},
                        'WINDOW_1': {'LOCATION': 1,
                                     'X': [deltas['BEAMS_ACROSS'], 'Diameter (beams)'],
                                     'Y': [deltas['Z0'],r'$\Delta$ Scaleheight (beams)']
                                     },
                        'WINDOW_2': {'LOCATION': 3,
                                     'X': [deltas['BEAMS_ACROSS'], 'Diameter (beams)'],
                                     'Y': [deltas['R_HI'],r'$\Delta$ R$_{\rm HI}$  (beams)']},
                        'WINDOW_3': {'LOCATION': 4,
                                     'X': [deltas['R_HI'], r'$\Delta$ R$_{\rm HI}$  (beams)'],
                                     'Y': [deltas['INCL'],r'$\Delta$ $i$ ($^{\circ}$)'],
                                      'NO_MEAN': True,
                                      'PATCH': Ellipse(xy=[np.nanmean(np.array([x[0] for x in deltas['R_HI']])),\
                                                         np.nanmean(np.array([x[0] for x in deltas['INCL']]))],\
                                                width=np.nanstd(np.array([x[0] for x in deltas['R_HI']])) ,\
                                                height=np.nanstd(np.array([x[0] for x in deltas['INCL']])), angle=0,\
                                             edgecolor='k', alpha=0., lw=6, facecolor=None, hatch = '//', zorder=-1)
                                     },
                        'WINDOW_4': {'LOCATION': 6,
                                     'X': [deltas['SNR'], 'SNR'],
                                     'Y': [deltas['TOTAL_FLUX'], r'$\Delta$ Tot Flux  (Jy/beam km/s)']},
                        'WINDOW_5': {'LOCATION': 7,
                                     'X': [deltas['SNR'], 'SNR'],
                                     'Y': [deltas['VROT'], r'$\Delta$ V$_{\rm rot}$  (channels)']},
                        'WINDOW_6': {'LOCATION': 8,
                                     'X': [deltas['SNR'], 'SNR'],
                                     'Y': [deltas['R_HI'],r'$\Delta$ R$_{\rm HI}$  (beams)']}}
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
        if LVHIS:
            coloring = [0 if x == 'ROTCUR' else 90. for x in deltas['INPUT_MODEL']]

        else:
            coloring = deltas['CENTRAL_INPUT_INCLINATION']
            colorbarlabel = 'Inclination'
        for i,key in enumerate(plot_assembly[plot]):
            if 'NO_MEAN' in plot_assembly[plot][key]:
                nomean = plot_assembly[plot][key]['NO_MEAN']
            else:
                nomean = False
            ax,legend_items = make_plot(plot_assembly[plot][key]['X'][0],\
                           plot_assembly[plot][key]['Y'][0],\
                           xlabel = plot_assembly[plot][key]['X'][1],\
                           ylabel = plot_assembly[plot][key]['Y'][1],\
                           location = gs[plot_assembly[plot][key]['LOCATION']],\
                           color=coloring,status=deltas['STATUS'],\
                           symbol=deltas['RCSHAPE'],No_Mean = nomean,corruption=deltas['CORRUPTION'])
            if 'PATCH' in plot_assembly[plot][key]:
                ax.add_patch( plot_assembly[plot][key]['PATCH'])
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
        a = np.array([[0,90.]])
        img = plt.imshow(a, cmap="rainbow")
        plt.gca().set_visible(False)
        cax = plt.axes([0.63, 0.4, 0.01, 0.45])
        barr = plt.colorbar(orientation="vertical", cax=cax)
        if LVHIS:
            barr.set_ticks([0, 90.])
            barr.ax.set_yticklabels(['ROTCUR', 'DISKFIT'],rotation=90,va = 'center')

        else:
            barr.set_label('Inclination', rotation=270, verticalalignment='bottom')


        labelfont= {'family':'Times New Roman',
                    'weight':'normal',
                    'size':37}
        plt.rc('font',**labelfont)
        galaxies = 0
        succes_galaxies=0
        for i,x in enumerate(deltas['INPUT_MODEL']):
            if x==deltas['INPUT_MODEL'][0]:
                galaxies+=1
                if deltas['STATUS'][i] > 0.1:
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



def make_plot(x_in,y_in, color= None, status= None, location = [0,1], symbol= None,
                    xlabel = '',ylabel = '', No_Mean = False, corruption = None):
        try:
            x = np.array([v[0] for v in x_in], dtype=float)
            x_err = np.array([v[1] for v in x_in], dtype=float)
        except TypeError:
            x = np.array([v for v in x_in], dtype=float)
            x_err = np.array([0. for v in x_in], dtype=float)
        except IndexError:
            x = np.array([v for v in x_in], dtype=float)
            x_err = np.array([0. for v in x_in], dtype=float)

        try:
            y = np.array([v[0] for v in y_in], dtype=float)
            y_err = np.array([v[1] for v in y_in], dtype=float)
        except TypeError:
            y = np.array([v for v in y_in], dtype=float)
            y_err = np.array([0. for v in y_in], dtype=float)
        except IndexError:
            y = np.array([v for v in y_in], dtype=float)
            y_err = np.array([0. for v in y_in], dtype=float)

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
            transparency = np.array(copy.deepcopy(status),dtype=float)
            for i,els in enumerate(stat_elements):
                transparency[status == els] = norm_elements[i]
            transparency = np.array(transparency,dtype = float)

        else:
            stat_elements = np.array([0.])
            norm_elements = [1]
            mean = np.nanmean(y[:])
            stdev = np.nanstd(y[:]-mean)
            transparency =np.ones(len(x[:]))

        if not corruption is None:
            corruption_types = np.unique(corruption)
            numerical_corruption= [x for x in range(len(corruption_types))]
            corruption_size = copy.deepcopy(corruption)
            for i,c_type in enumerate(corruption_types):
                corruption_size = [numerical_corruption[i] if x == c_type else x for x in corruption_size]
            corruption_size = np.array(corruption_size,dtype=float)
            corruption = np.array(corruption)
            size_legend_items= []
            if len(numerical_corruption) > 1:
                sizes = (np.array(numerical_corruption,dtype=float)+4)**3.
                for i,siz in enumerate(sizes):
                    lab_string = f'Corrupted with {corruption_types[i]}'
                    tmp_fig = plt.figure(1,figsize=(1,1),dpi=30,facecolor = 'w', edgecolor = 'k')
                    tmp_plot = plt.scatter([0,1],[0,1], c = 'k', s=siz, marker = 'o',label = lab_string)
                    size_legend_items.append(tmp_plot)
                    plt.close(tmp_fig)

        else:
            corruption_size = np.zeros(len(x[:]))

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
        shape_legend_items = []
        proc_lab_string = []
        for i,shaped in enumerate(req_no_elements):
            proc = np.where(symbol == shaped)[0]

            for j,transparency_val in enumerate(norm_elements):

                add = np.where(transparency[proc] == transparency_val)[0]

                if len(add) > 0:
                    lab_string = f'RC Shape {shaped}'
                    if lab_string not in proc_lab_string:
                        tmp_fig = plt.figure(1,figsize=(2,2),dpi=30,facecolor = 'w', edgecolor = 'k')
                        tmp_plot = plt.scatter(x[proc[add]],y[proc[add]],cmap= 'rainbow', c = 'k', s=4**3, marker = symbol_use[i],alpha = 1,label = lab_string)
                        shape_legend_items.append(tmp_plot)
                        proc_lab_string.append(lab_string)
                        plt.close(tmp_fig)
                    siz = (corruption_size[proc[add]]+4)**3.

                    #try:
                    #ax.scatter(x[proc[add]],y[proc[add]],cmap= 'rainbow', c = rgba_cols[0][proc[add]][:], s=2**4, marker = symbol_use[i],alpha = 1,label = lab_string)
                    ax.scatter(x[proc[add]],y[proc[add]],cmap= 'rainbow', c = rgba_cols[0][proc[add]][:], s=siz, marker = symbol_use[i],alpha = norm_elements[j])

                    plt.errorbar(x[proc[add]],y[proc[add]],xerr=x_err[proc[add]],yerr=y_err[proc[add]], linestyle="None", ecolor = rgba_cols[0][proc[add]][:],alpha = norm_elements[j])
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
def get_flux_values(directory):
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
    delta_max= (model_max-fat_max)/(bmaj*3600.)
    return [delta_max, np.max([model_last_ring,fat_last_ring])/(float(bmaj)*3600.)]

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
    print(file_extension)
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

    else:
        cfg = OmegaConf.structured(defaults)

        # read command line arguments anything list input should be set in '' e.g. pyROTMOD 'rotmass.MD=[1.4,True,True]'
        yaml_config = OmegaConf.load(filename)
#merge yml file with defaults
        cfg = OmegaConf.merge(cfg,yaml_config)
        # translate into our dictionary
        Configuration = setup_configuration(cfg)

    while not os.path.exists(Configuration['INPUT_CATALOGUE']):
        Configuration['INPUT_CATALOGUE'] = input('''
                    Your input catalogue ({}) does not exist.
                    Please provide the correct file name.
                    '''.format(Configuration['INPUT_CATALOGUE']))
    #The output catalogue only needs to be in a valid directory as we create it
    output_catalogue_dir = Configuration['OUTPUT_CATALOGUE'].split('/')
    if len(output_catalogue_dir) > 1:
        check_dir = '/'.join(output_catalogue_dir[:-1])
        while not os.path.isdir(check_dir):
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
def load_tirific(filename,Variables = ['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT','VROT_ERR',
                 'Z0', 'Z0_ERR', 'SBR','SBR_ERR', 'INCL','INCL_ERR','PA','PA_ERR','XPOS','XPOS_ERR','YPOS','YPOS_ERR','VSYS','VSYS_ERR','SDIS','SDIS_ERR'
                 ,'VROT_2','VROT_2_ERR',  'Z0_2','Z0_2_ERR','SBR_2','SBR_2_ERR',
                 'INCL_2', 'INCL_2_ERR','PA_2','PA_2_ERR','XPOS_2','XPOS_2_ERR','YPOS_2','YPOS_2_ERR','VSYS_2','VSYS_2_ERR','SDIS_2','SDIS_2_ERR','CONDISP','CFLUX','CFLUX_2'],
                LVHIS=False):
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
                print(line)
                output[var_concerned] = [float(x) for x in line.split('=')[1].rsplit()]
    for input in Variables:
        if input not in output:
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



def retrieve_deltas_and_RCs(database_config, database_inp_catalogue, database_out_catalogue,binary=False,LVHIS= False):
    '''for every galaxy read the Finalmodel and retrieve RCs and deltas'''
    if LVHIS:
        LVHIS_Names = load_LVHIS_Name(database_config["MAIN_DIRECTORY"])
    RCs = {}
    #deltas = {'NAME': [], 'BEAMS_ACROSS': [], 'SNR':[], 'MAX_EXTEND': [], 'R_HI':[]}
    deltas= {}
    for i, galaxy in enumerate(database_out_catalogue['DIRECTORY_NAME']):
        print(f'Processing directory {galaxy}')

        status = 0
        if binary:
            if int(database_out_catalogue['AC1'][i]) ==1:
                status =1
            if int(database_out_catalogue['AC2'][i]) ==1:
                status =2
        else:
            if database_out_catalogue['OS'][i]:
                status = 2

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
            total_flux,error_flux,average_flux = get_flux_values(f'{database_config["MAIN_DIRECTORY"]}/{galaxy}/')

        output_parameters = remove_too_faint(output_parameters)
        if 'DISTANCE' in model_parameters:
            distance=model_parameters['DISTANCE']
        else:
            distance=[model_parameters['VSYS'][0]/69.7] #km/s/Mpc
        if len(distance) == 0.:
            distance=[model_parameters['VSYS'][0]/69.7]

        #write our RCs to a variable
        # First read the the input model Class
        if LVHIS:
            print(model_parameters['RADI'][-1],model_parameters['BMAJ'][0])
            diameter_in_beams = [float(model_parameters['RADI'][-1])/float(model_parameters['BMAJ'][0])*2.]
            if 'RADI_2' in model_parameters:
                print('What happened')
                diameter_in_beams.append(float(model_parameters['RADI_2'][-1])/float(model_parameters['BMAJ'][0])*2.)
            else:
                diameter_in_beams.append(float(model_parameters['RADI'][-1])/float(model_parameters['BMAJ'][0])*2.)

            SNR = average_flux[0]/output_parameters['RMS']
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
        cubename = f'{database_inp_catalogue["CUBENAME"][database_inp_catalogue["DIRECTORYNAME"].index(galaxy)]}_preprocessed.fits'
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
            deltas['MAX_EXTEND'].append(
                get_diff_rmax(model_parameters[f'RADI{ext[j]}'],output_parameters['RADI'],hdr['BMAJ']))
            #for now leave R_HI
            deltas['TOTAL_FLUX'].append([total_flux[1]-total_flux[0], error_flux[1]+error_flux[0]])
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
                    print(model)
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

    sbr_msolar = columndensity({'OUTPUTLOG': None, 'DISTANCE': distance},np.array(sbr_profile[0],dtype=float)*1000.\
                        ,systemic=systemic,arcsquare=True,solar_mass_output=True)
    sbr_2_msolar = columndensity({'OUTPUTLOG': None, 'DISTANCE': distance},np.array(sbr_profile[1],dtype=float)*1000.\
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

    #if np.sum(radi) == 0. or np.sum(sbr_profile) == 0.:
    #    return  -1.
    #First read the model radius from the info file
    #for file in os.listdir(directory):
    #    if file.endswith("-Info.txt"):
    #        with open(f'{directory}/{file}') as Infofile:
    #            info_lines = Infofile.readlines()
    #        for line in info_lines:
    #            values = [x.strip() for x in line.split()]
    #            if values[0].lower() == 'distance':
    #                distance = float(values[1])
    #            if f'{values[0].lower()} {values[1].lower()}' == 'hi radius':
    #                kpc_radius = float(values[2])
    #model_kpc_arcsec = convertskyangle({'OUTPUTLOG': None, 'DISTANCE': distance},\
    #                    kpc_radius, physical=True)
#
#    av_sbr = np.array([np.mean([x,y]) for x,y in zip(sbr_profile[0],sbr_profile[1])],dtype=float)
#    M_solar_column=columndensity({'OUTPUTLOG': None, 'DISTANCE': distance},av_sbr,\
#                systemic=systemic,arcsquare=True,solar_mass_output=True)
#    start = np.where(M_solar_column < 1.)[0]
#    if len(start) == 0.:
#        radius =float('NaN')
#    else:
#        radius = radi[start[0]]
#    return radius-model_kpc_arcsec
