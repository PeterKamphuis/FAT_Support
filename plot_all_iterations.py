# -*- coding: future_fstrings -*-

from omegaconf import OmegaConf
import omegaconf
import os
import sys

from dataclasses import dataclass, field
from datetime import datetime
from omegaconf import MISSING
#from multiprocessing import cpu_count
from psutil import cpu_count
from typing import List, Optional
import support_functions as sf
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    from matplotlib.patches import Ellipse

@dataclass
class defaults:
    parameter: str = 'INCL'
    directory: str = './'
    stop_iteration: int = -1

def obtain_iteration(name):
    number = int(name.split('_')[-1].split('.')[0])
    return number

def iteration_sort(iterations):
    iterations = sorted(iterations , key = obtain_iteration)
    return iterations
def plot_parameter(parameter,color,parameter_name,label):
    plt.plot(parameter['RADI'],parameter[parameter_name],'o', ms = 10.,
        c=color,label = label)
    plt.plot(parameter['RADI'],parameter[parameter_name],'-',lw=3,c=color)
    plt.plot(parameter['RADI'],parameter[f'{parameter_name}_2'],'o', ms = 10.,c=color)
    plt.plot(parameter['RADI'],parameter[f'{parameter_name}_2'],'--',lw=3,c=color)

def main(argv):
    cfg = OmegaConf.structured(defaults)
    # read command line arguments anything list input should be set in '' e.g. pyROTMOD 'rotmass.MD=[1.4,True,True]'
    inputconf = OmegaConf.from_cli(argv)
    cfg_input = OmegaConf.merge(cfg,inputconf)
    iterations = []
    for file in os.listdir(f'{cfg_input.directory}/Fit_Tirific_OSC'):
        if '_Iteration_' in file:
            iterations.append(file)
    if cfg_input.stop_iteration == -1:
        cfg_input.stop_iteration = len(iterations)

    iterations = iteration_sort(iterations)
    iterations = iterations[:cfg_input.stop_iteration]
    plt.figure(89,figsize=(16,16),dpi=300,facecolor = 'w', edgecolor = 'k')
    cmap = plt.get_cmap('jet')
    for file in iterations:
        print(f'Processing {file}')
        parameter = sf.load_tirific(f"{cfg_input.directory}/Fit_Tirific_OSC/{file}",\
            Variables = ['RADI',cfg_input.parameter,f'{cfg_input.parameter}_2'])
        color = cmap(float(obtain_iteration(file))/len(iterations))
        plot_parameter(parameter,color,cfg_input.parameter,f'Iteration_{obtain_iteration(file)}')


    if os.path.exists(f"{cfg_input.directory}/ModelInput.def"):
        print(f'Processing the Model')
        parameter = sf.load_tirific(f"{cfg_input.directory}/ModelInput.def",\
            Variables = ['RADI',cfg_input.parameter,f'{cfg_input.parameter}_2'])
        plot_parameter(parameter,'k',cfg_input.parameter,'Model')

    plt.legend()
    plt.savefig(f"{cfg_input.directory}/All_Iterations_Of_{cfg_input.parameter}.png", bbox_inches='tight')
    plt.close()



if __name__ == '__main__':
    #with VizTracer(output_file="FAT_Run_Viztracer.json",min_duration=1000) as tracer:
    main(sys.argv[1:])
