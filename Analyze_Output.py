#!/usr/local/bin/ python3


# Is to analyze the output from both FAT BBarolo on the Artificial galaxies
#FATInput='FAT_INPUT.config'
FATInput='/data/users/kamphuis/Artificial/FAT_INPUT.config'


import os
import numpy as np
import multiprocessing
import subprocess
from collections import OrderedDict
import time
import copy
import sys
import re
from astropy.io import fits
import matplotlib.pyplot as plt
from scipy import interpolate
import matplotlib
from matplotlib.ticker import AutoMinorLocator
import matplotlib.patches as patches
from matplotlib.pyplot import cm
import matplotlib.gridspec as gridspec
#A function for converting RA and DEC
def convertRADEC(RA,DEC,invert=False, colon=False):
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
            if (sign < 0.): DEC[i]='-'+DEC[i]
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
            tmp = re.split(r"[a-z,:,',\"]+",ypos[i])
            DEC[i]=float(np.absolute(float(tmp[0]))+((float(tmp[1])+(float(tmp[2])/60.))/60.))*float(tmp[0])/np.absolute(float(tmp[0]))

        if len(RA) == 1:
            RA= float(RA[0])
            DEC = float(DEC[0])             
    return RA,DEC
# Code for creating  proper dictionary instead of this python unordered nonsense.
# This also includes an insert class.
class MyOrderedDict(OrderedDict):
     def insert(self,existing_key,new_key,key_value):
        done = False 
        if new_key in self:
            self[new_key]=key_value
            done = True
        else:
            new_orderded_dict=self.__class__()
            for key, value in self.items():
                new_orderded_dict[key]=value
                if key==existing_key:
                    new_orderded_dict[new_key]=key_value
                    done = True
            if not done:
                new_orderded_dict[new_key]=key_value
                done = True
                print("----!!!!!!!! YOUR new key was appended at the end as you provided a non-existing key to add it after!!!!!!---------")
            self.clear()
            self.update(new_orderded_dict)
        if not done:
            print("----!!!!!!!!We were unable to add your key!!!!!!---------")

start_time = time.time()
# First we get the catalogue of fitted galaxies
tmp = open(FATInput,'r')
Template_in = MyOrderedDict({})
unarranged = tmp.readlines()
# Separate the keyword names
for tmp in unarranged:
    # python is really annoying with needing endlines. Let's strip them here and add them when writing
    Template_in[tmp.split('=',1)[0].strip().upper()]=tmp.rstrip()
    if  tmp.split('=')[0] == 'catalogue':
        filefound=  tmp.split('=')[1].strip()
    if  tmp.split('=')[0] == 'startgalaxy':
        start= int(tmp.split('=')[1].strip())
    if  tmp.split('=')[0] == 'endgalaxy':
        end= int(tmp.split('=')[1].strip())

#replace space
filefound = "\ ".join(filefound.split())
#print(filefound)            
if end == -1:
    proc = subprocess.Popen(["wc", filefound], stdout=subprocess.PIPE)
    (out, err) = proc.communicate()
    end = int(out.split()[0]) -1
tmp = open(filefound,'r')
unarranged = tmp.readlines()
# Get all the dir names
dirname=[]
counter = 0
for tmp in unarranged:
	if counter != 0.:
		dirname.append(tmp.split('|')[2].rstrip())
	counter+= 1
dirname=dirname[0:2]
#print(dirname)
#create lists for the various parameters
#incval=[[[],[]],[[],[]]]
incval =np.zeros((len(dirname),2,2))
#paval = [[[],[]],[[],[]]]
paval =np.zeros((len(dirname),2,2))
#vrotval = [[[],[]],[[],[]]]
vrotval =np.zeros((len(dirname),2,2)) 
#sbrval =[[[],[]],[[],[]]]
sbrval =np.zeros((len(dirname),2,2)) 
#dispval =[[[],[]],[[],[]]]
dispval =np.zeros((len(dirname),2,2))
#scaleheightval =[[[],[]],[[],[]]]
scaleheightval =np.zeros((len(dirname),2,2)) 
#vsysval =[[[],[]],[[],[]]]
vsysval =np.zeros((len(dirname),2,2))
RAval =np.zeros((len(dirname),2,2))
DECval =np.zeros((len(dirname),2,2))
#[[[],[]],[[],[]]] 
#DECval =[[[],[]],[[],[]]] 
#R_maxval =[[[],[]],[[],[]]]
R_maxval =np.zeros((len(dirname),2,2))
incarrange = np.zeros(len(dirname))
beamarrange = np.zeros(len(dirname))
#vradval = [[[],[]],[[],[]]] 
# Some conversions
conv_pc_arcsec=605.7383*1.823E18*(2.*np.pi/(np.log(256.)))*1000./1.24756e+20


counter = 0
proc=dict()
for i in range(len(dirname)):
#for i in range(2):
    #os.chdir('/data/users/kamphuis/Artificial/'+dirname[i])
    os.chdir(dirname[i])
	# First we read the Model def file
    #print(dirname[i].split('bm')[1].split('-')[0])
   
    beamarrange[i]=float( dirname[i].split('ba')[1].split('SNR')[0])
    
    modelfile = open('ModelInput.def')
    modelline = modelfile.readlines()
    for line in modelline:
		#print(line.split('=')[0].strip())
        if line.split('=')[0].strip()== 'BMAJ':
            beamsize = float(line.split('=')[1].strip())
        if line.split('=')[0].strip()== 'RADI':
            radii =np.array( line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'PA':
            pa = np.array( line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'INCL':
            inc =np.array( line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'SBR':
            sbr = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'VROT':
            vrot = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'SDIS':
            disp = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'Z0':
            scaleheight = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'XPOS':
            RA = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'YPOS':
            DEC = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'VSYS':
            vsys = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'PA_2':
            pa2 = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'INCL_2':
            inc2 = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'SBR_2':
            sbr2 = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'VROT_2':
            vrot2 = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'SDIS_2':
            disp2 = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'Z0_2':
            scaleheight2 = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'XPOS_2':
            RA2 = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'YPOS_2':
            DEC2 = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'VSYS_2':
            vsys2 = np.array(line.split('=')[1].rsplit(),dtype=float)
        #Next up we read th FAT output
    #print(radii,pa,inc,vrot,sbr,scaleheight,pa2,inc2,disp,vsys)
    
    modelfile = open('Finalmodel/Finalmodel.def')
    modelline = modelfile.readlines()
    for line in modelline:
        #print(line.split('=')[0].strip())

        if line.split('=')[0].strip()== 'RADI':
            radiifat = np.array(line.split('=')[1].rsplit(),dtype=float)		
        if line.split('=')[0].strip()== 'PA':
            pafat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'INCL':
            incfat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'SBR':
            sbrfat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'VROT':
            vrotfat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'SDIS':
            dispfat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'Z0':
            scaleheightfat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'XPOS':
            RAfat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'YPOS':
            DECfat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'VSYS':
            vsysfat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'PA_2':
            pa2fat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'INCL_2':
            inc2fat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'SBR_2':
            sbr2fat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'VROT_2':
            vrot2fat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'SDIS_2':
            disp2fat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'Z0_2':
            scaleheight2fat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'XPOS_2':
            RA2fat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'YPOS_2':
            DEC2fat = np.array(line.split('=')[1].rsplit(),dtype=float)
        if line.split('=')[0].strip()== 'VSYS_2':
            vsys2fat = np.array(line.split('=')[1].rsplit(),dtype=float)
    #print(radiifat,pafat,incfat,vrotfat,sbrfat,scaleheightfat,pa2fat,inc2fat,dispfat,vsysfat)
    # then we need to pick upt the barolo parameters. This file has a rotcur outlay
    
    tmp1,radiiBB,vrotBB,dispBB,incBB,paBB,tmp2,scaleheightBB,sbrBB,xposBB,yposBB,vsysBB,vradBB = np.loadtxt('output/NONE/ringlog2.txt' , unpack =True, skiprows =1 )
    # in order to get the central position we need to relate back to the image
    cube_in = fits.open('Convolved_Cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    crval2rad=(cube_in[0].header["CRVAL2"]/180.)*np.pi
    cdelt1=(cube_in[0].header["CDELT1"]/(np.cos(crval2rad)))
    RABB=cube_in[0].header["CRVAL1"]+(xposBB-cube_in[0].header["CRPIX1"]+1)*(cube_in[0].header["CDELT1"])
    DECBB = cube_in[0].header["CRVAL2"]+(yposBB-cube_in[0].header["CRPIX2"]+1)*(cube_in[0].header["CDELT2"])
    #print(radiiBB,vrotBB,dispBB,incBB,paBB,scaleheightBB,sbrBB,RABB,DECBB,vsysBB,vradBB) 		

# Now let's first make a figure with the different profiles as fitted
    labelfont= {'family':'Times New Roman',
                'weight':'normal',
                'size':18}
    plt.rc('font',**labelfont)    
    plt.figure(45,figsize=(8,12),dpi=300,facecolor = 'w', edgecolor = 'k')
    gs = gridspec.GridSpec(6,1 )
    gs.update(wspace=0.0, hspace=0.0) 
    ax = plt.subplot(gs[5])
    #overview = plt.subplot(6,1,6)
    plt.plot(radii,sbr,'k')

    plt.plot(radiifat,sbrfat/1.5,'b')
    plt.plot(radiifat,sbrfat,'bo')
    plt.plot(radiifat,sbr2fat/1.5,'--g')
    plt.plot(radiifat,sbr2fat,'go')
    
    if np.sum(sbrBB)/len(sbrBB) != 1:
        plt.plot(radiiBB,sbrBB,'r')	
        plt.plot(radiiBB,sbrBB,'ro')
    ax.set_xlabel(r'Radius (arcsec)',**labelfont)
    ax.set_ylabel(r'SBR (Jy km s$^{-1}$ arcsec$^{-2}$)')
    plt.text
    ax.yaxis.set_minor_locator(AutoMinorLocator(4))
    ax.xaxis.set_minor_locator(AutoMinorLocator(4))                    
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(1.5)
    plt.tick_params(axis='both', which='minor', bottom='on',left='on',length=3)
    plt.tick_params(axis='both', which='major', labelsize=17, length=6)
    plt.tick_params(axis='both', which='both', direction = 'in', width=1.5 , bottom='on',left='on' ,right ='on', top='on')
    plt.tight_layout
    plt.margins(x=0., y=0.1)
   
   
    axe = plt.subplot(gs[4])
    plt.plot(radii,vrot,'k')
    plt.ylabel(r'V$_{{\rm rot}}$ (km s$^{-1}$)')
    plt.plot(radiifat,vrotfat,'b')
    plt.plot(radiifat,vrotfat,'bo')
    plt.plot(radiiBB,vrotBB,'r')
    plt.plot(radiiBB,vrotBB,'ro')
    axe.yaxis.set_minor_locator(AutoMinorLocator(4))
    axe.xaxis.set_minor_locator(AutoMinorLocator(4))                    
    for axis in ['top','bottom','left','right']:
        axe.spines[axis].set_linewidth(1.5)
    plt.tick_params(axis='both', which='minor', bottom='on',left='on',length=3)
    plt.tick_params(axis='both', which='major', labelsize=17, length=6)
    plt.tick_params(axis='both', which='both', direction = 'in', width=1.5 , bottom='on',left='on' ,right ='on', top='on')
    plt.tight_layout
    plt.margins(x=0., y=0.1)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom='on',      # ticks along the bottom edge are off
        top='on',         # ticks along the top edge are off
        labelbottom='off') #
    axe = plt.subplot(gs[3])
    plt.plot(radii,inc,'k')
    plt.plot(radii,inc2,'--',c='k')
    plt.ylabel(r'INCL ($\degree$)')
    plt.plot(radiifat,incfat,'b')
    plt.plot(radiifat,incfat,'bo')
    plt.plot(radiifat,inc2fat,'--g')
    plt.plot(radiifat,inc2fat,'go')
    plt.plot(radiiBB,incBB,'r')
    plt.plot(radiiBB,incBB,'ro')
    axe.yaxis.set_minor_locator(AutoMinorLocator(4))
    axe.xaxis.set_minor_locator(AutoMinorLocator(4))                    
    for axis in ['top','bottom','left','right']:
        axe.spines[axis].set_linewidth(1.5)
    plt.tick_params(axis='both', which='minor', bottom='on',left='on',length=3)
    plt.tick_params(axis='both', which='major', labelsize=17, length=6)
    plt.tick_params(axis='both', which='both', direction = 'in', width=1.5 , bottom='on',left='on' ,right ='on', top='on')
    plt.tight_layout
    plt.margins(x=0., y=0.1)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom='on',      # ticks along the bottom edge are off
        top='on',         # ticks along the top edge are off
        labelbottom='off') #
    axe = plt.subplot(gs[2])
    plt.plot(radii,pa,'k')
    plt.plot(radii,pa2,'--',c='k')
    plt.ylabel(r'PA ($\degree$)')
    plt.plot(radiifat,pafat,'b')
    plt.plot(radiifat,pafat,'bo')
    plt.plot(radiifat,pa2fat,'--g')
    plt.plot(radiifat,pa2fat,'go')
    plt.plot(radiiBB,paBB,'r')
    plt.plot(radiiBB,paBB,'ro')
    axe.yaxis.set_minor_locator(AutoMinorLocator(4))
    axe.xaxis.set_minor_locator(AutoMinorLocator(4))                    
    for axis in ['top','bottom','left','right']:
        axe.spines[axis].set_linewidth(1.5)
    plt.tick_params(axis='both', which='minor', bottom='on',left='on',length=3)
    plt.tick_params(axis='both', which='major', labelsize=17, length=6)
    plt.tick_params(axis='both', which='both', direction = 'in', width=1.5 , bottom='on',left='on' ,right ='on', top='on')
    plt.tight_layout
    plt.margins(x=0., y=0.1)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom='on',      # ticks along the bottom edge are off
        top='on',         # ticks along the top edge are off
        labelbottom='off') #
    axe = plt.subplot(gs[1])
    plt.plot(radii,disp,'k')
    plt.plot(radii,disp2,'--',c='k')
    plt.ylabel(r'Disp. (km s$^{-1}$)')
    plt.plot(radiifat,dispfat,'b')
    plt.plot(radiifat,dispfat,'bo')
    plt.plot(radiiBB,dispBB,'r')
    plt.plot(radiiBB,dispBB,'ro')
    axe.yaxis.set_minor_locator(AutoMinorLocator(4))
    axe.xaxis.set_minor_locator(AutoMinorLocator(4))                    
    for axis in ['top','bottom','left','right']:
        axe.spines[axis].set_linewidth(1.5)
    plt.tick_params(axis='both', which='minor', bottom='on',left='on',length=3)
    plt.tick_params(axis='both', which='major', labelsize=17, length=6)
    plt.tick_params(axis='both', which='both', direction = 'in', width=1.5 , bottom='on',left='on' ,right ='on', top='on')
    plt.tight_layout
    plt.margins(x=0., y=0.1)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom='on',      # ticks along the bottom edge are off
        top='on',         # ticks along the top edge are off
        labelbottom='off') #

    axe = plt.subplot(gs[0])
    plt.plot(radii,scaleheight,'k')
    plt.plot(radii,scaleheight2,'--',c='k')
    plt.ylabel(r'Scale height (arcsec)')
    plt.plot(radiifat,scaleheightfat,'b')
    plt.plot(radiifat,scaleheightfat,'bo')
    plt.plot(radiiBB,scaleheightBB,'r')
    plt.plot(radiiBB,scaleheightBB,'ro')
    axe.yaxis.set_minor_locator(AutoMinorLocator(4))
    axe.xaxis.set_minor_locator(AutoMinorLocator(4))                    
    for axis in ['top','bottom','left','right']:
        axe.spines[axis].set_linewidth(1.5)
    plt.tick_params(axis='both', which='minor', bottom='on',left='on',length=3)
    plt.tick_params(axis='both', which='major', labelsize=17, length=6)
    plt.tick_params(axis='both', which='both', direction = 'in', width=1.5 , bottom='on',left='on' ,right ='on', top='on')
    plt.tight_layout
    plt.margins(x=0., y=0.1)
    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='both',      # both major and minor ticks are affected
        direction = 'in',
        bottom='on',      # ticks along the bottom edge are off
        top='on',         # ticks along the top edge are off
        labelbottom='off') #
    # we also need to compare the central positions
   
    RAh,DECh = convertRADEC(RA[0],DEC[0])
    labelfont= {'family':'Times New Roman',
                'weight':'normal',
                'size':14}
    plt.text(0.0,1.5,'Central Positions',transform = axe.transAxes,**labelfont)
    plt.text(0.3,1.5,'Model RA={} Dec = {} Vsys = {}'.format(RAh,DECh,vsys[0]),transform = axe.transAxes, **labelfont)
    RAhfat,DEChfat = convertRADEC(RAfat[0],DECfat[0])
    plt.text(0.3,1.4,'FAT RA={} Dec = {} Vsys = {}'.format(RAhfat,DEChfat,vsysfat[0]),transform = axe.transAxes,color= 'b', **labelfont)
    # for barolo this is trickier as we need to convert from the pixel values
    RAhBB,DEChBB = convertRADEC(RABB[0],DECBB[0])
    plt.text(0.3,1.3,'BBarolo RA={} Dec = {} Vsys = {}'.format(RAhBB,DEChBB,vsysBB[0]),transform = axe.transAxes,color= 'r', **labelfont)
    plt.savefig('All_Fits_Overview.png',bbox_inches='tight',**labelfont)
    plt.close()


    # Then we need to make a single point for each parameter
    # First the sbr at the max radii
    # We will say that the maximum extend of the disc is where it drops below 5% of the maximum
    rmaxmod= radii[np.where(sbr > 1./conv_pc_arcsec)[0][-1]+1]*1.25
    #print(rmaxmod,np.where(sbr > np.max(sbr)/20.)[0][-1])
    sbrmod =  interpolate.interpolate.interp1d(radii,sbr*conv_pc_arcsec,fill_value="extrapolate")
    stepsize=int(beamsize/(radii[1]))
    minrad=radiiBB[0]
    min_index=np.where(radii < minrad)[0][0]
    #if len(R_maxval) == 0:
   
    R_maxval[i][0][:] = [(rmaxmod-radiifat[-1])/(cube_in[0].header['BMAJ']*3600.),float(sbrmod(radiifat[-1]))]
    R_maxval[i][1][:] = [(rmaxmod-radiiBB[-1])/(cube_in[0].header['BMAJ']*3600.),float(sbrmod(radiiBB[-1]))]

    # next up the difference in integratefdflux
    # which we get from the produced cubes not the fitted profiles
    # First calculate the flux within the mask
    mask =  fits.open('mask.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    totalflux = np.sum(cube_in[0].data[mask[0].data > 0.5])*cube_in[0].header["CDELT3"]/1000.
    # for fat and barolo we take the total of the moment 0
    mom0_fat= fits.open('Moments/Finalmodel_mom0.fits')
    totalflux_fat = np.sum(mom0_fat[0].data)
    mom0_BB= fits.open('output/NONE/NONE_mom0th.fits')
    beamarea=(np.pi*abs(mom0_BB[0].header['BMAJ']*3600.*mom0_BB[0].header['BMIN']*3600))/(4.*np.log(2.))
    # Convert arcsec to pixels            
    pixperbeam=beamarea/(abs(mom0_BB[0].header['CDELT1']*3600.)*abs(mom0_BB[0].header['CDELT2']*3600.))
    totalflux_BB = np.sum(mom0_BB[0].data)*pixperbeam


    
    sbrval[i][0][:] = [totalflux-totalflux_fat,totalflux/100.]
    sbrval[i][1][:] = [totalflux-totalflux_BB,totalflux/100.]
    
    # next up evaluate the rotation curve we do this from 1 beam out to the maximum radius of the shortest fit
    evalrad=np.arange(np.max([radiifat[0],radiiBB[0]]),np.min([radiifat[-1],radiiBB[-1]]),cube_in[0].header["BMAJ"]*3600.)
        
    #print(R_maxval)
    vrotint = interpolate.interpolate.interp1d(radii,vrot,fill_value="extrapolate")
    vrotfatint = interpolate.interpolate.interp1d(radiifat,vrotfat,fill_value="extrapolate")
    vrotBBint = interpolate.interpolate.interp1d(radiiBB,vrotBB,fill_value="extrapolate")
    # the error is the total difference with an error os std
   
    vrotval[i][0][:] = [np.mean(abs(vrotint(evalrad)-vrotfatint(evalrad))),np.std((vrotint(evalrad)-vrotfatint(evalrad)))]
    vrotval[i][1][:] = [np.mean(abs(vrotint(evalrad)-vrotBBint(evalrad))),np.std((vrotint(evalrad)-vrotBBint(evalrad)))]
   
   # then possible the inclination
    incint = interpolate.interpolate.interp1d(radii,inc,fill_value="extrapolate")
    incfatint = interpolate.interpolate.interp1d(radiifat,incfat,fill_value="extrapolate")
    inc2int = interpolate.interpolate.interp1d(radii,inc2,fill_value="extrapolate")
    inc2fatint = interpolate.interpolate.interp1d(radiifat,inc2fat,fill_value="extrapolate")
    incBBint = interpolate.interpolate.interp1d(radiiBB,incBB,fill_value="extrapolate")
    # the error is the total difference with an error os std
    incval[i][0][:] = [np.mean([abs(incint(evalrad)-incfatint(evalrad)),abs(inc2int(evalrad)-inc2fatint(evalrad))]),np.std([(incint(evalrad)-incfatint(evalrad)),(inc2int(evalrad)-inc2fatint(evalrad))])]
    incval[i][1][:] = [np.mean([abs(incint(evalrad)-incBBint(evalrad)),abs(inc2int(evalrad)-incBBint(evalrad))]),np.std([(incint(evalrad)-incBBint(evalrad)),(inc2int(evalrad)-incBBint(evalrad))])]
    #print(incval)
    incarrange[i] = inc[0]    
    # then possible the palination
    paint = interpolate.interpolate.interp1d(radii,pa,fill_value="extrapolate")
    pafatint = interpolate.interpolate.interp1d(radiifat,pafat,fill_value="extrapolate")
    pa2int = interpolate.interpolate.interp1d(radii,pa2,fill_value="extrapolate")
    pa2fatint = interpolate.interpolate.interp1d(radiifat,pa2fat,fill_value="extrapolate")
    paBBint = interpolate.interpolate.interp1d(radiiBB,paBB,fill_value="extrapolate")
    # the error is the total difference with an error os std
    paval[i][0][:] =[np.mean([abs(paint(evalrad)-pafatint(evalrad)),abs(pa2int(evalrad)-pa2fatint(evalrad))]),np.std([(paint(evalrad)-pafatint(evalrad)),(pa2int(evalrad)-pa2fatint(evalrad))])]
    paval[i][1][:] =[np.mean([abs(paint(evalrad)-paBBint(evalrad)),abs(pa2int(evalrad)-paBBint(evalrad))]),np.std([(paint(evalrad)-paBBint(evalrad)),(pa2int(evalrad)-paBBint(evalrad))])]
    #print(paval)
 

    # then possible the dispersion
    dispint = interpolate.interpolate.interp1d(radii,disp,fill_value="extrapolate")
    dispfatint = interpolate.interpolate.interp1d(radiifat,dispfat,fill_value="extrapolate")
    dispBBint = interpolate.interpolate.interp1d(radiiBB,dispBB,fill_value="extrapolate")
    # the error is the total difference with an error os std
    dispval[i][0][:] = [np.mean(abs(dispint(evalrad)-dispfatint(evalrad))),np.std((dispint(evalrad)-dispfatint(evalrad)))]
    dispval[i][1][:] = [np.mean(abs(dispint(evalrad)-dispBBint(evalrad))),np.std((dispint(evalrad)-dispBBint(evalrad)))]
    # then possible the scaleheight
    scaleheightint = interpolate.interpolate.interp1d(radii,scaleheight,fill_value="extrapolate")
    scaleheightfatint = interpolate.interpolate.interp1d(radiifat,scaleheightfat,fill_value="extrapolate")
    scaleheightBBint = interpolate.interpolate.interp1d(radiiBB,scaleheightBB,fill_value="extrapolate")
    # the error is the total difference with an error os std
    scaleheightval[i][0][:] = [np.mean(abs(scaleheightint(evalrad)-scaleheightfatint(evalrad))),np.std((scaleheightint(evalrad)-scaleheightfatint(evalrad)))]
    scaleheightval[i][1][:] = [np.mean(abs(scaleheightint(evalrad)-scaleheightBBint(evalrad))),np.std((scaleheightint(evalrad)-scaleheightBBint(evalrad)))]


    RAval[i][0][:]=[(RA[0]-RAfat[0])/cube_in[0].header["BMAJ"],0.05]
    RAval[i][1][:]=[(RA[0]-RABB[0])/cube_in[0].header["BMAJ"],0.05]
    DECval[i][0][:]=[(DEC[0]-DECfat[0])/cube_in[0].header["BMAJ"],0.05]
    DECval[i][1][:]=[(DEC[0]-DECBB[0])/cube_in[0].header["BMAJ"],0.05]
    vsysval[i][0][:]= [(vsys[0]-vsysfat[0]),0.1]
    vsysval[i][1][:]= [(vsys[0]-vsysBB[0]),0.1]
    print(cube_in[0].header["CDELT3"])
        
    os.chdir("../")

labelfont= {'family':'Times New Roman',
            'weight':'normal',
            'size':18}
plt.rc('font',**labelfont)    
plt.figure(89,figsize=(24,24),dpi=300,facecolor = 'w', edgecolor = 'k')
gs = gridspec.GridSpec(3,3 )
gs.update(wspace=0.2, hspace=0.2)
#First the RA and DEC 
ax = plt.subplot(gs[0])
ax.plot([RAval[i][0][0] for i in range(len(dirname))],[DECval[i][0][0] for i in range(len(dirname))],'bo')
ax.plot([RAval[i][1][0] for i in range(len(dirname))],[DECval[i][1][0] for i in range(len(dirname))],'ro')
plt.errorbar([RAval[i][0][0] for i in range(len(dirname))],[DECval[i][0][0] for i in range(len(dirname))],xerr=[RAval[i][0][1] for i in range(len(dirname))],yerr=[DECval[i][0][1] for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([RAval[i][1][0] for i in range(len(dirname))],[DECval[i][1][0] for i in range(len(dirname))],xerr=[RAval[i][1][1] for i in range(len(dirname))],yerr=[DECval[i][1][1] for i in range(len(dirname))], linestyle="None",c='r')
ax.set_xlabel('$\Delta$ RA (beams)')
ax.set_ylabel('$\Delta$ DEC (beams)')
# Then the vsys vs. sqrt(RA^2+DEC^2)
ax = plt.subplot(gs[1])
ax.plot([vsysval[i][0][0] for i in range(len(dirname))],[np.sqrt(DECval[i][0][0]**2+RAval[i][0][0]**2)  for i in range(len(dirname))],'bo')
ax.plot([vsysval[i][1][0] for i in range(len(dirname))],[np.sqrt(DECval[i][1][0]**2+RAval[i][1][0]**2)  for i in range(len(dirname))],'ro')
plt.errorbar([vsysval[i][0][0] for i in range(len(dirname))],[np.sqrt(DECval[i][0][0]**2+RAval[i][0][0]**2)  for i in range(len(dirname))],xerr=[vsysval[i][0][1] for i in range(len(dirname))],yerr=[np.sqrt(DECval[i][0][1]**2+RAval[i][0][1]**2)  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([vsysval[i][1][0] for i in range(len(dirname))],[np.sqrt(DECval[i][1][0]**2+RAval[i][1][0]**2)  for i in range(len(dirname))],xerr=[vsysval[i][1][1] for i in range(len(dirname))],yerr=[np.sqrt(DECval[i][1][1]**2+RAval[i][1][1]**2)  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_xlabel(r'$\Delta$ V$_{{\rm sys}}$ (km s$^{-1}$)')
ax.set_ylabel(r'$\Delta$ Central (beams)')
# inclination vs delta inclination
ax = plt.subplot(gs[2])
ax.plot([incarrange[i] for i in range(len(dirname))],[incval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([incarrange[i] for i in range(len(dirname))],[incval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[incval[i][0][0]  for i in range(len(dirname))],yerr=[incval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[incval[i][1][0]  for i in range(len(dirname))],yerr=[incval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ i ($^{\circ}$)')
ax.set_xlabel(r'i ($^{\circ}$)')
# inclination vs delta PA
ax = plt.subplot(gs[3])
ax.plot([incarrange[i] for i in range(len(dirname))],[paval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([incarrange[i] for i in range(len(dirname))],[paval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[paval[i][0][0]  for i in range(len(dirname))],yerr=[paval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[paval[i][1][0]  for i in range(len(dirname))],yerr=[paval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ PA ($^{\circ}$)')
ax.set_xlabel(r'i ($^{\circ}$)')
# inclination vs delta vrot
ax = plt.subplot(gs[4])
ax.plot([incarrange[i] for i in range(len(dirname))],[vrotval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([incarrange[i] for i in range(len(dirname))],[vrotval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[vrotval[i][0][0]  for i in range(len(dirname))],yerr=[vrotval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[vrotval[i][1][0]  for i in range(len(dirname))],yerr=[vrotval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ V$_{\rm rot}$  (km s$^{-1}$)')
ax.set_xlabel(r'i ($^{\circ}$)')
# inclination vs sbr
ax = plt.subplot(gs[5])
ax.plot([incarrange[i] for i in range(len(dirname))],[sbrval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([incarrange[i] for i in range(len(dirname))],[sbrval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[sbrval[i][0][0]  for i in range(len(dirname))],yerr=[sbrval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[sbrval[i][1][0]  for i in range(len(dirname))],yerr=[sbrval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ Tot Flux  (Jy/beam km/s)')
ax.set_xlabel(r'i ($^{\circ}$)')
# inclination vs dispersion
ax = plt.subplot(gs[6])
ax.plot([incarrange[i] for i in range(len(dirname))],[dispval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([incarrange[i] for i in range(len(dirname))],[dispval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[dispval[i][0][0]  for i in range(len(dirname))],yerr=[dispval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[dispval[i][1][0]  for i in range(len(dirname))],yerr=[dispval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ Dispersion  (km/s)')
ax.set_xlabel(r'i ($^{\circ}$)')

# inclination vs scaleheight
ax = plt.subplot(gs[7])
ax.plot([incarrange[i] for i in range(len(dirname))],[scaleheightval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([incarrange[i] for i in range(len(dirname))],[scaleheightval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[scaleheightval[i][0][0]  for i in range(len(dirname))],yerr=[scaleheightval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[scaleheightval[i][1][0]  for i in range(len(dirname))],yerr=[scaleheightval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ Scaleheight  (arcsec)')
ax.set_xlabel(r'i ($^{\circ}$)')
# inclination vs R_max
ax = plt.subplot(gs[8])
ax.plot([incarrange[i] for i in range(len(dirname))],[R_maxval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([incarrange[i] for i in range(len(dirname))],[R_maxval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[R_maxval[i][0][0]  for i in range(len(dirname))],yerr=[R_maxval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([incarrange[i]  for i in range(len(dirname))],[R_maxval[i][1][0]  for i in range(len(dirname))],yerr=[R_maxval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ R$_{\rm max}$  (beams)')
ax.set_xlabel(r'i ($^{\circ}$)')
#ax.set_ylim(ymin=-0.1, ymax=0.1)
#ax.set_xlim(xmin=-0.1, xmax=0.1)
print(beamarrange)
print([RAval[i][0][0] for i in range(len(dirname))])
print([DECval[i][0][0] for i in range(len(dirname))])
#print([b[0] for item in RAval for b in item[0]])
#RAval = np.array(RAval)
print(RAval.shape)

plt.savefig('All_Differences_Combined_Inclination.png', bbox_inches='tight')
plt.close()
# Then against beam across major
plt.figure(89,figsize=(24,24),dpi=300,facecolor = 'w', edgecolor = 'k')
gs = gridspec.GridSpec(3,3 )
gs.update(wspace=0.2, hspace=0.2)
#First the RA and DEC 
ax = plt.subplot(gs[0])
ax.plot([beamarrange[i] for i in range(len(dirname))],[vsysval[i][0][0] for i in range(len(dirname))],'bo')
ax.plot([beamarrange[i] for i in range(len(dirname))],[vsysval[i][1][0] for i in range(len(dirname))],'ro')
plt.errorbar([beamarrange[i] for i in range(len(dirname))],[vsysval[i][0][0] for i in range(len(dirname))],yerr=[vsysval[i][0][1] for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([beamarrange[i] for i in range(len(dirname))],[vsysval[i][1][0] for i in range(len(dirname))],yerr=[vsysval[i][1][1] for i in range(len(dirname))], linestyle="None",c='r')
ax.set_xlabel('beams across')
ax.set_ylabel(r'$\Delta$ V$_{\rm sys}$ (km s$^{-1}$) ')
# Then the vsys vs. sqrt(RA^2+DEC^2)
ax = plt.subplot(gs[1])
ax.plot([beamarrange[i] for i in range(len(dirname))],[np.sqrt(DECval[i][0][0]**2+RAval[i][0][0]**2)  for i in range(len(dirname))],'bo')
ax.plot([beamarrange[i] for i in range(len(dirname))],[np.sqrt(DECval[i][1][0]**2+RAval[i][1][0]**2)  for i in range(len(dirname))],'ro')
plt.errorbar([beamarrange[i] for i in range(len(dirname))],[np.sqrt(DECval[i][0][0]**2+RAval[i][0][0]**2)  for i in range(len(dirname))],yerr=[np.sqrt(DECval[i][0][1]**2+RAval[i][0][1]**2)  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([beamarrange[i] for i in range(len(dirname))],[np.sqrt(DECval[i][1][0]**2+RAval[i][1][0]**2)  for i in range(len(dirname))],yerr=[np.sqrt(DECval[i][1][1]**2+RAval[i][1][1]**2)  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_xlabel(r'beams across')
ax.set_ylabel(r'$\Delta$ Central (beams)')
# inclination vs delta inclination
ax = plt.subplot(gs[2])
ax.plot([beamarrange[i] for i in range(len(dirname))],[incval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([beamarrange[i] for i in range(len(dirname))],[incval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[incval[i][0][0]  for i in range(len(dirname))],yerr=[incval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[incval[i][1][0]  for i in range(len(dirname))],yerr=[incval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ beams across')
ax.set_xlabel(r'beams across')
# inclination vs delta PA
ax = plt.subplot(gs[3])
ax.plot([beamarrange[i] for i in range(len(dirname))],[paval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([beamarrange[i] for i in range(len(dirname))],[paval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[paval[i][0][0]  for i in range(len(dirname))],yerr=[paval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[paval[i][1][0]  for i in range(len(dirname))],yerr=[paval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ PA ($^{\circ}$)')
ax.set_xlabel(r'beams across')
# inclination vs delta vrot
ax = plt.subplot(gs[4])
ax.plot([beamarrange[i] for i in range(len(dirname))],[vrotval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([beamarrange[i] for i in range(len(dirname))],[vrotval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[vrotval[i][0][0]  for i in range(len(dirname))],yerr=[vrotval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[vrotval[i][1][0]  for i in range(len(dirname))],yerr=[vrotval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ V$_{\rm rot}$  (km s$^{-1}$)')
ax.set_xlabel(r'beams across')
# inclination vs sbr
ax = plt.subplot(gs[5])
ax.plot([beamarrange[i] for i in range(len(dirname))],[sbrval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([beamarrange[i] for i in range(len(dirname))],[sbrval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[sbrval[i][0][0]  for i in range(len(dirname))],yerr=[sbrval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[sbrval[i][1][0]  for i in range(len(dirname))],yerr=[sbrval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ Tot Flux  (Jy/beam km/s)')
ax.set_xlabel(r'beams across')
# inclination vs dispersion
ax = plt.subplot(gs[6])
ax.plot([beamarrange[i] for i in range(len(dirname))],[dispval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([beamarrange[i] for i in range(len(dirname))],[dispval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[dispval[i][0][0]  for i in range(len(dirname))],yerr=[dispval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[dispval[i][1][0]  for i in range(len(dirname))],yerr=[dispval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ Dispersion  (km/s)')
ax.set_xlabel(r'beams across')

# inclination vs scaleheight
ax = plt.subplot(gs[7])
ax.plot([beamarrange[i] for i in range(len(dirname))],[scaleheightval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([beamarrange[i] for i in range(len(dirname))],[scaleheightval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[scaleheightval[i][0][0]  for i in range(len(dirname))],yerr=[scaleheightval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[scaleheightval[i][1][0]  for i in range(len(dirname))],yerr=[scaleheightval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ Scaleheight  (arcsec)')
ax.set_xlabel(r'beams across')
# inclination vs R_max
ax = plt.subplot(gs[8])
ax.plot([beamarrange[i] for i in range(len(dirname))],[R_maxval[i][0][0]  for i in range(len(dirname))],'bo')
ax.plot([beamarrange[i] for i in range(len(dirname))],[R_maxval[i][1][0]  for i in range(len(dirname))],'ro')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[R_maxval[i][0][0]  for i in range(len(dirname))],yerr=[R_maxval[i][0][1]  for i in range(len(dirname))], linestyle="None",c='b')
plt.errorbar([beamarrange[i]  for i in range(len(dirname))],[R_maxval[i][1][0]  for i in range(len(dirname))],yerr=[R_maxval[i][1][1]  for i in range(len(dirname))], linestyle="None",c='r')
ax.set_ylabel(r'$\Delta$ R$_{\rm max}$  (beams)')
ax.set_xlabel(r'beams across')
#ax.set_ylim(ymin=-0.1, ymax=0.1)
#ax.set_xlim(xmin=-0.1, xmax=0.1)
print(beamarrange)
print([RAval[i][0][0] for i in range(len(dirname))])
print([DECval[i][0][0] for i in range(len(dirname))])
#print([b[0] for item in RAval for b in item[0]])
#RAval = np.array(RAval)
print(RAval.shape)

plt.savefig('All_Differences_Combined_Beams.png', bbox_inches='tight')
plt.close()

#It took 2454.8806738853455 seconds to run this    
