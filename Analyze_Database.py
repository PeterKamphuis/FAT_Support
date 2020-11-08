#!/usr/local/bin/ python3

# This script is to compare the any FAT fitting directory and produce the plots of the Current_Status report.

#The config file should be the only input required. The script assumes that the model/high resolution fit is in ModelInput.def
FATInput='/home/peter/FAT/Database/FAT_INPUT.config'


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
from scipy import interpolate
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
import support_functions as sf

start_time = time.time()
# First we get the catalogue of fitted galaxies


Template_in = sf.load_config_file(FATInput)
filefound = Template_in['CATALOGUE']
catalogue =sf.load_catalogue(filefound)
start = Template_in['STARTGALAXY']
if start < 0.:
    start = 0
if Template_in['ENDGALAXY'] == -1:
    end = len(catalogue['NUMBER'])-1
    if end == 0:
        end = 1
end = Template_in['ENDGALAXY']
dirname=catalogue['DIRECTORYNAME']

results = sf.load_output_catalogue(Template_in['OUTPUTCATALOGUE'])
fitted = np.zeros((len(dirname)))
incval =np.zeros((len(dirname),2))
#paval = [[[],[]],[[],[]]]
paval =np.zeros((len(dirname),2))
#vrotval = [[[],[]],[[],[]]]
vrotval =np.zeros((len(dirname),2))
#sbrval =[[[],[]],[[],[]]]
sbrval =np.zeros((len(dirname),2))
#dispval =[[[],[]],[[],[]]]
dispval =np.zeros((len(dirname),2))
#scaleheightval =[[[],[]],[[],[]]]
scaleheightval =np.zeros((len(dirname),2))
#vsysval =[[[],[]],[[],[]]]
vsysval =np.zeros((len(dirname),2))
RAval =np.zeros((len(dirname),2))
DECval =np.zeros((len(dirname),2))
rcshapes = ['Unspecified' for x in dirname]
RCs = {}
#[[[],[]],[[],[]]]
#DECval =[[[],[]],[[],[]]]
#R_maxval =[[[],[]],[[],[]]]
R_maxval =np.zeros((len(dirname),2))
incarrange = np.zeros(len(dirname))
beamarrange = np.zeros(len(dirname))
SNRarrange = np.zeros(len(dirname))
#vradval = [[[],[]],[[],[]]]
# Some conversions
conv_pc_arcsec=605.7383*1.823E18*(2.*np.pi/(np.log(256.)))*1000./1.24756e+20


counter = 0
proc=dict()
for i in range(len(dirname)):
#for i in range(2):
    print(f'Processing directory {dirname[i]}')
    if float(results['AC1'][results['DIRECTORY_NAME'].index(dirname[i])]) == 0.:
        diameter_in_beams, SNR, RCshape = sf.get_name_info(dirname[i])
        rcshapes[i] =RCshape
        if RCshape not in RCs:
            RCs[RCshape] = {dirname[i]: {'RADIUS':[], 'RC':[], 'DISTANCE': [], 'STATUS': 0.}}
        else:
            RCs[RCshape][dirname[i]] = {'RADIUS':[], 'RC':[], 'DISTANCE': [], 'STATUS': 0.}
        beamarrange[i]= diameter_in_beams
        SNRarrange[i] = float(SNR)
        continue

    #os.chdir('/data/users/kamphuis/Artificial/'+dirname[i])
    os.chdir(f"{Template_in['MAINDIR']}{dirname[i]}")
	# First we read the Model def file
    #print(dirname[i].split('bm')[1].split('-')[0])

    minbeamsize,beamsize,beamangle,noise,distance,numrings,radii,vrot,vrot_err,scaleheight,scaleheight_err,sbr, sbr_err,\
    inc,inc_err,pa,pa_err, RA,DEC,vsys,disp,disp_err,vrot2, vrot2_err,scaleheight2,scaleheight2_err,sbr2,sbr2_err,inc2,inc2_err,pa2,pa2_err,RA2,DEC2,vsys2,\
    disp2,disp2_err,instdisp,cflux,cflux2 = sf.load_tirific('ModelInput.def')
    beamsize= beamsize[0] #This is 0. for the Models

    diameter_in_beams, SNR, RCshape = sf.get_name_info(dirname[i])

    if RCshape not in RCs:
        RCs[RCshape] = {'MODEL':{'RADIUS':radii, 'RC':vrot, 'DISTANCE': distance[0], 'STATUS': -1}}
    else:
        if 'MODEL' not in RCs[RCshape]:
            RCs[RCshape]['MODEL'] = {'RADIUS':radii, 'RC':vrot, 'DISTANCE': distance[0], 'STATUS': -1}
    beamarrange[i]= diameter_in_beams
    SNRarrange[i] = float(SNR)
    if float(results['AC2'][results['DIRECTORY_NAME'].index(dirname[i])]) == 0.:
        output_name = 'Finalmodel/Finalmodel.def'
        fitted[i] = 1
    else:
        output_name = 'Finalmodel/Finalmodel.def'
        fitted[i] = 2
    try:
        minbeamsizefat,beamsizefat,beamanglefat,noisefat,distancefat,numringsfat,radiifat,vrotfat,vrotfat_err,scaleheightfat,scaleheightfat_err,\
        sbrfat, sbrfat_err, incfat, incfat_err,pafat,pafat_err, RAfat,DECfat,vsysfat,dispfat,dispfat_err,vrot2fat,vrot2fat_err,scaleheight2fat \
        ,scaleheight2fat_err,sbr2fat, sbr2fat_err,inc2fat,inc2fat_err,pa2fat, pa2fat_err, \
        RA2fat,DEC2fat,vsys2fat,disp2fat,disp2fat_err,instdispfat,cfluxfat,cflux2fat = sf.load_tirific(output_name)

        illegit =np.where(sbrfat < 1e-8)[0]
        if len(illegit) > 0:
            sbrfat[illegit] = 0.
            vrotfat[illegit] = 0.
            vrotfat_err[illegit] = 0.
            scaleheightfat[illegit] = 0.
            scaleheightfat_err[illegit] = 0.
            sbrfat[illegit] = 0.
            sbrfat_err[illegit] = 0.
            incfat[illegit] = 0.
            incfat_err[illegit] = 0.
            pafat[illegit] = 0.
            pafat_err[illegit] = 0.
            RAfat[illegit] = 0.
            DECfat[illegit] = 0.
            vsysfat[illegit] = 0.
            dispfat[illegit] = 0.
            dispfat_err[illegit] = 0.
        illegit2 =np.where(sbr2fat < 1e-8)[0]
        if len(illegit2) > 0:
            sbr2fat[illegit2] = 0.
            vrot2fat[illegit2] = 0.
            vrot2fat_err[illegit2] = 0.
            scaleheight2fat[illegit2] = 0.
            scaleheight2fat_err[illegit2] = 0.
            sbr2fat[illegit2] = 0.
            sbr2fat_err[illegit2] = 0.
            inc2fat[illegit2] = 0.
            inc2fat_err[illegit2] = 0.
            pa2fat[illegit2] = 0.
            pa2fat_err[illegit2] = 0.
            RA2fat[illegit2] = 0.
            DEC2fat[illegit2] = 0.
            vsys2fat[illegit2] = 0.
            disp2fat[illegit2] = 0.
            disp2fat_err[illegit2] = 0.

            # temporary sol for mix up in

    except:
        fitted[i]= 0
        RCs[RCshape][dirname[i]] = {'RADIUS':[], 'RC':[], 'DISTANCE': [], 'STATUS': fitted[i]}
        beamarrange[i]= diameter_in_beams
        rcshapes[i] =RCshape
        continue

    # temporary sol for mix up in ModelInput
    if RCshape == 'NGC_3198':
        vrot2 = vrot
        vrot2_err = vrot_err
        scaleheight2 = scaleheight
        scaleheight2_err = scaleheight_err
        sbr2 = sbr
        sbr2_err = sbr_err
        inc2 = inc
        inc2_err = inc_err
        pa2 = pa
        pa2_err =pa_err
        RA2 = RA
        DEC2 = DEC
        vsys2 = vsys
        disp2 = disp
        disp2_err = disp_err


    beamsizefat= beamsizefat[0]
    #RCs[RCshape] = {dirname[i]: [radii,vrot], 'DISTANCE': distance, 'STATUS': -1 }
    RCs[RCshape][dirname[i]] = {'RADIUS':radiifat, 'RC':vrotfat, 'DISTANCE': distancefat[0], 'STATUS': fitted[i]}
    #RCs[RCshape][dirname[i]] =  [radiifat,vrotfat,distancefat]
    rcshapes[i] =RCshape
    # in order to get the central position we need to relate back to the image
    cube_in = fits.open('Convolved_Cube.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
    crval2rad=(cube_in[0].header["CRVAL2"]/180.)*np.pi
    cdelt1=(cube_in[0].header["CDELT1"]/(np.cos(crval2rad)))



    fatmax = np.max(np.hstack((radiifat[sbrfat> 1e-8],radiifat[sbr2fat> 1e-8])))
    # Then we need to make a single point for each parameter
    # First the sbr at the max radii
    # We will say that the maximum extend of the disc is where it drops below 5% of the maximum
    if np.sum(sbr) == 0.:
        rmaxmod = radii[-1]
        sbrvaluse =5e-6
    else:
        rmaxmod= np.max([radii[np.where(sbr > np.max(sbr)*0.05)[0][-1]+1],radii[np.where(sbr2 > np.max(sbr2)*0.05)[0][-1]+1]])
        #as the outer ring in the database can go crazy we cut the profiles at this radius
        sbr[radii > rmaxmod] = 0.
        vrot[radii > rmaxmod] = 0.
        vrot_err[radii > rmaxmod] = 0.
        scaleheight[radii > rmaxmod] = 0.
        scaleheight_err[radii > rmaxmod] = 0.
        sbr[radii > rmaxmod] = 0.
        sbr_err[radii > rmaxmod] = 0.
        inc[radii > rmaxmod] = 0.
        inc_err[radii > rmaxmod] = 0.
        pa[radii > rmaxmod] = 0.
        pa_err[radii > rmaxmod] = 0.
        RA[radii > rmaxmod] = 0.
        DEC[radii > rmaxmod] = 0.
        vsys[radii > rmaxmod] = 0.
        disp[radii > rmaxmod] = 0.
        disp_err[radii > rmaxmod] = 0.
        sbr2[radii > rmaxmod] = 0.
        vrot2[radii > rmaxmod] = 0.
        vrot2_err[radii > rmaxmod] = 0.
        scaleheight2[radii > rmaxmod] = 0.
        scaleheight2_err[radii > rmaxmod] = 0.
        sbr2[radii > rmaxmod] = 0.
        sbr2_err[radii > rmaxmod] = 0.
        inc2[radii > rmaxmod] = 0.
        inc2_err[radii > rmaxmod] = 0.
        pa2[radii > rmaxmod] = 0.
        pa2_err[radii > rmaxmod] = 0.
        RA2[radii > rmaxmod] = 0.
        DEC2[radii > rmaxmod] = 0.
        vsys2[radii > rmaxmod] = 0.
        disp2[radii > rmaxmod] = 0.
        disp2_err[radii > rmaxmod] = 0.

    #print(rmaxmod,np.where(sbr > np.max(sbr)/20.)[0][-1])
        sbrmod =  interpolate.interpolate.interp1d(radii,sbr,fill_value="extrapolate")
        sbrvaluse = sbrmod(radiifat[radiifat == fatmax])

    # We need to check whether the profiles have brightness

    #if len(R_maxval) == 0:

    R_maxval[i,:] = [(rmaxmod-fatmax)/(cube_in[0].header['BMAJ']*3600.),sbrvaluse/1e-5]
    ch_width = cube_in[0].header["CDELT3"]
    if ch_width > 500:
        ch_width = ch_width/1000.
    # next up the difference in integratefdflux
    # which we get from the produced cubes not the fitted profiles
    # First calculate the flux within the mask
    beamarea=(np.pi*abs(cube_in[0].header['BMAJ']*cube_in[0].header['BMIN']))/(4.*np.log(2.))
    pixperbeam = beamarea/(abs(cube_in[0].header['CDELT1'])*abs(cube_in[0].header['CDELT2']))
    try:
        mask =  fits.open('mask.fits',uint = False, do_not_scale_image_data=True,ignore_blank = True)
        totalflux = np.sum(cube_in[0].data[mask[0].data > 0.5])/pixperbeam*cube_in[0].header["CDELT3"]/1000.
    except:
        #totalflux = sf.get_flux_from_info(dirname[i])
        try:
            mom0= fits.open('Finalmodel/Finalmodel_mom0.fits')
        except:
            try:
                mom0= fits.open('Moments/Convolved_Cube_preprocessed_mom0.fits')
            except:
                try:
                    mom0= fits.open('Moments/Convolved_Cube_preprocessed_mom0_small.fits')
                except:
                    try:
                        mom0= fits.open('Moments/Convolved_Cube_mom0_small.fits')
                    except:
                        mom0= fits.open('Moments/Convolved_Cube_mom0.fits')
        totalflux = np.sum(mom0[0].data)/pixperbeam
    # for fat and barolo we take the total of the moment 0
    try:
        mom0_fat= fits.open('Finalmodel/Finalmodel_mom0.fits')
    except:
        mom0_fat= fits.open('Moments/Finalmodel_mom0.fits')
    totalflux_fat = np.sum(mom0_fat[0].data)/pixperbeam
    sbrval[i,:] = [totalflux-totalflux_fat,totalflux/100.]

    # next up evaluate the rotation curve we do this from 1 beam out to the maximum radius of the shortest fit
    evalrad=copy.deepcopy(radiifat)
    #Apparently there are 0's among the errors
    zeros = np.where(vrotfat_err == 0.)[0]
    if len(zeros) > 0:
        vrotfat_err[zeros] = ch_width*0.5

    vrotval[i,:] = sf.get_diff(vrotfat,vrot, radii= radiifat, model_radii= radii, errors =vrotfat_err ,norm = ch_width)
    if not np.isfinite(vrotval[i,0]):
        print(vrotfat,vrot,vrotfat_err)
        exit()


    incval[i,:] = sf.get_diff(incfat,inc, radii= radiifat, model_radii= radii, errors =incfat_err \
                                ,second = inc2fat, second_model = inc2, second_errors = inc2fat_err)
    incarrange[i] = inc[0]
    # then possible the palination
    paval[i,:] = sf.get_diff(pafat,pa, radii= radiifat, model_radii= radii, errors =pafat_err \
                                ,second = pa2fat, second_model = pa2, second_errors = pa2fat_err)
    # then possible the dispersion
    if np.sum(disp) > 0:
        dispval[i,:] = sf.get_diff(dispfat,disp, radii= radiifat, model_radii= radii, errors =dispfat_err \
                                ,second = disp2fat, second_model = disp2, second_errors = disp2fat_err,norm = ch_width)
    else:
        dispval[i,:] = float('NaN')
# then possible the Z0
    if np.sum(scaleheight) > 0:
        scaleheightval[i,:] = sf.get_diff(scaleheightfat,scaleheight, radii= radiifat, model_radii= radii, errors =scaleheightfat_err \
                                ,second = scaleheight2fat, second_model = scaleheight2, second_errors = scaleheight2fat_err,norm =cube_in[0].header["BMAJ"]*3600. )
    else:
        scaleheightval[i,:] = float('NaN')


    #Error of 10 percent of the beam and channel on the central coordinates
    RAval[i,:]=[(RA[0]-RAfat[0])/cube_in[0].header["BMAJ"],0.1]
    DECval[i,:]=[(DEC[0]-DECfat[0])/cube_in[0].header["BMAJ"],0.1]
    vsysval[i,:]= [(vsys[0]-vsysfat[0])/ch_width,0.1]
    #print(cube_in[0].header["CDELT3"])

    os.chdir(f"{Template_in['MAINDIR']}")

labelfont= {'family':'Times New Roman',
            'weight':'normal',
            'size':28}
plt.rc('font',**labelfont)
plt.figure(89,figsize=(24,24),dpi=300,facecolor = 'w', edgecolor = 'k')
#this runs counter to figsize. How can a coding language be this illogical?
gs = gridspec.GridSpec(3,3 )
gs.update(wspace=0.2, hspace=0.2)
#First the RA and DEC
RAval = np.array(RAval,dtype=float)
DECval = np.array(DECval,dtype=float)

ax = sf.make_plot(RAval,DECval, color=incarrange, status= fitted,location = gs[0], symbol=rcshapes,\
                    xlabel = '$\Delta$ RA (beams)',ylabel = '$\Delta$ DEC (beams)', colorbarlabel = 'Inclination', No_Mean = True)


# Then the beam vs Delt Central
ax = sf.make_plot(beamarrange, np.sqrt(RAval**2+DECval**2), color=incarrange, status= fitted,location = gs[1], symbol=rcshapes,\
                    xlabel = 'Diameter (beams)',ylabel = '$\Delta$ Central (beams)', colorbarlabel = 'Inclination')

# make a legend
labelfont= {'family':'Times New Roman',
            'weight':'normal',
            'size':18}
plt.rc('font',**labelfont)
chartBox = ax.get_position()
ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*1.0, chartBox.height])
ax.legend(loc='upper left', bbox_to_anchor=(1.25, 1.0), shadow=True, ncol=1)
# Beams vs delta inclination
labelfont= {'family':'Times New Roman',
            'weight':'normal',
            'size':28}
plt.rc('font',**labelfont)


# Beam vs vsys
ax = sf.make_plot(beamarrange, vsysval, color=incarrange, status= fitted,location = gs[3], symbol=rcshapes,\
                    xlabel = 'Diameter (beams)',ylabel = '$\Delta$ vsys (channels)', colorbarlabel = 'Inclination')

ax = sf.make_plot(beamarrange, incval, color=incarrange, status= fitted,location = gs[4], symbol=rcshapes,\
                    xlabel = 'Diameter (beams)',ylabel = '$\Delta$ $i$ ($^{\circ}$)', colorbarlabel = 'Inclination')
# Beams vs Delta PA
ax = sf.make_plot(beamarrange, paval, color=incarrange, status= fitted,location = gs[6], symbol=rcshapes,\
                    xlabel = 'Diameter (beams)',ylabel = '$\Delta$ PA ($^{\circ}$)', colorbarlabel = 'Inclination')

# Beams vs Delta VROT

ax = sf.make_plot(beamarrange, vrotval, color=incarrange, status= fitted,location = gs[7], symbol=rcshapes,\
                    xlabel = 'Diameter (beams)',ylabel = r'$\Delta$ V$_{\rm rot}$  (channels)', colorbarlabel = 'Inclination')
# Beams vs Delta SBR
ax = sf.make_plot(beamarrange, sbrval, color=incarrange, status= fitted,location = gs[8], symbol=rcshapes,\
                    xlabel = 'Diameter (beams)',ylabel = r'$\Delta$ Tot Flux  (Jy/beam km/s)', colorbarlabel = 'Inclination')


#Make a color bar for the inlination

ax = plt.subplot(gs[5])
#Make a color bar for the inlination
a = np.array([[0,90.]])
img = plt.imshow(a, cmap="rainbow")
plt.gca().set_visible(False)
cax = plt.axes([0.63, 0.4, 0.01, 0.45])
barr = plt.colorbar(orientation="vertical", cax=cax)
barr.set_label('Inclination', rotation=270, verticalalignment='bottom')


labelfont= {'family':'Times New Roman',
            'weight':'normal',
            'size':37}
plt.rc('font',**labelfont)
plt.figtext(0.5,0.91,f'Out of {len(incarrange)} galaxies, {len(np.where(fitted > 0.1)[0])} were succesfully fitted', horizontalalignment='center')

#print([b[0] for item in RAval for b in item[0]])
#RAval = np.array(RAval)
print(RAval.shape)

plt.savefig('Release_v2.0_All_1.png', bbox_inches='tight')
plt.close()
labelfont= {'family':'Times New Roman',
            'weight':'normal',
            'size':28}
plt.rc('font',**labelfont)
plt.figure(89,figsize=(24,24),dpi=300,facecolor = 'w', edgecolor = 'k')
#this runs counter to figsize. How can a coding language be this illogical?
gs = gridspec.GridSpec(3,3 )
gs.update(wspace=0.2, hspace=0.2)



# Beams vs Delta Dispersion
ax = sf.make_plot(beamarrange, dispval, color=incarrange, status= fitted,location = gs[0], symbol=rcshapes,\
                    xlabel = 'Diameter (beams)',ylabel = r'$\Delta$ Dispersion  (channels)', colorbarlabel = 'Inclination')

# Beams vs Delta scale height
ax = sf.make_plot(beamarrange, scaleheightval, color=incarrange, status= fitted,location = gs[1], symbol=rcshapes,\
                    xlabel = 'Diameter (beams)',ylabel = r'$\Delta$ Scaleheight (beams)', colorbarlabel = 'Inclination')

# make a legend
labelfont= {'family':'Times New Roman',
            'weight':'normal',
            'size':18}
plt.rc('font',**labelfont)
chartBox = ax.get_position()
ax.set_position([chartBox.x0, chartBox.y0, chartBox.width*1.0, chartBox.height])
ax.legend(loc='upper left', bbox_to_anchor=(1.25, 1.0), shadow=True, ncol=1)
# Beams vs delta inclination
labelfont= {'family':'Times New Roman',
            'weight':'normal',
            'size':28}
plt.rc('font',**labelfont)



# Beams vs R_max
ax = sf.make_plot(beamarrange, R_maxval, color=incarrange, status= fitted,location = gs[3], symbol=rcshapes,\
                    xlabel = 'Diameter (beams)',ylabel = r'$\Delta$ R$_{\rm max}$  (beams)', colorbarlabel = 'Inclination')
# Beams vs R_max

fittoplot = np.array([[x,0.] for x in fitted])

ax = sf.make_plot(beamarrange, fittoplot, color=incarrange, status= fitted,location = gs[4], symbol=rcshapes,\
                    xlabel = 'Diameter (beams)',ylabel = r'Fit Result', colorbarlabel = 'Inclination')

# SNR vs Flux
ax = sf.make_plot(SNRarrange, sbrval, color=incarrange, status= fitted,location = gs[6], symbol=rcshapes,\
                    xlabel = 'SNR',ylabel = r'$\Delta$ Tot Flux  (Jy/beam km/s)', colorbarlabel = 'Inclination')
#

#SNR vs vrot
ax = sf.make_plot(SNRarrange, vrotval, color=incarrange, status= fitted,location = gs[7], symbol=rcshapes,\
                    xlabel = 'SNR',ylabel = r'$\Delta$ V$_{\rm rot}$  (channels)', colorbarlabel = 'Inclination')
# B
#SNR vs Rmax
ax = sf.make_plot(SNRarrange, R_maxval, color=incarrange, status= fitted,location = gs[8], symbol=rcshapes,\
                    xlabel = 'SNR',ylabel = r'$\Delta$ R$_{\rm max}$  (beams)', colorbarlabel = 'Inclination')
# SNR vs fitresult
#ax = sf.make_plot(SNRarrange, fittoplot, color=incarrange, status= fitted,location = gs[15], symbol=rcshapes,\
#                    xlabel = 'SNR',ylabel = r'Fit Result', colorbarlabel = 'Inclination')

#Make a color bar for the inlination

ax = plt.subplot(gs[5])
#Make a color bar for the inlination
a = np.array([[0,90.]])
img = plt.imshow(a, cmap="rainbow")
plt.gca().set_visible(False)
cax = plt.axes([0.63, 0.4, 0.01, 0.45])
barr = plt.colorbar(orientation="vertical", cax=cax)
barr.set_label('Inclination', rotation=270, verticalalignment='bottom')


labelfont= {'family':'Times New Roman',
            'weight':'normal',
            'size':37}
plt.rc('font',**labelfont)
plt.figtext(0.5,0.91,f'Out of {len(incarrange)} galaxies, {len(np.where(fitted > 0.1)[0])} were succesfully fitted', horizontalalignment='center')

#print([b[0] for item in RAval for b in item[0]])
#RAval = np.array(RAval)
print(RAval.shape)

plt.savefig('Release_v2.0_All_2.png', bbox_inches='tight')
plt.close()

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

for i,key in enumerate(RCs):
    ax = plt.subplot(gs[i])
    failed = 0
    tot = 0
    for indi in RCs[key]:
        print(indi)
        tot += 1
        kpcradii = sf.convertskyangle(RCs[key][indi]['RADIUS'],distance= RCs[key][indi]['DISTANCE'])

        if indi == 'MODEL':
            ax.plot(kpcradii[RCs[key][indi]['RC'] > 0.], RCs[key][indi]['RC'][RCs[key][indi]['RC'] > 0.], 'r',zorder= 2)
            ax.plot(kpcradii[RCs[key][indi]['RC'] > 0.], RCs[key][indi]['RC'][RCs[key][indi]['RC'] > 0.], 'ro',zorder=2)
        else:
            if RCs[key][indi]['STATUS'] == 0:
                failed += 1
            elif RCs[key][indi]['STATUS'] == 1:
                #ymin, ymax = ax.get_ylim()
                ax.plot(kpcradii[RCs[key][indi]['RC'] > 0.], RCs[key][indi]['RC'][RCs[key][indi]['RC'] > 0.], 'k--',zorder= 1, alpha =0.5)
                #ax.set_ylim(ymin,ymax)
            else:
                ax.plot(kpcradii[RCs[key][indi]['RC'] > 0.], RCs[key][indi]['RC'][RCs[key][indi]['RC'] > 0.], 'k',zorder= 1, alpha =0.75)
    ax.set_xlabel('Radius (kpc)', **labelfont)
    ax.set_ylabel('V$_{rot}$ (km s$^{-1}$)', **labelfont)
    ax.set_title(key)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin,ymax+(ymax-ymin)/10.)
    ax.text(0.95,0.95,f'Out of {tot} galaxies, {failed} failed to fit. ', transform=ax.transAxes,horizontalalignment= 'right', verticalalignment='top')

plt.savefig('Release_v2.0_RC.png', bbox_inches='tight')
plt.close()

#It took 2454.8806738853455 seconds to run this
