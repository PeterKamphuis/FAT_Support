Pro Compare

compile_opt idl2

main_dir='/home/peter/FAT_Main/Test_Sets/From_Bochum/LVHIS-26/'


catname=main_dir+'Galaxies.txt'

spacecheck=strtrim(str_sep(main_dir,' '),2)
main_dirsl=STRJOIN(spacecheck,'\ ')


RESOLVE_ROUTINE,'gettirific'
RESOLVE_ROUTINE,'getrc'
RESOLVE_ROUTINE,'getdiskfit'
RESOLVE_ROUTINE,'plotparameters'
RESOLVE_ROUTINE,'plotdifflinev3'
RESOLVE_ROUTINE,'plotdiffcircv3'
RESOLVE_ROUTINE,'interpolate'
RESOLVE_ROUTINE,'calcavv2'
RESOLVE_ROUTINE,'fat_ploterror'
;RESOLVE_ROUTINE,'calcav'
RESOLVE_ROUTINE,'calcav_vobs'
RESOLVE_ROUTINE,'convert_ps'

readcol,main_dir+'Gal_Names_Sub.txt',format='(A,A,A,A)',HPASSNameDF,JNameDF,JNumberDF,OthernameDF,skipline=1
;readcol,'../Galaxy_Names.txt',format='(A,A,A,A)',HPASSNameDF,JNameDF,JNumberDF,OthernameDF,skipline=1
readcol,main_dir+'Galaxy_Names.txt',format='(A,A,A,A)',HPASSName,JName,JNumber,Othername,skipline=1
;readcol,'/DATA/SERPENS_2/kam036/Warps/CatalogueData/Data_LVHIS.txt',format='(F,A,A,A,F,F,F,F,F,F,F,F,F,F,A,A,A,A)',number,RAHIPASS,DECHIPASS,epoch,PAdata,inclinationdata,inclinationdevdata,Systemicvelocitydata,MaxRotdata,Distance,noisedata,BMMajdata,BMMinndata,BeamPAdata,Dirnamedata,Cubenamedata,Moment0data,Moment1data,skipline=1,delimiter='|'
Distance=[3.7,2.6,1.8,4.8,3.4,6.1,0.0,4.5,6.5,5.3,4.9,5.1,18.8,8.2,4.9,5.9,5.3,0.0,6.0,6.8,4.0,4.1,3.2,4.2,6.8]
litmajax=[1.90000,3.00000,6.00000,0.600000,2.40000,1.00000,0.00000,2.10000,3.20000,3.50000,3.50000,1.60000,12.0000,3.80000,2.60000,4.00000,2.00000,3.00000,4.00000,9.00000,12.0000,26.0000,32.0000,7.00000,18.0000]


first=2
makeplots=0
;Ok so first we collect the Tirific info from Basic
HPASSNAMEin1=strarr(n_elements(HPASSname))
HPASSNAMEin2=strarr(n_elements(HPASSname))
RA=strarr(n_elements(HPASSname))
DEC=strarr(n_elements(HPASSname))
RAdeg=dblarr(n_elements(HPASSname))
DECdeg=dblarr(n_elements(HPASSname))
vsys=strarr(n_elements(HPASSname))
vrot=strarr(n_elements(HPASSname))
pa=strarr(n_elements(HPASSname))
incl=strarr(n_elements(HPASSname))
Z0=dblarr(n_elements(HPASSname))
channelwidth=dblarr(n_elements(HPASSname))
SDIS=dblarr(n_elements(HPASSname))
pahalf=dblarr(n_elements(HPASSname))
inclhalf=dblarr(n_elements(HPASSname))
RAerr=strarr(n_elements(HPASSname))
DECerr=strarr(n_elements(HPASSname))
vsyserr=strarr(n_elements(HPASSname))
vroterr=strarr(n_elements(HPASSname))
paerr=strarr(n_elements(HPASSname))
inclerr=strarr(n_elements(HPASSname))
radius=dblarr(n_elements(HPASSname))
radiuskpc=dblarr(n_elements(HPASSname))
beam=dblarr(n_elements(HPASSname))
vrotout=dblarr(n_elements(HPASSname))
flux=dblarr(n_elements(HPASSname))
pixperbeam=dblarr(n_elements(HPASSname))
Mdyn=dblarr(n_elements(HPASSname))
MHI=dblarr(n_elements(HPASSname))
;next up are our rotcur solutions
RCHPASSNAMEin=strarr(n_elements(HPASSname))

RCRA=strarr(n_elements(HPASSname))
RCDEC=strarr(n_elements(HPASSname))
RCRAdeg=strarr(n_elements(HPASSname))
RCDECdeg=strarr(n_elements(HPASSname))


RCvsys=strarr(n_elements(HPASSname))
RCvrot=strarr(n_elements(HPASSname))
RCpa=strarr(n_elements(HPASSname))
RCincl=strarr(n_elements(HPASSname))
RCpahalf=strarr(n_elements(HPASSname))
RCinclhalf=strarr(n_elements(HPASSname))
RCRAerr=strarr(n_elements(HPASSname))
RCDECerr=strarr(n_elements(HPASSname))
RCvsyserr=strarr(n_elements(HPASSname))
RCvroterr=strarr(n_elements(HPASSname))
RCpaerr=strarr(n_elements(HPASSname))
RCinclerr=strarr(n_elements(HPASSname))
RCradius=dblarr(n_elements(HPASSname))
SHOrotcurdir=main_dir+'rotcur_outputs_LVHIS_subsample'
SHOvelfields=main_dir+'Hermite\ velocity\ fields\ LVHIS\ rotcur\ fits/'

;next up are our rotcur solutions
DFHPASSNAMEin=strarr(n_elements(HPASSname))
DFRA=strarr(n_elements(HPASSname))
DFDEC=strarr(n_elements(HPASSname))
DFRAdeg=strarr(n_elements(HPASSname))
DFDECdeg=strarr(n_elements(HPASSname))
DFvsys=strarr(n_elements(HPASSname))
DFvrot=strarr(n_elements(HPASSname))
DFpa=strarr(n_elements(HPASSname))
DFincl=strarr(n_elements(HPASSname))
DFpahalf=strarr(n_elements(HPASSname))
DFinclhalf=strarr(n_elements(HPASSname))
DFRAerr=dblarr(n_elements(HPASSname))
DFDECerr=dblarr(n_elements(HPASSname))
DFvsyserr=strarr(n_elements(HPASSname))
DFvroterr=strarr(n_elements(HPASSname))
DFpaerr=strarr(n_elements(HPASSname))
DFinclerr=strarr(n_elements(HPASSname))
RCRA=strarr(n_elements(HPASSname))
RCDEC=strarr(n_elements(HPASSname))
RCRAdeg=strarr(n_elements(HPASSname))
RCDECdeg=strarr(n_elements(HPASSname))
DFradius=dblarr(n_elements(HPASSname))
DFrotcurdir=main_dir+'Diskfit/'
DFinclrms=dblarr(n_elements(HPASSname))
DFinclrmsmean=dblarr(n_elements(HPASSname))
DFinclrmsbm=dblarr(n_elements(HPASSname))
DFinclrmsavbm=dblarr(n_elements(HPASSname))
DFinclrmsav=dblarr(n_elements(HPASSname))
DFparms=dblarr(n_elements(HPASSname))
DFparmsmean=dblarr(n_elements(HPASSname))
DFparmsav=dblarr(n_elements(HPASSname))
DFvrotrms=dblarr(n_elements(HPASSname))
DFvrotrmsmean=dblarr(n_elements(HPASSname))
DFvrotrmsav=dblarr(n_elements(HPASSname))
DFvobsrms=dblarr(n_elements(HPASSname))
DFvobsrmsmean=dblarr(n_elements(HPASSname))
DFvobsrmsav=dblarr(n_elements(HPASSname))
DFRCvrotrmsav=dblarr(n_elements(HPASSname))
DFRCinclrms=dblarr(n_elements(HPASSname))
DFRCinclrmsmean=dblarr(n_elements(HPASSname))
DFRCinclrmsav=dblarr(n_elements(HPASSname))
DFRCinclrmsbm=dblarr(n_elements(HPASSname))
DFRCinclrmsavbm=dblarr(n_elements(HPASSname))
DFRCparms=dblarr(n_elements(HPASSname))
DFRCparmsmean=dblarr(n_elements(HPASSname))
DFRCparmsav=dblarr(n_elements(HPASSname))
DFRCvrotrms=dblarr(n_elements(HPASSname))
DFRCvrotrmsmean=dblarr(n_elements(HPASSname))
DFRCvrotrmsav=dblarr(n_elements(HPASSname))
DFRCvobsrms=dblarr(n_elements(HPASSname))
DFRCvobsrmsmean=dblarr(n_elements(HPASSname))
DFRCvobsrmsav=dblarr(n_elements(HPASSname))
RCinclrms=dblarr(n_elements(HPASSname))
RCinclrmsmean=dblarr(n_elements(HPASSname))
RCinclrmsav=dblarr(n_elements(HPASSname))
RCinclrmsbm=dblarr(n_elements(HPASSname))
RCinclrmsavbm=dblarr(n_elements(HPASSname))
RCparms=dblarr(n_elements(HPASSname))
RCparmsmean=dblarr(n_elements(HPASSname))
RCparmsav=dblarr(n_elements(HPASSname))
RCvrotrms=dblarr(n_elements(HPASSname))
RCvrotrmsmean=dblarr(n_elements(HPASSname))
RCvrotrmsav=dblarr(n_elements(HPASSname))
RCvobsrms=dblarr(n_elements(HPASSname))
RCvobsrmsmean=dblarr(n_elements(HPASSname))
RCvobsrmsav=dblarr(n_elements(HPASSname))

openw,45,main_dir+'info_onRC.txt'
printf,45,format='(A10,8A10)','Name','RApos','DECpos','3 PA','End PA','3 INC','END INC','3 Radius','End radius'
close,45
openw,45,main_dir+'info_onDF.txt'
printf,45,format='(A10,8A10)','Name','RApos','DECpos','3 PA','End PA','3 INC','END INC','3 Radius','End radius'
close,45
for i=0,n_elements(HPASSname)-1 do begin
;for i=21,24 do begin
   ;first we'll read in the parameters from the 1stfit.
   print,'this is the galaxy we are going to look at'
   print,HPASSName[i]
 print,'this is the galaxy we are going to look at'
   print,HPASSName[i]
   ;let's obtain the channelwidth
   Cube=readfits(main_dir+'/'+HPASSname[i]+'/Cube.fits',hed)
   case strtrim(sxpar(hed,'CUNIT3'),2) OF
      'M/S': begin
         channelwidth[i]=sxpar(hed,'CDELT3')/1000.
      end
      'KM/S': channelwidth[i]=sxpar(hed,'CDELT3')
      ELSE:BEGIN
         IF sxpar(hed,'CDELT3') GT 500. then channelwidth[i]=sxpar(hed,'CDELT3')/1000. else channelwidth[i]=sxpar(hed,'CDELT3')
      END
   ENDCASE

   IF FILE_TEST(main_dir+'/'+HPASSname[i]+'/No_Warp/No_Warp.def') then begin
      Tirparameters=['RADI','XPOS','YPOS','VSYS','INCL','PA','VROT','SBR','SBR_2','BMAJ']
      gettirific,main_dir+'/'+HPASSname[i]+'/No_Warp/No_Warp.def',Tirparameters,Tirresult
      HPASSNAMEin1[i]=HPASSname[i]
   ENDIF ELSE BEGIN
      print,"We couldn't find the No Warp model Model"
      HPASSNAMEin1[i]='-1'
      Tirresult=dblarr(2,10)
   ENDELSE
   IF FILE_TEST(main_dir+'/'+HPASSname[i]+'/Finalmodel/Finalmodel.def') then begin
      Tirparameters2=['RADI','XPOS','YPOS','VSYS','INCL','PA','VROT','SBR','INCL_2','PA_2','VROT_2','SBR_2','INCL_ERR','INCL_2_ERR','PA_ERR','PA_2_ERR','VROT_ERR','VROT_2_ERR','BMAJ','Z0','SDIS']
      gettirific,main_dir+'/'+HPASSname[i]+'/Finalmodel/Finalmodel.def',Tirparameters2,Tirresult2,/errors
      print,'are we actually doing this?'
      IF FILE_TEST(main_dir+'/'+HPASSname[i]+'/Intermediate/Finalmodel_unsmoothed.def') then begin
         gettirific,main_dir+'/'+HPASSname[i]+'/Intermediate/Finalmodel_unsmoothed.def',Tirparameters2,Tirunsmooth
      ENDIF ELSE BEGIN
         IF FILE_TEST(main_dir+'/'+HPASSname[i]+'/2ndfituncor.def') then begin
            gettirific,main_dir+'/'+HPASSname[i]+'/2ndfituncor.def',Tirparameters2,Tirunsmooth
         ENDIF
      ENDELSE
      tmp=readfits(main_dir+'/'+HPASSname[i]+'/Moments/Finalmodel_mom0.fits',hed)
      flux[i]=TOTAL(tmp)
      print,flux[i]

      beamarea=(!pi*ABS(double(sxpar(hed,'BMAJ'))*double(sxpar(hed,'BMIN'))))/(4.*ALOG(2.))
      pixperbeam[i]=beamarea/(ABS(sxpar(hed,'CDELT1'))*ABS(sxpar(hed,'CDELT2')))
      HPASSNAMEin2[i]=HPASSname[i]
 
   ENDIF ELSE begin
      HPASSNAMEin2[i]='-1'
      Tirresult2=dblarr(2,19)
      if HPASSNAMEin1[i] EQ'-1' then begin
         print,'We found no models, aborting'
         stop
      endif
      first=1
   ENDELSE

   print,Tirresult2[*,0]

   name=StrJoin(StrSplit(othername[i], '-', /Regex, /Extract, $
                          /Preserve_Null), '_')
   IF STRMID(name,0 , 1 ) NE 'e' then name=STRUPCASE(name)
   print,name

   IF name eq 'eso154_g23' then name='J0256_54'
   IF name eq 'NGC1313' then name='J0317_66'
   IF name eq 'J1305_49F1' then name='J1305_49'
   IF name eq 'NGC253' then name='J0047_25'
   IF name eq 'CIRC_MAP_1' then name='J1413_65'
   IF name eq 'M83MAP_1' then name='J1337_29'

 ;  print,name
   CD,SHOrotcurdir,current=old_dir
   spawn,'ls '+name+'*',output
   if output then begin
     ; print,output
      tmp=output
      output=tmp[0]
      getrc,output,RCParameters,RCresult,main_dir+'/'+HPASSname[i]+'/Cube.fits'
      tmpp=WHERE(RCParameters EQ 'PA')
      IF RCResult[0,tmpp] LT 0. then  RCResult[*,tmpp]=RCResult[*,tmpp]+360.
      RCHPASSNAMEin[i]=HPASSname[i]
      velfielddir=JName[i]+'-'+StrJoin(StrSplit(Jnumber[i], '--', /Regex, /Extract, $
                          /Preserve_Null), '-')
      CD,SHOvelfields+'/'+velfielddir,current=backhere
      spawn,'ls *filt.fits',outputfits

      if outputfits then begin
         print,'cp '+outputfits+' '+main_dirsl+HPASSname[i]+'/OriginalVF.fits'
         spawn,'cp '+outputfits+' '+main_dirsl+HPASSname[i]+'/OriginalVF.fits'
      ENDIF
      spawn,'ls *model.fits',outputfits
      if outputfits then begin
         spawn,'cp '+outputfits+' '+main_dirsl+HPASSname[i]+'/OriginalRC.fits'
         spawn,'cp '+SHOrotcurdir+'/'+output+' '+main_dirsl+HPASSname[i]+'/Parameters.RC'
      ENDIF
      CD,backhere
      print,'We have a rotcur file'
   ENDIF ELSE BEGIN
      RCHPASSNAMEin[i]='-1'
      RCParameters=0.
      RCResult=dblarr(2,13)
      print,'We do not have a rotcur file'
   ENDELSE


   CD,old_dir
   checkDF=WHERE(HPASSNAME[i] EQ HPASSNAMEDF)
   name=StrSplit(othername[i], '-', /Regex, /Extract, $
                          /Preserve_Null)
   IF STRMID(name[0],0 , 1 ) NE 'e' then name=STRUPCASE(name[0]) else name=name[0]+'g'
   IF name EQ 'NGC5237' then name='N5237'
   IF name eq 'eso154g' then name='J0256'
   IF name eq 'NGC1313' then name='J0317'
   IF name eq 'J1305_49F1' then name='J1305'
   IF name eq 'NGC253' then name='J0047'
   IF name eq 'CIRC_MAP_1' then name='J1413'
   IF name eq 'M83MAP_1' then name='J1337'
   print,DFrotcurdir+name,checkDF[0],FILE_TEST(DFrotcurdir+name)
   IF FILE_TEST(DFrotcurdir+name) AND checkDF[0] NE -1 then begin
     ; print,name
      CD,DFrotcurdir+name,current=old_dir
      spawn,'ls *.out',output
      print,'Mysterious'
      print,output
      if output then begin
         getdiskfit,output,DFParameters,DFresult,main_dir+'/'+HPASSname[i]+'/Cube.fits'
         DFHPASSNAMEin[i]=HPASSname[i]
      ENDIF ELSE begin
         DFHPASSNAMEin[i]='-1'
         DFParameters=0.
         DFResult=dblarr(2,13)
      ENDELSE
      spawn,'ls *.mod.fits',outputfits
      if outputfits then begin
         spawn,'cp '+outputfits+' '+main_dirsl+HPASSname[i]+'/OriginalDF.fits'
         spawn,'cp '+output+' '+main_dirsl+HPASSname[i]+'/Parameters.DF'

      ENDIF
      CD,old_dir
   ENDIF ELSE begin
      DFHPASSNamein[i]='-1'
      DFParameters=0.
      DFResult=dblarr(2,13)
   ENDELSE

   IF makeplots then begin
      IF FILE_TEST(main_dir+'/'+HPASSname[i]+'/No_Warp/No_Warp.def') then begin
         spawn,'mv '+main_dirsl+'/'+HPASSname[i]+'/Combinedplot1.ps '+main_dirsl+'/'+HPASSname[i]+'/Combinedplot1old.ps'
         plotparameters,main_dir+'/'+HPASSname[i]+'/Combinedplot1.ps',Tirresult,RCResult,DFresult,1,HPASSNAME[i]
      ENDIF
      IF FILE_TEST(main_dir+'/'+HPASSname[i]+'/Finalmodel/Finalmodel.def') then begin
         spawn,'mv '+main_dirsl+'/'+HPASSname[i]+'/Combinedplot2.ps '+main_dirsl+'/'+HPASSname[i]+'/Combinedplot2old.ps'
         plotparameters,main_dir+'/'+HPASSname[i]+'/Combinedplot2.ps',Tirresult2,RCResult,DFresult,2,HPASSNAME[i]
         plotparameters,main_dir+'/'+HPASSname[i]+'/Combinedplot2uns.ps',Tirresult2,RCResult,DFresult,2,HPASSNAME[i],unsmooth=Tirunsmooth
      ENDIF
   ENDIF
                                ;now we need to make some general
                                ;overview parameters for the whole
                                ;thing
   ;Let's weight them by SBR
   ;First average SBR
   IF first EQ 1 then begin
      SBR=(Tirresult[*,7]+TirResult[*,8])/2.
     ; print,Tirresult[*,7],TirResult[*,8]
      IF TOTAL(SBR) EQ 0 then SBR=1
   ;First Central coordinates
      radius[i]=TirResult[n_elements(TirResult[*,0])-1,0]+ABS((TirResult[n_elements(TirResult[*,0])-1,0]-TirResult[n_elements(TirResult[*,0])-2,0])/2.)
      beam[i]=TirResult[0,9]
      RAdeg[i]=TirResult[0,1]
      DECdeg[i]=TirResult[0,2]
      VSYS[i]=TirResult[0,3]
      PA[i]=TirResult[0,5]
      INCL[i]=TirResult[0,4]
;      print,TirResult[0,6],TirResult[0,4],TirResult[1,6],SBR
      VROT[i]=MEDIAN(TirResult[*,6])
   ENDIF else begin
      help,Tirresult2
      SBR=(Tirresult2[*,7]+TirResult2[*,11])/2.
      IF TOTAL(SBR) EQ 0 then begin
         SBR=(Tirresult[*,7]+TirResult[*,8])/2.
         beam[i]=TirResult2[0,18]
         radius[i]=0.
         Tirresult2[*,7]=1.
         Tirresult2[*,11]=1.
         PA[i]=0.
         INCL[i]=0.
         PAhalf[i]=0.
         INCLhalf[i]=0.
         VROT[i]=0.
         RAdeg[i]=0
         DECdeg[i]=0
         VSYS[i]=0.
      ENDIF ELSE BEGIN
         help,Tirresult2,SBR
   ;First Central coordinates
         radius[i]=TirResult2[n_elements(TirResult2[*,0])-1,0]+ABS((TirResult2[n_elements(TirResult2[*,0])-1,0]-TirResult2[n_elements(TirResult2[*,0])-2,0])/2.)
         beam[i]=TirResult2[0,18]
         RAdeg[i]=TirResult2[0,1]
         DECdeg[i]=TirResult2[0,2]
         VSYS[i]=TirResult2[0,3]
         ;let identify the optical size of a galaxy.
         tmp=TirResult2[0,19]
         res=convertskyanglefunction(tmp,Distance[i])
         Z0[i]=res
         SDIS[i]=double(TirResult2[0,20])
         PA[i]=MEDIAN([TirResult2[*,5],Tirresult2[*,9]])
         INCL[i]=MEDIAN([TirResult2[*,4],Tirresult2[*,8]])
         IF litmajax[i] NE 0 then begin
            tmphalf1=1
            tmphalf2=(litmajax[i]*60.)/2.
            PAav=(TirResult2[*,5]+Tirresult2[*,9])/2.
            print,PAav
            interpolate,PAav,TirResult2[*,0],output=tmphalf1,newradii=tmphalf2
            print,tmphalf2,TirResult2[n_elements(TirResult2[*,0])-1,0]
            print,'The Reul',tmphalf1

            PAhalf[i]=tmphalf1
            INCLav=(TirResult2[*,4]+Tirresult2[*,8])/2.
            interpolate,INCLav,TirResult2[*,0],output=tmphalf1,newradii=tmphalf2
            INCLhalf[i]=tmphalf1
         ENDIF ELSE BEgin
            PAhalf[i]=MEDIAN([TirResult2[0:fix((n_elements(TirResult2[*,5])-1)/2.),5],Tirresult2[0:fix((n_elements(TirResult2[*,9])-1)/2.),9]])
            INCLhalf[i]=MEDIAN([TirResult2[0:fix((n_elements(TirResult2[*,4])-1)/2.),4],Tirresult2[0:fix((n_elements(TirResult2[*,8])-1)/2.),8]])
         ENDELSE
;         PA[i]=TOTAL(Tirresult2[1:n_elements(SBR)-1,7]*TirResult2[1:n_elements(SBR)-1,5]+Tirresult2[1:n_elements(SBR)-1,11]*TirResult2[1:n_elements(SBR)-1,9])/(TOTAL(SBR[1:n_elements(SBR)-1])*2.)
;         INCL[i]=TOTAL(Tirresult2[1:n_elements(SBR)-1,7]*TirResult2[1:n_elements(SBR)-1,4]+Tirresult2[1:n_elements(SBR)-1,11]*TirResult2[1:n_elements(SBR)-1,8])/(TOTAL(SBR[1:n_elements(SBR)-1])*2.)
         print,([TirResult2[1:n_elements(SBR)-1,6],Tirresult2[1:n_elements(SBR)-1,10]])
         VROT[i]=MEDIAN([TirResult2[1:n_elements(SBR)-1,6],Tirresult2[1:n_elements(SBR)-1,10]])
         vrotout[i]=TOTAL(TirResult2[n_elements(SBR)-3:n_elements(SBR)-1,6])/3.
         radiuskpc[i]=convertskyanglefunction(radius[i],distance[i])
         parsec=3.08568025E16 ;m
         Solarmass=1.98892E30 ;kg
         G=6.67300E-11        ;m^3kg^-1 s^-2
         G2=(G*Solarmass)     ;G m^3 M_sol^-1 s^-22
;print,G2

         Mdyn[i]=((radiuskpc[i]*1E3*parsec)*(vrotout[i]*1e3)^2)/G2
         MHI[i]=2.36E5*Distance[i]^2*flux[i]/pixperbeam[i]
      ENDELSE
   ENDELSE

   ;then the others
 ;  print,first, HPASSNamein1[i]
   IF first eq 1 OR HPASSNamein2[i] EQ '-1' then radii=TirResult[*,0] else radii=Tirresult2[*,0]
   ;first make the arrays the same distribution
   IF RCHPASSNamein[i] NE '-1' then begin
    ;  print,'And the total result'
    ;  print,radii,'what'
    ;  print,RCresult[*,7],RCresult[*,0]
      tmp=1
     ; interpolate,RCresult[*,7],RCresult[*,0],output=tmp,newradii=radii
    ;  print,tmp,radii

      RCradius[i]=RCresult[n_elements(RCResult[*,0])-1,0]+ABS((RCResult[n_elements(RCResult[*,0])-1,0]-RCResult[n_elements(RCResult[*,0])-2,0])/2.)
      RCincl[i]=MEDIAN(RCresult[*,7])



   ;   RCincL[i]=TOTAL(tmp[1:n_elements(SBR)-1]*SBR[1:n_elements(SBR)-1])/TOTAL(SBR[1:n_elements(SBR)-1])
      RCincLerr[i]=TOTAL(RCresult[*,8])/n_elements(RCresult[*,8])
      tmp=1
      tmp2=1
      tmperr=1
      tmperr2=1
      tmpinc=1
      tmpinc2=1.
      tmpincerr=1
      tmpincerr2=1.
      rctmpinc=1
;INCL

      print,radius,Tirresult2[*,0]
      print,'weird'
      calcavv2,[[TirResult2[*,4]],[TirResult2[*,8]]],RCresult[*,7],radii,RCresult[*,0],average=av,mean=mean,rms=rms,error1=[[TirResult2[*,12]],[TirResult2[*,13]]],error2=RCResult[*,8]

      RCinclrms[i]=rms
;      RCinclrms[i]=rms
      RCinclrmsmean[i]=mean
      RCinclrmsav[i]=av
      print,'Firt the new ones'
      print,rms,mean,av
      print,'Then the old'
      print,RCinclrms[i],RCinclrmsmean[i],RCinclrmsav[i]

      diff=[RCresult[*,0]*COS(RCresult[*,7]*!DtoR)/beam[i]-RCresult[*,0]*COS(tmpinc*!DtoR)/beam[i],RCresult[*,0]*COS(RCresult[*,7]*!DtoR)/beam[i]-RCresult[*,0]*COS(tmpinc2*!DtoR)/beam[i]]
      print,diff

      RCinclrmsbm[i]=STDDEV(diff)
      RCinclrmsavbm[i]=TOTAL(diff)/n_elements(diff)


                                ; RCpa[i]=TOTAL(tmp[1:n_elements(SBR)-1]*SBR[1:n_elements(SBR)-1])/TOTAL(SBR[1:n_elements(SBR)-1])
      RCpa[i]=MEDIAN(RCresult[*,5])
      IF litmajax[i] NE 0 then begin
         tmphalf1=1
         tmphalf2=(litmajax[i]*60.)/2.
         PAav=RCresult[*,5]
         interpolate,PAav,RCresult[*,0],output=tmphalf1,newradii=tmphalf2
         RCpahalf[i]=tmphalf1
         INCLav=RCresult[*,7]
         interpolate,INCLav,RCresult[*,0],output=tmphalf1,newradii=tmphalf2
         RCinclhalf[i]=tmphalf1
      ENDIF ELSE BEgin
         RCinclhalf[i]=MEDIAN(RCresult[0:fix((n_elements(RCResult[*,7])-1)/2.),7])
         RCpahalf[i]=MEDIAN(RCresult[0:fix((n_elements(RCResult[*,5])-1)/2.),5])
      ENDELSE
                                ;PA
      RCpaerr[i]=TOTAL(RCresult[*,6])/n_elements(RCresult[*,6])

      calcavv2,[[TirResult2[*,5]],[TirResult2[*,9]]],RCresult[*,5],radii,RCresult[*,0],average=av,mean=mean,rms=rms,error1=[[TirResult2[*,14]],[TirResult2[*,15]]],error2=RCResult[*,6]
      IF rms GT 1e-10 then RCparms[i]=rms else  RCparms[i]=ABS((TirResult2[0,5]-RCresult[0,5])+(TirResult2[0,9]-RCresult[0,5]))/2.

      RCparmsmean[i]=mean
      RCparmsav[i]=av

      RCvsys[i]=Median(RCresult[*,1])

;      RCvsys[i]=TOTAL(tmp[1:n_elements(SBR)-1]*SBR[1:n_elements(SBR)-1])/TOTAL(SBR[1:n_elements(SBR)-1])
      RCvsyserr[i]=TOTAL(RCresult[*,2])/n_elements(RCresult[*,2])
     ;;;; vrot
      RCvrot[i]=MEDIAN(RCResult[*,3])
      RCvroterr[i]=TOTAL(RCResult[*,4])/n_elements(RCResult[*,4])
;VROT

      calcavv2,[[TirResult2[*,6]],[TirResult2[*,10]]],RCresult[*,3],radii,RCresult[*,0],average=av,mean=mean,rms=rms,error1=[[TirResult2[*,16]],[TirResult2[*,17]]],error2=RCResult[*,4]
      RCvrotrms[i]=rms
      RCvrotrmsmean[i]=mean
      RCvrotrmsav[i]=av
;;;;vobs

      Errorontir=[[SQRT(SIN(Tirresult2[*,4]*!DtoR)^2*TirResult2[*,16]^2+TirResult2[*,6]^2*COS(Tirresult2[*,4]*!DtoR)^2*TirResult2[*,12]^2)],[SQRT(SIN(Tirresult2[*,8]*!DtoR)^2*TirResult2[*,17]^2+TirResult2[*,10]^2*COS(Tirresult2[*,8]*!DtoR)^2*TirResult2[*,13]^2)]]
      ErrorRC=[SQRT(SIN(RCresult[*,7]*!DtoR)^2*RCResult[*,4]^2+RCResult[*,3]^2*COS(RCresult[*,7]*!DtoR)^2*RCResult[*,8]^2)]

     ; calcavv2,[[TirResult2[*,6]*SIN(Tirresult2[*,4]*!DtoR)],[TirResult2[*,10]*SIN(Tirresult2[*,8]*!DtoR)]],RCresult[*,3]*SIN(RCresult[*,7]*!DtoR),radii,RCresult[*,0],average=av,mean=mean,rms=rms,error1=Errorontir,error2=ErrorRC
      calcav_vobs,[[TirResult2[*,6]],[TirResult2[*,10]]],[[Tirresult2[*,4]],[Tirresult2[*,8]]],RCresult[*,3],RCresult[*,7],radii,RCresult[*,0],average=av,mean=mean,rms=rms,errorvr1=[[TirResult2[*,16]],[TirResult2[*,17]]],errorinc1=[[TirResult2[*,12]],[TirResult2[*,13]]],errorvr2=RCResult[*,4],errorinc2=RCResult[*,8],/rem_cent
      RCvobsrms[i]=rms
      RCvobsrmsmean[i]=mean
      RCvobsrmsav[i]=av
 ;     RCvobsrms[i]=SQRT(SIN(RCinclrmsav[i]*!DtoR)^2*RCvrotrms[i]^2+RCvrotrmsav[i]^2*COS(RCinclrmsav[i]*!DtoR)^2*RCinclrms[i]^2)
  ;    RCvobsrmsmean[i]=RCvrotrmsmean[i]*SIN(RCinclrmsmean[i]*!DtoR)
   ;   RCvobsrmsav[i]=RCvrotrmsav[i]*SIN(RCinclrmsav[i]*!DtoR)
      RCRAdeg[i]=RCResult[0,9]
      RCRAerr[i]=RCResult[0,10]
      RCDECdeg[i]=RCResult[0,11]
      RCDECerr[i]=RCResult[0,12]

      openu,45,main_dir+'info_onRC.txt',/APPEND
      IF n_elements(RCresult[*,7])-1 GE 3 then begin
         printf,45,format='(A10,5F10.5,3F10.2)',RCHPASSNamein[i],RCresult[0,9],RCresult[0,11],RCresult[2,5],RCresult[n_elements(RCresult[*,7])-1,5],RCresult[2,7],RCresult[n_elements(RCresult[*,7])-1,7],RCresult[2,0],RCresult[n_elements(RCresult[*,7])-1,0]+(RCresult[n_elements(RCresult[*,7])-1,0]-RCresult[n_elements(RCresult[*,7])-2,0])/2.
      Endif ELSE BEGIN
         printf,45,format='(A10,5F10.5,3F10.2)',RCHPASSNamein[i],RCresult[0,9],RCresult[0,11],0.,RCresult[n_elements(RCresult[*,7])-1,5],0.,RCresult[n_elements(RCresult[*,7])-1,7],0.,RCresult[n_elements(RCresult[*,7])-1,0]+(RCresult[n_elements(RCresult[*,7])-1,0]-RCresult[n_elements(RCresult[*,7])-2,0])/2.
      ENDELSE
      close,45


   ENDIF
   print,'I have no idea what is happening'
   IF DFHPASSNamein[i] NE '-1' then begin
      DFradius[i]=DFresult[n_elements(DFResult[*,0])-1,0]+ABS((DFResult[n_elements(DFResult[*,0])-1,0]-DFResult[n_elements(DFResult[*,0])-2,0])/2.)
      calcavv2,[[TirResult2[*,4]],[TirResult2[*,8]]],DFresult[*,7],radii,DFresult[*,0],average=av,mean=mean,rms=rms,error1=[[TirResult2[*,12]],[TirResult2[*,13]]],error2=DFResult[*,8]
      IF rms GT 1e-10 then  DFinclrms[i]=rms else   DFinclrms[i]=ABS((TirResult2[0,4]-DFresult[0,7])+(TirResult2[0,8]-DFresult[0,7]))/2.

      DFinclrmsmean[i]=mean
      DFinclrmsav[i]=av

      DFincl[i]=MEDIAN(DFresult[*,7])


 ;     DFinclrmsav[i]=TOTAL(diff)/n_elements(diff)

      diff=[DFresult[*,0]*COS(DFresult[*,7]*!DtoR)/beam[i]-DFresult[*,0]*COS(tmpinc*!DtoR)/beam[i],DFresult[*,0]*COS(DFresult[*,7]*!DtoR)/beam[i]-DFresult[*,0]*COS(tmpinc2*!DtoR)/beam[i]]
      DFinclrmsbm[i]=STDDEV(diff)
      DFinclrmsavbm[i]=TOTAL(diff)/n_elements(diff)



      interpolate,RCresult[*,7],RCresult[*,0],output=RCtmpinc,newradii=DFresult[*,0]
      diff=[RCtmpinc-DFresult[*,7]]
   ;   print,i,ROBUST_STDDEV(diff)
      help,DFRCinclrms
;      DFRCinclrms[i]=ROBUST_STDDEV(diff)
      DFRCinclrms[i]=STDDEV(diff)
      IF  DFRCinclrms[i] LT 1e-10 then   DFRCinclrms[i]=ABS([RCtmpinc[0]-DFresult[0,7]])
      DFRCinclrmsmean[i]=MEDIAN(diff)
      DFRCinclrmsav[i]=TOTAL(ABS(diff))/n_elements(diff)
      diff=[DFresult[*,0]*COS(RCtmpinc*!DtoR)/beam[i]-DFresult[*,0]*COS(DFresult[*,7]*!DtoR)/beam[i]]
      DFRCinclrmsbm[i]=STDDEV(diff)
      DFRCinclrmsavbm[i]=TOTAL(diff)/n_elements(diff)

      IF litmajax[i] NE 0 then begin
         tmphalf1=1
         tmphalf2=(litmajax[i]*60.)/2.
         PAav=DFresult[*,5]
         interpolate,PAav,DFresult[*,0],output=tmphalf1,newradii=tmphalf2
         DFpahalf[i]=tmphalf1
         INCLav=DFresult[*,7]
         interpolate,INCLav,DFresult[*,0],output=tmphalf1,newradii=tmphalf2
         DFinclhalf[i]=tmphalf1
      ENDIF ELSE BEgin
         DFinclhalf[i]=MEDIAN(DFresult[0:fix((n_elements(DFResult[*,7])-1)/2.),7])
         DFpahalf[i]=MEDIAN(DFresult[0:fix((n_elements(DFResult[*,5])-1)/2.),5])
      ENDELSE
      DFpa[i]=MEDIAN(DFresult[*,5])
    calcavv2,[[TirResult2[*,5]],[TirResult2[*,9]]],DFresult[*,5],radii,DFresult[*,0],average=av,mean=mean,rms=rms,error1=[[TirResult2[*,14]],[TirResult2[*,15]]],error2=DFResult[*,6]

      IF rms GT 1e-10 then  DFparms[i]=rms else DFparms[i]=ABS((TirResult2[0,5]-DFresult[0,5])+(TirResult2[0,9]-DFresult[0,5]))/2.
      DFparmsmean[i]=mean
      DFparmsav[i]=av

      DFvsys[i]=Median(DFresult[*,1])

;      DFvsys[i]=TOTAL(tmp[1:n_elements(SBR)-1]*SBR[1:n_elements(SBR)-1])/TOTAL(SBR[1:n_elements(SBR)-1])
      DFvsyserr[i]=TOTAL(DFresult[*,2])/n_elements(DFresult[*,2])
     ;;;; vrot
      DFvrot[i]=MEDIAN(DFResult[*,3])
      DFvroterr[i]=TOTAL(DFResult[*,4])/n_elements(DFResult[*,4])
;VROT

      calcavv2,[[TirResult2[*,6]],[TirResult2[*,10]]],DFresult[*,3],radii,DFresult[*,0],average=av,mean=mean,rms=rms,error1=[[TirResult2[*,16]],[TirResult2[*,17]]],error2=DFResult[*,4]
      DFvrotrms[i]=rms
      DFvrotrmsmean[i]=mean
      DFvrotrmsav[i]=av
;;;;vobs

      Errorontir=[[SQRT(SIN(Tirresult2[*,4]*!DtoR)^2*TirResult2[*,16]^2+TirResult2[*,6]^2*COS(Tirresult2[*,4]*!DtoR)^2*TirResult2[*,12]^2)],[SQRT(SIN(Tirresult2[*,8]*!DtoR)^2*TirResult2[*,17]^2+TirResult2[*,10]^2*COS(Tirresult2[*,8]*!DtoR)^2*TirResult2[*,13]^2)]]

      ErrorDF=[SQRT(SIN(DFresult[*,7]*!DtoR)^2*DFResult[*,4]^2+DFResult[*,3]^2*COS(DFresult[*,7]*!DtoR)^2*DFResult[*,8]^2)]
   ;   print,radii
 ;     calcavv2,[[TirResult2[*,6]*SIN(Tirresult2[*,4]*!DtoR)],[TirResult2[*,10]*SIN(Tirresult2[*,8]*!DtoR)]],DFresult[*,3]*SIN(DFresult[*,7]*!DtoR),radii,DFresult[*,0],average=av,mean=mean,rms=rms,error1=Errorontir,error2=ErrorDF

      calcav_vobs,[[TirResult2[*,6]],[TirResult2[*,10]]],[[Tirresult2[*,4]],[Tirresult2[*,8]]],DFresult[*,3],DFresult[*,7],radii,DFresult[*,0],average=av,mean=mean,rms=rms,errorvr1=[[TirResult2[*,16]],[TirResult2[*,17]]],errorinc1=[[TirResult2[*,12]],[TirResult2[*,13]]],errorvr2=DFResult[*,4],errorinc2=DFResult[*,8],/rem_cent
     ; set_plot,'x'

     ; plot,radii,[TirResult2[*,6]*SIN(Tirresult2[*,4]*!DtoR)]
     ; oplot,radii,[TirResult2[*,6]],color=200
     ; oplot,DFresult[*,0],DFresult[*,3]*SIN(DFresult[*,7]*!DtoR),linestyle=2
      DFvobsrms[i]=rms
      DFvobsrmsmean[i]=mean
      DFvobsrmsav[i]=av
      print,rms,mean,av
      print,DFResult[*,4],DFResult[*,8]

   ;   DFvobsrms[i]=SQRT(SIN(DFinclrmsav[i]*!DtoR)^2*DFvrotrms[i]^2+DFvrotrmsav[i]^2*COS(DFinclrmsav[i]*!DtoR)^2*DFinclrms[i]^2)
   ;   DFvobsrmsmean[i]=DFvrotrmsmean[i]*SIN(DFinclrmsmean[i]*!DtoR)
    ;  DFvobsrmsav[i]=DFvrotrmsav[i]*SIN(DFinclrmsav[i]*!DtoR)
    ;  DFparmsav[i]=TOTAL(diff)/n_elements(diff)
      interpolate,RCresult[*,5],RCresult[*,0],output=tmp,newradii=DFresult[*,0]
      diff=[tmp-DFresult[*,5]]
    ;  DFRCparms[i]=ROBUST_STDDEV(diff)
      DFRCparms[i]=STDDEV(diff)
      IF  DFRCparms[i] LT 1e-10 then  DFRCparms[i]=ABS(tmp[0]-DFresult[0,5])
      DFRCparmsmean[i]=MEDIAN(diff)
      DFRCparmsav[i]=TOTAL(ABS(diff))/n_elements(diff)


      DFvsys[i]=MEDIAN(DFresult[*,1])
;  interpolate,DFresult[*,7],DFresult[*,0],output=tmp,newradii=radii
   ;   DFincL[i]=TOTAL(tmp[1:n_elements(SBR)-1]*SBR[1:n_elements(SBR)-1])/TOTAL(SBR[1:n_elements(SBR)-1])
      DFincLerr[i]=TOTAL(DFresult[*,8])/n_elements(DFresult[*,8])
    ;  interpolate,DFresult[*,5],DFresult[*,0],output=tmp,newradii=radii
    ;  DFpa[i]=TOTAL(tmp[1:n_elements(SBR)-1]*SBR[1:n_elements(SBR)-1])/TOTAL(SBR[1:n_elements(SBR)-1])
      DFpaerr[i]=TOTAL(DFresult[*,6])/n_elements(DFresult[*,6])
     ; interpolate,DFresult[*,1],DFresult[*,0],output=tmp,newradii=radii
     ; DFvsys[i]=TOTAL(tmp[1:n_elements(SBR)-1]*SBR[1:n_elements(SBR)-1])/TOTAL(SBR[1:n_elements(SBR)-1])
      DFvsyserr[i]=TOTAL(DFresult[*,2])/n_elements(DFresult[*,2])
      DFvrot[i]=MEDIAN(DFResult[*,3])
      DFvroterr[i]=TOTAL(DFResult[*,4])/n_elements(DFResult[*,4])


      interpolate,RCresult[*,3],RCresult[*,0],output=tmp,newradii=DFresult[*,0]
      diff=[tmp-DFresult[*,3]]
  ;    DFRCvrotrms[i]=ROBUST_STDDEV(diff)
      DFRCvrotrms[i]=STDDEV(diff)
      DFRCvrotrmsmean[i]=MEDIAN(diff)
      DFRCvrotrmsav[i]=TOTAL(ABS(diff))/n_elements(diff)
      diff=[(tmp*SIN(RCtmpinc*!DtoR))-(DFresult[*,3]*SIN(DFresult[*,7]*!DtoR))]
    ;  DFRCvobsrms[i]=ROBUST_STDDEV(diff)
      DFRCvobsrms[i]=STDDEV(diff)
      DFRCvobsrmsmean[i]=MEDIAN(diff)
      DFRCvobsrmsav[i]=TOTAL(ABS(diff))/n_elements(diff)

      print,'I have no idea what is happening'
      DFRAdeg[i]=DFResult[0,9]
      DFRAerr[i]=DFResult[0,10]
      DFDECdeg[i]=DFResult[0,11]
      DFDECerr[i]=DFResult[0,12]
      openu,45,main_dir+'info_onDF.txt',/APPEND
      printf,45,format='(A10,5F10.5,3F10.2)',DFHPASSNamein[i],DFresult[0,9],DFresult[0,11],DFresult[2,5],DFresult[n_elements(DFresult[*,7])-1,5],DFresult[2,7],DFresult[n_elements(DFresult[*,7])-1,7],DFresult[2,0],DFresult[n_elements(DFresult[*,7])-1,0]+(DFresult[n_elements(DFresult[*,7])-1,0]-DFresult[n_elements(DFresult[*,7])-2,0])/2.
      close,45
   ENDIF

ENDFOR


;We should also read the pipeline result
close,1
print,catname
h=''
openr,1,catname
readf,1,h
fitresult=dblarr(n_elements(HPASSName))

WHILE ~EOF(1) do begin
   readf,1,h
;   print,h
   values=str_sep(strtrim(strcompress(h)),' ')
;   print,values
   tmp=WHERE(values[0] EQ HPASSName)
;   print,HPASSName[tmp],values[1]
   case 1 OF
      values[1] EQ 'True': fitresult[tmp] = 1
      values[1] EQ 'False': fitresult[tmp] = 0
      else: fitresult[tmp] = 0;
   ENDCASE
ENDWHILE
close,1

;fitresult=dblarr(n_elements(HPASSName))

;fitresult[*]=1
;fitresult[11]=0
;fitresult[14]=0
;fitresult[20]=0
;clean away the galaxies not fitted with rotcur



VOBS=double(VROT)*SIN(double(INCL)*!pi/180.)
RCVOBS=double(RCVROT)*SIN(double(RCINCL)*!pi/180.)
DFVOBS=double(DFVROT)*SIN(double(DFINCL)*!pi/180.)



!p.font=1
!p.charsize=1.1
A = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a unit circle with 16 points,
; and set the filled flag:
USERSYM, COS(A), SIN(A), /FILL
SET_PLOT, 'PS'



openw,1,main_dir+'FittedAverages'+string(first,format='(I1)')+'.txt'
printf,1,format='(19A10)','Name','PA_'+string(first,format='(I1)'),'PA_RC','PA_DF','INC_'+string(first,format='(I1)'),'INC_RC','INC_DF','RA_'+string(first,format='(I1)'),'RA_RC','RA_DF','DEC_'+string(first,format='(I1)'),'DEC_RC','DEC_DF','VSYS_'+string(first,format='(I1)'),'VSYS_RC','VSYS_DF','VROT_'+string(first,format='(I1)'),'VROT_RC','VROT_DF'
for i=0,n_elements(HPASSName)-1 do begin
   printf,1,format='(A10,18F10.4)',HPASSName[i],PA[i],RCPA[i],DFPA[i],INCL[i],RCINCL[i],DFINCL[i],RAdeg[i],RCRAdeg[i],DFRAdeg[i],DECdeg[i],RCDECdeg[i],DFDECdeg[i],VSYS[i],RCVSYS[i],DFVSYS[i],VROT[i],RCVROT[i],DFVROT[i]
endfor
close,1
tmpsort=SORT(double(INCL))
openw,1,main_dir+'FittedAverages_Sorted'+string(first,format='(I1)')+'.txt'
printf,1,format='(19A10)','Name','PA_'+string(first,format='(I1)'),'PA_RC','PA_DF','INC_'+string(first,format='(I1)'),'INC_RC','INC_DF','RA_'+string(first,format='(I1)'),'RA_RC','RA_DF','DEC_'+string(first,format='(I1)'),'DEC_RC','DEC_DF','VSYS_'+string(first,format='(I1)'),'VSYS_RC','VSYS_DF','VROT_'+string(first,format='(I1)'),'VROT_RC','VROT_DF'
for i=0,n_elements(HPASSName)-1 do begin
   printf,1,format='(A10,18F10.4)',HPASSName[tmpsort[i]],PA[tmpsort[i]],RCPA[tmpsort[i]],DFPA[tmpsort[i]],INCL[tmpsort[i]],RCINCL[tmpsort[i]],DFINCL[tmpsort[i]],RAdeg[tmpsort[i]],RCRAdeg[tmpsort[i]],DFRAdeg[tmpsort[i]],DECdeg[tmpsort[i]],RCDECdeg[tmpsort[i]],DFDECdeg[tmpsort[i]],VSYS[tmpsort[i]],RCVSYS[tmpsort[i]],DFVSYS[tmpsort[i]],VROT[tmpsort[i]],RCVROT[tmpsort[i]],DFVROT[tmpsort[i]]
endfor
close,1

print,'Dow e get here?'

print,first
first=2
openw,1,main_dir+'FittedRingDifferences'+string(first,format='(I1)')+'.txt'
printf,1,format='(10A10)','Name','PA_'+string(first,format='(I1)'),'PA_RC','PA_DF','INC_'+string(first,format='(I1)'),'INC_RC','INC_DF','VROT_'+string(first,format='(I1)'),'VROT_RC','VROT_DF'
for i=0,n_elements(HPASSName)-1 do begin
   printf,1,format='(A10,9F10.4)',HPASSName[i],DFRCPArmsmean[i],RCPArmsmean[i],DFPArmsmean[i],DFRCINCLrmsmean[i],RCINCLrmsmean[i],DFINCLrmsmean[i],DFRCVROTrmsmean[i],RCVROTrmsmean[i],DFVROTrmsmean[i]
endfor
close,1

openw,1,main_dir+'AvFittedRingDifferences'+string(first,format='(I1)')+'.txt'
printf,1,format='(13A10)','Name','PA_'+string(first,format='(I1)'),'PA_RC','PA_DF','INC_'+string(first,format='(I1)'),'INC_RC','INC_DF','VROT_'+string(first,format='(I1)'),'VROT_RC','VROT_DF','VOBS_'+string(first,format='(I1)'),'VOBS_RC','VOBS_DF'
for i=0,n_elements(HPASSName)-1 do begin
   printf,1,format='(A10,12F10.4)',HPASSName[i],DFRCPArmsav[i],RCPArmsav[i],DFPArmsav[i],DFRCINCLrmsav[i],RCINCLrmsav[i],DFINCLrmsav[i],DFRCVROTrmsav[i],RCVROTrmsav[i],DFVROTrmsav[i],DFRCvobsrmsav[i],RCvobsrmsav[i],DFvobsrmsav[i]
endfor
close,1

IF first EQ 2 then begin
   incllim=[10.,90.]
   beamlow=1.9
  !p.charsize=1.1
  !p.thick=7
  fitresultsmall=fitresult
  tmp=WHERE(fitresult EQ 1 AND incl LT incllim[0] )
  IF tmp[0] NE -1 then fitresult[tmp]=1.5
  tmp=WHERE(fitresult EQ 1 AND incl GT incllim[1] )
  IF tmp[0] NE -1 then fitresult[tmp]=1.5
  tmp=WHERE(fitresult EQ 1 AND radius/beam LT beamlow )
  IF tmp[0] NE -1 then fitresult[tmp]=1.5
  tmp=WHERE(radius/beam LT beamlow )
  IF tmp[0] NE -1 then fitresultsmall[tmp]=1.5

  Jstripped=dblarr(n_elements(Jname))
  for j=0,n_elements(Jstripped)-1 do begin
     print,Jnumber[j]
     tmp=STRMID(Jnumber[j],1,4)
     tmp2=STRMID(Jnumber[j],7,2)
     print,tmp,tmp2
     Jstripped[j]=double(tmp+tmp2)
  endfor
  Jstrippedin=SORT(jstripped)
  openw,1,main_dir+'Table_For_LVHIS_Paper.txt'
  for j=0,n_elements(Jstrippedin)-1 do begin
     print,Jstrippedin[j],fitresult[Jstrippedin[j]]
     IF Distance[Jstrippedin[j]] EQ 0 then begin
        IF fitresult[Jstrippedin[j]] EQ 1. then begin
           print,'HIPASS '+Jnumber[Jstrippedin[j]]
           print, Mdyn[Jstrippedin[j]],radiuskpc[Jstrippedin[j]],vrotout[Jstrippedin[j]]
           print, MHI[Jstrippedin[j]],Distance[Jstrippedin[j]],flux[Jstrippedin[j]],channelwidth[Jstrippedin[j]],pixperbeam[Jstrippedin[j]]
           printf,1,STRJOIN(['HIPASS '+Jnumber[Jstrippedin[j]],othername[Jstrippedin[j]],strtrim(string(vrotout[Jstrippedin[j]],format='(F10.1)'),2),'-',strtrim(string(inclhalf[Jstrippedin[j]],format='(F10.1)'),2),strtrim(string(PAhalf[Jstrippedin[j]],format='(F10.1)'),2),'-','-'],'&')+'\\'
        ENDIF
        IF fitresult[Jstrippedin[j]] EQ 1.5 then begin
           print,'HIPASS '+Jnumber[Jstrippedin[j]]
           print, Mdyn[Jstrippedin[j]],radiuskpc[Jstrippedin[j]],vrotout[Jstrippedin[j]]
           print, MHI[Jstrippedin[j]],Distance[Jstrippedin[j]],flux[Jstrippedin[j]],channelwidth[Jstrippedin[j]],pixperbeam[Jstrippedin[j]]
           printf,1,STRJOIN(['HIPASS '+Jnumber[Jstrippedin[j]]+'$^*$',othername[Jstrippedin[j]],strtrim(string(vrotout[Jstrippedin[j]],format='(F10.1)'),2),'-',strtrim(string(inclhalf[Jstrippedin[j]],format='(F10.1)'),2),strtrim(string(PAhalf[Jstrippedin[j]],format='(F10.1)'),2),'-','-'],'&')+'\\'
        ENDIF
        IF fitresult[Jstrippedin[j]] EQ 0 then begin
           printf,1,STRJOIN(['HIPASS '+Jnumber[Jstrippedin[j]]+'$^**$',othername[Jstrippedin[j]],'-','-','-','-','-','-'],'&')+'\\'
        ENDIF
     endif else begin
        IF fitresult[Jstrippedin[j]] EQ 1. then begin
           print,'HIPASS '+Jnumber[Jstrippedin[j]]
           print, Mdyn[Jstrippedin[j]],radiuskpc[Jstrippedin[j]],vrotout[Jstrippedin[j]]
           print, MHI[Jstrippedin[j]],Distance[Jstrippedin[j]],flux[Jstrippedin[j]],channelwidth[Jstrippedin[j]],pixperbeam[Jstrippedin[j]]
           printf,1,STRJOIN(['HIPASS '+Jnumber[Jstrippedin[j]],othername[Jstrippedin[j]],strtrim(string(vrotout[Jstrippedin[j]],format='(F10.1)'),2),strtrim(string(Radiuskpc[Jstrippedin[j]],format='(F10.1)'),2),strtrim(string(inclhalf[Jstrippedin[j]],format='(F10.1)'),2),strtrim(string(PAhalf[Jstrippedin[j]],format='(F10.1)'),2),strtrim(string(Mdyn[Jstrippedin[j]],format='(E10.3)'),2),strtrim(string(MHI[Jstrippedin[j]]/Mdyn[Jstrippedin[j]],format='(E10.3)'),2)],'&')+'\\'
        ENDIF
        IF fitresult[Jstrippedin[j]] EQ 1.5 then begin
           print,'HIPASS '+Jnumber[Jstrippedin[j]]
           print, Mdyn[Jstrippedin[j]],radiuskpc[Jstrippedin[j]],vrotout[Jstrippedin[j]]
           print, MHI[Jstrippedin[j]],Distance[Jstrippedin[j]],flux[Jstrippedin[j]],channelwidth[Jstrippedin[j]],pixperbeam[Jstrippedin[j]]
           printf,1,STRJOIN(['HIPASS '+Jnumber[Jstrippedin[j]]+'$^*$',othername[Jstrippedin[j]],strtrim(string(vrotout[Jstrippedin[j]],format='(F10.1)'),2),strtrim(string(Radiuskpc[Jstrippedin[j]],format='(F10.1)'),2),strtrim(string(inclhalf[Jstrippedin[j]],format='(F10.1)'),2),strtrim(string(PAhalf[Jstrippedin[j]],format='(F10.1)'),2),strtrim(string(Mdyn[Jstrippedin[j]],format='(E10.3)'),2),strtrim(string(MHI[Jstrippedin[j]]/Mdyn[Jstrippedin[j]],format='(E10.3)'),2)],'&')+'\\'
        ENDIF
        IF fitresult[Jstrippedin[j]] EQ 0 then begin
           printf,1,STRJOIN(['HIPASS '+Jnumber[Jstrippedin[j]]+'$^**$',othername[Jstrippedin[j]],'-','-','-','-','-','-'],'&')+'\\'
        ENDIF
     ENDELSE
  endfor
  close,1

  litincl=[35.396801,59.239250,51.528217,36.344627,73.109802,60.472172,0.,20.131441,90.000000,79.980278,65.143723,38.511372,70.095177,88.152992,55.560135,73.109802,64.114426,76.206825,44.460140,90.000000,36.344627,85.252037,83.193771,64.114426,2.0000000]
   litinclerr=[0.0000000,6.7687607,4.0319176,0.0000000,0.0000000,0.0000000,0.0000000,18.680689,0.0000000,1.6487579,0.0000000,4.1923065,6.9918747,4.9592209,8.5542908,2.3108139,27.769806,11.168510,19.347107,0.073509216,7.3697472,4.7479630,0.27157593,4.6811829,27.543171]
   litpa=[128.000,45.0000,100.000,0.00000,0.00000,0.00000,0.00000,0.00000,75.0000,40.0000,0.00000,0.00000,48.0000,147.000,62.0000,30.0000,72.0000,170.000,135.000,39.0000,0.00000,43.0000,52.0000,40.0000,0.00000]

   print,radius
   inputradius=(radius)/beam

   xaxisinp=inputradius*COS(incl*!DtoR)

   xaxisinphalf=(litmajax*60.)/beam*COS(inclhalf*!DtoR)
   RCxaxisinphalf=(litmajax*60.)/beam*COS(RCinclhalf*!DtoR)
   DFxaxisinphalf=(litmajax*60.)/beam*COS(DFinclhalf*!DtoR)
   litxaxisinp=(litmajax*60.)/beam*COS(litincl*!DtoR)
 ;  print, "let's do some simple tests like what is the size of the arrays we want to plot. All should be 25"
 ;  help,litincl,litinclerr,litpa,HPASSNAME,pahalf,inclhalf,RCpahalf,RCinclhalf,DFpahalf,DFinclhalf
 ;  print,"whoohoo they are"
;   litpa=[128.,90.,0.,0.,75.,36.,0.,48.,330.-180.,340.-180.,0.,325.-180.,39.]
;   litincl=[38.,31.,74.,40.,90.,66.,66.,77.,68.,66.,0.,55.,63.]
  ;Start the plotting
  ;Define the hex colors
  black = '000000'x
  grey = '808080'x
  blue = '0000FF'x
  light_blue = '7CACFF'x
  red = 'FF0000'x
  light_red = 'FFCCCC'x
   !p.charsize=1.1
   SET_PLOT, 'PS'
   DEVICE, FILENAME=main_dir+'Overviewlit.ps',/color,/PORTRAIT,/ENCAPSULATED,xsize=21,ysize=10.5,/DECOMPOSED
   ;loadct,40

   thick=3
   x=findgen(100)*10
;PA
   for i=0,n_elements(pahalf)-1 do begin
      print,HPASSNAME[i],ABS(litpa[i]-(pahalf[i]-180)),ABS(litpa[i]-pahalf[i]),pahalf[i]
      if pahalf[i] GT 180. AND (ABS(litpa[i]-(pahalf[i]-180)) LT ABS(litpa[i]-pahalf[i]))  then begin

         pahalf[i]=pahalf[i]-180.
      ENDIF
      if RCpahalf[i] GT 180. AND (ABS(litpa[i]-(RCpahalf[i]-180)) LT ABS(litpa[i]-RCpahalf[i])) then RCpahalf[i]=RCpahalf[i]-180.
      if DFpahalf[i] GT 180. AND (ABS(litpa[i]-(DFpahalf[i]-180)) LT ABS(litpa[i]-DFpahalf[i])) then DFpahalf[i]=DFpahalf[i]-180.
   endfor

;   print,inputhname,tmp,gven
 ;  print,'Need to print some otherr shit'
   tmp=WHERE(double(fitresult) EQ 1 OR fitresult EQ 1.5 AND double(litpa) NE 0.)
   tmp2=WHERE(double(RCpahalf) NE 0. AND double(litpa) NE 0.)
   tmp3=WHERE(DFHPASSNamein NE '-1' AND double(litpa) NE 0.)
   for i=0,n_elements(litpa)-1 do print,DFHPASSNamein[i],litpa[i]
;   print,PAhalf
;   print,inclhalf
;   print,xaxisinp

   print,tmp3
   maxpay=MAX(double([double(litPA[tmp]-PAhalf[tmp]),double(litPA[tmp2])-double(RCPAhalf[tmp2]),double(litPA[tmp3]-DFPAhalf[tmp3])]),min=minpay)

   maxpa=MAX(double([(double(INCL[tmp])),double(INCL[tmp3]),double(INCL[tmp2])]),min=minpa)
   maxpay=maxpay+5.
   minpay=minpay-5
   maxpa=maxpa+1.
   minpa=minpa-1
;   print,minpay,maxpay,minpa,maxpa
;   print,'we printed the pa max min which are Nan'
   plot,INCL,PA,position=[0.1,0.2,0.4,0.8],xtitle='INCL FAT (Deg.)',yrange=[minpay,maxpay],xthick=thick,ythick=thick,charthick=thick,thick=thick,psym=8,xrange=[minpa,maxpa],YSTYLE=1,xstyle=1,/NODATA,ytitle='Diff. PA (Deg.)',/NOERASE
   tmp=WHERE(double(fitresult) EQ 1 AND double(litpa) NE 0.)
   av1=TOTAL([double(litPA[tmp]-PAhalf[tmp])])/n_elements([double(litPA[tmp]-PAhalf[tmp])])
   error1=STDDEV([double(litPA[tmp]-PAhalf[tmp])])


   oplot,[minpa,maxpa],[av1-error1,av1-error1],thick=thick,color=grey,linestyle=2
   oplot,[minpa,maxpa],[av1,av1],thick=thick,color=grey
   oplot,[minpa,maxpa],[av1+error1,av1+error1],thick=thick,color=grey,linestyle=2

   filename=  main_dir+'Scatterliterature.txt' ;that first
   openw,1,filename

   av=TOTAL(double(litPA[tmp])-double(PAhalf[tmp]))/n_elements(tmp)
   print,"This should result in 16 galaxies"
   printf,1,'\sigma_{PA}= '+string(STDDEV(double(litPA[tmp])-double(PAhalf[tmp])),format='(F6.2)')+', av='+string(av,format='(F6.2)')+'No gal = '+string(n_elements(tmp))
   print,'test this'
   print,'\sigma_{PA}= '
   print,string(STDDEV(double(litPA[tmp])-double(PAhalf[tmp])),format='(F6.2)')
   print,', av=',av
   print,string(MEAN(double(litPA[tmp])-double(PAhalf[tmp])),format='(F6.2)')
   print,'No gal = '
   print,string(n_elements(tmp))
   print,'here we plot our info on PA literature in fat range'
   for i=0,n_elements(tmp)-1 do begin
      print,HPASSNAME[tmp[i]],double(INCL[tmp[i]]),double(litPA[tmp[i]])-double(PAhalf[tmp[i]]),double(litPA[tmp[i]]),double(PAhalf[tmp[i]])
   endfor
   tmp=WHERE(double(fitresult) EQ 1 AND double(litpa) NE 0.)
   print,tmp,fitresult,litpa
   PLOTSYM, 0 , /FILL
   ;oplot,double(INCL[tmp]),double(litPA[tmp])-double(PAhalf[tmp]),color=black,psym=8
   xerr = dblarr(n_elements(tmp))
   yerr= xerr
   fat_ploterror,double(INCL[tmp]),double(litPA[tmp])-double(PAhalf[tmp]),xerr,yerr,psym = 8, $
                     color=black,ERRCOLOR = black, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
   tmp=WHERE(double(fitresult) EQ 1.5 AND double(litpa) NE 0.)

   if tmp[0] NE -1 then begin
      print,'here we plot our info on PA literature out fat range'
      for i=0,n_elements(tmp)-1 do begin
         print,HPASSNAME[tmp[i]],double(INCL[tmp[i]]),double(litPA[tmp[i]])-double(PAhalf[tmp[i]]),double(litPA[tmp[i]]),double(PAhalf[tmp[i]])
      endfor
      PLOTSYM, 3 , /FILL
      xerr = dblarr(n_elements(tmp))
      yerr= xerr
      fat_ploterror,double(INCL[tmp]),double(litPA[tmp])-double(PAhalf[tmp]),xerr,yerr,psym = 8, $
                    color=grey,ERRCOLOR = grey, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
   ENDIF
   ;oplot,double(INCL[tmp]),double(litPA[tmp])-double(PAhalf[tmp]),color=black,psym=8,symsize=0.7
   ;Then we want where there is a RC result and a PA
    tmp2=WHERE(double(RCpahalf) NE 0. AND double(litpa) NE 0. AND double(fitresult) EQ 1)
    if tmp2[0] NE -1 then begin
       av=TOTAL(double(litPA[tmp2])-double(RCPAhalf[tmp2]))/n_elements(tmp2)
       printf,1,'\sigma_{PA-RC}= '+string(STDDEV(double(litPA[tmp2])-double(RCPAhalf[tmp2])),format='(F6.2)')+', av-RC='+string(av,format='(F6.2)')+'No gal = '+string(n_elements(tmp2))

       PLOTSYM, 0 , /FILL
       print,tmp2
       print,'Do we ahve any ble?'
       print,double(INCL[tmp2]),double(litPA[tmp2])-double(RCPAhalf[tmp2])

       oplot,double(INCL[tmp2]),double(litPA[tmp2])-double(RCPAhalf[tmp2]),color= blue, psym=8
       xerr = dblarr(n_elements(tmp2))
       yerr= xerr
       fat_ploterror,double(INCL[tmp2]),double(litPA[tmp2])-double(RCPAhalf[tmp2]),xerr,yerr,psym = 8, $
                     color=blue,ERRCOLOR = blue, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
       xxx=WHERE(ABS(double(litPA[tmp2])-double(RCPAhalf[tmp2])) GT 40)
       Print,"here we go"
    endif
    tmp2=WHERE(double(RCpahalf) NE 0. AND double(litpa) NE 0. AND double(fitresult) EQ 1.5)

    PLOTSYM, 3 , /FILL
    light_blue = '0092FF'x
    if tmp2[0] NE -1 then begin
       xerr = dblarr(n_elements(tmp2))
       yerr= xerr
       fat_ploterror,double(INCL[tmp2]),double(litPA[tmp2])-double(RCPAhalf[tmp2]),xerr,yerr,psym = 8, $
                     color=light_blue,ERRCOLOR = light_blue, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
    endif
   ; for nj=0,n_elements(xxx)-1 do begin
   ;    print,HPASSNAme[tmp2[xxx[nj]]],double(litPA[tmp2[xxx[nj]]])-double(RCPAhalf[tmp2[xxx[nj]]]),double(INCL[tmp2[xxx[nj]]]),litPA[tmp2[xxx[nj]]],double(litPA[tmp2[xxx[nj]]])-double(PAhalf[tmp2[xxx[nj]]]),double(PAhalf[tmp2[xxx[nj]]]),double(RCPAhalf[tmp2[xxx[nj]]])

   ; endfor


    ; and then DF (DFHPASSNamein NE '-1'
    print,DFHPASSNamein,litPA
    print,'Dunno'
    tmp3=WHERE(DFHPASSNamein NE '-1' AND double(litpa) NE 0.  AND double(fitresult) EQ 1)
    if tmp3[0] NE -1 then begin
       av=TOTAL(double(litPA[tmp3])-double(DFPAhalf[tmp3]))/n_elements(tmp3)
       printf,1,'\sigma_{PA-DF}= '+string(STDDEV(double(litPA[tmp3])-double(DFPAhalf[tmp3])),format='(F6.2)')+', av-DF='+string(av,format='(F6.2)')+'No gal = '+string(n_elements(tmp3))


       PLOTSYM, 0 , /FILL
       print,double(INCL[tmp3]),double(litPA[tmp3])-double(DFPAhalf[tmp3])
       xerr = dblarr(n_elements(tmp3))
       yerr= xerr
       fat_ploterror,double(INCL[tmp3]),double(litPA[tmp3])-double(DFPAhalf[tmp3]),xerr,yerr,psym = 8, $
                     color=red,ERRCOLOR = red, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
    endif
    ;oplot,double(INCL[tmp3]),double(litPA[tmp3])-double(DFPAhalf[tmp3]),color='FF0000',psym=8,thick=thick
    tmp3=WHERE(DFHPASSNamein NE '-1' AND double(litpa) NE 0.  AND double(fitresult) EQ 1.5)
    PLOTSYM, 3 , /FILL
    if tmp3[0] NE -1 then begin
       xerr = dblarr(n_elements(tmp3))
       yerr= xerr
       fat_ploterror,double(INCL[tmp3]),double(litPA[tmp3])-double(DFPAhalf[tmp3]),xerr,yerr,psym = 8, $
                     color=light_red,ERRCOLOR = light_red, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
    endif
    ;oplot,double(INCL[tmp3]),double(litPA[tmp3])-double(DFPAhalf[tmp3]),color='ffb6c1',psym=8,thick=thick,symsize=0.7
    close,1



;now the inclination

  tmp=WHERE((fitresult EQ 1 or fitresult EQ 1.5) AND litincl NE 0. and INCL GT 0 AND INCL LT 90)
  tmp2=WHERE(RCinclhalf NE 0. AND litincl NE 0. and INCL GT 0 AND INCL LT 90)
  tmp3=WHERE(DFHPASSNamein NE '-1' AND double(litincl) NE 0.  and INCL GT 0 AND INCL LT 90)
  print,n_elements(tmp),n_elements(tmp2),n_elements(tmp3)

  maxpay=MAX(double([double(litINCL[tmp]-INCLhalf[tmp]),double(litINCL[tmp2])-double(RCINCLhalf[tmp2]),double(litINCL[tmp3]-DFINCLhalf[tmp3])]),min=minpay)

  maxpa=MAX(double([(double(INCL[tmp])),double(INCL[tmp3]),double(INCL[tmp2])]),min=minpa)
  maxpay=maxpay+5.
  minpay=minpay-5
  maxpa=maxpa+5.
  minpa=minpa-5.
;  print,minpay,maxpay,minpa,maxpa
 ;  print,'we printed the pa max min which are Nan'
   plot,INCL,litINCL,position=[0.55,0.2,0.85,0.8],xtitle='INCL FAT (Deg.)',yrange=[minpay,maxpay],xthick=thick,ythick=thick,charthick=thick,thick=thick,psym=8,xrange=[minpa,maxpa],YSTYLE=1,xstyle=1,/NODATA,ytitle='Diff. INCL (Deg.)',/NOERASE
   tmp=WHERE(fitresult EQ 1 AND litincl NE 0. and INCL GT 0 AND INCL LT 90)
   tmp2=WHERE(RCinclhalf NE 0. AND fitresult EQ 1. AND litincl NE 0. and INCL GT 0 AND INCL LT 90)
   tmp3=WHERE(DFHPASSNamein NE '-1'  AND fitresult EQ 1. AND double(litincl) NE 0.  and INCL GT 0 AND INCL LT 90)

   av1=TOTAL([double(litINCL[tmp]-INCLhalf[tmp])])/n_elements([double(litINCL[tmp]-INCLhalf[tmp])])
   error1=STDDEV([double(litINCL[tmp]-INCLhalf[tmp])])
   av2=TOTAL([double(litINCL[tmp2]-RCINCLhalf[tmp2])])/n_elements([double(litINCL[tmp2]-RCINCLhalf[tmp2])])
   error2=STDDEV([double(litINCL[tmp2]-RCINCLhalf[tmp2])])
   av3=TOTAL([double(litINCL[tmp3]-DFINCLhalf[tmp3])])/n_elements([double(litINCL[tmp3]-DFINCLhalf[tmp3])])
   error3=STDDEV([double(litINCL[tmp3]-DFINCLhalf[tmp3])])



   oplot,[minpa,maxpa],[av1-error1,av1-error1],thick=thick,color=grey,linestyle=2
   oplot,[minpa,maxpa],[av1,av1],thick=thick,color=grey
   oplot,[minpa,maxpa],[av1+error1,av1+error1],thick=thick,color=grey,linestyle=2



   print,'The averages'
   print,av1,av2,av3
   print,'the sigmas'
   print,error1,error2,error3

  openu,1,filename,/append
  print,'here we plot our info on INCL literature in FAT range'
   for i=0,n_elements(tmp)-1 do begin
      print,HPASSNAME[tmp[i]],double(INCL[tmp[i]]),double(litINCL[tmp[i]])-double(INCLhalf[tmp[i]]),double(litINCL[tmp[i]]),double(INCLhalf[tmp[i]])
   endfor
   av=TOTAL(double(litINCL[tmp])-double(INCLhalf[tmp]))/n_elements(tmp)
   printf,1,'\sigma_{INCL}= '+string(STDDEV(double(litINCL[tmp])-double(INCLhalf[tmp])),format='(F6.2)')+', av='+string(av,format='(F6.2)')+'No gal = '+string(n_elements(tmp))
   av=TOTAL(double(litINCL[tmp2])-double(RCINCLhalf[tmp2]))/n_elements(tmp2)
  printf,1,'\sigma_{INCL-RC}= '+string(STDDEV(double(litINCL[tmp2])-double(RCINCLhalf[tmp2])),format='(F6.2)')+', av-RC='+string(av,format='(F6.2)')+'No gal = '+string(n_elements(tmp2))
  av=TOTAL(double(litINCL[tmp3])-double(DFINCLHALF[tmp3]))/n_elements(tmp3)
  printf,1,'\sigma_{INCL-DF}= '+string(STDDEV(double(litINCL[tmp3])-double(DFINCLHALF[tmp3])),format='(F6.2)')+', av-DF='+string(av,format='(F6.2)')+'No gal = '+string(n_elements(tmp3))
  close,1
  tmp=WHERE(fitresult EQ 1 AND litincl NE 0. and INCL GT 0 AND INCL LT 90 )
  PLOTSYM, 0 , /FILL
  ;oplot,double(INCL[tmp]),double(litINCL[tmp])-double(INCLhalf[tmp]),color=0,psym=8
  ;errplot,double(INCL[tmp]),double(litINCL[tmp])-double(INCLhalf[tmp])-litinclerr[tmp],double(litINCL[tmp])-double(INCLhalf[tmp])+litinclerr[tmp],
  ;thick=thick,color=0
  xerr = dblarr(n_elements(tmp))
  fat_ploterror,double(INCL[tmp]),double(litINCL[tmp])-double(INCLhalf[tmp]),xerr,litinclerr[tmp],psym = 8, $
                   color=black,ERRCOLOR = black, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot


  tmp=WHERE(fitresult EQ 1.5 AND litincl NE 0. and INCL GT 0 AND INCL LT 90 )
  print,'here we plot our info on INCL literature out FAT range'
   for i=0,n_elements(tmp)-1 do begin
      print,HPASSNAME[tmp[i]],double(INCL[tmp[i]]),double(litINCL[tmp[i]])-double(INCLhalf[tmp[i]]),double(litINCL[tmp[i]]),double(INCLhalf[tmp[i]])
   endfor

  PLOTSYM, 3 , /FILL
  ;oplot,double(INCL[tmp]),double(litINCL[tmp])-double(INCLhalf[tmp]),color=100,psym=8,symsize=0.7
  ;errplot,double(INCL[tmp]),double(litINCL[tmp])-double(INCLhalf[tmp])-litinclerr[tmp],double(litINCL[tmp])-double(INCLhalf[tmp])+litinclerr[tmp],thick=thick,color=100

  if tmp[0] NE -1 then begin
     xerr = dblarr(n_elements(tmp))
     fat_ploterror,double(INCL[tmp]),double(litINCL[tmp])-double(INCLhalf[tmp]),xerr,litinclerr[tmp],psym = 8, $
                   color=grey,ERRCOLOR = grey, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
  endif


  tmp2=WHERE(RCinclhalf NE 0. AND litincl NE 0. and INCL GT 0 AND INCL LT 90 AND fitresult EQ 1)

  PLOTSYM, 0 , /FILL
                                ;oplot,double(INCL[tmp2]),double(litINCL[tmp2])-double(RCINCLhalf[tmp2]),color=50,psym=8,thick=thick
                                ;errplot,double(INCL[tmp2]),double(litINCL[tmp2])-double(RCINCLhalf[tmp2])-litinclerr[tmp2],double(litINCL[tmp2])-double(RCINCLhalf[tmp2])+litinclerr[tmp2],thick=thick,color=50
  if tmp2[0] NE -1 then begin
     xerr = dblarr(n_elements(tmp2))
     fat_ploterror,double(INCL[tmp2]),double(litINCL[tmp2])-double(RCINCLhalf[tmp2]),xerr,litinclerr[tmp2],psym = 8, $
                   color=blue,ERRCOLOR = blue, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
  endif

  tmp2=WHERE(RCinclhalf NE 0. AND litincl NE 0. and INCL GT 0 AND INCL LT 90 AND fitresult EQ 1.5)

  PLOTSYM, 3 , /FILL
                                ;oplot,double(INCL[tmp2]),double(litINCL[tmp2])-double(RCINCLhalf[tmp2]),color=65,psym=8,thick=thick,symsize=0.7
                                ;errplot,double(INCL[tmp2]),double(litINCL[tmp2])-double(RCINCLhalf[tmp2])-litinclerr[tmp2],double(litINCL[tmp2])-double(RCINCLhalf[tmp2])+litinclerr[tmp2],thick=thick,color=65
  if tmp2[0] NE -1 then begin 
     xerr = dblarr(n_elements(tmp2))
     fat_ploterror,double(INCL[tmp2]),double(litINCL[tmp2])-double(RCINCLhalf[tmp2]),xerr,litinclerr[tmp2],psym = 8, $
                   color=light_blue,ERRCOLOR = light_blue, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
  endif

  tmp3=WHERE(DFHPASSNamein NE '-1' AND double(litincl) NE 0.  and INCL GT 0 AND INCL LT 90 AND fitresult EQ 1)
  PLOTSYM, 0 , /FILL
                                ;oplot,double(INCL[tmp3]),double(litINCL[tmp3])-double(DFINCLhalf[tmp3]),color=254,psym=8,thick=thick
                                ;errplot,double(INCL[tmp3]),double(litINCL[tmp3])-double(DFINCLhalf[tmp3])-litinclerr[tmp3],double(litINCL[tmp3])-double(DFINCLhalf[tmp3])+litinclerr[tmp3],thick=thick,color=254
  if tmp3[0] NE -1 then begin 
     xerr = dblarr(n_elements(tmp3))
     fat_ploterror,double(INCL[tmp3]),double(litINCL[tmp3])-double(DFINCLhalf[tmp3]),xerr,litinclerr[tmp3],psym = 8, $
                   color=red,ERRCOLOR = red ,ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot

  endif
  tmp3=WHERE(DFHPASSNamein NE '-1' AND double(litincl) NE 0.  and INCL GT 0 AND INCL LT 90 AND fitresult EQ 1.5)
  PLOTSYM, 3 , /FILL
                                ;oplot,double(INCL[tmp3]),double(litINCL[tmp3])-double(DFINCLhalf[tmp3]),color=240,psym=8,thick=thick,symsize=0.7
                                ;errplot,double(INCL[tmp3]),double(litINCL[tmp3])-double(DFINCLhalf[tmp3])-litinclerr[tmp3],double(litINCL[tmp3])-double(DFINCLhalf[tmp3])+litinclerr[tmp3],thick=thick,color=240
  if tmp3[0] NE -1 then begin 
     xerr = dblarr(n_elements(tmp3))
     fat_ploterror,double(INCL[tmp3]),double(litINCL[tmp3])-double(DFINCLhalf[tmp3]),xerr,litinclerr[tmp3],psym = 8, $
                   color=light_red,ERRCOLOR =light_red, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
   endif
   tmpx=WHERE(INCL GT 0. AND INCL LT 40.)
   print,'these galxies ar above 70 degrees.'
   if tmpx[0] NE -1 then begin 
      for i=0,n_elements(tmpx)-1 do begin
         print,STRJOIN([JNAME[tmpx[i]],JNumber[tmpx[i]],string(INCL[tmpx[i]]),string(INCLhalf[tmpx[i]]),string(litINCL[tmpx[i]]),string(RCinclhalf[tmpx[i]]),string(DFINCLhalf[tmpx[i]]),string(DFRCinclrmsav[tmpx[i]]),string(RCinclrmsav[tmpx[i]]),string(DFinclrmsav[tmpx[i]]),string(fitresult[tmpx[i]])],'**')
      endfor
   endif
   tmp=WHERE(fitresult EQ 1 AND litincl NE 0. and INCL GT 0. AND INCL LT 40.)
   tmp2=WHERE(RCinclhalf NE 0. AND fitresult EQ 1. AND litincl NE 0. and INCL GT 0. AND INCL LT 40.)
   tmp3=WHERE(DFHPASSNamein NE '-1'  AND fitresult EQ 1. AND double(litincl) NE 0.  and INCL GT 0. AND INCL LT 40.)

   av1=TOTAL([double(litINCL[tmp]-INCLhalf[tmp])])/n_elements([double(litINCL[tmp]-INCLhalf[tmp])])
   error1=STDDEV([double(litINCL[tmp]-INCLhalf[tmp])])
   print,'The tirific average offset and scatter',av1,error1,n_elements(tmp)
   av2=TOTAL([double(litINCL[tmp2]-RCINCLhalf[tmp2])])/n_elements([double(litINCL[tmp2]-RCINCLhalf[tmp2])])
   error2=STDDEV([double(litINCL[tmp2]-RCINCLhalf[tmp2])])
   print,'The RC average offset and scatter',av2,error2,n_elements(tmp2)
   av3=TOTAL([double(litINCL[tmp3]-DFINCLhalf[tmp3])])/n_elements([double(litINCL[tmp3]-DFINCLhalf[tmp3])])
   error3=STDDEV([double(litINCL[tmp3]-DFINCLhalf[tmp3])])
   print,'The DF average offset and scatter',av3,error3,n_elements(tmp3)

 DEVICE,/CLOSE
  convert_ps,main_dir+'Overviewlit.ps',/trim,/png,/delete
  print,'What happended to my plots'
 
;Some alternative plots
  !p.charsize=1.1
  !p.thick=thick

  openw,1,main_dir+'SingleTirIncl.txt'
  for i=0,n_elements(INCL)-1 do begin
     printf,1,format='(A10,F10.3)',HPASSNAME[i],INCL[i]
  endfor
  close,1

   SET_PLOT, 'PS'
   DEVICE, FILENAME=main_dir+'Overviewmean.ps',/color,/PORTRAIT,/ENCAPSULATED,xsize=21,ysize=21*0.175/0.125,/DECOMPOSED
        ;Next up is the PA which we can plot
                                ;either against radius but probably
                                ;inclination is mor usefull so we do
   openw,1,main_dir+'Proc_results.txt'
   printf,1,format='(A20,2F10.3)','Name','AC1','AC2'

   for i=0,n_elements(fitresult)-1 do begin
      printf,1,format='(A20,2F10.3)',HPASSNAME[i],fitresultsmall[i],fitresult[i]
   endfor
   close,1
   filename=  main_dir+'Scattermean.txt'                           ;that first
   openw,1,filename
   xpos1=[0.1,0.45]
   xpos2=[0.55,0.9]

   ypos1=[0.795,0.92]
   ypos2=[0.59,0.715]
   ypos3=[0.385,0.51]
   ypos4=[0.18,0.305]


   openw,88,main_dir+'averageandsigma.txt'
   printf,88,format='(2A10)','Average','Sigma'
   close,88

   tmp=WHERE(fitresult EQ 1 OR fitresult EQ 1.5)
   printf,1,'Out of '+strtrim(string(n_elements(RCRadius)),2)+' Galaxies '+strtrim(string(n_elements(tmp)),2)+' were accepted'
   printf,1,'Out of '+strtrim(string(n_elements(RCRadius)),2)+' Galaxies '+strtrim(string(n_elements(RCRadius)-n_elements(tmp)),2)+' galaxies failed'
   close,1
   openw,1,filename+'_outliers.txt'
   printf,1,'This file contains the individual galaxies that lie outside the 3 sigma from the average'
   close,1
   ;INCL against inclination
   Errorsinput=[[DFRCinclrms],[RCinclrms],[DFinclrms]]
   PAinput=[[DFRCinclrmsav],[RCinclrmsav],[DFinclrmsav]]
   print,'Here we printhe inclination and their and their diffs'
   for j=0,n_elements(PAinput[*,0])-1 do begin
      print,HPASSNAMe[j],PAinput[j,0],Errorsinput[j,0],PAinput[j,1],Errorsinput[j,1],PAinput[j,2],Errorsinput[j,2],INCL[j],fitresult[j]
   endfor

   plotdifflinev3,fitresult,PAinput,INCL,position=[xpos2[0],ypos2[0],xpos2[1],ypos2[1]],xtitle='INCL FAT (Deg.)',ytitle='Diff INCL (Deg.)',sigmanames='Incl',filename=filename,errors=Errorsinput,anno='d',name=HPASSNAME

 ;Vro against inlincation
   Errorsinput=[[DFRCvrotrms/channelwidth],[RCvrotrms/channelwidth],[DFvrotrms/channelwidth]]
   PAinput=[[DFRCvrotrmsav/channelwidth],[RCvrotrmsav/channelwidth],[DFvrotrmsav/channelwidth]]
   print,'Here we print the vrot diffrences  inclination'
   for j=0,n_elements(PAinput[*,0])-1 do begin
      print,HPASSNAMe[j],PAinput[j,0],PAinput[j,1],PAinput[j,2],INCL[j],fitresult[j]
   endfor
   plotdifflinev3,fitresult,PAinput,INCL,position=[xpos1[0],ypos3[0],xpos1[1],ypos3[1]],xtitle='INCL FAT (deg.)',ytitle='Diff V_rot (channels)',sigmanames='Vrot',/NOERASE,filename=filename,errors=Errorsinput,anno='e',name=HPASSNAME

 ;Vobs against inlincation
   Errorsinput=[[DFRCvobsrms/channelwidth],[RCvobsrms/channelwidth],[DFvobsrms/channelwidth]]
   PAinput=[[DFRCvobsrmsav/channelwidth],[RCvobsrmsav/channelwidth],[DFvobsrmsav/channelwidth]]
;   print,HPASSNAME[WHERE(DFvobsrmsmean/channelwidth GT 6)]

  ; print,'Here we print the vobs diffrences  inclination'
 ;  for j=0,n_elements(PAinput[*,0])-1 do begin
 ;     print,HPASSNAMe[j],PAinput[j,0],PAinput[j,1],PAinput[j,2],INCL[j],fitresult[j]
 ;  endfor
   plotdifflinev3,fitresult,PAinput,INCL,position=[xpos2[0],ypos3[0],xpos2[1],ypos3[1]],xtitle='INCL_{FAT} (Deg.)',ytitle='Diff. V_{los} (channels)',sigmanames='Vobs',/NOERASE,filename=filename,errors=Errorsinput,anno='f',name=HPASSNAME

   notsitmp=WHERE(fitresult EQ 1)
   print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
   tmp=WHERE(fitresult EQ 1)
   check1=WHERE(double(PAinput[tmp,0]) NE 0.)
   check2=WHERE(double(PAinput[tmp,1]) NE 0.)
   check3=WHERE(double(PAinput[tmp,2]) NE 0.)
   print,'the vobs means are'
   print,'RC-DF'
   print,mean(double(DFRCvobsrmsmean[check1]))
    print,'RC-TF'
   print,mean(double(RCvobsrmsmean[check2]))
   print,'DF-TF'
   print,mean(double(DFvobsrmsmean[check3]))
    print,'!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print,DFvobsrmsmean[check3],RCvobsrmsmean[check2],DFRCvobsrmsmean[check1]
    print,'Positions angle against inclination'
   Errorsinput=[[DFRCparms],[RCparms],[DFparms]]
   PAinput=[[DFRCparmsav],[RCparmsav],[DFparmsav]]
   plotdifflinev3,fitresult,PAinput,INCL,position=[xpos1[0],ypos2[0],xpos1[1],ypos2[1]],xtitle='INCL FAT (Deg.)',ytitle='Diff PA (Deg.)',sigmanames='PA',/NOERASE,filename=filename,errors=Errorsinput,anno='c',name=HPASSNAME
;  print,'Here we printhe position angles diffrences and the inclination'
;   for j=0,n_elements(PAinput[*,0])-1 do begin
 ;     print,HPASSNAMe[j],PAinput[j,0],PAinput[j,1],PAinput[j,2],INCL[j],fitresult[j]
 ;  endfor

   check1=WHERE(double(RCRAdeg) EQ 0.)
   check2=WHERE(double(DFRAdeg) EQ 0.)
   check3=WHERE(double(RAdeg) EQ 0.)
   ;The RA against DEC
   RAinput=[[(double(RCRAdeg)-double(DFRAdeg))*(3600./beam)*COS(RCDECdeg*!DtoR)],[(double(RCRAdeg)-double(RAdeg))*(3600./beam)*COS(DECdeg*!DtoR)],[(double(DFRAdeg)-double(RAdeg))*(3600./beam)*COS(DECdeg*!DtoR)]]
   DECinput=[[(double(RCDECdeg)-double(DFDECdeg))*3600./beam],[(double(RCDECdeg)-double(DECdeg))*3600./beam],[(double(DFDECdeg)-double(DECdeg))*3600./beam]]
   IF check1[0] NE -1 then RAinput[check1,0]=0.
   IF check2[0] NE -1 then RAinput[check2,0]=0.
   IF check2[0] NE -1 then RAinput[check2,2]=0.
   IF check3[0] NE -1 then RAinput[check3,2]=0.
   IF check1[0] NE -1 then RAinput[check1,1]=0.
   IF check3[0] NE -1 then RAinput[check3,1]=0.
   IF check1[0] NE -1 then DECinput[check1,0]=0.
   IF check2[0] NE -1 then DECinput[check2,0]=0.
   IF check2[0] NE -1 then DECinput[check2,2]=0.
   IF check3[0] NE -1 then DECinput[check3,2]=0.
   IF check1[0] NE -1 then DECinput[check1,1]=0.
   IF check3[0] NE -1 then DECinput[check3,1]=0.
 ;  Print,'Yea yea'
  PRint,RAinput
   print,'What is the problem'
 ;  print,(double(RCDECdeg)-double(DFDECdeg))*3600.
   plotdiffcircv3,fitresult,RAinput,DECinput,position=[xpos2[0],ypos1[0],xpos2[1],ypos1[1]],xtitle='Diff RA (beams)',ytitle='Diff. DEC (beams)',sigmanames=['RA','DEC'],filename=filename,xticks=2,yticks=2,xminor=5,yminor=5,anno='b',/NOERASE,xrange=[-3.0,3.0],yrange=[-3.0,3.0],xtickv=[-1.5,0,1.5],ytickv=[-1.5,0,1.5],name=HPASSNAME,radius=[0.5,0.5]
   ;The vsys against the squared error in the position
   tmp=DECinput
   DECinput=dblarr(n_elements(RAinput[*,1]))
   DECinput=SQRT(RAinput[*,1]^2+tmp[*,1]^2)
   tmp=WHERE(ABS(DECinput) GT 30)
 ;  print,'Well here we go'
 ;  help,DECinput
 ;  for h=0,n_elements(tmp)-1 do begin
  ;    print,tmp[h],hpassname[tmp[h]],DECinput[tmp[h]]
  ; endfor
   RAinput=[[(double(RCvsys)-double(DFvsys))/channelwidth],[(double(RCvsys)-double(vsys))/channelwidth],[(double(DFvsys)-double(vsys))/channelwidth]]

   IF check1[0] NE -1 then RAinput[check1,0]=0.
   IF check2[0] NE -1 then RAinput[check2,0]=0.
   IF check2[0] NE -1 then RAinput[check2,2]=0.
   IF check3[0] NE -1 then RAinput[check3,2]=0.
   IF check1[0] NE -1 then RAinput[check1,1]=0.
   IF check3[0] NE -1 then RAinput[check3,1]=0.
   plotdifflinev3,fitresult,RAinput,DECinput,position=[xpos1[0],ypos1[0],xpos1[1],ypos1[1]],xtitle='Diff. Central (beams)',ytitle='Diff. Vsys (channels)',sigmanames='vsys',/NOERASE,filename=filename,anno='a',xticks=2,xminor=5,yminor=5,name=HPASSNAME,xtickv=[0.5,1.5,2.5]

   ;and lastly we plot the vrot differece against the total size
;   Errorsinput=[[DFRCvrotrms/channelwidth],[RCvrotrms/channelwidth],[DFvrotrms/channelwidth]]
;   PAinput=[[DFRCvrotrmsav/channelwidth],[RCvrotrmsav/channelwidth],[DFvrotrmsav/channelwidth]]

  Errorsinput=[[DFRCinclrms],[RCinclrms],[DFinclrms]]
   PAinput=[[DFRCinclrmsav],[RCinclrmsav],[DFinclrmsav]]

   IF check1[0] NE -1 then PAinput[check1,0]=0.
   IF check2[0] NE -1 then PAinput[check2,0]=0.
   IF check2[0] NE -1 then PAinput[check2,2]=0.
   IF check3[0] NE -1 then PAinput[check3,2]=0.
   IF check1[0] NE -1 then PAinput[check1,1]=0.
   IF check3[0] NE -1 then PAinput[check3,1]=0.

   inputradius=(radius)/beam
   tmp=WHERE(fitresult EQ 1)
 ;  print,inputradius[tmp],HPASSNAME[tmp],radius[tmp],beam[tmp]
 ;  print,"Man o man"
   plotdifflinev3,fitresult,PAinput,inputradius*2.,position=[xpos1[0],ypos4[0],xpos1[1],ypos4[1]],xtitle='Diameter FAT (beams)',ytitle='Diff. INCL (Deg.)',sigmanames='Rad',/NOERASE,filename=filename,anno='g',errors=Errorsinput,name=HPASSNAME
 ;and the difference in radius as a function of diameter
   PAinput=[[(double(RCradius)-double(DFradius))/beam],[(double(RCradius)-double(radius))/beam],[(double(DFradius)-double(radius))/beam]]
   tmpx=WHERE(inputradius LT 5.)
   print,'These galaxies are too small'
   for i=0,n_elements(tmpx)-1 do begin
      print,STRJOIN([JNAME[tmpx[i]],JNumber[tmpx[i]]],' '),inputradius[tmpx[i]]
   endfor
   IF check1[0] NE -1 then PAinput[check1,0]=0.
   IF check2[0] NE -1 then PAinput[check2,0]=0.
   IF check2[0] NE -1 then PAinput[check2,2]=0.
   IF check3[0] NE -1 then PAinput[check3,2]=0.
   IF check1[0] NE -1 then PAinput[check1,1]=0.
   IF check3[0] NE -1 then PAinput[check3,1]=0.
  plotdifflinev3,fitresult,PAinput*2.,inputradius*2.,position=[xpos2[0],ypos4[0],xpos2[1],ypos4[1]],xtitle='Diameter FAT (beams)',ytitle='Diff Diameter (beams)',sigmanames='Dia',/NOERASE,filename=filename,anno='h',name=HPASSNAME



   DEVICE,/CLOSE
   convert_ps,main_dir+'Overviewmean.ps',/trim,/png

   tmp=WHERE(fitresult EQ 1 AND Z0 NE 0.)
   maxpay=MAX(double(Z0[tmp]),min=minpay)

   maxpa=MAX(double(INCL[tmp]),min=minpa)
   maxpay=maxpay+0.1
   minpay=minpay-0.1
   maxpa=maxpa+5.
   minpa=minpa-5.
   SET_PLOT, 'PS'
   DEVICE, FILENAME=main_dir+'OverviewdispZ0.ps',/color,/PORTRAIT,/ENCAPSULATED,xsize=21,ysize=10.5,/DECOMPOSED
   plot,INCL[tmp],Z0[tmp],position=[0.1,0.2,0.4,0.8],xtitle='INCL FAT (Deg.)',yrange=[minpay,maxpay],xthick=thick,ythick=thick,charthick=thick,thick=thick,psym=8,xrange=[minpa,maxpa],YSTYLE=1,xstyle=1,ytitle='Scale height (kpc)'

   avZ0=TOTAL(Z0[tmp])/n_elements(tmp)
   oplot,[minpa,maxpa],[avZ0-1*STDDEV(double(z0[tmp])),avZ0-1*STDDEV(double(z0[tmp]))],thick=thick,color=grey,linestyle=2
   oplot,[minpa,maxpa],[avZ0,avZ0],thick=thick,color=grey
   oplot,[minpa,maxpa],[avZ0+STDDEV(double(z0[tmp])),avZ0+STDDEV(double(z0[tmp]))],thick=thick,color=grey,linestyle=2

   tmp=WHERE(fitresult EQ 1)
   maxpay=MAX(double(SDIS[tmp]),min=minpay)

   maxpa=MAX(double(INCL[tmp]),min=minpa)
   maxpay=maxpay+1
   minpay=minpay-1
   maxpa=maxpa+5.
   minpa=minpa-5.
   avSDIS=TOTAL(SDIS[tmp])/n_elements(tmp)
   plot,INCL[tmp],SDIS[tmp],position=[0.55,0.2,0.85,0.8],xtitle='INCL FAT  (Deg.)',yrange=[minpay,maxpay],xthick=thick,ythick=thick,charthick=thick,thick=thick,psym=8,xrange=[minpa,maxpa],YSTYLE=1,xstyle=1,ytitle='Dispersion (km s!E-1!N)',/NOERASE


   oplot,[minpa,maxpa],[avSDIS-1*STDDEV(double(SDIS[tmp])),avSDIS-1*STDDEV(double(SDIS[tmp]))],thick=thick,color=grey,linestyle=2
   oplot,[minpa,maxpa],[avSDIS,avSDIS],thick=thick,color=grey
   oplot,[minpa,maxpa],[avSDIS+STDDEV(double(SDIS[tmp])),avSDIS+STDDEV(double(SDIS[tmp]))],thick=thick,color=grey,linestyle=2


   DEVICE,/CLOSE
   convert_ps,main_dir+'OverviewdispZ0.ps',/trim,/png
   tmp=WHERE(SDIS GT 20)
   if tmp[0] NE -1 then print,HPASSNAME[tmp]



endif
end
