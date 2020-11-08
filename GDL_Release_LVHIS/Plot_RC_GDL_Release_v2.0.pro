Pro Plotallrc
  compile_opt idl2

  main_dir='LVHIS-26_3/'


  catname=main_dir+'Galaxies.txt'
;catname='/Users/kam036/WALLABY/LVHIS-26/LVHIS-26_In Article/LVHIS_3_results.txt'
;This file is produced by compare.
readcol,main_dir+'FittedAverages2.txt',format='(A,F)',HPASSName,PA,RCPA,DFPA,INCL,RCINCL,DFINCL,RAdeg,RCRAdeg,DFRAdeg,DECdeg,RCDECdeg,DFDECdeg,VSYS,RCVSYS,DFVSYS,VROT,RCVROT,DFVROT,skipline=1
first=2
alphabet=string(bindgen(1,26)+(byte('a'))[0])
black = '000000'x
grey = '808080'x
red = '0000FF'x
light_blue = '7CACFF'x
blue = 'FF0000'x
light_red = 'FFCCCC'x
 PLOTSYM, 0 , /FILL
print,INCL

spacecheck=strtrim(str_sep(main_dir,' '),2)
main_dirsl=STRJOIN(spacecheck,'\ ')
print,main_dir
acount=0
RESOLVE_ROUTINE,'gettirific'
RESOLVE_ROUTINE,'getrc'
RESOLVE_ROUTINE,'getdiskfit'
RESOLVE_ROUTINE,'interpolate'
RESOLVE_ROUTINE,'fat_ploterror'
RESOLVE_ROUTINE,'convert_ps'
readcol,main_dirsl+'Gal_Names_Sub.txt',format='(A,A,A,A)',HPASSNameDF,JNameDF,JNumberDF,OthernameDF,skipline=1
readcol,main_dirsl+'Galaxy_Names.txt',format='(A,A,A,A)',HPASSNameAll,JName,JNumber,Othername,skipline=1

for i=0,n_elements(Jnumber)-1 do begin
   tmp=strsplit(JNumber[i],'--',/extract,/regex)
   Jnumber[i]=STRJOIN(tmp,'-')
endfor
h=' '
openr,1,catname
readf,1,h
fitresult=dblarr(n_elements(HPASSName))

WHILE ~EOF(1) do begin
   readf,1,h
   values=str_sep(strtrim(strcompress(h)),' ')
   tmp=WHERE(values[0] EQ HPASSName)
;   print,HPASSName[tmp],values[1],values[2]
   case 1 OF
        float(values[1]) EQ 1 and float(values[2]) EQ 1: fitresult[tmp] = 1
        float(values[1]) EQ 1 and float(values[2]) EQ 0: fitresult[tmp] = 1.5
        float(values[1]) EQ 0 and float(values[2]) EQ 0: fitresult[tmp] = 0;
        else: fitresult[tmp] = 0;
   ENDCASE
ENDWHILE
close,1

fitresult=fitresult[SORT(INCL)]
Othername=Othername[SORT(INCL)]
JName=JName[SORT(INCL)]
JNumber=Jnumber[SORT(INCL)]
HPASSName=HPASSName[SORT(INCL)]
print,HPASSName,INCL,SORT(INCL)
INCL=INCL[SORT(INCL)]

;Ok so first we collect the Tirific info from Basic
HPASSNAMEin=strarr(n_elements(HPASSname))
beam=dblarr(n_elements(HPASSname))
Tirparameters2=['RADI','INCL','INCL_ERR','INCL_2','INCL_2_ERR','VROT','VROT_ERR','SBR','SBR_2']
Tirific=dblarr(2,n_elements(Tirparameters2),n_elements(HPASSNAME))
;next up are our rotcur solutions
RCHPASSNAMEin=strarr(n_elements(HPASSname))
SHOrotcurdir=main_dirsl+'rotcur_outputs_LVHIS_subsample'
RCParameters=['RADI','INCL','INCL_ERR','VROT','VROT_ERR']
Rotcur=dblarr(2,n_elements(RCparameters),n_elements(HPASSNAME))
;next up are our diskfit solutions
DFHPASSNAMEin=strarr(n_elements(HPASSname))
DFrotcurdir=main_dirsl+'Diskfit/'
DFParameters=['RADI','INCL','INCL_ERR','VROT','VROT_ERR']
Diskfit=dblarr(2,n_elements(DFparameters),n_elements(HPASSNAME))
for i=0,n_elements(HPASSname)-1 do begin
   ;first we'll read in the parameters from the 1stfit.
   print,'this is the galaxy we are going to look at'
   print,HPASSName[i]
   IF FILE_TEST(main_dir+'/'+HPASSname[i]+'/Finalmodel/Finalmodel.def') then begin
      gettirific,main_dir+'/'+HPASSname[i]+'/Finalmodel/Finalmodel.def',Tirparameters2,Tirresult2,/errors
      tmp=WHERE(Tirparameters2 EQ 'SBR_2')

      HPASSNAMEin[i]=HPASSname[i]
      IF n_elements(Tirific[*,0,0]) LT n_elements(Tirresult2[*,0]) then begin
         tmp=Tirific
         Tirific=dblarr(n_elements(Tirresult2[*,0]),n_elements(Tirparameters2),n_elements(HPASSNAME))
         Tirific[0:n_elements(tmp[*,0,0])-1,*,*]=tmp[0:n_elements(tmp[*,0,0])-1,*,*]
         Tirific[*,*,i]=Tirresult2[*,*]
      ENDIF ELSE BEGIN
         Tirific[0:n_elements(Tirresult2[*,0])-1,*,i]=Tirresult2[*,*]
      ENDELSE
      tmp=WHERE(Tirparameters2 EQ 'SBR')


      tmp2=WHERE(Tirific[0:n_elements(Tirresult2[*,0])-1,tmp,i] LT 1E-15)
      print,tmp2,Tirific[0:n_elements(Tirresult2[*,0])-1,tmp,i]
      tmp=WHERE(Tirparameters2 EQ 'SBR_2')
      tmp3=WHERE(Tirific[0:n_elements(Tirresult2[*,0])-1,tmp,i] LT 1E-15)
      print,tmp3,Tirific[0:n_elements(Tirresult2[*,0])-1,tmp,i],Tirresult2[*,tmp]

      IF tmp2[0] NE -1 and tmp3[0] NE -1 then MAXin=MAX([tmp2,tmp3]) else maxin=n_elements(Tirresult2[*,0])
      print,MAXin, n_elements(Tirresult2[*,0])-1
      print,'whatever'
      h=' '
    ;  read,h

      IF MAXin LE n_elements(Tirresult2[*,0])-1  AND maxin GT -1 then begin
        ; print,'we"re here'
        ; print,maxin,tmp2,tmp3,Tirific[0:n_elements(Tirresult2[*,0])-1,tmp,i]
        ; wait,10.
         for j=Maxin,n_elements(Tirresult2[*,0])-1 do begin
            Tirific[j,*,i]=0.
         endfor
      endif
   ENDIF ELSE begin
      HPASSNAMEin[i]='-1'
   ENDELSE
   name=StrJoin(StrSplit(othername[i], '-', /Regex, /Extract, $
                          /Preserve_Null), '_')
   IF STRMID(name,0 , 1 ) NE 'e' then name=STRUPCASE(name)
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
      RCHPASSNAMEin[i]=HPASSname[i]
      IF n_elements(Rotcur[*,0,0]) LT n_elements(RCResult[*,0]) then begin
         tmp=Rotcur
         Rotcur=dblarr(n_elements(RCResult[*,0]),n_elements(RCparameters),n_elements(HPASSNAME))
         help,Rotcur,tmp
         print,n_elements(tmp[*,0,0])-1
         Rotcur[0:n_elements(tmp[*,0,0])-1,*,*]=tmp[0:n_elements(tmp[*,0,0])-1,*,*]
         Rotcur[*,*,i]=RCResult[*,*]
      ENDIF ELSE BEGIN
         Rotcur[0:n_elements(RCResult[*,0])-1,*,i]=RCResult[*,*]
      ENDELSE
   ENDIF ELSE BEGIN
      RCHPASSNAMEin[i]='-1'
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
   IF FILE_TEST(DFrotcurdir+name) AND checkDF[0] NE -1 then begin
     ; print,name
      CD,DFrotcurdir+name,current=old_dir
      spawn,'ls *.out',output
      if output then begin
         getdiskfit,output,DFParameters,DFresult,main_dir+'/'+HPASSname[i]+'/Cube.fits'
         DFHPASSNAMEin[i]=HPASSname[i]
         IF n_elements(Diskfit[*,0,0]) LT n_elements(DFResult[*,0]) then begin
            tmp=Diskfit
            Diskfit=dblarr(n_elements(DFResult[*,0]),n_elements(DFparameters),n_elements(HPASSNAME))
            Diskfit[0:n_elements(tmp[*,0,0])-1,*,*]=tmp[0:n_elements(tmp[*,0,0])-1,*,*]
            Diskfit[*,*,i]=DFResult[*,*]
         ENDIF ELSE BEGIN
            Diskfit[0:n_elements(DFResult[*,0])-1,*,i]=DFResult[*,*]
         ENDELSE
      ENDIF ELSE begin
         DFHPASSNAMEin[i]='-1'
      ENDELSE
      CD,old_dir
   ENDIF ELSE begin
      DFHPASSNamein[i]='-1'
   ENDELSE
ENDFOR

;Now first let's make a simple plot where all Rotationcurves of
;one galaxy get the same color but a differen linestyle
count=0
for i=0,n_elements(HPASSNAme)-1 do begin
   IF TOTAL(Tirific[*,5,i]) NE 0. then count++
   print,HPASSName[i],fitresult[i]
endfor

;colorjump=254/n_elements(HPASSNAME)

xsize=32.
ysize=xsize/5.*ceil(count/5.)+(xsize/2.)*0.15
;ysize= 32.
ixsize=0.16
iysize=0.8/ceil(count/5.)
!p.font=1
!p.charsize=1.
!p.thick=3.
!P.FONT=1
!p.background=0
A = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a unit circle with 16 points,
; and set the filled flag:
USERSYM, COS(A), SIN(A), /FILL
SET_PLOT, 'PS'

DEVICE, FILENAME=main_dirsl+'AllRC.ps',/color,/PORTRAIT,/ENCAPSULATED,xsize=xsize,ysize=ysize,/DECOMPOSED
maxRAdius=MAX([Tirific[*,0,*],Rotcur[*,0,*],Diskfit[*,0,*]])
MaxVel=MAX([Tirific[*,5,*],Rotcur[*,3,*],Diskfit[*,3,*]])
Maxvel=210.
MaxRadius=500.
;loadct,40,/SILENT
;print,Tirific[*,6,*]
;plot,Tirific[*,0,0],Tirific[*,5,0],position=[0.1,0.1,0.9,0.9],xtitle=textoidl('Radius (arcsec)'),yrange=[0,maxVel],xrange=[0,maxRadius],xthick=!p.thick,ythick=!p.thick,charthick=!p.thick,thick=!p.thick,YSTYLE=1,xstyle=1,/NODATA,ytitle=textoidl('Velocity (km s^{-1})')
xnum=1.
ynum=0.
print,HPASSNAme,fitresult

for i=0,n_elements(HPASSNAme)-1 do begin
   IF TOTAL(Tirific[*,5,i]) NE 0.  AND fitresult[i] EQ 1 OR fitresult[i] EQ 1.5 then begin
      maxRAdius=MAX([Tirific[*,0,i],Rotcur[*,0,i],Diskfit[*,0,i]])
      MaxVel=MAX([Tirific[*,5,i],Rotcur[*,3,i],Diskfit[*,3,i]])+10.
      position=[0.1+xnum*ixsize,0.1+ynum*iysize,0.1+(xnum+1)*ixsize-0.02,0.1+(ynum+1)*iysize-0.1/ceil(count/5.)]
      position=[0.9-(xnum+1)*ixsize+0.02,0.1+ynum*iysize,0.9-xnum*ixsize,0.1+(ynum+1)*iysize-0.15/ceil(count/5.)]
      print,position
      PLOTSYM, 0 , /FILL
      tmp=WHERE(Tirific[*,5,i] NE 0.)


      IF xnum EQ 1 and ynum EQ 0 then begin
         plot,Tirific[0:tmp[n_elements(tmp)-1],0,i],Tirific[0:tmp[n_elements(tmp)-1],5,i],position=position,yrange=[0,maxVel],xrange=[0,maxRadius],xthick=!p.thick,ythick=!p.thick,charthick=!p.thick,thick=!p.thick,YSTYLE=1,xstyle=1
      ENDIF ELSE BEGIN
         plot,Tirific[0:tmp[n_elements(tmp)-1],0,i],Tirific[0:tmp[n_elements(tmp)-1],5,i],position=position,yrange=[0,maxVel],xrange=[0,maxRadius],xthick=!p.thick,ythick=!p.thick,charthick=!p.thick,thick=!p.thick,YSTYLE=1,xstyle=1,/NOERASE
      ENDELSE
      ssize=0.75
      xerr= dblarr(n_elements(Tirific[0:tmp[n_elements(tmp)-1],0,i]))
      fat_ploterror,Tirific[0:tmp[n_elements(tmp)-1],0,i],Tirific[0:tmp[n_elements(tmp)-1],5,i],xerr,Tirific[0:tmp[n_elements(tmp)-1],6,i],psym = 8, $
                     color=black,ERRCOLOR = black, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
      ;errplot,Tirific[0:tmp[n_elements(tmp)-1],0,i],Tirific[0:tmp[n_elements(tmp)-1],5,i]-Tirific[0:tmp[n_elements(tmp)-1],6,i],Tirific[0:tmp[n_elements(tmp)-1],5,i]+Tirific[0:tmp[n_elements(tmp)-1],6,i]
      IF TOTAL(Rotcur[*,3,i]) NE 0. then begin
         tmp=WHERE(Rotcur[*,3,i] NE 0.)
         print,blue
         oplot,Rotcur[0:tmp[n_elements(tmp)-1],0,i],Rotcur[0:tmp[n_elements(tmp)-1],3,i], color=blue,linestyle=2
         xerr= dblarr(n_elements(Rotcur[0:tmp[n_elements(tmp)-1],0,i]))
         fat_ploterror,Rotcur[0:tmp[n_elements(tmp)-1],0,i],Rotcur[0:tmp[n_elements(tmp)-1],3,i],xerr,Rotcur[0:tmp[n_elements(tmp)-1],4,i],psym = 8, $
                     color=blue,ERRCOLOR = blue, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
         ;errplot,Rotcur[0:tmp[n_elements(tmp)-1],0,i],Rotcur[0:tmp[n_elements(tmp)-1],3,i]-Rotcur[0:tmp[n_elements(tmp)-1],4,i],Rotcur[0:tmp[n_elements(tmp)-1],3,i]+Rotcur[0:tmp[n_elements(tmp)-1],4,i],color=50
      ENDIF
      IF TOTAL(Diskfit[*,3,i]) NE 0. then begin
         tmp=WHERE(Diskfit[*,3,i] NE 0.)
         oplot,Diskfit[0:tmp[n_elements(tmp)-1],0,i],Diskfit[0:tmp[n_elements(tmp)-1],3,i], color=red,linestyle=3
         xerr= dblarr(n_elements(Diskfit[0:tmp[n_elements(tmp)-1],0,i]))
         fat_ploterror,Diskfit[0:tmp[n_elements(tmp)-1],0,i],Diskfit[0:tmp[n_elements(tmp)-1],3,i],xerr,Diskfit[0:tmp[n_elements(tmp)-1],4,i],psym = 8, $
                     color=red,ERRCOLOR = red, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot

         ;errplot,Diskfit[0:tmp[n_elements(tmp)-1],0,i],Diskfit[0:tmp[n_elements(tmp)-1],3,i]-Diskfit[0:tmp[n_elements(tmp)-1],4,i],Diskfit[0:tmp[n_elements(tmp)-1],3,i]+Diskfit[0:tmp[n_elements(tmp)-1],4,i],color=254
      ENDIF
  ;    XYOUTS,0.1+(xnum+1)*ixsize-0.025,0.1+(ynum)*iysize+0.02*iysize,JName[i]+' '+JNumber[i],/normal,alignment=1.0
  ;    XYOUTS,0.1+(xnum+1)*ixsize-0.025,0.1+(ynum)*iysize+0.8*iysize,alphabet[n_elements(HPASSNAme)-1-acount],/normal,alignment=1.0
      XYOUTS,position[2]-0.01,position[1]+0.02*iysize,JName[i]+' '+JNumber[i],/normal,alignment=1.0
      XYOUTS,position[2]-0.01,position[1]+0.75*iysize,alphabet[n_elements(HPASSNAme)-2-acount],/normal,alignment=1.0
      plotsym, 3 ,/fill
      if fitresult[i] EQ 1.5 then oplot,[maxRadius/16.],[maxVel-MaxVel/16.],psym=8
;XYOUTS,position[0]+0.01,position[1]+0.8*iysize,'Star',/normal,alignment=1.0
     acount++
   ;   XYOUTS,0.1+(xnum+1)*ixsize-0.025,0.1+(ynum)*iysize+0.075*iysize,HPASSName[i],/normal,alignment=1.0
      IF xnum LT 4 then xnum++ else begin
         xnum=0
         ynum++
      ENDELSE

   ENDIF
ENDfor
XYOUTS,0.5,0.05,'Radius (arcsec)',alignment=0.5,/normal,charsize=2.
XYOUTS,0.05,0.5,'Velocity (km s!E-1!N)',alignment=0.5,ORIENTATION=90  ,/normal,charsize=2.

DEVICE,/CLOSE
convert_ps,main_dirsl+'AllRC.ps',/png,/trim

end
