Pro plotparameters,fileoutput,TirResult,RCResult,DFResult,map,name,unsmooth=TirExtra

;This is version 1.3 for general purposes with all plots estimated from
;the values  and the possibility to give the values, filenames plot
;multiple disk and other things 
  

tryinputagain:
thick=4
xthick=4
ythick=4
lthick=5
ssize=1.
charsize=3
charthick=0.1
IF map EQ 1 then seconddisk='N' else seconddisk='Y'
disp='SBR VROT PA INCL'
Variables=str_sep(strtrim(strcompress(disp),2),' ')
Variables=strupcase(Variables)
varunits=strarr(n_elements(Variables))

;attach Units to the chosen varible
for i=0,n_elements(Variables)-1 do begin
   tmp=str_sep(strtrim(strcompress(Variables[i]),2),'_')
   case tmp[0] of
      'VROT':Varunits[i]=textoidl('(km s^{-1})')
      'SDIS':Varunits[i]=textoidl('(km s^{-1})')
      'PA' :Varunits[i]=textoidl('(Degrees)')
      'INCL':Varunits[i]=textoidl('(Degrees)')
      'SBR':Varunits[i]=textoidl('(Jy km s^{-1} arcsec^{-2})')
      'Z0':Varunits[i]='(Arcsec)'
      'DVRO':Varunits[i]=textoidl('(km s^{-1} kpc^{-1})')
      'GA1A':Varunits[i]=textoidl('(Jy km s^{-1} arcsec^{-2})')
      'GA1P':Varunits[i]=textoidl('(Degrees)')
      'GA1D':Varunits[i]=textoidl('(Arcsec)')
      'VRAD':Varunits[i]=textoidl('(km s^{-1})')
      else:Varunits[i]=''
   endcase
endfor
rings=[n_elements(TirResult[*,0]),n_elements(RCResult[*,0]),n_elements(DFResult[*,0])]
size=MAX(Rings)
IF size eq 1 then begin
   openw,1,'NoFitsFound.txt'
   printf,1,'We could not find any type of radial fit for this galaxy'
   close,1
   goto,nodataavailable
ENDIF
;make a string in hh:mm:ss for central RA and DEC
IF TOTAL(TirResult) EQ 0. then begin
   IF map EQ 1 then TirResult=dblarr(2,10) else begin
     openw,1,'NoSecondFound.txt'
     printf,1,'You are trying to plot the 2 disk solution but it is not found'
     printf,1,'Try plotting the one disk solution if you want only rotcur and Diskfit' 
     close,1
     goto,nodataavailable  
  ENDELSE
ENDIF
IF TOTAL(DFResult) EQ 0. then DFResult=dblarr(2,13)
IF TOTAL(RCResult) EQ 0. then RCResult=dblarr(2,13)
RA=strarr(3)
DEC=strarr(3)
vsys=strarr(3)

tmpRA=TirResult[0,1]
tmpDEC=TirResult[0,2]
convertRADEC,tmpRA,tmpDEC
RA[0]=tmpRA
DEC[0]=tmpDEC
vsys[0]=strtrim(strcompress(TirResult[0,3]),2)
tmpRA=RCResult[0,9]
tmpDEC=RCResult[0,11]
convertRADEC,tmpRA,tmpDEC
RA[1]=tmpRA
DEC[1]=tmpDEC
vsys[1]=strtrim(strcompress(RCResult[0,1]),2)
tmpRA=DFResult[0,9]
tmpDEC=DFResult[0,11]
convertRADEC,tmpRA,tmpDEC
RA[2]=tmpRA
DEC[2]=tmpDEC
vsys[2]=strtrim(strcompress(DFResult[0,1]),2)


;then obtain the maximums and minimums for the variables to be
;plotted.

maxvar=dblarr(n_elements(Variables))
minvar=dblarr(n_elements(Variables))
buffer=dblarr(n_elements(Variables))
minvar[*]=100000
maxvar[*]=-10000

if map EQ 1 then begin
   tmp1=MAX([Tirresult[*,7],Tirresult[*,8]],min=tmp2)
   maxvar[0]=tmp1
   minvar[0]=tmp2
   case (1) of
      TOTAL(TirResult) EQ 0 AND TOTAL(RCResult) EQ 0 AND TOTAL(DFResult) NE 0:begin          
          tmp1=MAX([DFresult[*,3]],min=tmp2)
          maxvar[1]=tmp1
          minvar[1]=tmp2
          tmp1=MAX([DFresult[*,5]],min=tmp2)
          maxvar[2]=tmp1
          minvar[2]=tmp2
          tmp1=MAX([DFresult[*,7]],min=tmp2)
          maxvar[3]=tmp1
          minvar[3]=tmp2
       end
      TOTAL(TirResult) EQ 0 AND TOTAL(DFResult) EQ 0 AND TOTAL(RCResult) NE 0:begin          
         tmp1=MAX([RCresult[*,3]],min=tmp2)
         maxvar[1]=tmp1
         minvar[1]=tmp2
         tmp1=MAX([RCresult[*,5]],min=tmp2)
         maxvar[2]=tmp1
         minvar[2]=tmp2
         tmp1=MAX([RCresult[*,7]],min=tmp2)
         maxvar[3]=tmp1
         minvar[3]=tmp2      
      end
      TOTAL(TirResult) EQ 0 AND TOTAL(DFResult) NE 0 AND TOTAL(RCResult) NE 0:begin
         tmp1=MAX([RCresult[*,3],DFresult[*,3]],min=tmp2)
         maxvar[1]=tmp1                                 
         minvar[1]=tmp2                                 
         tmp1=MAX([RCresult[*,5],DFresult[*,5]],min=tmp2)
         maxvar[2]=tmp1                                 
         minvar[2]=tmp2                                 
         tmp1=MAX([RCresult[*,7],DFresult[*,7]],min=tmp2)          
         maxvar[3]=tmp1
         minvar[3]=tmp2  
      end
      TOTAL(TirResult) NE 0. AND TOTAL(DFResult) NE 0. AND TOTAL(RCResult) EQ 0.:begin
         tmp1=MAX([Tirresult[*,6],DFresult[*,3]],min=tmp2)
         maxvar[1]=tmp1                                       
         minvar[1]=tmp2                                       
         tmp1=MAX([Tirresult[*,5],DFresult[*,5]],min=tmp2)
         maxvar[2]=tmp1                                       
         minvar[2]=tmp2                                       
         tmp1=MAX([Tirresult[*,4],DFresult[*,7]],min=tmp2)
         maxvar[3]=tmp1
         minvar[3]=tmp2
      end

      TOTAL(TirResult) NE 0. AND TOTAL(DFResult) EQ 0. AND TOTAL(RCResult) NE 0.:begin
         tmp1=MAX([Tirresult[*,6],RCresult[*,3]],min=tmp2)
         maxvar[1]=tmp1                                       
         minvar[1]=tmp2                                       
         tmp1=MAX([Tirresult[*,5],RCresult[*,5]],min=tmp2)
         maxvar[2]=tmp1                                       
         minvar[2]=tmp2                                       
         tmp1=MAX([Tirresult[*,4],RCresult[*,7]],min=tmp2)
         maxvar[3]=tmp1
         minvar[3]=tmp2
      end
      else:begin
         tmp1=MAX([Tirresult[*,6],RCresult[*,3],DFresult[*,3]],min=tmp2)
         maxvar[1]=tmp1                                                     
         minvar[1]=tmp2                                                     
         tmp1=MAX([Tirresult[*,5],RCresult[*,5],DFresult[*,5]],min=tmp2)
         maxvar[2]=tmp1                                                     
         minvar[2]=tmp2                                                     
         tmp1=MAX([Tirresult[*,4],RCresult[*,7],DFresult[*,7]],min=tmp2)
         maxvar[3]=tmp1
         minvar[3]=tmp2      
      end
   endcase
endif else begin
   tmp1=MAX([Tirresult[*,7],Tirresult[*,11]],min=tmp2)
   maxvar[0]=tmp1
   minvar[0]=tmp2
   case (1) of
      TOTAL(TirResult) NE 0. AND TOTAL(DFResult) NE 0. AND TOTAL(RCResult) EQ 0.:begin
         tmp1=MAX([Tirresult[*,6],Tirresult[*,10],DFresult[*,3]],min=tmp2)
         maxvar[1]=tmp1                                                       
         minvar[1]=tmp2                                                       
         tmp1=MAX([Tirresult[*,5],Tirresult[*,9],DFresult[*,5]],min= tmp2)
         maxvar[2]=tmp1                                                       
         minvar[2]=tmp2                                                       
         tmp1=MAX([Tirresult[*,4],Tirresult[*,8],DFresult[*,7]],min= tmp2)
         maxvar[3]=tmp1
         minvar[3]=tmp2   
      end
      TOTAL(TirResult) NE 0. AND TOTAL(DFResult) EQ 0. AND TOTAL(RCResult) NE 0.:begin
         tmp1=MAX([Tirresult[*,6],Tirresult[*,10],RCresult[*,3]],min=tmp2)
         maxvar[1]=tmp1                                                       
         minvar[1]=tmp2                                                       
         tmp1=MAX([Tirresult[*,5],Tirresult[*,9],RCresult[*,5]],min= tmp2)
         maxvar[2]=tmp1
         minvar[2]=tmp2 
         tmp1=MAX([Tirresult[*,4],Tirresult[*,8],RCresult[*,7]],min= tmp2    )
         maxvar[3]=tmp1
         minvar[3]=tmp2    
      end                                                                     
      else:begin                                                         
         tmp1=MAX([Tirresult[*,6],Tirresult[*,10],RCresult[*,3],DFresult[*,3]],min=tmp2)
         maxvar[1]=tmp1
         minvar[1]=tmp2 
         tmp1=MAX([Tirresult[*,5],Tirresult[*,9],RCresult[*,5],DFresult[*,5]],min= tmp2)
         maxvar[2]=tmp1
         minvar[2]=tmp2 
         tmp1=MAX([Tirresult[*,4],Tirresult[*,8],RCresult[*,7],DFresult[*,7]],min= tmp2)
         maxvar[3]=tmp1
         minvar[3]=tmp2 
      end
   endcase
   
endelse

;get the maximum radius
Maxradii=MAX([Tirresult[*,0],RCresult[*,0],DFresult[*,0]],Min=Minradii)
MaxRadii=MAxRadii+Tirresult[1,0]

;Set Some plotting values
!x.style=1.5
!y.style=1.5
!p.charsize=1.7

;estimate the ysize of the plots bases on the amount of variables
ysize=0.75/n_elements(Variables)
xsize=0.68
;/n_elements(existing)
print,'xsize',xsize
!P.FONT=1
mydevice = !D.NAME
SET_PLOT, 'PS', /COPY
;open plot

; Make a vector of 16 points, A[i] = 2pi/16:  
A = FINDGEN(17) * (!PI*2/16.)  
; Define the symbol to be a unit circle with 16 points,   
; and set the filled flag:  
USERSYM, COS(A), SIN(A), /FILL  

DEVICE, FILENAME=fileoutput,/color,/PORTRAIT,/ENCAPSULATED,xsize=40,ysize=48,SET_FONT='Times', /TT_FONT  
loadct,40,/SILENT

;first plot SBR
plot,Tirresult[*,0] ,Tirresult[*,7],position=[0.17,0.1,0.17+xsize,0.1+ysize],xtitle='Radius (arcmin)',$
                         xrange=[0.,maxradii],yrange=[minvar[0],maxvar[0]+1e-5],xthick=xthick,ythick=ythick,charthick=charthick,thick=thick,linestyle=2,charsize=charsize,/NODATA
XYOUTs,0.07,0.1+ysize/2.,'SBR',/NORMAL,alignment=0.5,ORIENTATION=90, CHARTHICK=charthick,charsize=charsize*1.25
XYOUTs,0.11,0.1+ysize/2.,'Jy/arcsec^2*KM/S' ,/NORMAL,alignment=0.5,ORIENTATION=90, CHARTHICK=charthick,charsize=charsize
oplot,Tirresult[*,0] ,Tirresult[*,7],linestyle=0,thick=lthick
oplot,Tirresult[*,0] ,Tirresult[*,7],psym=8,symsize=ssize
IF map eq 1 then begin
   oplot,Tirresult[*,0] ,Tirresult[*,8],linestyle=2,thick=lthick,color=50
   oplot,Tirresult[*,0] ,Tirresult[*,8],psym=8,symsize=ssize,color=50
ENDIF else begin
   oplot,Tirresult[*,0] ,Tirresult[*,11],linestyle=2,thick=lthick,color=50
   oplot,Tirresult[*,0] ,Tirresult[*,11],psym=8,symsize=ssize,color=50
   IF n_elements(TirExtra) GT 0 then begin
      oplot,TirExtra[*,0] ,TirExtra[*,11],linestyle=2,thick=lthick,color=200
      oplot,TirExtra[*,0] ,TirExtra[*,11],psym=8,symsize=ssize,color=200
      oplot,TirExtra[*,0] ,TirExtra[*,7],linestyle=0,thick=lthick,color=100
      oplot,TirExtra[*,0] ,TirExtra[*,7],psym=8,symsize=ssize,color=100
   ENDIF
ENDELSE
i=1
plot,Tirresult[*,0] ,Tirresult[*,6],position=[0.17,0.1+ysize*i,0.17+xsize,0.1+ysize*(i+1)],$
                    xrange=[0.,maxradii],yrange=[minvar[1]-5,maxvar[1]+5],xthick=xthick,ythick=ythick,charthick=charthick,thick=thick,charsize=charsize,$
                    /noerase,linestyle=0,xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],/NODATA          
oplot,Tirresult[*,0] ,Tirresult[*,6],psym=8,symsize=ssize
oplot,Tirresult[*,0] ,Tirresult[*,6],linestyle=0,thick=lthick

IF map ne 1 then begin
   print,'What the f is going on her?'
   print,Tirresult[*,16]
   errplot,Tirresult[*,0] ,Tirresult[*,6]-Tirresult[*,16],Tirresult[*,6]+Tirresult[*,16],thick=lthick/2.
   oplot,Tirresult[*,0] ,Tirresult[*,10],linestyle=2,thick=lthick,color=50
   oplot,Tirresult[*,0] ,Tirresult[*,10],psym=8,symsize=ssize,color=50
   errplot,Tirresult[*,0] ,Tirresult[*,10]-Tirresult[*,17],Tirresult[*,10]+Tirresult[*,17],color=50,thick=lthick/2.
   IF n_elements(TirExtra) GT 0 then begin
      oplot,TirExtra[*,0] ,TirExtra[*,10],linestyle=2,thick=lthick,color=200
      oplot,TirExtra[*,0] ,TirExtra[*,10],psym=8,symsize=ssize,color=200
      oplot,TirExtra[*,0] ,TirExtra[*,6],linestyle=0,thick=lthick,color=100
      oplot,TirExtra[*,0] ,TirExtra[*,6],psym=8,symsize=ssize,color=100
   ENDIF
ENDIF
oplot,RCresult[*,0] ,RCresult[*,3],linestyle=0,thick=lthick,color=150
oplot,RCresult[*,0] ,RCresult[*,3],psym=8,symsize=ssize,color=150
oplot,DFresult[*,0] ,DFresult[*,3],linestyle=2,thick=lthick,color=250
oplot,DFresult[*,0] ,DFresult[*,3],psym=8,symsize=ssize,color=250
XYOUTs,0.07,0.1+ysize*i+ysize/2.,'Velocity',/NORMAL,alignment=0.5,ORIENTATION=90, CHARTHICK=charthick,charsize=charsize*1.25
XYOUTs,0.11,0.1+ysize*i+ysize/2.,'KM/S',/NORMAL,alignment=0.5,ORIENTATION=90, CHARTHICK=charthick,charsize=charsize*1.25
i=2
plot,Tirresult[*,0] ,Tirresult[*,5],position=[0.17,0.1+ysize*i,0.17+xsize,0.1+ysize*(i+1)],$
                    xrange=[0.,maxradii],yrange=[minvar[i]-5,maxvar[i]+5],xthick=xthick,ythick=ythick,charthick=charthick,thick=thick,charsize=charsize,$
                    /noerase,linestyle=0,xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],/NODATA          
oplot,Tirresult[*,0] ,Tirresult[*,5],psym=8,symsize=ssize
oplot,Tirresult[*,0] ,Tirresult[*,5],linestyle=0,thick=lthick

IF map ne 1 then begin
   errplot,Tirresult[*,0] ,Tirresult[*,5]-Tirresult[*,14],Tirresult[*,5]+Tirresult[*,14],thick=lthick/2.
   oplot,Tirresult[*,0] ,Tirresult[*,9],linestyle=2,thick=lthick,color=50
   oplot,Tirresult[*,0] ,Tirresult[*,9],psym=8,symsize=ssize,color=50
   errplot,Tirresult[*,0] ,Tirresult[*,9]-Tirresult[*,15],Tirresult[*,9]+Tirresult[*,15],color=50,thick=lthick/2.
   IF n_elements(TirExtra) GT 0 then begin
      oplot,TirExtra[*,0] ,TirExtra[*,9],linestyle=2,thick=lthick,color=200
      oplot,TirExtra[*,0] ,TirExtra[*,9],psym=8,symsize=ssize,color=200
      oplot,TirExtra[*,0] ,TirExtra[*,5],linestyle=0,thick=lthick,color=100
      oplot,TirExtra[*,0] ,TirExtra[*,5],psym=8,symsize=ssize,color=100
   ENDIF
ENDIF
oplot,RCresult[*,0] ,RCresult[*,5],linestyle=0,thick=lthick,color=150
oplot,RCresult[*,0] ,RCresult[*,5],psym=8,symsize=ssize,color=150
oplot,DFresult[*,0] ,DFresult[*,5],linestyle=2,thick=lthick,color=250
oplot,DFresult[*,0] ,DFresult[*,5],psym=8,symsize=ssize,color=250
XYOUTs,0.07,0.1+ysize*i+ysize/2.,'PA',/NORMAL,alignment=0.5,ORIENTATION=90, CHARTHICK=charthick,charsize=charsize*1.25
XYOUTs,0.11,0.1+ysize*i+ysize/2.,'degree',/NORMAL,alignment=0.5,ORIENTATION=90, CHARTHICK=charthick,charsize=charsize*1.25

i=3
plot,Tirresult[*,0] ,Tirresult[*,4],position=[0.17,0.1+ysize*i,0.17+xsize,0.1+ysize*(i+1)],$
                    xrange=[0.,maxradii],yrange=[minvar[i]-5,maxvar[i]+5],xthick=xthick,ythick=ythick,charthick=charthick,thick=thick,charsize=charsize,$
                    /noerase,linestyle=0,xtickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],/NODATA          
oplot,Tirresult[*,0] ,Tirresult[*,4],psym=8,symsize=ssize
oplot,Tirresult[*,0] ,Tirresult[*,4],linestyle=0,thick=lthick
IF map ne 1 then begin
   errplot,Tirresult[*,0] ,Tirresult[*,4]-Tirresult[*,12],Tirresult[*,4]+Tirresult[*,12],thick=lthick/2.
   oplot,Tirresult[*,0] ,Tirresult[*,8],linestyle=2,thick=lthick,color=50
   oplot,Tirresult[*,0] ,Tirresult[*,8],psym=8,symsize=ssize,color=50
   errplot,Tirresult[*,0] ,Tirresult[*,8]-Tirresult[*,13],Tirresult[*,8]+Tirresult[*,13],color=50,thick=lthick/2.
   IF n_elements(TirExtra) GT 0 then begin
      oplot,TirExtra[*,0] ,TirExtra[*,8],linestyle=2,thick=lthick,color=200
      oplot,TirExtra[*,0] ,TirExtra[*,8],psym=8,symsize=ssize,color=200
      oplot,TirExtra[*,0] ,TirExtra[*,4],linestyle=0,thick=lthick,color=100
      oplot,TirExtra[*,0] ,TirExtra[*,4],psym=8,symsize=ssize,color=100
   ENDIF
ENDIF
oplot,RCresult[*,0] ,RCresult[*,7],linestyle=0,thick=lthick,color=150
oplot,RCresult[*,0] ,RCresult[*,7],psym=8,symsize=ssize,color=150
oplot,DFresult[*,0] ,DFresult[*,7],linestyle=2,thick=lthick,color=250
oplot,DFresult[*,0] ,DFresult[*,7],psym=8,symsize=ssize,color=250
XYOUTs,0.07,0.1+ysize*i+ysize/2.,'INCL',/NORMAL,alignment=0.5,ORIENTATION=90, CHARTHICK=charthick,charsize=charsize*1.25
XYOUTs,0.11,0.1+ysize*i+ysize/2.,'degree',/NORMAL,alignment=0.5,ORIENTATION=90, CHARTHICK=charthick,charsize=charsize*1.25






colordes = [0,150,250]
names=['TiRiFiC','Rotcur','DiskFit']
for j=0,2 do begin
   XYOUTS,0.17+(0.2)*xsize,0.1+ysize*(i+1)+0.025+j*0.01,'Systemic='+strtrim(string(vsys[j],format='(F10.1)'),2)+' R.A.='+RA[j]+' DEC.='+DEC[j],charthick=7,/normal,alignment=0.5,color=colordes[j]
    XYOUTS,0.17+(0.6)*xsize,0.1+ysize*(i+1)+0.025+j*0.01,names[j],charthick=7,/normal,alignment=0.5,color=colordes[j]
endfor
XYOUTS,0.17+(0.4)*xsize,0.1+ysize*(i+1)+0.025+(j+1)*0.01,name,charthick=7,/normal,alignment=0.5,color=0


 

DEVICE,/CLOSE
convert_ps,fileoutput,/DELETE,/trim
NoDataavailable:
end
