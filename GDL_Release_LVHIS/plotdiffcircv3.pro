Pro plotdiffcircv3,fitresult,RAinput,DECinput,position=position,xtitle=xtitle,ytitle=ytitle,radius=lim,sigmanames=sigma,noerase=noerase,filename=filename,anno=anno,_EXTRA=ex,name=name

  black = '000000'x
  grey = '808080'x
  blue = '0000FF'x
  light_blue = '7CACFF'x
  red = 'FF0000'x
  light_red = 'FFCCCC'x

  A = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a unit circle with 16 points,
; and set the filled flag:
  USERSYM, COS(A), SIN(A), /FILL
      ;Next up is the PA which we can plot
  IF n_elements(DECinput[0,*]) GT 1  then begin
     IF n_elements(sigma) EQ 2 then begin
        tmp=sigma
        sigma=strarr(3,2)
        sigma[0,0]='\sigma_{'+tmp[0]+'}= '
        sigma[1,0]='\sigma_{'+tmp[0]+'-RC}= '
        sigma[2,0]='\sigma_{'+tmp[0]+'-DF}= '
        sigma[0,1]='\sigma_{'+tmp[1]+'}= '
        sigma[1,1]='\sigma_{'+tmp[1]+'-RC}= '
        sigma[2,1]='\sigma_{'+tmp[1]+'-DF}= '
     ENDIF

  ENDIF

                                ;either against radius but probably
                                ;inclination is mor usefull so we do
                                ;that first
  print,fitresult
  tmp=WHERE(fitresult EQ 1. OR fitresult EQ 1.5)
  print,tmp
  print,n_elements(tmp),n_elements(RAinput[*,0])
   check1=WHERE(double(RAinput[tmp,0]) NE 0.)
   check2=WHERE(double(RAinput[tmp,1]) NE 0.)
   check3=WHERE(double(RAinput[tmp,2]) NE 0.)
   maxpay=MAX(double([double(DECinput[tmp[check3],2]),double(DECinput[tmp[check2],1]),double(DECinput[tmp[check1],0])]),min=minpay)
   maxpay=maxpay+5.
   minpay=minpay-5
   maxpa=MAX(double([double(RAinput[tmp[check3],2]),double(RAinput[tmp[check2],1]),double(RAinput[tmp[check1],0])]),min=minpa)
   maxpa=maxpa+5.
   minpa=minpa-5
   trigger=0

   tmp=WHERE(fitresult EQ 1.)
   check1=WHERE(double(RAinput[tmp,0]) NE 0.)
   check2=WHERE(double(RAinput[tmp,1]) NE 0.)
   check3=WHERE(double(RAinput[tmp,2]) NE 0.)
  if n_elements(lim) EQ 0 then begin
     lim=[(SIGMA(double([RAinput[tmp[check2],1],RAinput[tmp[check3],2]]))),SIGMA(double([DECinput[tmp[check2],1],DECinput[tmp[check3],2]]))]
     trigger=1
  endif
;print,minpay
  IF keyword_set(noerase) then begin
     plot,RAinput[*,0],DECinput[*,0],position=[position[0]+(position[2]-position[0])/2.,position[1],position[2],position[3]],yrange=[minpay,maxpay],psym=8,xrange=[minpa,maxpa],YSTYLE=1,xstyle=1,/NODATA,/NOERASE,ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],_STRICT_EXTRA=ex
  ENDIF ELSE BEGIN
     plot,RAinput[*,0],DECinput[*,0],position=position,yrange=[minpay,maxpay],psym=8,xrange=[minpa,maxpa],YSTYLE=1,xstyle=1,/NODATA,ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],_STRICT_EXTRA=ex
  ENDELSE
;   loadct,0,/SILENT
   oplot,[minpa,maxpa],[0,0],thick=thick,color=grey
   oplot,[0,0],[minpay,maxpay],thick=thick,color=grey
   A = FINDGEN(170) * (!PI*2/16.)
   IF n_elements(lim) EQ 2 then begin
      x=COS(A)*lim[0]
      y=SIN(A)*lim[1]
   ENDIF else begin
      x=COS(A)*lim
      y=SIN(A)*lim
   ENDELSE

   tmpx=WHERE(RAinput[*,1] NE 0. AND (RAinput[*,1] LT -3.*lim[0] OR RAinput[*,1] GT 3.*lim[0]))
   IF tmpx[0] NE -1 then begin
      openu,1,filename+'_outliers.txt',/APPEND
      printf,1,sigma[1,0]
      for hj=0,n_elements(tmpx)-1 do begin
         printf,1,string(name[tmpx[hj]],format='(A80)')+string(RAinput[tmpx[hj],1],format='(F20.10)')
      endfor
      close,1
   ENDIF
   tmpx=WHERE(RAinput[*,2] NE 0. AND (RAinput[*,2] LT -3.*lim[0] OR RAinput[*,2] GT 3.*lim[0]))
   IF tmpx[0] NE -1 then begin
      openu,1,filename+'_outliers.txt',/APPEND
      printf,1,sigma[2,0]
      for hj=0,n_elements(tmpx)-1 do begin
         printf,1,string(name[tmpx[hj]],format='(A80)')+string(RAinput[tmpx[hj],2],format='(F20.10)')
      endfor
      close,1
   ENDIF
   tmpx=WHERE(DECinput[*,1] NE 0. AND (DECinput[*,1] LT -3.*lim[1] OR DECinput[*,1] GT 3.*lim[1]))
   IF tmpx[0] NE -1 then begin
      openu,1,filename+'_outliers.txt',/APPEND
      printf,1,sigma[1,1]
      for hj=0,n_elements(tmpx)-1 do begin
         printf,1,string(name[tmpx[hj]],format='(A80)')+string(DECinput[tmpx[hj],1],format='(F20.10)')
      endfor
      close,1
   ENDIF
   tmpx=WHERE(DECinput[*,2] NE 0. AND (DECinput[*,2] LT -3.*lim[1] OR DECinput[*,2] GT 3.*lim[1]))
   IF tmpx[0] NE -1 then begin
      openu,1,filename+'_outliers.txt',/APPEND
      printf,1,sigma[2,1]
      for hj=0,n_elements(tmpx)-1 do begin
         printf,1,string(name[tmpx[hj]],format='(A80)')+string(DECinput[tmpx[hj],2],format='(F20.10)')
      endfor
      close,1
   ENDIF


  openu,88,'../averageandsigma.txt',/APPEND
  printf,88,'RA all'
  printf,88,format='(2F10.5)',TOTAL(double([RAinput[tmp[check2],1],RAinput[tmp[check3],2]]))/n_elements(double([RAinput[tmp[check2],1],RAinput[tmp[check3],2]])),lim[0]
   printf,88,'DEC all'
  printf,88,format='(2F10.5)',TOTAL(double([DECinput[tmp[check2],1],DECinput[tmp[check3],2]]))/n_elements(double([DECinput[tmp[check2],1],DECinput[tmp[check3],2]])),lim[1]
  close,88
   oplot,x,y,color=grey
   ;loadct,40,/SILENT
   tmp=WHERE(fitresult EQ 1.)
   check1=WHERE(double(RAinput[tmp,0]) NE 0.)
   check2=WHERE(double(RAinput[tmp,1]) NE 0.)
   check3=WHERE(double(RAinput[tmp,2]) NE 0.)
   IF n_elements(filename) GT 0 then begin
      openu,1,filename,/APPEND
      printf,1,sigma[0,1]+string(SIGMA(double(DECinput[tmp[check1],0])),format='(F6.2)')+', av='+string(TOTAL(double(DEcinput[tmp[check1],0]))/n_elements(DECinput[tmp[check1],0]),format='(F6.2)')
      tmp=WHERE(fitresult EQ 1)
      printf,1,sigma[1,1]+string(SIGMA(double(DECinput[tmp[check2],1])),format='(F6.2)')+', av='+string(TOTAL(double(DECinput[tmp[check2],1]))/n_elements(DECinput[tmp[check2],1]),format='(F6.2)')
      printf,1,sigma[2,1]+string(SIGMA(double(DECinput[tmp[check3],2])),format='(F6.2)')+', av='+string(TOTAL(double(DECinput[tmp[check3],2]))/n_elements(DECinput[tmp[check3],2]),format='(F6.2)')

      printf,1,sigma[0,0]+string(SIGMA(double(RAinput[tmp[check1],0])),format='(F6.2)')+', av='+string(TOTAL(double(DEcinput[tmp[check1],0]))/n_elements(RAinput[tmp[check1],0]),format='(F6.2)')
      tmp=WHERE(fitresult EQ 1)
      printf,1,sigma[1,0]+string(SIGMA(double(RAinput[tmp[check2],1])),format='(F6.2)')+', av='+string(TOTAL(double(RAinput[tmp[check2],1]))/n_elements(RAinput[tmp[check2],1]),format='(F6.2)')
      printf,1,sigma[2,0]+string(SIGMA(double(RAinput[tmp[check3],2])),format='(F6.2)')+', av='+string(TOTAL(double(RAinput[tmp[check3],2]))/n_elements(RAinput[tmp[check3],2]),format='(F6.2)')
      close,1
   ENDIF

   tmp=WHERE(fitresult EQ 1)
   PLOTSYM, 0 , /FILL
   if tmp[0] NE -1 then begin
      oplot,double(RAinput[tmp[check2],1]),double(DECinput[tmp[check2],1]),color=blue,psym=8
      oplot,double(RAinput[tmp[check3],2]),double(DECinput[tmp[check3],2]),color=red,psym=8,thick=thick
   ENDIF
   tmp=WHERE(fitresult EQ 1.5)
   check1=WHERE(double(RAinput[tmp,0]) NE 0.)
   check2=WHERE(double(RAinput[tmp,1]) NE 0.)
   check3=WHERE(double(RAinput[tmp,2]) NE 0.)
   PLOTSYM, 3 , /FILL
   ;colormaps,'sauron_colormap'
   if tmp[0] NE -1 then begin
      oplot,double(RAinput[tmp[check2],1]),double(DECinput[tmp[check2],1]),color=light_blue,psym=8,symsize=0.7
      oplot,double(RAinput[tmp[check3],2]),double(DECinput[tmp[check3],2]),color=light_red,psym=8,thick=thick,symsize=0.7
   ENDIF


   plot,RAinput[*,0],DECinput[*,0],position=[position[0],position[1],position[0]+(position[2]-position[0])/2,position[3]],yrange=[minpay,maxpay],psym=8,xrange=[minpa,maxpa],YSTYLE=1,xstyle=1,/NODATA,/NOERASE,_STRICT_EXTRA=ex
;   loadct,0,/SILENT
   oplot,[minpa,maxpa],[0,0],thick=thick,color=grey
   oplot,[0,0],[minpay,maxpay],thick=thick,color=grey
   tmp=WHERE(fitresult EQ 1.)
   check1=WHERE(double(RAinput[tmp,0]) NE 0.)
   check2=WHERE(double(RAinput[tmp,1]) NE 0.)
   check3=WHERE(double(RAinput[tmp,2]) NE 0.)
   A = FINDGEN(170) * (!PI*2/16.)
   if trigger then begin
      lim=[(SIGMA(double(RAinput[tmp[check1],0]))),SIGMA(double(DECinput[tmp[check1],0]))]
      trigger=1
   endif

   IF n_elements(lim) EQ 2 then begin
      x=COS(A)*lim[0]
      y=SIN(A)*lim[1]
   ENDIF else begin
      x=COS(A)*lim
      y=SIN(A)*lim
   ENDELSE
   oplot,x,y,color=grey

   tmpx=WHERE(RAinput[*,0] NE 0. AND (RAinput[*,0] LT -3.*lim[0] OR RAinput[*,0] GT 3.*lim[0]))
   IF tmpx[0] NE -1 then begin
      openu,1,filename+'_outliers.txt',/APPEND
      printf,1,sigma[0,0]
      for hj=0,n_elements(tmpx)-1 do begin
         printf,1,string(name[tmpx[hj]],format='(A80)')+string(RAinput[tmpx[hj],0],format='(F20.10)')
      endfor
      close,1
   ENDIF
   tmpx=WHERE(DECinput[*,0] NE 0. AND (DECinput[*,0] LT -3.*lim[1] OR DECinput[*,0] GT 3.*lim[1]))
   IF tmpx[0] NE -1 then begin
      openu,1,filename+'_outliers.txt',/APPEND
      printf,1,sigma[0,1]
      for hj=0,n_elements(tmpx)-1 do begin
         printf,1,string(name[tmpx[hj]],format='(A80)')+string(DECinput[tmpx[hj],0],format='(F20.10)')
      endfor
      close,1
   ENDIF



   ;loadct,40,/SILENT
   PLOTSYM, 0 , /FILL
   oplot,double(RAinput[tmp[check1],0]),double(DECinput[tmp[check1],0]),color=black,psym=8
   ;IF n_elements(errors) GT 0 then begin
  ;    ERRPLOT,double(RAinput[tmp[check1],0]),double(DECinput[tmp[check1],0])-errors[tmp[check1],0],double(DECinput[tmp[check1],0])+errors[tmp[check1],0],color=0,thick=!p.thick/2.
   ;ENDIF

   tmp=WHERE(fitresult EQ 1.5)
   check1=WHERE(double(RAinput[tmp,0]) NE 0.)
   check2=WHERE(double(RAinput[tmp,1]) NE 0.)
   check3=WHERE(double(RAinput[tmp,2]) NE 0.)
   PLOTSYM, 3 , /FILL
   ;loadct,0
   oplot,double(RAinput[tmp[check1],0]),double(DECinput[tmp[check1],0]),color=grey,psym=8,symsize=0.7
   ;IF n_elements(errors) GT 0 then begin
    ;  ERRPLOT,double(RAinput[tmp[check1],0]),double(DECinput[tmp[check1],0])-errors[tmp[check1],0],double(DECinput[tmp[check1],0])+errors[tmp[check1],0],color=100,thick=!p.thick/2.
   ;ENDIF


    XYOUTS,position[0]-0.05,position[1]+(position[3]-position[1])/2,ytitle,/normal,alignment=0.5,orientation=90
    XYOUTS,position[0]+(position[2]-position[0])/2,position[1]-0.05,xtitle,/normal,alignment=0.5
    IF n_elements(anno) NE 0 then XYOUTS,position[0]+0.01,position[3]-0.014,anno,/normal,alignment=0.5

END
