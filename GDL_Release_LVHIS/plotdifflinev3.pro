Pro plotdifflinev3,fitresult,PAinput,RCINCL,position=position,xtitle=xtitle,ytitle=ytitle,limit=lim,sigmanames=sigma,noerase=noerase,filename=filename,errors=errors,_EXTRA=ex,anno=anno,name=name

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
  IF n_elements(PAinput[0,*]) GT 1  then begin
     IF n_elements(RCINCL[0,*]) LT 2 then begin
        tmp=RCINCL
        RCINCL=dblarr(n_elements(tmp),3)
        RCINCL[*,0]=tmp
        RCINCL[*,1]=tmp
        RCINCL[*,2]=tmp
     ENDIF
     IF n_elements(errors) GT 0 then begin
        IF n_elements(errors[0,*]) LT 2 then begin
           tmp=errors
           RCINCL=dblarr(n_elements(tmp),3)
           errors[*,0]=tmp
           errors[*,1]=tmp
           errors[*,2]=tmp
        ENDIF
     ENDIF
     IF n_elements(sigma) EQ 1 then begin
        tmp=sigma
        sigma=strarr(3)
        sigma[0]='\sigma_{'+tmp+'}= '
        sigma[1]='\sigma_{'+tmp+'-RC}= '
        sigma[2]='\sigma_{'+tmp+'-DF}= '
     ENDIF

                                ;either against radius but probably
                                ;inclination is mor usefull so we do
                                ;that first
     tmp=WHERE(fitresult EQ 1. OR fitresult EQ 1.5)
     check1=WHERE(double(PAinput[tmp,0]) NE 0.)
     check2=WHERE(double(PAinput[tmp,1]) NE 0.)
     check3=WHERE(double(PAinput[tmp,2]) NE 0.)
     maxpay=MAX(double([double(PAinput[tmp[check3],2]),double(PAinput[tmp[check2],1]),double(PAinput[tmp[check1],0])]),min=minpay)
     maxpay=maxpay+maxpay/10.
     minpay=minpay-maxpay/10.
     check4=WHERE(double(RCINCL[tmp,0]) NE 0.)
     maxpa=MAX(double([double(RCINCL[tmp[check3],2]),double(RCINCL[tmp[check2],1]),double(RCINCL[tmp[check1],0])]),min=minpa)
     maxpa=maxpa+maxpa/10.
     minpa=minpa-maxpa/10.
      trigger=0
     tmp=WHERE(fitresult EQ 1.)
     check1=WHERE(double(PAinput[tmp,0]) NE 0.)
     check2=WHERE(double(PAinput[tmp,1]) NE 0.)
     check3=WHERE(double(PAinput[tmp,2]) NE 0.)
;  if n_elements(lim) EQ 0 then begin
;     lim=[(STDDEV(double([RAinput[tmp[check2],1],RAinput[tmp[check3],2]]))),STDDEV(double([DECinput[tmp[check2],1],DECinput[tmp[check3],2]]))]
;     trigger=1
;  endif
   ;  if n_elements(lim) EQ 0 then lim=STDDEV(double(PAinput[tmp[check1],0]))
;print,minpay
     IF keyword_set(noerase) then begin
        plot,RCINCL[*,0],PAinput[*,0],position=[position[0]+(position[2]-position[0])/2,position[1],position[2],position[3]],yrange=[minpay,maxpay],psym=8,xrange=[minpa,maxpa],YSTYLE=1,xstyle=1,/NODATA,/NOERASE,ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],_STRICT_EXTRA=ex
     ENDIF ELSE BEGIN
        plot,RCINCL[*,0],PAinput[*,0],position=[position[0]+(position[2]-position[0])/2,position[1],position[2],position[3]],yrange=[minpay,maxpay],psym=8,xrange=[minpa,maxpa],YSTYLE=1,xstyle=1,/NODATA,ytickname=[' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '],_STRICT_EXTRA=ex
     ENDELSE
      IF n_elements(errors) EQ 0 then BEGIN
          errors = PAinput
          errors[*] =0
        ENDIF
     ;loadct,40,/SILENT
     IF n_elements(PAinput[0,*]) GT 1 then begin
        IF n_elements(filename) GT 0  then begin
           IF n_elements(errors) GT 0 then begin
              openu,1,filename,/APPEND
              printf,1,sigma[0]+string(STDDEV(double(PAinput[tmp[check1],0])),format='(F6.2)')+', av='+string(TOTAL(double(PAinput[tmp[check1],0])/ABS(errors[tmp[check1],0]))/TOTAL(1./ABS([errors[tmp[check1],0]])),format='(F6.2)')
              tmp=WHERE(fitresult EQ 1)
              printf,1,sigma[1]+string(STDDEV(double(PAinput[tmp[check2],1])),format='(F6.2)')+', av='+string(TOTAL(double(PAinput[tmp[check2],1])/ABS(errors[tmp[check2],1]))/TOTAL(1./ABS([errors[tmp[check2],1]])),format='(F6.2)')
              printf,1,sigma[2]+string(STDDEV(double(PAinput[tmp[check3],2])),format='(F6.2)')+', av='+string(TOTAL(double(PAinput[tmp[check3],2])/ABS(errors[tmp[check3],2]))/TOTAL(1./ABS([errors[tmp[check3],2]])),format='(F6.2)')
              close,1
           ENDIF ELSE BEGIN
              openu,1,filename,/APPEND
              printf,1,sigma[0]+string(STDDEV(double(PAinput[tmp[check1],0])),format='(F6.2)')+', av='+string(TOTAL(double(PAinput[tmp[check1],0]))/n_elements(PAinput[*,0]),format='(F6.2)')
              printf,1,sigma[1]+string(STDDEV(double(PAinput[tmp[check2],1])),format='(F6.2)')+', av='+string(TOTAL(double(PAinput[tmp[check2],1]))/n_elements(PAinput[tmp[check2],1]),format='(F6.2)')
              printf,1,sigma[2]+string(STDDEV(double(PAinput[tmp[check3],2])),format='(F6.2)')+', av='+string(TOTAL(double(PAinput[tmp[check3],2]))/n_elements(PAinput[tmp[check3],2]),format='(F6.2)')
           close,1
           ENDELSE

        ENDIF
      ;  loadct,40

        tmp=WHERE(fitresult EQ 1)
        PLOTSYM, 0 , /FILL
        if tmp[0] NE -1 then begin

          xerr= dblarr(n_elements(tmp[check2]))
          fat_ploterror,double(RCINCL[tmp[check2],1]),double(PAinput[tmp[check2],1]),xerr,errors[tmp[check2],1],psym = 8, $
                     color=blue,ERRCOLOR = blue, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
           ;oplot,double(RCINCL[tmp[check2],1]),double(PAinput[tmp[check2],1]),color=blue,psym=8
           ;IF n_elements(errors) GT 0 then begin
          ;    ERRPLOT,double(RCINCL[tmp[check2],1]),double(PAinput[tmp[check2],1])-errors[tmp[check2],1],double(PAinput[tmp[check2],1])+errors[tmp[check2],1],color=50,thick=!p.thick/2.
          ; ENDIF
          xerr= dblarr(n_elements(tmp[check3]))
          fat_ploterror,double(RCINCL[tmp[check3],2]),double(PAinput[tmp[check3],2]),xerr,errors[tmp[check3],2],psym = 8, $
                     color=red,ERRCOLOR = red, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
           ;oplot,double(RCINCL[tmp[check3],2]),double(PAinput[tmp[check3],2]),color=254,psym=8
           ;IF n_elements(errors) GT 0 then begin
            ;  ERRPLOT,double(RCINCL[tmp[check3],2]),double(PAinput[tmp[check3],2])-errors[tmp[check3],2],double(PAinput[tmp[check3],2])+errors[tmp[check3],2],color=254,thick=!p.thick/2.
          ; ENDIF
        ENDIF
        IF TOTAL(errors) GT 0 then begin
           av=TOTAL([PAinput[tmp[check2],1]/ABS(errors[tmp[check2],1]),PAinput[tmp[check3],2]/ABS(errors[tmp[check3],2])])/TOTAL(1./ABS([errors[tmp[check2],1],errors[tmp[check3],2]]))
        ENDIF ELSE av=TOTAL([PAinput[tmp[check2],1],PAinput[tmp[check3],2]])/N_elements([PAinput[tmp[check2],1],PAinput[tmp[check3],2]])


        print,'getting an av from these values'
        print,PAinput[tmp[check2],1]
        print,PAinput[tmp[check3],2]
        print,'getting an av from these values'
        IF n_elements(errors) GT 0 then print,errors[tmp[check2],1]
        IF n_elements(errors) GT 0 then print,errors[tmp[check3],2]
        print,'av is here'
        openu,88,'../averageandsigma.txt',/APPEND
        printf,88,sigma[0],'Tir'
        printf,88,format='(2F10.5)',av,STDDEV(double([PAinput[tmp[check2],1],PAinput[tmp[check3],2]]))
        close,88

        print,av,STDDEV(double([PAinput[tmp[check2],1],PAinput[tmp[check3],2]]))
        ;loadct,0,/SILENT
        trigger=0
     ;   if n_elements(lim) EQ 0 then begin
        lim=STDDEV(double([PAinput[tmp[check2],1],PAinput[tmp[check3],2]]))
        trigger=1
      ;  endif
        oplot,[minpa,maxpa],[-1*lim+av,-1*lim+av],color=grey,linestyle=2
        oplot,[minpa,maxpa],[av,av],color=grey
        oplot,[minpa,maxpa],[lim+av,lim+av],color=grey,linestyle=2
        ;loadct,40,/SILENT

        tmp=WHERE(fitresult EQ 1.5)
        check1=WHERE(double(PAinput[tmp,0]) NE 0.)
        check2=WHERE(double(PAinput[tmp,1]) NE 0.)
        check3=WHERE(double(PAinput[tmp,2]) NE 0.)
        PLOTSYM, 3 , /FILL
        ;colormaps,'sauron_colormap'
        if tmp[0] NE -1 then begin
          xerr= dblarr(n_elements(tmp[check2]))
          fat_ploterror,double(RCINCL[tmp[check2],1]),double(PAinput[tmp[check2],1]),xerr,errors[tmp[check2],1],psym = 8, $
                     color=light_blue,ERRCOLOR = light_blue, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
           ;oplot,double(RCINCL[tmp[check2],1]),double(PAinput[tmp[check2],1]),color=blue,psym=8
           ;IF n_elements(errors) GT 0 then begin
          ;    ERRPLOT,double(RCINCL[tmp[check2],1]),double(PAinput[tmp[check2],1])-errors[tmp[check2],1],double(PAinput[tmp[check2],1])+errors[tmp[check2],1],color=50,thick=!p.thick/2.
                                ; ENDIF
          if check3[0] NE -1 then begin
             xerr=dblarr( n_elements(tmp[check3]))          
             fat_ploterror,double(RCINCL[tmp[check3],2]),double(PAinput[tmp[check3],2]),xerr,errors[tmp[check3],2],psym = 8, $
                           color=light_red,ERRCOLOR = light_red, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
          endif
           ;oplot,double(RCINCL[tmp[check2],1]),double(PAinput[tmp[check2],1]),color=65,psym=8,symsize=0.7
           ;IF n_elements(errors) GT 0 then begin
          ;    ERRPLOT,double(RCINCL[tmp[check2],1]),double(PAinput[tmp[check2],1])-errors[tmp[check2],1],double(PAinput[tmp[check2],1])+errors[tmp[check2],1],color=65,thick=!p.thick/2.
           ;ENDIF

           ;oplot,double(RCINCL[tmp[check3],2]),double(PAinput[tmp[check3],2]),color=240,psym=8,symsize=0.7
           ;IF n_elements(errors) GT 0 then begin
            ;  ERRPLOT,double(RCINCL[tmp[check3],2]),double(PAinput[tmp[check3],2])-errors[tmp[check3],2],double(PAinput[tmp[check3],2])+errors[tmp[check3],2],color=240,thick=!p.thick/2.
           ;ENDIF
        ENDIF
        IF n_elements(errors) GT 0 then begin
           IF n_elements(errors[*,1]) GT 1 then addi=errors[*,1] else begin
              if n_elements(errors[*,1]) EQ 1 then addi=replicate(errors[0,1],n_elements(PAinput[*,1])) else  addi=replicate(0.,n_elements(PAinput[*,1]))
           ENDELSE
        ENDIF ELSE addi=replicate(0.,n_elements(PAinput[*,1]))
        tmpx=WHERE(PAinput[*,1] NE 0. AND (PAinput[*,1]+addi LT -3.*lim+av OR PAinput[*,1]-addi GT 3.*lim+av))
        IF tmpx[0] NE -1 then begin
           openu,1,filename+'_outliers.txt',/APPEND
           printf,1,sigma[1]
           for hj=0,n_elements(tmpx)-1 do begin
              IF N_ELEMENTS(ERRORS) GT 1 THEN BEGIN
                 print,n_elements(name),n_elements(errors),n_elements(PAinput),tmpx[hj]
                 printf,1,string(name[tmpx[hj]],format='(A80)')+string(PAinput[tmpx[hj],1],format='(F20.10)')+string(errors[tmpx[hj],1],format='(F20.10)')
              ENDIF ELSE printf,1,string(name[tmpx[hj]],format='(A80)')+string(PAinput[tmpx[hj],1],format='(F20.10)')
           endfor
           close,1
        ENDIF
        IF n_elements(errors) GT 0 then begin
           IF n_elements(errors[*,2]) GT 1 then addi=errors[*,2] else begin
              if n_elements(errors[*,2]) EQ 1 then addi=replicate(errors[0,2],n_elements(PAinput[*,2])) else  addi=replicate(0.,n_elements(PAinput[*,2]))
           ENDELSE
        ENDIF ELSE addi=replicate(0.,n_elements(PAinput[*,2]))
        tmpx=WHERE(PAinput[*,2] NE 0. AND (PAinput[*,2]+addi LT -3.*lim+av OR PAinput[*,2]-addi GT 3.*lim+av))
        IF tmpx[0] NE -1 then begin
           openu,1,filename+'_outliers.txt',/APPEND
           printf,1,sigma[2]
           for hj=0,n_elements(tmpx)-1 do begin
              IF N_ELEMENTS(ERRORS) GT 1 THEN BEGIN
                 print,n_elements(name),n_elements(errors),n_elements(PAinput),tmpx[hj]
                 printf,1,string(name[tmpx[hj]],format='(A80)')+string(PAinput[tmpx[hj],2],format='(F20.10)')+string(errors[tmpx[hj],2],format='(F20.10)')
              ENDIF ELSE printf,1,string(name[tmpx[hj]],format='(A80)')+string(PAinput[tmpx[hj],2],format='(F20.10)')
           endfor
           close,1
        ENDIF


        plot,RCINCL[*,0],PAinput[*,0],position=[position[0],position[1],position[0]+(position[2]-position[0])/2,position[3]],yrange=[minpay,maxpay],psym=8,xrange=[minpa,maxpa],YSTYLE=1,xstyle=1,/NODATA,/NOERASE,_STRICT_EXTRA=ex
        ;loadct,0,/SILENT

        tmp=WHERE(fitresult EQ 1.)
        check1=WHERE(double(PAinput[tmp,0]) NE 0.)
        check2=WHERE(double(PAinput[tmp,1]) NE 0.)
        check3=WHERE(double(PAinput[tmp,2]) NE 0.)
        IF TOTAL(errors) GT 0 then begin
           av=TOTAL(PAinput[tmp[check1],0]/ABS(errors[tmp[check1],0]))/TOTAL(1./ABS([errors[tmp[check1],0]]))
        ENDIF ELSE av=TOTAL([PAinput[tmp[check1],0]])/N_elements([PAinput[tmp[check1],0]])
        print,'getting an av from these values'
        print,PAinput[tmp[check1],0]

        print,'getting an av from these values'
        IF n_elements(errors) GT 0 then print,errors[tmp[check1],0]

        print,'av is here'
        print,av


        ;loadct,0,/SILENT
        if trigger then lim=STDDEV(PAinput[tmp[check1],0])
        openu,88,'../averageandsigma.txt',/APPEND
        printf,88,sigma[0]+'   RC-DF'
        printf,88,format='(2F10.5)',av,lim
        close,88
        oplot,[minpa,maxpa],[-1*lim+av,-1*lim+av],color=grey,linestyle=2
        oplot,[minpa,maxpa],[av,av],color=grey
        oplot,[minpa,maxpa],[lim+av,lim+av],color=grey,linestyle=2
        ;loadct,40,/SILENT
        PLOTSYM, 0 , /FILL
        IF tmp[0] NE -1 then begin

          xerr= dblarr(n_elements(tmp[check1]))
          fat_ploterror,double(RCINCL[tmp[check1],0]),double(PAinput[tmp[check1],0]),xerr,errors[tmp[check1],0],psym = 8, $
                     color=black,ERRCOLOR = black, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot

          ; oplot,double(RCINCL[tmp[check1],0]),double(PAinput[tmp[check1],0]),color=0,psym=8
           ;IF n_elements(errors) GT 0 then begin
          ;    ERRPLOT,double(RCINCL[tmp[check1],0]),double(PAinput[tmp[check1],0])-errors[tmp[check1],0],double(PAinput[tmp[check1],0])+errors[tmp[check1],0],color=0,thick=!p.thick/2.
          ; ENDIF
        ENDIF

        tmp=WHERE(fitresult EQ 1.5)
        check1=WHERE(double(PAinput[tmp,0]) NE 0.)
        check2=WHERE(double(PAinput[tmp,1]) NE 0.)
        check3=WHERE(double(PAinput[tmp,2]) NE 0.)

        PLOTSYM, 3 , /FILL
        ;loadct,0
        IF check1[0] NE -1 and tmp[0] NE -1 then begin
           
          xerr= dblarr(n_elements(tmp[check1]))
          fat_ploterror,double(RCINCL[tmp[check1],0]),double(PAinput[tmp[check1],0]),xerr,errors[tmp[check1],0],psym = 8, $
                     color=grey,ERRCOLOR = grey, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot
           ;oplot,double(RCINCL[tmp[check1],0]),double(PAinput[tmp[check1],0]),color=100,psym=8,symsize=0.7
           ;IF n_elements(errors) GT 0 then begin
          ;    ERRPLOT,double(RCINCL[tmp[check1],0]),double(PAinput[tmp[check1],0])-errors[tmp[check1],0],double(PAinput[tmp[check1],0])+errors[tmp[check1],0],color=100,thick=!p.thick/2.
          ; ENDIF
        ENDIF
        XYOUTS,position[0]-0.05,position[1]+(position[3]-position[1])/2,ytitle,/normal,alignment=0.5,orientation=90
        XYOUTS,position[0]+(position[2]-position[0])/2,position[1]-0.05,xtitle,/normal,alignment=0.5
        IF n_elements(errors) GT 0 then begin
           IF n_elements(errors[*,0]) GT 1 then addi=errors[*,0] else begin
              if n_elements(errors[*,0]) EQ 1 then addi=replicate(errors[0,0],n_elements(PAinput[*,0])) else  addi=replicate(0.,n_elements(PAinput[*,0]))
           ENDELSE
        ENDIF ELSE addi=replicate(0.,n_elements(PAinput[*,0]))
        tmpx=WHERE(PAinput[*,0] NE 0. AND (PAinput[*,0]+addi LT -3.*lim+av OR PAinput[*,0]-addi GT 3.*lim+av))
        IF tmpx[0] NE -1 then begin
           openu,1,filename+'_outliers.txt',/APPEND
           printf,1,sigma[0]
           for hj=0,n_elements(tmpx)-1 do begin
              IF N_ELEMENTS(ERRORS) GT 1 THEN BEGIN
                 print,n_elements(name),n_elements(errors),n_elements(PAinput),tmpx[hj]
                 printf,1,string(name[tmpx[hj]],format='(A80)')+string(PAinput[tmpx[hj],0],format='(F20.10)')+string(errors[tmpx[hj],0],format='(F20.10)')
              ENDIF ELSE printf,1,string(name[tmpx[hj]],format='(A80)')+string(PAinput[tmpx[hj],0],format='(F20.10)')
           endfor
           close,1
        ENDIF

        IF n_elements(anno) NE 0 then XYOUTS,position[0]+0.01,position[3]-0.014,anno,/normal,alignment=0.5
     ENDIF ELSE BEGIN
        IF n_elements(filename) GT 0 then begin
           openu,1,filename,/APPEND
           printf,1,sigma[0]+string(STDDEV(double(PAinput[*])),format='(F6.2)')+', av='+string(TOTAL(double(PAinput[*]))/n_elements(PAinput[*]),format='(F6.2)')
           close,1
        ENDIF
        xerr= dblarr(n_elements(RCINCL))
        fat_ploterror,double(RCINCL),double(PAinput),xerr,errors,psym = 8, $
                   color=black,ERRCOLOR = black, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot

        ;oplot,double(RCINCL),double(PAinput),color=0,psym=8
        ;IF n_elements(errors) GT 0 then begin
        ;   ERRPLOT,double(RCINCL),double(PAinput)-errors,double(PAinput)+errors,color=0,thick=!p.thick/2.
        ;ENDIF
     ENDELSE
  ENDIF ELSE BEGIN
     tmp=WHERE(fitresult EQ 1 AND PAinput NE 0.)
     maxpay=MAX(double(PAinput[tmp]),min=minpay)
     maxpa=MAX(double(RCINCL[tmp]),min=minpa)
     maxpay=maxpay+maxpay/10.
     minpay=minpay-maxpay/10.
     if n_elements(lim) EQ 0 then lim=STDDEV(double(PAinput[tmp,0]))
     maxpa=maxpa+maxpa/10.
     minpa=minpa-maxpa/10.

     plot,RCINCL[tmp],PAinput[tmp],position=[position[0],position[1],position[2],position[3]],yrange=[minpay,maxpay],psym=8,xrange=[minpa,maxpa],YSTYLE=1,xstyle=1,/NODATA,/NOERASE,_STRICT_EXTRA=ex

         IF TOTAL(errors) GT 0 then begin
           av=TOTAL(PAinput[tmp]/errors[tmp])/TOTAL(1./errors[tmp])
        ENDIF ELSE av=TOTAL(PAinput[tmp])/N_elements(PAinput[tmp])

      ;  loadct,0,/SILENT

        oplot,[minpa,maxpa],[-1*lim+av,-1*lim+av],color=grey,linestyle=2
        oplot,[minpa,maxpa],[av,av],color=grey
        oplot,[minpa,maxpa],[lim+av,lim+av],color=grey,linestyle=2

    ; loadct,40,/SILENT


     IF n_elements(filename) GT 0 then begin
        openu,1,filename,/APPEND
        printf,1,sigma[0]+string(STDDEV(double(PAinput[*])),format='(F6.2)')+', av='+string(TOTAL(double(PAinput[*]))/n_elements(PAinput[*]),format='(F6.2)')
        close,1
     ENDIF
     xerr= dblarr(n_elements(tmp))
     fat_ploterror,double(RCINCL[tmp]),double(PAinput[tmp]),xerr,errors[tmp],psym = 8, $
                color=black,ERRCOLOR = black, ERRTHICK=!p.thick*0.4,symsize=ssize,/over_plot

     ;oplot,double(RCINCL[tmp]),double(PAinput[tmp]),color=0,psym=8
     ;IF n_elements(errors) GT 0 then begin
    ;    ERRPLOT,double(RCINCL[tmp]),double(PAinput[tmp])-errors[tmp],double(PAinput[tmp])+errors[tmp],color=0,thick=!p.thick/2.
     ;ENDIF
     IF n_elements(anno) NE 0 then XYOUTS,position[0]+0.01,position[3]-0.02,anno,/normal,alignment=0.5
     XYOUTS,position[0]-0.05,position[1]+(position[3]-position[1])/2,ytitle,/normal,alignment=0.5,orientation=90
     XYOUTS,position[0]+(position[2]-position[0])/2,position[1]-0.063,xtitle,/normal,alignment=0.5
     IF n_elements(anno) NE 0 then XYOUTS,position[0]+0.01,position[3]-0.014,anno,/normal,alignment=0.5
  ENDELSE



END
