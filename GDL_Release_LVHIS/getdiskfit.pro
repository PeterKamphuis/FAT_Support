Pro getdiskfit,FileName,Parameters,Arrays,CubeName,noconv=noconv

                                ;Template = Array with the Template That the Ne file needs to be written to
                                ;NewFileName = Is the File Name Of the
                                ;Tirific file containing the new parameters
                                ;Arrays is a 2D array which will contain the
                                ;arrays [*,#changed variables] the ordering is
                                ;the same as the order of the variables to be
                                ;changed. However string values will be
                                ;excluded even if they are switched provide an empty array there to get
                                ;the default ordering of variableChange and their names
                           
                                ;Variables= gives an array indicating the
                                ;variables in the Template  it is not
                                ;necessary to provide but saves time
                                ;/EXTRACT will only get the arrays for the
                                ;keywords in Variable Change
  IF NOT keyword_set(Parameters) then begin
     Parameters=['RADI','VSYS','VSYS_ERR','VROT','VROT_ERR','PA','PA_ERR','INCL','INCL_ERR','XPOS','XPOS_ERR','YPOS','YPOS_ERR']
  ENDIF
  openr,lun,FileName,/GET_LUN
  h=''
  trigger=0
  WHILE ~ EOF(lun) do begin
     readf,lun,h
     values=strtrim(strcompress(str_sep(strtrim(strcompress(h),2),' ')),2)
     IF STRUPCASE(values[0]) EQ 'BEST' then trigger=1
     IF n_elements(values) GT 1 then begin
        IF values[1] EQ 'PA,' then begin
           DFpa=values[4]
           IF trigger then DFpaerr=values[6] else DFpaerr=1.
        ENDIF
        IF values[1] EQ 'incl' then begin
           DFincl=values[3]
           IF trigger then DFinclerr=values[5] else DFinclerr=1.
        ENDIF
        IF values[1] EQ 'eps:' then begin
           DFincl=double(acos(SQRT(((1-double(values[2]))^2-0.2^2)/0.96))*180./!pi+2.)  
           IF trigger then DFinclerr=values[4] else DFinclerr=1.
        ENDIF
        IF values[0] EQ 'r_w' then begin
           rw=double(values[3])  
           print,rw
           print,'What is going on'
           IF trigger then rwerr=values[5] else rwerr=1.
        ENDIF
        IF n_elements(values) GT 3 then begin
           
           IF values[2] EQ 'wphim:' then begin
            
              phim=double(values[3])  
            
              IF trigger then phimerr=values[5] else phimerr=1.
           ENDIF
           IF values[2] EQ 'welm:' then begin
              elm=double(values[3])  
              IF trigger then elmerr=values[5] else elmerr=1.
           ENDIF
        ENDIF
        IF values[0] EQ 'x,y' then begin
           DFRA=values[4]
           IF trigger then begin
              tmp=strtrim(strcompress(str_sep(values[6],',')))
              DFRAerr=double(tmp[0])
              DFDEC=values[7]
              DFDECerr=double(values[9])
           ENDIF else begin 
              DFDEC=values[5]
              DFRAerr=1.
              DFDECerr=1.
           ENDELSE
        ENDIF
        IF STRUPCASE(values[0]) EQ 'VSYS' then begin
           DFvsys=values[2]
           IF trigger then DFvsyserr=values[4] else DFvsyserr=1.
        ENDIF
     ENDIF
     
  ENDWHILE
 
  free_lun,lun
  inp=string(Filename[0])
  readcol,inp,DFradius,DFnpts,DFVt,DFeVt,DFVr,DFeVr,DFVmt,DFeVmt,DFVmr,eVmr,skipline=63
 ; print,FileName
  tmp=WHERE(Parameters EQ 'XPOS' OR Parameters EQ 'XPOS_ERR' OR Parameters EQ 'YPOS' OR Parameters EQ 'YPOS_ERR' OR Parameters EQ 'RADI')
  IF tmp[0] NE -1 AND NOT keyword_set(noconv) then Cube=readfits(Cubename,hed,/SILENT)
  original=n_elements(DFRadius)-1
  
  Arrays=dblarr(n_elements(DFRadius),n_elements(Parameters))
 FOR i=0,n_elements(Parameters)-1 do begin
     CASE Parameters[i] of
        'RADI': Arrays[*,i]=DFradius*(ABS(sxpar(hed,'CDELT1')*3600.)+ABS(sxpar(hed,'CDELT2')*3600.))/2.
        'VSYS':  Arrays[*,i]=DFvsys
        'VSYS_ERR': Arrays[*,i]=DFvsyserr
        'VROT':Arrays[0:original,i]=DFVt
        'VROT_ERR':Arrays[0:original,i]=DFeVt
        'PA':begin
           IF n_elements(phim) NE 0. then begin
              print,'yip'
              rmax=DFRadius[n_elements(DFRadius)-1]
              phi_0=phim/(rmax-rw)^2. 
              for x=0,n_elements(DFRadius)-1 do begin
                 IF DFRadius[x] LT rw then Arrays[x,i]=DFpa else  Arrays[x,i]=DFpa+phi_0*(DFRadius[x]-rw)^2
              endfor
           ENDIF ELSE begin
              Arrays[*,i]=DFpa
           ENDELSE
           IF Arrays[0,i] LT 0 then Arrays[*,i]=Arrays[*,i]+360.
        end
        'PA_ERR':Arrays[*,i]=DFpaerr
        'INCL':begin
           IF n_elements(elm) NE 0. then begin
              print,'yip'
              rmax=DFRadius[n_elements(DFRadius)-1]
              el_0=elm/(rmax-rw)^2. 
              for x=0,n_elements(DFRadius)-1 do begin
                 IF DFRadius[x] LT rw then Arrays[x,i]=DFincl else  Arrays[x,i]=DFincl+el_0*(DFRadius[x]-rw)^2
              endfor
           ENDIF ELSE begin
              Arrays[*,i]=DFincl
           ENDELSE
        end
        'INCL_ERR': Arrays[*,i]=DFinclerr
        'XPOS':Begin
           IF NOT keyword_set(noconv) then begin
              tmp=sxpar(hed,'CRVAL2')+sxpar(hed,'CDELT2')*(double(DFDEC)-sxpar(hed,'CRPIX2'))
              Arrays[*,i]=sxpar(hed,'CRVAL1')+sxpar(hed,'CDELT1')*(double(DFRA)-sxpar(hed,'CRPIX1'))/COS(tmp*!DtoR)
           ENDIF ELSE Arrays[*,i]=DFRA
        end
        'XPOS_ERR':Begin
           IF NOT keyword_set(noconv) then begin
              Arrays[*,i]=ABS(DFRAerr*sxpar(hed,'CDELT1')*3600.)
           ENDIF ELSE Arrays[*,i]=DFRAerr
        end
        'YPOS':Begin
            IF NOT keyword_set(noconv) then begin
              Arrays[*,i]=sxpar(hed,'CRVAL2')+sxpar(hed,'CDELT2')*(double(DFDEC)-sxpar(hed,'CRPIX2'))
           ENDIF ELSE Arrays[*,i]=DFDEC
        end
        'YPOS_ERR':Begin
           IF NOT keyword_set(noconv) then begin
              Arrays[*,i]=ABS(DFDECerr*sxpar(hed,'CDELT2')*3600.)
           ENDIF ELSE Arrays[*,i]=DFDECerr
        end
        else: Arrays[*,i]=-1
     ENDCASE
  ENDFOR
 



END
