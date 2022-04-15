Pro getrc,FileName,Parameters,Arrays,Cubename,noconv=noconv

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
  inp=string(Filename[0])
  readcol,inp,radiusv,width,vsysrot,vsysroterr,rotationv,errorv,vexp,vexperr,pav,paverr,incv,incverr,xposv,xposverr,yposv,yposverr,point,sigv
  IF NOT keyword_set(Parameters) then begin
     Parameters=['RADI','VSYS','VSYS_ERR','VROT','VROT_ERR','PA','PA_ERR','INCL','INCL_ERR','XPOS','XPOS_ERR','YPOS','YPOS_ERR']
  ENDIF
  Arrays=dblarr(n_elements(radiusv),n_elements(Parameters))
  tmp=WHERE(Parameters EQ 'XPOS' OR Parameters EQ 'XPOS_ERR' OR Parameters EQ 'YPOS' OR Parameters EQ 'YPOS_ERR')
  IF tmp[0] NE -1 then Cube=readfits(Cubename,hed,/SILENT)
  FOR i=0,n_elements(Parameters)-1 do begin
     CASE Parameters[i] of
        'RADI': Arrays[*,i]=radiusv
        'VSYS':  Arrays[*,i]=vsysrot
        'VSYS_ERR': Arrays[*,i]=vsysroterr
        'VROT':Arrays[*,i]=rotationv
        'VROT_ERR':Arrays[*,i]=errorv
        'PA':Arrays[*,i]=pav
        'PA_ERR':Arrays[*,i]=paverr
        'INCL':Arrays[*,i]=incv
        'INCL_ERR': Arrays[*,i]=incverr
        'XPOS':Begin
           IF NOT keyword_set(noconv) then begin
              tmp=sxpar(hed,'CRVAL2')+sxpar(hed,'CDELT2')*double(yposv)
              Arrays[*,i]=sxpar(hed,'CRVAL1')+sxpar(hed,'CDELT1')*double(xposv)/COS(tmp*!DtoR)
           ENDIF ELSE Arrays[*,i]=xposv
        end
        'XPOS_ERR':Begin
           IF NOT keyword_set(noconv) then begin
              Arrays[*,i]=ABS(xposverr*sxpar(hed,'CDELT1')*3600.)
           ENDIF ELSE Arrays[*,i]=xposverr
        end
        'YPOS':Begin
            IF NOT keyword_set(noconv) then begin
              Arrays[*,i]=sxpar(hed,'CRVAL2')+sxpar(hed,'CDELT2')*double(yposv)
           ENDIF ELSE Arrays[*,i]=yposv
        end
        'YPOS_ERR':Begin
           IF NOT keyword_set(noconv) then begin
              Arrays[*,i]=ABS(yposverr*sxpar(hed,'CDELT2')*3600.)
           ENDIF ELSE Arrays[*,i]=yposverr
        end
        else: Arrays[*,i]=-1
     ENDCASE
  ENDFOR

END
