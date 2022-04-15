Pro gettirific,NewFileName,VariableChange,Arrays,Errors=Errors


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
                               
 
;Lets make an default array with variables we want to transfer if it is not given
  IF NOT keyword_set(VariableChange) then begin
     IF keyword_set(errors) then begin
        VariableChange=['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT','VROT_ERR' , 'Z0',$
                        'SBR', 'INCL','INCL_ERR','PA','PA_ERR','XPOS','YPOS','VSYS','SDIS','VROT_2', 'VROT_2_ERR'  , 'Z0_2','SBR_2', 'INCL_2','INCL_2_ERR','PA_2','PA_2_ERR','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP']
     ENDIF ELSE BEGIN
        VariableChange=['BMIN','BMAJ','BPA','RMS','DISTANCE','NUR','RADI','VROT',  'Z0',$
                        'SBR', 'INCL','PA','XPOS','YPOS','VSYS','SDIS','VROT_2',  'Z0_2','SBR_2', 'INCL_2','PA_2','XPOS_2','YPOS_2','VSYS_2','SDIS_2','CONDISP']
     ENDELSE

  ENDIF
   
                                ;then open the previous fit ;or when looping back open the previous fit

  close,1
  
  h=' '
  openr,1,NewFileName
  rings=0
  WHILE rings EQ 0 Do begin
     readf,1,h
     tmp=str_sep(strtrim(strcompress(h),2),'=')
     IF tmp[0] EQ 'NUR' then begin
        Arrays=dblarr(fix(tmp[1]),n_elements(VariableChange))
        rings++
     ENDIF
  ENDWHILE
  close,1
  h=' '
  openr,1,NewFileName

  WHILE (NOT EOF(1)) DO BEGIN
     readf,1,h
     tmp=str_sep(strtrim(strcompress(h),2),'=')
     varpos = WHERE(tmp[0] EQ VariableChange)
     IF varpos[0] NE -1 then begin
      ;  print,tmp,n_elements(tmp)
        arr=str_sep(strtrim(strcompress(tmp[1]),2),' ')
        IF isnumeric(arr[0]) then begin
           Arrays[0:n_elements(arr)-1,varpos]=double(arr)
        ENDIF
     endif ELSE BEGIN
        IF keyword_set(Errors) then begin
           tmp1=str_sep(strtrim(strcompress(tmp[0]),2),' ')
           IF tmp1[0] EQ '#' then begin
            ;  print,tmp1,n_elements(tmp1)
              tmp3=WHERE(tmp1[1] EQ VariableChange)
              IF tmp3[0] NE -1 then begin
                 arr=str_sep(strtrim(strcompress(tmp[1]),2),' ')
                 IF isnumeric(arr[0]) then begin
              ;      print,tmp1[1],tmp3,arr
             ;       help,Arrays,n_elements(arr)-1,arr
                    IF n_elements(arr) GT n_elements(Arrays[*,0]) then arr=[0,0]
                    Arrays[0:n_elements(arr)-1,tmp3]=double(arr)
                 ENDIF
              ENDIF
           ENDIF
        ENDIF
     ENDELSE
  ENDWHILE
  close,1  
END
