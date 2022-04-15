Pro calcav_vobs,Resultvr1ch,Resultinc1ch,Resultvr2ch,Resultinc2ch,Radius1ch,Radius2ch,errorvr1=Errorsvr1ch,errorvr2=Errorsvr2ch,errorinc1=Errorsinc1ch,errorinc2=Errorsinc2ch,rms=rms,average=av,mean=mean,rem_cent=rem_cent

Radius1=Radius1ch
Radius2=Radius2ch
Resultvr1=Resultvr1ch
Resultvr2=Resultvr2ch
Errorsvr1=Errorsvr1ch
Errorsvr2=Errorsvr2ch
Resultinc1=Resultinc1ch
Resultinc2=Resultinc2ch
Errorsinc1=Errorsinc1ch
Errorsinc2=Errorsinc2ch

print,Resultvr1[*,0]
print,'some check'
print,Radius1
print,'more'
print,Radius2
print,'and more'
print,Resultvr1[*,1]
tmp=1
tmpinc=1
tmpincerr=1
tmp2=1
tmperr=1
tmperr2=1
trigerr2=0
trigerr1=0
IF n_elements(Errorsvr2) GT 0 then trigerr2=1
IF n_elements(Errorsvr1) GT 0 then trigerr1=1


;IF (Resultvr1[0,0] EQ 0. AND Resultvr1[0,1] EQ 0.) OR keyword_set(rem_cent) then begin
 ;  IF  keyword_set(rem_cent) then stop
   tmp=Resultvr1
   Resultvr1=dblarr(n_elements(Resultvr1[*,0])-1,2)
   Resultvr1[*,0]=tmp[1:n_elements(tmp[*,0])-1,0]
   Resultvr1[*,1]=tmp[1:n_elements(tmp[*,1])-1,1]
   tmp=Resultinc1
   Resultinc1=dblarr(n_elements(Resultinc1[*,0])-1,2)
   Resultinc1[*,0]=tmp[1:n_elements(tmp[*,0])-1,0]
   Resultinc1[*,1]=tmp[1:n_elements(tmp[*,1])-1,1]
   tmp=Radius1
   Radius1=dblarr(n_elements(Radius1)-1)
   Radius1[*]=tmp[1:n_elements(tmp)-1]
   tmp=Errorsvr1
   Errorsvr1=dblarr(n_elements(Errorsvr1[*,0])-1,2)
   Errorsvr1[*,0]=tmp[1:n_elements(tmp[*,0])-1,0]
   Errorsvr1[*,1]=tmp[1:n_elements(tmp[*,1])-1,1]
   tmp=Errorsinc1
   Errorsinc1=dblarr(n_elements(Errorsinc1[*,0])-1,2)
   Errorsinc1[*,0]=tmp[1:n_elements(tmp[*,0])-1,0]
   Errorsinc1[*,1]=tmp[1:n_elements(tmp[*,1])-1,1]

;ENDIF
  
print,trigerr1,trigerr2
IF trigerr1 then begin
   IF trigerr2 then begin
      case 1 of
         n_elements(Resultvr1[0,*]) GT 1 AND n_elements(Errorsvr1[0,*]) GT 1: begin
            print,Radius1
            IF (Radius1[0]-Radius1[1]) NE (Radius2[0]-Radius2[1])   OR n_elements(Radius1) NE n_elements(Radius2) then begin
               interpolate,Resultvr2[*],Radius2,output=tmp,newradii=Radius1	       
               interpolate,Errorsvr2[*],Radius2,output=tmperr,newradii=Radius1	 
               interpolate,Resultinc2[*],Radius2,output=tmpinc,newradii=Radius1	       
               interpolate,Errorsinc2[*],Radius2,output=tmpincerr,newradii=Radius1	     
            ENDIF ELSE BEGIN
               tmp=Resultvr2[*]
               tmperr=Errorsvr2[*]
               tmpin=Resultinc2[*]
               tmpincerr=Errorsinc2[*]
            ENDELSE
            Errorsontir=[[SQRT(SIN(Resultinc1[*,0]*!DtoR)^2*Errorsvr1[*,0]^2+Resultvr1[*,0]^2*COS(Resultinc1[*,0]*!DtoR)^2*Errorsinc1[*,0]^2)],[SQRT(SIN(Resultinc1[*,1]*!DtoR)^2*Errorsvr1[*,1]^2+Resultvr1[*,1]^2*COS(Resultinc1[*,1]*!DtoR)^2*Errorsinc1[*,1]^2)]]
            ErrorRC=[SQRT(SIN(tmpinc*!DtoR)^2*tmperr^2+tmp^2*COS(tmpinc*!DtoR)^2*tmpincerr^2)]

    ;  calcavv2,[[TirResult2[*,6]*SIN(Tirresult2[*,4]*!DtoR)],[TirResult2[*,10]*SIN(Tirresult2[*,8]*!DtoR)]],RCresult[*,3]*SIN(RCresult[*,7]*!DtoR),radii,RCresult[*,0],average=av,mean=mean,rms=rms,error1=Errorontir,error2=ErrorRC,/rem_cent            



            diff=[((tmp[*]*SIN(tmpinc*!DtoR)-Resultvr1[*,0]*SIN(RESULTinc1[*,0]*!DtoR))/SQRT(Errorsontir[*,0]^2+ErrorRC[*]^2)),((tmp[*]*SIN(tmpinc*!DtoR)-Resultvr1[*,1]*SIN(RESULTinc1[*,1]*!DtoR))/SQRT(Errorsontir[*,1]^2+ErrorRC[*]^2))]   
            diff2=[tmp[*]*SIN(tmpinc*!DtoR)-Resultvr1[*,0]*SIN(RESULTinc1[*,0]*!DtoR),tmp[*]*SIN(tmpinc*!DtoR)-Resultvr1[*,1]*SIN(RESULTinc1[*,1]*!DtoR)]
         ;   diff=diff2
            av=Total(ABS(diff))/TOTAL(1./ABS([SQRT(Errorsontir[*,0]^2+ErrorRC[*]^2),SQRT(Errorsontir[*,1]^2+ErrorRC[*]^2)]))/SQRT(n_elements(diff))
                                ; rms=SIGMA(diff2)
                                ;           rms=TOTAL(ABS([tmperr,tmperr2]))/n_elements([tmperr,tmperr2])/SQRT(n_elements([tmperr,tmperr2]))
        ;    rms=SQRT((TOTAL(Errors1)/n_elements(Errors1))^2+(TOTAL(tmperr)/n_elements(tmperr))^2)/SQRT(n_elements([Errors1[*,0],Errors1[*,1]]))
            rms=TOTAL([SQRT(Errorsontir[*,0]^2+ErrorRC^2),SQRT(Errorsontir[*,1]^2+ErrorRC^2)])/n_elements(Errorsontir)/SQRT(n_elements(Errorsontir));            IF rms LT 0.05 then rms=ABS(ABS((TOTAL([Errors1,tmperr])/n_elements([Errors1,tmperr]))+ABS(diff2[0]))/2.)
         end
         n_elements(Resultvr1[0,*]) EQ 1 AND n_elements(Errorsvr1[0,*]) EQ 1: begin
            stop
             IF (Radius1[0]-Radius1[1]) NE   (Radius2[0]-Radius2[1])   OR n_elements(Radius1) NE n_elements(Radius2) then begin
               interpolate,Result2[*],Radius2,output=tmp,newradii=Radius1	       
               interpolate,Errors2[*],Radius2,output=tmperr,newradii=Radius1	     
            ENDIF ELSE BEGIN
               tmp=Result2[*]
               tmperr=Errors2[*]
            ENDELSE

            diff=[(tmp[*]-Result1[*])/SQRT(Errors1[*]^2+tmperr[*]^2)]   
            diff2=[(tmp[*]-Result1[*])]
            av=TOTAL(ABS(diff))/TOTAL(1./ABS([SQRT(Errors1[*]^2+tmperr[*]^2)]))
                                ; rms=SIGMA(diff2)
                                ;           rms=TOTAL(ABS([tmperr,tmperr2]))/n_elements([tmperr,tmperr2])/SQRT(n_elements([tmperr,tmperr2]))
          ;  rms=SQRT((TOTAL(Errors1)/n_elements(Errors1))^2+(TOTAL(tmperr)/n_elements(tmperr))^2)/SQRT(n_elements([Errors1[*,0],Errors1[*,1]]))
               rms=TOTAL([SQRT(Errors1^2+tmperr^2)])/n_elements([SQRT(Errors1[*]^2+tmperr^2)])/SQRT(n_elements([SQRT(Errors1[*]^2+tmperr^2)]))
         end
         ELSE:begin
            print,'yeah need to code this up'
         end
      endcase
   endif else begin
      stop
      case 1 of
         n_elements(Result1[0,*]) GT 1 AND n_elements(Errors1[0,*]) GT 1: begin
            IF (Radius1[0]-Radius1[1]) NE   (Radius2[0]-Radius2[1])   OR n_elements(Radius1) NE n_elements(Radius2) then begin
               interpolate,Result2[*],Radius2,output=tmp,newradii=Radius1	         
            ENDIF ELSE BEGIN
               tmp=Result2[*]
            ENDELSE

            diff=[(tmp[*]-Result1[*,0])/SQRT(Errors1[*,0]^2),(tmp[*]-Result1[*,1])/SQRT(Errors1[*,0]^2)]   
            diff2=[(tmp[*]-Result1[*,0]),(tmp[*]-Result1[*,1])]
            av=TOTAL(ABS(diff))/TOTAL(1./ABS([SQRT(Errors1[*,0]^2),SQRT(Errors1[*,1]^2)]))
                                ; rms=SIGMA(diff2)
                                ;           rms=TOTAL(ABS([tmperr,tmperr2]))/n_elements([tmperr,tmperr2])/SQRT(n_elements([tmperr,tmperr2]))
            rms=TOTAL(Errors1)/n_elements(Errors1)/SQRT(n_elements(Errors1))
;            IF rms LT 0.05 then rms=ABS(ABS((TOTAL([Errors1,tmperr])/n_elements([Errors1,tmperr]))+ABS(diff2[0]))/2.)
         end
         n_elements(Result1[0,*]) EQ 1 AND n_elements(Errors1[0,*]) EQ 1: begin
             IF (Radius1[0]-Radius1[1]) NE   (Radius2[0]-Radius2[1])   OR n_elements(Radius1) NE n_elements(Radius2) then begin
               interpolate,Result2[*],Radius2,output=tmp,newradii=Radius1	       
            ENDIF ELSE BEGIN
               tmp=Result2[*]
             
            ENDELSE

            diff=[(tmp[*]-Result1[*])/SQRT(Errors1[*]^2)]   
            diff2=[(tmp[*]-Result1[*])]
            av=TOTAL(ABS(diff))/TOTAL(1./ABS([SQRT(Errors1[*]^2)]))
                                ; rms=SIGMA(diff2)
                                ;           rms=TOTAL(ABS([tmperr,tmperr2]))/n_elements([tmperr,tmperr2])/SQRT(n_elements([tmperr,tmperr2]))
            rms=TOTAL(Errors1)/n_elements(Errors1)/SQRT(n_elements(Errors1))
               
         end
         ELSE:begin
            print,'yeah need to code this up'
         end
      endcase
   endelse
endif




mean=MEDIAN(diff2)
end
    
