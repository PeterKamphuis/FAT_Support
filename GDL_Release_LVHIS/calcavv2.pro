Pro calcavv2,Result1ch,Result2ch,Radius1ch,Radius2ch,error1=Errors1ch,error2=Errors2ch,rms=rms,average=av,mean=mean,rem_cent=rem_cent

Radius1=Radius1ch
print,Radius1ch,Radius2ch
Radius2=Radius2ch
Result1=Result1ch
Result2=Result2ch
Errors1=Errors1ch
Errors2=Errors2ch

print,Result1[*,0]
print,'some check'
print,Radius1
print,'more'
print,Radius2
print,'and more'
print,Result1[*,1]
tmp=1
tmp2=1
tmperr=1
tmperr2=1
trigerr2=0
trigerr1=0
IF n_elements(Errors2) GT 0 then trigerr2=1
IF n_elements(Errors1) GT 0 then trigerr1=1


IF (Result1[0,0] EQ 0. AND Result1[0,1] EQ 0.) OR keyword_set(rem_cent) then begin
 ;  IF  keyword_set(rem_cent) then stop
   tmp=Result1
   Result1=dblarr(n_elements(Result1[*,0])-1,2)
   Result1[*,0]=tmp[1:n_elements(tmp[*,0])-1,0]
   Result1[*,1]=tmp[1:n_elements(tmp[*,1])-1,1]
   tmp2=Radius1
   Radius1=dblarr(n_elements(Radius1)-1)
   Radius1[*]=tmp2[1:n_elements(tmp2)-1]
   tmp=Errors1
   Errors1=dblarr(n_elements(Errors1[*,0])-1,2)
   Errors1[*,0]=tmp[1:n_elements(tmp[*,0])-1,0]
   Errors1[*,1]=tmp[1:n_elements(tmp[*,1])-1,1]

ENDIF
tmp=1
tmp2=1
print,trigerr1,trigerr2
print,"So confusing"
print,Radius1
print,"WTf"
print,Radius2
IF trigerr1 then begin
   IF trigerr2 then begin
      case 1 of
         n_elements(Result1[0,*]) GT 1 AND n_elements(Errors1[0,*]) GT 1: begin
            IF (Radius1[0]-Radius1[1]) NE   (Radius2[0]-Radius2[1]) OR n_elements(Radius1) NE n_elements(Radius2) then begin
               interpolate,Result2[*],Radius2,output=tmp,newradii=Radius1
               interpolate,Errors2[*],Radius2,output=tmperr,newradii=Radius1
            ENDIF ELSE BEGIN
               tmp=Result2[*]
               tmperr=Errors2[*]
            ENDELSE

            diff=[(tmp[*]-Result1[*,0])/SQRT(Errors1[*,0]^2+tmperr[*]^2),(tmp[*]-Result1[*,1])/SQRT(Errors1[*,0]^2+tmperr[*]^2)]
            diff2=[(tmp[*]-Result1[*,0]),(tmp[*]-Result1[*,1])]
            av=TOTAL(ABS(diff))/TOTAL(1./ABS([SQRT(Errors1[*,0]^2+tmperr[*]^2),SQRT(Errors1[*,1]^2+tmperr[*]^2)]))
                                ; rms=SIGMA(diff2)
                                ;           rms=TOTAL(ABS([tmperr,tmperr2]))/n_elements([tmperr,tmperr2])/SQRT(n_elements([tmperr,tmperr2]))
            rms=SQRT((TOTAL(Errors1)/n_elements(Errors1))^2+(TOTAL(tmperr)/n_elements(tmperr))^2)/SQRT(n_elements([Errors1[*,0],Errors1[*,1]]))
            rms=TOTAL([SQRT(Errors1[*,0]^2+tmperr^2),SQRT(Errors1[*,1]^2+tmperr^2)])/n_elements([SQRT(Errors1[*,0]^2+tmperr^2),SQRT(Errors1[*,1]^2+tmperr^2)])/SQRT(n_elements([SQRT(Errors1[*,0]^2+tmperr^2),SQRT(Errors1[*,1]^2+tmperr^2)]))
;            IF rms LT 0.05 then rms=ABS(ABS((TOTAL([Errors1,tmperr])/n_elements([Errors1,tmperr]))+ABS(diff2[0]))/2.)
         end
         n_elements(Result1[0,*]) EQ 1 AND n_elements(Errors1[0,*]) EQ 1: begin
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
