Pro convert_ps,filename,resolution=resolution,delete=delete,trim=trim,png=png,transparent=transparent,negate=negate

;Program to convert large ps files to smaller files without a significant loss of quality
;note that it uses imconvert which is actually convert but because of
;a clash with kvis it is renamed to imconvert
compile_opt idl2
 CATCH,Error_status
  IF  Error_status NE 0. THEN BEGIN
     print, 'Oops the following went wrong:'
     print, !ERROR_STATE.MSG
     print,'Use convert_ps in this way:filename,directory,resolution=resolution,/DELETE'
     print, "CALLING SEQUENCE: "
     print,  '/DELETE deletes the original file, if not set the new file will have the ending filename_converted.ps'
     print, 'filename= name of the file to be converted'
     print, 'resolution= the resolution of the outpu image [300]'
     print, '/TRIM the edges to make the image as smaal as possible'
     print, '/png to get a png image instead of PS'
     print, '/transparent to get a transparent backgroun in a png image'
     goto,ending
  endif
directory=''
IF NOT n_elements(resolution) then resolution=300.

dircheck=str_sep(filename,'/')
IF n_elements(dircheck) GT 1. then begin
   directory=dircheck[0]+'/'
   for i=1,n_elements(dircheck)-2 do begin
      directory=directory+dircheck[i]+'/'
   endfor
   filename=dircheck[n_elements(dircheck)-1]
   CD,directory,CURRENT=Old_Dir
endif

checkfile=FILE_TEST(filename)
IF not checkfile then begin
   print,'The file '+directory+filename+' cannot be found'
   print,'Aborting'
   IF n_elements(dircheck) GT 1. then begin
      CD,Old_Dir
   endif
   goto,ending
ENDIF
split=str_sep(filename,'.')
IF STRUPCASE(split[n_elements(split)-1]) NE 'PS' AND  STRUPCASE(split[n_elements(split)-1]) NE 'EPS' AND STRUPCASE(split[n_elements(split)-1]) NE 'PDF' then begin
   print,'The file '+directory+filename+' is not a ps or pdf file'
   print,'Aborting'
   IF n_elements(dircheck) GT 1. then begin
      CD,Old_Dir
   endif
   goto,ending
ENDIF
IF n_elements(split) GT 2. then begin
   base=split[0]
    for i=1,n_elements(split)-2 do begin
      base=base+'.'+split[i]
   endfor
 ENDIF else base=split[0]

spawn,'gs -r'+strtrim(strcompress(string(resolution,format='(i10)')),2)+' -dEPSCrop -dTextAlphaBits=4 -sDEVICE=png16m -sOutputFile="'+base+'.png" -dBATCH -dNOPAUSE "'+filename+'"'


converted=0
spawn,'convert --help',result,notfound
IF n_elements(result) GT 1 then begin
   gf=strsplit(result[0],' ',/extract)
   IF n_elements(gf) GT 1 then begin
      IF strtrim(gf[1],2) EQ 'ImageMagick' then begin
         converted=1
         conv_prog = 'convert'
         IF keyword_set(trim) then spawn,'convert "'+base+'.png" -trim "'+base+'.png"'
      ENDIF
   ENDIF
ENDIF
IF converted LT 1 then begin
   spawn,'/usr/bin/convert --help',result,notfound
   IF n_elements(result) GT 1 then begin
      gf=strsplit(result[0],' ',/extract)
      IF n_elements(gf) GT 1 then begin
         IF strtrim(gf[1],2) EQ 'ImageMagick' then begin
            converted=2
            conv_prog = '/usr/bin/convert'
            IF keyword_set(trim) then spawn,'/usr/bin/convert "'+base+'.png" -trim "'+base+'.png"'
         ENDIF
      ENDIF
   ENDIF
ENDIF
IF converted LT 1 then begin
   spawn,'imconvert --help',result,notfound
   IF n_elements(result) GT 1 then begin
      gf=strsplit(result[0],' ',/extract)
      IF n_elements(gf) GT 1 then begin
         IF strtrim(gf[1],2) EQ 'ImageMagick' then begin
            converted=3
            conv_prog = 'imconvert'
            IF keyword_set(trim) then spawn,'imconvert "'+base+'.png" -trim "'+base+'.png"'
         ENDIF
      ENDIF
   ENDIF
ENDIF
IF converted LT 1 then begin
   spawn,'convert-im6 --help',result,notfound
   IF n_elements(result) GT 1 then begin
      gf=strsplit(result[0],' ',/extract)
      IF n_elements(gf) GT 1 then begin
         IF strtrim(gf[1],2) EQ 'ImageMagick' then begin
            converted=4
            conv_prog = 'convert-im6'
            IF keyword_set(trim) then spawn,'convert-im6 "'+base+'.png" -trim "'+base+'.png"'
         ENDIF
      ENDIF
   ENDIF
ENDIF



if keyword_set(negate) then begin
   IF keyword_set(delete) then begin
      spawn,'rm -f '+filename
      if not keyword_set(png) then spawn,conv_prog+' "'+base+'.png" -negate eps3:"'+base+'.ps"'
   ENDIF else begin
      if not keyword_set(png) then spawn,conv_prog+' "'+base+'.png" -negate eps3:"'+base+'_converted.ps"' else spawn,'magick "'+base+'.png" -negate eps3:"'+base+'.ps"'
                                ;Sometimes the colorsget inverted due
                                ;to this command. If this haappens add
                                ;-negate before eps3:
   ENDELSE
ENDIF ELSE BEGIN
   IF keyword_set(delete) then begin
      spawn,'rm -f '+filename
      if not keyword_set(png) then spawn,conv_prog+' "'+base+'.png" eps3:"'+base+'.ps"'
   ENDIF else begin
      if not keyword_set(png) then spawn,conv_prog+' "'+base+'.png" eps3:"'+base+'_converted.ps"' else spawn,'magick "'+base+'.png" eps3:"'+base+'.ps"'
                                ;Sometimes the colorsget inverted due
                                ;to this command. If this haappens add
                                ;-negate before eps3:
   ENDELSE
ENDELSE
if not keyword_set(png) then spawn,'rm -f  "'+base+'.png"' else  begin
   IF keyword_set(delete) then spawn,'rm -f  "'+base+'.ps"'
ENDELSE
if keyword_set(png) AND keyword_set(transparent) then spawn,conv_prog+' "'+base+'.png" -transparent white "'+base+'.png"'


IF n_elements(dircheck) GT 1. then begin
 CD,Old_Dir
endif



ending:
end
