pro SB_findTaylor



; RESTORE RM_SUBTRACTED.SAV
restore, 'rm_subtracted.sav', /verbose
help, rm, /struct

RA = []
Dec = []
;READ IN TABLE OF SOURCES
sed_data = READ_CSV('SB_SOURCES.csv', HEADER=SedHeader)
print, SedHeader

;OPEN TEXT FILE WITH SOURCES + COORDINATES
openr, lun, '12B-336_all.txt', /get_lun
array2 = ''
line= '' 
;CREATE ARRAY OUT OF SOURCES
WHILE NOT EOF(lun) DO BEGIN
readf, lun, line
array2 = [array2, line]
ENDWHILE 
free_lun, lun

; SPLIT EACH ENTRY IN SOURCES ARRAY AND MATCH THE TABLE OF SOURCES WITH THEIR COORESPONDING ENTRY
j=0
FOR j=0,n_elements(sed_data.field1)-1 DO BEGIN
IF sed_data.field1[j] LT 1000 THEN BEGIN
source_name = '0' + strtrim(sed_data.field1[j],2)
ENDIF ELSE BEGIN
source_name = strtrim(sed_data.field1[j],2)
ENDELSE

; MATCH UP SOURCES AND COORDINATES
information = array2[where(strmatch(array2, 'SMITH_0'+source_name+'*', /FOLD_CASE) EQ 1)]

;ISOLATE RA AND DEC INTO INDIVIDUAL LIST
coordinates2 = strmid(information, 41, 27)
RA2 = strmid(coordinates2, 0,13)
Dec2 = strmid(coordinates2, 14,13)
RA = [RA, RA2]
Dec = [Dec, Dec2]
ENDFOR


; MANUALLY ENTER THE RA AND DEC AS GIVEN IN .SAV FILE 
;value = ''
;read, value, prompt = 'enter RA and Dec: '
;SEPARATES THE RA AND DECLINATION
;RA = strmid(value,0,13)
;Dec = strmid(value,14,13)

; DO NOT USE FOR LOOP IF ENTERING COORDINATES MANUALLY.  
w=0
FOR w=0, 228 DO BEGIN

; DELETE [W] IN RA AND DEC IF ENTERING COORDINATES MANUALLY
; CONVERT RA AND DEC TO GL AND GB
glactc, ten(RA[w]), ten(Dec[w]), 2000, gl, gb, 1

;takes the angular distance betwen points on a sphere
d = sphdist(GL, GB, rm.l, rm.b, /degree)

;find the smallest distance between the points
mn = min(d, location)

;finds the index of that point
ind = array_indices(d, location)

;print the rotation measure of the source 
IF sed_data.field1[w] LT 1000 THEN BEGIN
source_name = '0' + strtrim(sed_data.field1[w],2)
ENDIF ELSE BEGIN
source_name = strtrim(sed_data.field1[w],2)
ENDELSE
print, 'SMITH_0'+source_name

print, RA[w], Dec[w]
print, rm[ind].rm, rm[ind].erm
print, sixty(rm[ind].ra*24/360)
print, sixty(rm[ind].dec)
print, GL, GB


ENDFOR

RETURN

END



