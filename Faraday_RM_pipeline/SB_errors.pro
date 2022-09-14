pro errors, source

spawn, 'ls /students/sbetti/12B-336/IDL_Scripts/Fits_Images/SMITH_0' + strtrim(source,2) + '.fits', fnlist

; reads in SB_SOURCES.csv
;sed_data = READ_CSV('SB_SOURCES.csv', HEADER=SedHeader)
;print, SedHeader
;help, sed_data

; opens the .fits file for each source
;j=0
;FOR j=0,n_elements(sed_data.field1)-1 DO BEGIN
;   print,j
;   IF sed_data.field1[j] LT 1000 THEN BEGIN
;  spawn, 'ls /students/sbetti/12B-336/IDL_Scripts/Fits_Images/SMITH_00' + strtrim(sed_data.field1[j],2) + '.fits', fnlist
;  ENDIF ELSE BEGIN
;  spawn, 'ls /students/sbetti/12B-336/IDL_Scripts/Fits_Images/SMITH_0' + strtrim(sed_data.field1[j],2) + '.fits', fnlist
;  ENDELSE


i=0
!p.multi = 0
FOR i=0, n_elements(fnlist) - 1 DO BEGIN
  image = readfits(fnlist[i], hdr) 
  print, fnlist[i]

  ; finds the name
  name = strmid(fnlist[i], 54, 11)
  phi = findgen(400)*10 - 2000

  ; reads in readsfitsheader.pro
  freq = findgen(readfitsheader(hdr, 'NAXIS3'))*readfitsheader(hdr,'CDELT3') + readfitsheader(hdr, 'CRVAL3') 
  lambda = 3e8/freq 
  ; pulls Stokes Q, U from .fits file
  nf = n_elements(freq) - 1 
  Q=image[*,*,0:nf,1] 
  U=image[*,*,0:nf,2] 

  ; SIGMA LAMBDA^2
  totchan = 126.0
  
  lambdanot2 = 0.045260040  ;from rmclean.pro
  t=0
  summation = []
  FOR t=0, 125 DO BEGIN   
    summation = [summation, (lambda[t] - sqrt(lambdanot2))^2]    
  ENDFOR 
  summation = TOTAL(summation)
  sigmalamb = sqrt((summation /(totchan-1.0)))
 
  ; SIGMA CHI
  P = sqrt(Q^2 + U^2)  
  
  sigmachi = STDDEV(U) / (2*mean(P))
  
  ; SIGMA RM
  sigmaRM = sigmachi / (sigmalamb * sqrt(totchan-2))
  print, sigmaRM 


ENDFOR
;ENDFOR

RETURN
END

