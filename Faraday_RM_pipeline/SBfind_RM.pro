pro SBfind_RM
  

;RESTORE RM_SUBTRACTED.SAV -- to find GB, GL
restore, 'rm_subtracted.sav', /verbose
;help, rm, /struct

; list of RA and Dec for each source -- to find GB, GL
RA = []
Dec = []

; reads in SB_SOURCES.csv -- list of sources
sed_data = READ_CSV('SB_SOURCES.csv', HEADER=SedHeader)

;OPEN TEXT FILE WITH SOURCES + COORDINATES -- to find GB, GL
openr, lun, '12B-336_all.txt', /get_lun
array2 = ''
line= '' 
;CREATE ARRAY OUT OF SOURCES -- to find GB, GL
WHILE NOT EOF(lun) DO BEGIN
readf, lun, line
array2 = [array2, line]
ENDWHILE 
free_lun, lun


; opens the .fits file for each source from SB_SOURCES.csv
j=0
FOR j=0,n_elements(sed_data.field1)-1 DO BEGIN
   print,j
   IF sed_data.field1[j] LT 1000 THEN BEGIN
  spawn, 'ls /students/sbetti/12B-336/IDL_Scripts/Fits_Images/SMITH_00' + strtrim(sed_data.field1[j],2) + '.fits', fnlist
  source_name = '0' + strtrim(sed_data.field1[j],2)
  ENDIF ELSE BEGIN
  spawn, 'ls /students/sbetti/12B-336/IDL_Scripts/Fits_Images/SMITH_0' + strtrim(sed_data.field1[j],2) + '.fits', fnlist
  source_name = strtrim(sed_data.field1[j],2)
  ENDELSE

; MATCH UP SOURCES AND COORDINATES -- to find GB,GL
information = array2[where(strmatch(array2, 'SMITH_0'+source_name+'*', /FOLD_CASE) EQ 1)]
;ISOLATE RA AND DEC INTO INDIVIDUAL LIST -- to find GB, GL
coordinates2 = strmid(information, 41, 27)
RA2 = strmid(coordinates2, 0,13)
Dec2 = strmid(coordinates2, 14,13)
RA = [RA, RA2]
Dec = [Dec, Dec2]


; START FINDING EVERYTHING...
i=0
!p.multi = 0
FOR i=0, n_elements(fnlist) - 1 DO BEGIN
  image = readfits(fnlist[i], hdr) 
  print, fnlist[i]
  ; finds the name
  name = strmid(fnlist[i], 49, 11)
  phi = findgen(400)*10 - 2000

  ; reads in readsfitsheader.pro
  freq = findgen(readfitsheader(hdr, 'NAXIS3'))*readfitsheader(hdr,'CDELT3') + readfitsheader(hdr, 'CRVAL3') 
 
  lambda = 3e8/freq 
  ; pulls Stokes Q, U from .fits file
  nf = n_elements(freq) - 1 
  Q=image[*,*,0:nf,1] 
  U=image[*,*,0:nf,2] 

  ; RM SYNTHESIS
  rmsynthesis, transpose(q), transpose(u), lambda[0:nf]^2, phi, fdf_cube 
  
  aa = sed_data.field2[j]
  bb = strsplit(aa, ':', /EXTRACT)
  cc = strsplit(bb, ',', /EXTRACT)

  IF n_elements(cc) eq 3 THEN BEGIN
  cube_to_clean=total(total(fdf_cube[*,cc[0]:cc[1,0],cc[1,1]:cc[2]],2),2)
  rmclean, cube_to_clean, phi, lambda[0:nf]^2, 0.0005, fdf_cleaned, fdf_residuals, clean_components 
  ENDIF ELSE BEGIN
  rmclean, fdf_cube[*,cc[0],cc[1]], phi, lambda[0:nf]^2, 0.0005, fdf_cleaned, fdf_residuals, clean_components 
  ENDELSE
;-------------------------------------------------------------------  
  ; PLOT 
  rmarr = float(sqrt(fdf_cleaned*conj(fdf_cleaned))) 
  ;plot, phi, rmarr, title=name, xtitle='phi'

  ; SAVE PLOT as .PS
  ;set_plot, 'PS'
  ;device, filename='RM_Plots/' + name + '.ps'
  ;plot, phi, rmarr, title=name, xtitle='phi'
  ;device, /close
  ;set_plot, 'X'

  ; FIND SIGNAL/NOISE + RM
  standard_dev=STDDEV(rmarr)
  signal_to_noise = MAX(rmarr) / standard_dev
  rotation_measure = phi(where(rmarr eq max(rmarr)))
  ;print, rotation_measure,signal_to_noise

;--------------------------------------------------------------------
  ; PLOT POLARIZATION(Wavelength)
  ; POLARIZATION FUNCTION
  ;X = (1.0/2.0) * atan(U,Q)
  ;greeklamb = '!4' + string("153B) + '!X'
  ;greekchi = '!4' + string("166B) + '!X'
  ; SAVE PLOT
  ;set_plot, 'PS'
  ;device, filename= '/students/sbetti/12B-336_2016/IDL_Scripts/PolWav/'+ name + '.ps'
  ; LAMBDA^2 FOR FITTING
  ;lsq = findgen(100) / 1000.0 +0.001
  ; USE MEAN OF BOX TO CALCULATE POLARIZATION
  ;IF n_elements(cc) eq 3 THEN BEGIN
  ;low = mean([fix(cc[0]), fix(cc[1,0])])
  ;high = mean([fix(cc[1,1]), fix(cc[2])]) 
  ; PLOT POLARIZATION VS WAVELENTH 
  ;plot, lambda^2, X[low,high,*], title=name+', Polarization vs '+greeklamb+'!E2'+'!N, RM = '+strtrim(rotation_measure,2), xtitle=greeklamb + '!E2!N [m!E2!N]', ytitle = greekchi + ' (rad)',  font = -2
  ;help, strtrim(rotation_measure[0],2)  
  ; MANUALLY ENTER AMPLITUDE CORRECTION
  ;B=''
  ;READ, B, prompt='enter amplitude correction: '
  ; OVERPLOT LINEAR FIT
  ;B= 1.0
  ;oplot, lsq, lsq * fix(strtrim(rotation_measure[0],2)) + B, psym = 2
  ;ENDIF ELSE BEGIN
  ;ENDELSE
  ;device, /close
  ;set_plot, 'X'
;-----------------------------------------------------------------


  ; ERROR ANALYSIS
  ; SIGMA LAMBDA^2
  totchan = 126.0
  ;r=0
  ;channelwidthlamb = [] 
  ;FOR r=0, 125 DO BEGIN
  ;  channelwidthlamb = [channelwidthlamb, lambda[r]^2 - lambda[r-1]^2]     
  ;ENDFOR
  ;den = TOTAL(channelwidthlamb[1:*]^2,/double)
  ;num = TOTAL(channelwidthlamb[1:*]^2 * lambda^2, /double)
  ;lambdanot2 = num/den

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

  ;CONVERT RA AND DEC TO GL AND GB
  glactc, ten(RA[j]), ten(Dec[j]), 2000, gl, gb, 1
  print, name
  print, RA[j], Dec[j]
  print, GL, GB

  ;SAVE VALUES TO .DAT FILE
  fname = 'RM_values.txt'
  openw,2,fname, /append
  printf, 2, string(name), rotation_measure[0], sigmaRM, signal_to_noise, float(GL), float(GB)
  close,2


ENDFOR
ENDFOR

RETURN
END
