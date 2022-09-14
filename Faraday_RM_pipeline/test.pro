pro test,source

; run each source individually
spawn, 'ls /students/sbetti/12B-336/IDL_Scripts/Fits_Images/SMITH_0' + strtrim(source,2) + '.fits', fnlist


; run all files in folder
;spawn, 'ls /students/sbetti/12B-336/IDL_Scripts/Fits_Images/*',fnlist


;----find the X and Y location of the brighest pixel in Stokes I-----

;xval = []
;yval = []
;array = readfits(fnlist[0])
;t=0
;FOR t=0, 125 DO BEGIN
;  mx = max(array[*,*,t,0], location)
;  ind = array_indices(array[*,*,t,0], location)
;  xval = [xval, ind[0]]
;  yval = [yval, ind[1]]
;ENDFOR
;print, xval
;no0x = xval[where(xval NE 0)]
;no0y = yval[where(yval NE 0)]
;voidx = max(histogram(no0x, OMIN=mnx), mxposx)
;voidy = max(histogram(no0y, OMIN=mny), mxposy)
;modex = mnx + mxposx
;modey = mny + mxposy
;print, modex, modey
;-------------------------------------------------------


i=0
;!p.multi = [0, 3, 4]
!p.multi = 0
FOR i=0, n_elements(fnlist) - 1 DO BEGIN
  image = readfits(fnlist[i], hdr) 
  print, fnlist[i]
  name = strmid(fnlist[i], 54, 11)
  phi = findgen(400)*10 - 2000 
  freq = findgen(readfitsheader(hdr, 'NAXIS3'))*readfitsheader(hdr,'CDELT3') + readfitsheader(hdr, 'CRVAL3') 
  lambda = 3e8/freq 
  nf = n_elements(freq) - 1 
  Q=image[*,*,0:nf,1] 
  U=image[*,*,0:nf,2] 

  ;--------RM SYNTHESIS----------------------
  rmsynthesis, transpose(q), transpose(u), lambda[0:nf]^2, phi, fdf_cube 

  ;rmclean, fdf_cube[*,modex,modey], phi, lambda[0:nf]^2, 0.0005, fdf_cleaned, fdf_residuals, clean_components 
  
  ;rmclean, fdf_cube[*,64,64], phi, lambda[0:nf]^2, 0.0005, fdf_cleaned, fdf_residuals, clean_components 
 
  cube_to_clean=total(total(fdf_cube[*,58:59,56:57],2),2)  
  ;rmclean, cube_to_clean, phi, lambda[0:nf]^2, 0.0005, fdf_cleaned, fdf_residuals, clean_components 
  
  ;---------PLOT------------------------------
  rmarr = float(sqrt(fdf_cleaned*conj(fdf_cleaned))) 
  plot, phi, rmarr, title=name, xtitle='phi'

  ; PLOT POLARIZATION(Wavelength)
  ; POLARIZATION FUNCTION
  X = (1.0/2.0) * atan(U,Q)
  greeklamb = '!4' + string("153B) + '!X'
  greekchi = '!4' + string("166B) + '!X'
  
  ; LAMBDA^2 FOR FITTING
  lsq = findgen(100) / 1000.0 +0.001
 
  ; PLOT POLARIZATION VS WAVELENTH 
  plot, lambda^2, X[64,64,*], title=name+', Polarization vs '+greeklamb+'!E2'+'!N, RM = '+strtrim(rotation_measure,2), xtitle=greeklamb + '!E2!N [m!E2!N]', ytitle = greekchi + ' (rad)',  font = -2
  
  ; MANUALLY ENTER AMPLITUDE CORRECTION
  B=''
  READ, B, prompt='enter amplitude correction: '
  ; OVERPLOT LINEAR FIT
  oplot, lsq, lsq * fix(strtrim(rotation_measure[0],2)) + B, psym = 2
  

  ;-------SAVE PLOTS---------------------------
  ; Rotation Measure

  ;set_plot, 'PS'
  ;device, filename='RM_Plots/' + name + '.ps'
  ;plot, phi, rmarr, title=name, xtitle='phi'
  ;device, /close
  ;set_plot, 'X'

  ;Polarization(wave)

  ;set_plot, 'PS'
  ;device, filename= '/students/sbetti/12B-336/IDL_Scripts/PolWav/'+ name + '.ps'
  ;plot, lambda^2, X[64,64,*], title=name+', Polarization vs '+greeklamb+'!E2'+'!N, RM = '+strtrim(rotation_measure,2), xtitle=greeklamb + '!E2!N [m!E2!N]', ytitle = greekchi + ' (rad)',  font = -2
  ;oplot, lsq, lsq * fix(strtrim(rotation_measure[0],2)) + B, psym = 2
  ;device, /close
  ;set_plot, 'X'

  ;------SIGNAL/NOISE + RM-------------------------------
  standard_dev=STDDEV(rmarr)
  signal_to_noise = MAX(rmarr) / standard_dev
  rotation_measure = phi(where(rmarr eq max(rmarr)))
  
  ;---------ERROR-----------------------
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
  
  print, rotation_measure, sigmaRM, signal_to_noise

ENDFOR
RETURN

PI = abs(transpose(fdf_cube))
print, PI
phi = (2*sqrt(3))/


;writefits, 'PI.fits', PI

END

