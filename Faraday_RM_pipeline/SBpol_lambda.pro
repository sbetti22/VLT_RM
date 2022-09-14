pro pol_lambda, source, pixels, amplitude, RotMeasure


; FROM THE SOURCE, FINDS THE CORRESPONDING .FITS IMAGE
IF source LT 1000 THEN BEGIN
  spawn, 'ls /students/sbetti/12B-336/IDL_Scripts/Fits_Images/SMITH_00' + strtrim(source,2) + '.fits', fnlist
ENDIF ELSE BEGIN
  spawn, 'ls /students/sbetti/12B-336/IDL_Scripts/Fits_Images/SMITH_0' + strtrim(source,2) + '.fits', fnlist
ENDELSE

;spawn, 'ls /students/sbetti/12B-336/IDL_Scripts/Fits_Images/SMITH_0' + strtrim(source,2) + '.fits', fnlist

i=0
!p.multi = 0
FOR i=0, n_elements(fnlist) - 1 DO BEGIN
  image = readfits(fnlist[i], hdr) 
  print, fnlist[i]
  name = strmid(fnlist[i], 49, 11)
  phi = findgen(400)*10 - 2000
  freq = findgen(readfitsheader(hdr, 'NAXIS3'))*readfitsheader(hdr,'CDELT3') + readfitsheader(hdr, 'CRVAL3')  
  lambda = 3e8/freq 
  nf = n_elements(freq) - 1 
  Q=image[*,*,0:nf,1] 
  U=image[*,*,0:nf,2] 
  X = (1.0/2.0) * atan(U,Q)

  lsq = findgen(100) / 1000.0 +0.001

  ; FINDS THE MEAN VALUE OF THE PIXELS TO USE FOR THE POLARIZATION OF THE SOURCE
  ;aa = ''
  ;READ, aa, prompt='enter pixels: '
  bb = strsplit(pixels, ':', /EXTRACT)
  cc = strsplit(pixels, ',', /EXTRACT)
  low = mean([fix(cc[0]), fix(cc[1,0])])
  high = mean([fix(cc[1,1]), fix(cc[2])]) 

  ;PLOT POLARIZATION VS WAVELENGTH
  greeklamb = '!4' + string("153B) + '!X'
  greekchi = '!4' + string("166B) + '!X'

  ;set_plot, 'PS'
  ;device, filename= '/students/sbetti/12B-336/IDL_Scripts/PolWav/'+ name + '.ps'
  plot, lambda^2, X[low,high,*], title= name + ', Polarization vs ' + greeklamb + '!E2' + '!N, RM = ' + RotMeasure, xtitle=greeklamb + '!E2!N [m!E2!N]', ytitle = greekchi + ' (rad)',  font = -2
 
  ;plot, lambda^2, X[65,65,*], title= name + ', Polarization vs ' + greeklamb + '!E2' + '!N, RM = -35.0', xtitle=greeklamb + '!E2!N [m!E2!N]', ytitle = greekchi + ' (rad)',  font = -2

  ; OVERPLOT A LINEAR FIT
  oplot, lsq, lsq * (RotMeasure) + amplitude, psym=2  
  ;device, /close
  ;set_plot, 'X'


 
ENDFOR
RETURN
END



 
