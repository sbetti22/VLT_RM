pro rmclean, fdf_dirty, $
             phi, $
             lambdasq, $
             cutoff, $
             ; OUTPUTS...
             fdf_cleaned, $
             fdf_residuals, $
             clean_components, $
             ; KEYWORDS...
             WEIGHT=weight_in, $
             GAIN=gain, $
             MAXITER=maxiter, $
             ITERATIONS=iterations, $
             RMSF_IN=rmsf_in, $
             RMSF_OUT=rmsf_out, $
             LAMBDA0SQ=lambda0sq,$
             WATCH=watch, $
             DEROTATE=derotate

a = findgen(400)
b = findgen(400)+300
c = dcomplex(a,b)
phi = findgen(400)*10 - 2000
freq = findgen(126)*8e6 + 9.975e8
l = 3e8/freq
.r rmclean
rmclean, c, phi, l^2, 0.05, fdf_cleaned


file = 'test_data.txt'
data=fltarr(2,400) ; create the data array
openr,lun,file,/get_lun 
readf,lun,data 
close,/all 
x=reform(data(0,*)) 
y=reform(data(1,*))

c = dcomplex(x,y)
print 
phi = findgen(400)*10 - 2000
freq = findgen(126)*8e6 + 9.975e8
l = 3e8/freq
.r rmclean
rmclean, c, phi, l^2, 0.0005, fdf_cleaned

;+
; NAME:
;       RMCLEAN
;
; PURPOSE:
;       Deconvolve the Rotation Measure Spread Function from a Faraday 
;       dispersion function.
;
; CALLING SEQUENCE:
;       RMCLEAN, fdf_dirty, phi, lambdasq, cutoff, fdf_cleaned[,
;                /WATCH][,plot WEIGHT=vector][,
;                RMSF_IN=variable][,RMSF_OUT=variable][,
;                GAIN=scalar][,MAXITER=scalar][,ITERATIONS=variable][,
;                LAMBDA0SQ=variable]
;
; INPUTS:
;       FDF_DIRTY - dirty Faraday dispersion function ; complex vector
;
;       PHI - Faraday depth vector ; a floating point vector that must have
;             the same size as fdf_dirty ; there must be a 1-to-1
;             correspondence between the channels in fdf and the Faraday
;             depths in this array ; phi must be monotonically increasing
;             and regularly spaced
;
;       LAMBDASQ - the lambda-squared sampling of the original Stokes
;                  measurements in inverse meters squared; a floating point
;                  vector; the size of this vector is independent of the
;                  size of fdf and phi.  !!NOTA BENE!!: If channels in the
;                  original Q or U spectra have been blanked out, then the
;                  user must be sure to NOT pass in the value(s) of
;                  lambda-squared for these channels, otherwise the RMSF
;                  will be incorrect and the cleaning will fail.  The user
;                  can also choose to send in all of the original values of
;                  lambda-squared and just set the weight(s) of the
;                  corresponding blanked channels to zero via the WEIGHT
;                  keyword.
;
;       CUTOFF - the polarized intensity level to clean down to; this
;                should be set to something reasonable like 5 times the
;                polarized RMS level.
;
; KEYWORD PARAMETERS:
;       /WATCH - a window is created and the cleaning process is animated
;
;       WEIGHT = the weight function (a.k.a. the sampling function) for the 
;                complex polarized surface brightness; floating point vector 
;                that must have the same length as the LAMBDASQ array.  If
;                not set, a uniform weighting is assumed.
;
;       RMSF_OUT = pass out a structure that contains the complex RM spread
;                  function and the corresponding Faraday depth values.
;                  The RMCLEAN algorithm requires an RMSF with twice the
;                  number of channels of FDF_DIRTY. Set this to a named
;                  variable and the RMSF will be passed out via a structure
;                  of the form:
;
;                  {rmsf:dcomplexarr(2*n_fdf),phi_rmsf:dblarr(2*n_fdf)}
;
;       RMSF_IN = pass in a structure that contains the complex RM spread
;                 function and the corresponding Faraday depth values.  By
;                 the nature of the RMCLEAN algorithm, we need to have an
;                 RMSF that has twice the number of channels of
;                 FDF_DIRTY. The input must be a structure of the form:
;
;                 {rmsf:dcomplexarr(2*n_fdf),phi_rmsf:dblarr(2*n_fdf)}
;
;                 This allows the RMSF to be calculated once and returned
;                 via the RMSF_OUT keyword, and passed in via RMSF_IN in
;                 subsequent calls so that the RMSF doesn't get
;                 recalculated needlessly.
;
;       GAIN = the gain factor for the clean. Default is 0.1.
;
;       MAXITER = set this to the maximum number of iterations for the
;                 cleaning. Default is 1000.
;
;       ITERATIONS = set this to a named variable that will hold the number
;                    of iterations for convergence.
;
;       LAMBDA0SQ - set this optional keyword to a named variable to pass
;                   back the weighted mean lambda-squared.
;
;       /DEROTATE - derotate clean components and FDF residuals back to
;                   lambda=0
;
; OUTPUTS:
;       FDF_CLEANED - the cleaned Faraday dispersion function; a
;                     complex vector with the same length as FDF_DIRTY.
;
;       FDF_RESIDUALS - the residuals of the cleaned FDF from the dirty FDF
;                       in the weighted-mean lambda-squared domain; a
;                       complex vector with the same length as FDF_DIRTY.
;
;       CLEAN_COMPONENTS - the model clean components for the deconvolution;
;                          a complex vector with the same length as 
;                          FDF_DIRTY.
;
; SIDE EFFECTS:
;       If /WATCH is set, a window will be created on the XWINDOWS device.
;
; NOTES:
;       The user must be careful to remove any blank channels from the
;       input Stokes Q and U spectrum
;
; CAVEAT:
;       The FWHM of the RMSF that is always quoted as 2*sqrt(3) /
;       (max(lambda^2)-min(lambda^2)) is only valid for uniform weighting!
;       This needs to be adjusted to take account of the weight vector.
;       Nobody has built this into their RMCLEAN algorithms yet.
;
; PROCEDURE:
;       Based on the PEAKP scheme outlined in Appendix A of Heald et
;       al. (2009).
;
; PROCEDURES USED:
;       Robishaw's SETCOLORS.
;
; MODIFICATION HISTORY:
;       Written by Tim Robishaw, USyd 02 Jun 2010
;       Based on Ann Mao's CLEAN_RMS and George Heald's
;       MIRIAD code RMCLEAN.
;       TR 22 Oct 2012 cleaned up passing in/out RMSF an cleaned up
;       documentation.
;-


; GET THE ATTRIBUTES OF THE DIRTY FARADAY DISPERSION FUNCTION...
sz = size(fdf_dirty)

; MAKE SURE INPUT FARADAY DISPERSION FUNCTION IS ONE DIMENSIONAL...
if (sz[0] ne 1) then $
   message, 'FDF_DIRTY must be 1 dimensional.'

; MAKE SURE INPUT FARADAY DISPERSION FUNCTION IS COMPLEX...
fdf_type = sz[sz[0]+1]
if (fdf_type ne 6) AND (fdf_type ne 9) then $
   message, 'FDF_DIRTY must be a complex array'

; WHAT IS THE LENGTH OF FARADAY DISPERSION FUNCTION...
n_fdf = sz[sz[0]+2]

; IF WE'VE PASSED IN THE RMSF, THEN WE HAVE TO MAKE SURE IT'S PROPERLY
; FORMATTED: IT MUST HAVE TWICE THE NUMBER OF CHANNELS AS FDF_DIRTY AND
; MUST BE COMPLEX...
n_rmsf_in = N_elements(rmsf_in)
if (n_rmsf_in gt 0) then begin

   ; MAKE SURE RMSF IS A STRUCTURE...
   if size(RMSF_IN,/TYPE) ne 8 then $
      message, 'RMSF_IN must be a structure'

   ; MAKE SURE RMSF STRUCTURE DEFINITION IS CORRECT...
   if not array_equal(tag_names(rmsf_in),['RMSF','PHI_RMSF']) then $
      message, 'RMSF_IN must be a structure with tags RMSF and PHI_RMSF'

   szRMSF = size(RMSF_IN.rmsf)
   n_RMSF = szRMSF[szRMSF[0]+2]

   ; MAKE SURE INPUT RMSF IS ONE DIMENSIONAL...
   if (szRMSF[0] ne 1) then $
      message, 'RMSF_IN.RMSF must be 1 dimensional.'

   ; MAKE SURE INPUT RMSF IS COMPLEX...
   rmsf_type = szRMSF[szRMSF[0]+1]
   if (rmsf_type ne 6) AND (rmsf_type ne 9) then $
      message, 'RMSF_IN.RMSF must be a complex array'

   ; MAKE SURE INPUT RMSF HAS TWICE THE NUMBER OF ELEMENTS AS FDF_DIRTY...
   if (n_RMSF ne 2*n_fdf) then $
      message, 'RMSF_IN.RMSF must have twice the number of elements as FDF_DIRTY.'

   ; MAKE SURE RMSF.PHI_RMSF HAS TWICE NUMBER OF CHANNELS AS FDF_DIRTY...
   if (N_elements(rmsf_in.phi_rmsf) ne 2*n_fdf) then $
      message, 'RMSF_IN.PHI_RMSF must have twice the number of elements as FDF_DIRTY.'

endif

; DEFAULT FOR MAXITER WILL BE 1000...
if (N_elements(maxiter) eq 0) then maxiter = 1000l

; DEFAULT FOR GAIN IS 0.1...
if (N_elements(gain) eq 0) then gain = 0.1

; GET THE SIZE OF THE FARADAY DEPTH AND LAMBDA-SQUARED ARRAYS...
nphi = N_elements(phi)
nlambdasq = N_elements(lambdasq)

; MAKE SURE PHI AND FDF ARRAYS ARE THE SAME SIZE...
if (n_fdf ne nphi) $
   then message, 'FDF_DIRTY and PHI vectors must have the same size.'

; CHECK FOR EXISTENCE OF WEIGHT KEYWORD...
nweight = N_elements(weight_in)
if (nweight eq 0) then begin
   ; IF NOT PASSED IN THEN JUST USE UNIFORM WEIGHTING...
   weight = dblarr(nlambdasq)+1.0
endif else begin
   weight = weight_in
   ; MAKE SURE LAMBDASQ AND WEIGHT ARRAYS ARE THE SAME SIZE...
   if (nweight ne nlambdasq) $
      then message, 'LAMBDASQ and WEIGHT vectors must have the same size.'
endelse

;!!!
; vvv this won't apply for non-uniform weighting
; GET THE THEORETICAL FWHM OF RMSF...
; THIS IS EQUATION (61) OF BdB05...
fwhm = 2d0*sqrt(3.0)/(max(lambdasq)-min(lambdasq))
print, fwhm

; BdB EQUATIONS (24) AND (38) GIVE THE INVERSE SUM OF THE WEIGHTS...
K = 1d0 / total(weight,/DOUBLE)

; GET THE MEAN OF THE LAMBDA-SQUARED DISTRIBUTION...
; THIS IS EQUATION (32) OF BdB05...
lambda0sq = K * total(weight * lambdasq,/DOUBLE)
;print, lambda0sq

icomp = dcomplex(0,1)

; CALCULATE THE RMSF IF IT HASN'T ALREADY BEEN PASSED IN...
if (n_RMSF_in gt 0) then begin
   RMSF = rmsf_in.rmsf
   phi2 = rmsf_in.phi_rmsf
endif else begin

   ; GEORGE HARDWIRES RMSF TO HAVE AN ODD NUMBER OF CHANNELS 2*NPHI+1...
   ; HE ALWAYS PUTS THE MAXIMUM OF THE RMSF AT PHI=0 IN CHANNEL
   ; NPHI (OUT OF 2*NPHI+1 CHANNELS)
   ; FORTRAN HAS INTEGER DIVISION, SO NUMPHI/2 = 500 FOR NUMPHI=1001

   ; THE RMSF WILL HAVE TWICE THE EXTENT OF OUR FARADAY DISPERSION FUNCTION...
   sample_phi = findgen(2*nphi)-nphi/2
   phi2 = interpol(phi,lindgen(nphi),sample_phi)

   ; THIS IS EQUATION (26) OF BdB05...
   RMSF = K * total(rebin(reform(weight,1,nlambdasq,/OVERWRITE),2*nphi,nlambdasq,/SAMPLE) * $
                    exp(-2d0 * icomp * phi2 # (lambdasq - lambda0sq)),2,/DOUBLE)

   if arg_present(RMSF_OUT) $
      then rmsf_out = {rmsf:rmsf,phi_rmsf:phi2}

endelse

; FIND THE POSITION OF THE PEAK OF THE RMSF AMPLITUDE...
useless = max(abs(RMSF), center_RMSF)


; IF WE WANT TO WATCH THE CLEAN PROCESS THEN DEFINE PIXMAP TEMPLATES...
if keyword_set(WATCH) then begin
   !p.charsize=0.8
   xsize = 1024
   ysize = 512
   ; ESTABLISH X WINDOW...
   window, /FREE, XS=xsize, YS=ysize
   xwin = !d.window
   ; ESTABLISH A PIXMAP WINDOW TO HOLD THE PLOTS THAT WON'T CHANGE...
   window, /FREE, XS=xsize, YS=ysize, /PIXMAP
   pixwin_static = !d.window
   ; SET UP 2-BY-2 PLOT WITH RMSF IN UPPER LEFT, FDF IN UPPER RIGHT,
   ; FDF RESIDUALS IN LOWER LEFT AND THE CLEAN MODEL COMPONENTS IN LOWER
   ; RIGHT...
   !p.multi=[0,2,2]
   plot, phi2, abs(RMSF), XSTYLE=3, YSTYLE=19, $
         XTIT=textoidl('Faraday Depth \phi (rad m^{-2})'), $
         YTIT='RMSF'
   !p.multi=[3,2,2]
   plot, phi, abs(fdf_dirty), XSTYLE=3, YSTYLE=1, $
         XTIT=textoidl('Faraday Depth \phi (rad m^{-2})'), $
         YTIT='Dirty Faraday Dispersion Function'
   oplot, phi, replicate(cutoff,nphi), LINES=1, COLOR=!red
   fdfpxy = {p:!p,x:!x,y:!y}
   yrange = !y.crange
   ; ESTABLISH PIXMAP WINDOW...
   window, /FREE, XS=xsize, YS=ysize, /PIXMAP
   pixwin = !d.window
endif

; SET REISDUALS TO INPUT DIRTY FARADAY DISPERSION FUNCTION...
fdf_residuals = fdf_dirty

; INITIALIZE STUFF...
iterations = 0
clean_components = dcomplexarr(nphi)
fdf_cleaned = dcomplexarr(nphi)

; KEEP GOING UNTIL THE PEAK OF THE AMPLITUDE OF THE RESIDUAL FDF GETS BELOW
; THE CUTOFF OR WE REACH THE MAXIMUM NUMBER OF ITERATIONS...
while ( (max(abs(fdf_residuals)) ge cutoff) AND (iterations lt maxiter) ) do begin

   ; WHERE IS THE PEAK OF THE AMPLITUDE OF THE RESIDUAL FDF...
   useless = max(abs(fdf_residuals), peak_chan)

   ; STORE THE VALUE OF THE FDF RESIDUALS AT THE CHANNEL WHERE THE ABSOLUTE
   ; AMPLITUDE OF THE RESIDUALS PEAKS...
   max_abs_fdf_resids = fdf_residuals[peak_chan]

   ; WHAT IS THE FARADAY DEPTH AT THIS CHANNEL...
   phi_peak = phi[peak_chan]

   ; MULTIPLY THE MAX ABS AMPLITUDE OF THE RESIDUALS BY THE GAIN AND THIS
   ; WILL BE OUR CLEAN COMPONENT...
   clean_comp = gain * max_abs_fdf_resids

   ; STORE THE CLEAN COMPONENT...
   clean_components[peak_chan] += clean_comp

   ; WHICH CHANNEL IS THIS CLEAN COMPONENT LOCATED AT IN THE RMSF ARRAY...
   useless = min(abs(phi2-phi_peak), peak_index_RMSF)

   ; SHIFT THE RMSF IN FARADAY DEPTH SO THAT ITS PEAK IS CENTERED ABOVE
   ; THIS CLEAN COMPONENT...
   shifted_RMSF = shift(RMSF, peak_index_RMSF-center_RMSF)

   index0 = center_RMSF - nphi/2
   index1 = center_RMSF + nphi/2 - (nphi mod 2 eq 0)

   ; NOW THAT WE'VE SHIFTED THE ENTIRE RMSF, TAKE ONLY THE CHANNELS
   ; THAT CORRESPOND TO OUR FARADAY DEPTH RANGE...
   subtract_RMSF = shifted_RMSF[index0:index1]

   ; SUBTRACT THE PRODUCT OF THE GAIN, THE SHIFTED RMSF, AND THE PEAK
   ; RESIDUAL FROM THE RESIDUALS...
   fdf_residuals -= clean_comp * subtract_RMSF

   ; DO WE WANT TO DEROTATE BACK TO LAMBDA=0...
   if keyword_set(DEROTATE) then $
      clean_comp = clean_comp * exp(-2d0 * icomp * phi_peak * lambda0sq)

   ; WHAT'S THE DIFFERENCE HERE...
   ; HEALD ET AL 2009 SECTION A.3 SAYS WE SET THE IMAGINARY PART TO ZERO...

   ; ADD TO THE CLEANED FDF THE CLEAN COMPONENT MULTIPLIED BY A GAUSSIAN
   ; WITH THE FWHM OF THE RMSF...
   fdf_cleaned += clean_comp * exp(-2.77258872224d0 * ((phi-phi_peak)/fwhm)^2)

   iterations +=1

   if not keyword_set(WATCH) then continue

   ; SHOW THE PROGRESS OF THE CLEANING...
   wait, 0.05
   wset, pixwin
   device, copy=[0,0,!d.x_vsize,!d.y_vsize,0,0,pixwin_static]
   !p.multi=[2,2,2]
   plot, phi, abs(fdf_residuals), XSTYLE=3, YSTYLE=1, $
         XTIT=textoidl('Faraday Depth \phi (rad m^{-2})'), $
         YTIT='Residuals'
   oplot, phi, replicate(cutoff,nphi), LINES=1, COLOR=!red
   !p.multi=[1,2,2]
   plot, phi, abs(fdf_cleaned), XSTYLE=3, YSTYLE=1, YRANGE=yrange, $
         XTIT=textoidl('Faraday Depth \phi (rad m^{-2})'), $
         YTIT='Clean Components'
   oplot, phi, clean_components, COLOR=!magenta
   wset, xwin
   device, copy=[0,0,!d.x_vsize,!d.y_vsize,0,0,pixwin]

endwhile

; ADD THE RESIDUALS TO THE CLEANED FARADAY DISPERSION FUNCTION...
fdf_cleaned += not keyword_set(DEROTATE) $
               ? fdf_residuals $
               ; DO WE WANT TO DEROTATE RESIDUALS BACK TO LAMBDA=0...
               ;??? THIS IS WHAT GEORGE DOES...
               : fdf_residuals * exp(-2d0 * icomp * phi * lambda0sq)

; OVERPLOT THE FINAL CLEANED SPECTRUM...
if keyword_set(WATCH) AND (iterations gt 0) then begin
   wdelete, pixwin_static, pixwin
   wset, xwin
   !p=fdfpxy.p & !x=fdfpxy.x & !y=fdfpxy.y
   oplot, phi, abs(fdf_cleaned), COLOR=!green
   !p.multi=0
   !p.charsize=0
endif

; IF THE USER PASSED IN THE DIRTY FDF AS SINGLE-PRECISION COMPLEX, THEN
; RETURN EVERYTHING AS SINGLE-PRECISION...
if (fdf_type eq 9) then return

fdf_cleaned = complex(fdf_cleaned)
fdf_residuals = complex(fdf_residuals)
clean_components = complex(clean_components)

end ; rmclean
