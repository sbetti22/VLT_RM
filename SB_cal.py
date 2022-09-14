#######################################CALIBRATION######################################
#set flux density of the absolute flux calibrator
#### ASH 2013-2-7: is this causing trouble because the model image varies across our band?
setjy(vis=split_ms_file, field=flux_calibrator_field, spw=bandspw, modimage=cal_model, standard=myflux_standard, usescratch=F, scalebychan=True)

if antpos:
    gencal(vis=split_ms_file, caltable=outputdir+split_ms_name+'.antpos', caltype='antpos')

#obtain parallel-hand delay solution, which is SPW based
without_edgechan='1~6'
gaincal(vis=split_ms_file, caltable=outputdir+split_ms_name+'_'+band+'_delay.cal', spw=bandspw+':'+without_edgechan, field=flux_calibrator_field, refant=refant, minsnr=5, selectdata=T, gaintype='K')

#Inspect delay solution
if doplots:
    plotcal(caltable=outputdir+split_ms_name+'_'+band+'_delay.cal', xaxis='antenna', yaxis='delay', spw=bandspw, iteration='spw', figfile=outputdir + 'delay_'+band+'.png')

#obtain per channel bandpass solution
#bandpass(vis=split_ms_file, caltable=outputdir+split_ms_name+'_'+band+'_bandpass.cal', field=flux_calibrator_field, spw=bandspw, solint='inf', refant=refant, solnorm=T, bandtype='B', gaintable=[outputdir+split_ms_name+'_'+band+'_delay.cal'], interp=['nearest'])

#bandpass calibration, note that here we introduce binning to obtain better signal to noise
if antpos:
    bandpass(vis=split_ms_file, caltable=outputdir+split_ms_name+'_'+band+'_bandpass.cal', field=flux_calibrator_field, spw=bandspw, solint='inf', refant=refant, solnorm=T, bandtype='B', gaintable=[outputdir+split_ms_name+'.antpos', outputdir+split_ms_name+'_'+band+'_delay.cal'], interp=['', 'nearest'])
else:
    bandpass(vis=split_ms_file, caltable=outputdir+split_ms_name+'_'+band+'_bandpass.cal', field=flux_calibrator_field, spw=bandspw, solint='inf', refant=refant, solnorm=T, bandtype='B', gaintable=[outputdir+split_ms_name+'_'+band+'_delay.cal'], interp=['nearest'])

#bandpass(vis=split_ms_file, caltable='temp.cal', field=flux_calibrator_field, spw='2', solint='inf, 1ch', refant=refant, solnorm=T, bandtype='B', gaintable=[outputdir+split_ms_name+'_'+band+'_delay.cal'], interp=['nearest'])

#inspect bandpass solution
if doplots:
    plotcal(caltable=outputdir+split_ms_name+'_'+band+'_bandpass.cal', xaxis='freq', yaxis='amp', spw=bandspw, iteration='antenna', subplot=331, figfile=outputdir + 'bandpass_'+band+'.png')

#obtain gain calibration
gaincal(vis=split_ms_file, caltable=outputdir+split_ms_name+'_'+band+'_gain.cal', field=gaincal_fields, spw=bandspw, gaintype='G', refant=refant, gaintable=[outputdir+split_ms_name+'_'+band+'_delay.cal', outputdir+split_ms_name+'_'+band+'_bandpass.cal'], gainfield=[flux_calibrator_field, flux_calibrator_field], interp=['nearest', 'nearest'], calmode='ap')

# Inspect gain solutions
if doplots:
    #Look at R and L separately
    plotcal(caltable=outputdir+split_ms_name+'_'+band+'_gain.cal', iteration='antenna', spw=bandspw, subplot=321, xaxis='time', yaxis='amp', poln='R', figfile=outputdir + 'gain_R_'+band+'.png')
    plotcal(caltable=outputdir+split_ms_name+'_'+band+'_gain.cal', iteration='antenna', spw=bandspw, subplot=321, xaxis='time', yaxis='amp', poln='L', figfile=outputdir + 'gain_L_'+band+'.png')
    #Look at amplitude and phase
    plotcal(caltable=outputdir+split_ms_name+'_'+band+'_gain.cal', iteration='antenna', spw=bandspw, subplot=321, xaxis='time', yaxis='phase', figfile=outputdir + 'gain_phase_'+band+'.png')
    plotcal(caltable=outputdir+split_ms_name+'_'+band+'_gain.cal', iteration='antenna', spw=bandspw, subplot=321, xaxis='time', yaxis='amp', figfile=outputdir + 'gain_'+band+'.png')

########################################POLARIZATION CALIBRATION choose either (i) or (ii) #################################################
## First step: To be safe, we delete the flux information in the header of the polarization calibrator to avoid confusion, since we are defining it in the model column in the next section
#delmod(vis=split_ms_file, otf=T, field=polang_calibrator_field)

if polcal_using_model == True:
######################################(i) Using a Polarization Model Image ################################################################################
# To construct a polarizaiton model can be difficult = > Making a large image cube is slow!
# to speed things up, we can bin to coarser channels and using smaller image size

# Resolution of data at different bands:
#at 1.5 GHz, synthesized beam size is 1.3 arcsec for A Config  (ranging from 1.69998626374", 1.13332417582", 0.849993131868") cellsize ~ 0.3 arcsec imsize=160  spw 0~15
#at 3.0 GHz, synthesized beam size is  0.65 arsrec for A Config (ranging from  0.849993131868", 0.566662087912", 0.424996565934") cellsize ~ 0.15 arcsec  imsize=320spw 16~31
#at 6.0 GHz, synthesized beam size is  0.33 arcsec for A Config: 
#C1 4-6GHz (ranging from  0.424996565934",0.339997252747",0.283331043956") cellsize ~ 0.1 arcsec  imsize=480  spw 32-47
#C2 6-8GHz (ranging from  0.283331043956, 0.242855180534",0.212498282967") cellsize ~ 0.07 arcsec imsize=800 spw 48~63
#also note that when we make a image cube the resolution is set to be the worst one within the frequency cube

#Beam has to be well sampled, or else fourier transform of the model image (what is entered into the MODEL_DATA column of the measurement set) will not be correct
#Also, need to specify interpolation=nearest, otherwise the frequency output may be strange: some channels may be skipped, resulting in problems when filling in values back into the ms
#Meanwhile, the image size has to be large enough too (cannot use a very small image size, or else get crazy numbers)
#    Bmax=36.4e3 #maximum baseline length for A config
#    nu = [1.0e9,1.5e9,2.0e9,3.0e9,4.0e9,5.0e9,6.0e9,7.0e9,8.0e9]
#    for i in nu:
#         wavelength=3e8/i
#         res=wavelength/Bmax*206265
#         print res


#Define properties of the modeled image depending on the band that we are working on:
     if band == 'L':
          mycellsize='0.3arcsec'
          nx=160 
          ny=160

     if band == 'S':
          mycellsize='0.15arcsec'
          nx=320
          ny=320

     if band =='C1':
          mycellsize='0.1arcsec'
          nx=480
          ny=480

     if band =='C2':
          mycellsize='0.07arcsec'
          nx=700
          ny=700

     imsize=[nx,ny]
     x0=nx/2
     y0=ny/2

    
# Note that this needs to be binned so that we have a reasonable number of images!!! Otherwise the cube will be LARGE!
     for myspw in bandspw_list:
          os.system('rm -rf zero_'+myspw+'*') #remove zero_spw.*
          clean(vis=split_ms_file, imagename='zero_'+myspw, field=polang_calibrator_field, spw= myspw, width=binned_channel_width, start=0, psfmode='hogbom', mode='channel', interpolation='nearest', niter=1, gain=0.1, threshold='10.0mJy', imsize=imsize, cell=mycellsize, stokes='IQU', weighting='natural')
          os.system('rm -rf zero.modim')
          immath(imagename='zero_'+myspw+'.model', mode='evalexpr', outfile='zero.modim', expr='IM0*0', varnames='', sigma='0.0mJy/beam')
          setqu(x0, y0, 'zero.modim', polang_calibrator_source) 
          ft(vis=split_ms_file, field=polang_calibrator_field, model='zero.modim', spw=myspw)
          os.system('rm -rf zero_'+myspw+'*') #remove zero_spw.* in the directory

#Examine how RL and LR phase and amp changes over the frequency bands in plotms


#####################################################################################################################################################


if polcal_using_model == False:
     ######################################(ii)  Use SETJY to set polang, only for 3C286 ##########################################################################
    #Note that this only works for 3C286, which has a constant PA pretty much across 1-50GHz (does change a little bit, but for our purpose it is 33 deg)

    #Fitted spectral index, curvature and higher order coefficients for 3C286 (Table 10 in Perley & Butler 2013),
    #Hold off using this for now, since Perley & Butler 2013 is not available for 3C138..and I don't know quite how to scale that yet 
    spectraldict ={'3C286':[1.2515,-0.4605,-0.1715,0.0336]} 

    # Perley & Butler 2010 density scale (AIPS MANUAL http://www.aips.nrao.edu/cgi-bin/ZXHLP2.PL?SETJY,verified by checking its flux at
     # 4.536 GHz is 7.81694 Jy, same as what tutorial on 3C391 online gives )
    #spectraldict={'3C286':[1.2361,-0.4127,-0.1864,0.0294]} #Perley & Butler 2010 denisty scale


    # Setjy in every single SPW, lookup channel zero value from listobs and output into freq.txt
    # Format:
    # 0          64 TOPO  994         1000          64000       RR  RL  LR  LLW output

     os.system('grep TOPO ' + outputdir+split_ms_name+'.ms.listobs' + ' > freq.txt')

     freq_file= open("freq.txt", 'r')

     spw_name=[];CH0freqMHz=[];totwidthkHz=[]
     for line in freq_file:
#          col1, col2, col3, col4, col5, col6, col7, col8, col9, col10= line.split()
          collist = line.split()
          spw_name.append((collist[1])) 
          CH0freqMHz.append(float(collist[4]))
          totwidthkHz.append(float(collist[6]))

     CH0freqGHz = [x/1.0e3 for x in CH0freqMHz]
     totwidthGHz = [x/1.0e6 for x in totwidthkHz]
     
     #lookup and interpoate the fractional polarization of 3C286
     pol_file = open("polcal_table.txt", 'r')
     freq_mhz_pol_file=[];pp_3C286=[]
     for line in pol_file:
          col1, col2, col3, col4, col5, col6, col7, col8= line.split()
          freq_mhz_pol_file.append(float(col1)) 
          pp_3C286.append(float(col8))

     freq_ghz_pol_file=[x/1.0e3 for x in freq_mhz_pol_file]

     
     bandspw_list = spw_name
     for bnum in bandspw_list:
          ch0=CH0freqGHz[spw_name.index(bnum)]
          ch1=ch0+totwidthGHz[spw_name.index(bnum)]
          #print ch0, ch1
        
          stokes_i = pow(10, (spectraldict['3C286'][0]+spectraldict['3C286'][1]*np.log10(ch0)+spectraldict['3C286'][2]*pow(np.log10(ch0), 2)+spectraldict['3C286'][3]*pow(np.log10(ch0), 3)))
            
          stokes_i_lastchan =  pow(10, (spectraldict['3C286'][0]+spectraldict['3C286'][1]*np.log10(ch1)+spectraldict['3C286'][2]*pow(np.log10(ch1), 2)+spectraldict['3C286'][3]*pow(np.log10(ch1), 3)))
        
          #print stokes_i, stokes_i_lastchan
          alpha=log(stokes_i_lastchan/stokes_i)/log(ch1/ch0) #spectral index
          pp = extrap(ch0, freq_ghz_pol_file, pp_3C286) #extrapolated percent polarization at ch0
          pol_int = 0.01*pp*stokes_i #polarized intensity
          stokes_q = pol_int*cos(66.0*pi/180.0) #assume that the RL phase difference is 66 degrees
          stokes_u = pol_int*sin(66.0*pi/180.0) 
          setjy(vis=split_ms_file, field=polang_calibrator_field, spw=str(bnum), fluxdensity=[stokes_i, stokes_q, stokes_u, 0], spix=alpha, reffreq=str(ch0)+'GHz', scalebychan=True, usescratch=False, standard='manual')
          #print ch0, stokes_i, stokes_q, stokes_u, alpha
          
#CHECKED RL phas look okay and continous
#RL amplitude looks discontinuous across SPW, this is due to the fact that we have applied a percent pol that changes across the band
#in particular the percent pol increases as frequency increases, causing what should be a shallower polarization spectral index
#however, the task setjy does not know and apply the spectral index alpha obtained from stokesi to do the scalebychan
#this results in steeper than expected spectral index that does not connect across spw can eliminate this using a constant percent pol across all SPW
#good that the only information needed is the angle... 

###################################### Leakage Calibration #########################################

# Solving for cross-hand delays, Kcross does not know about the solint channel component
gaincal(vis=split_ms_file, caltable=outputdir+split_ms_name+'_'+band+'_kcross.cal', field=polang_calibrator_field, spw=bandspw, refant=refant, gaintype='KCROSS', calmode='ap', append=False, gaintable=[outputdir+split_ms_name+'_'+band+'_delay.cal', outputdir+split_ms_name+'_'+band+'_bandpass.cal', outputdir+split_ms_name+'_'+band+'_gain.cal'], gainfield=[flux_calibrator_field, flux_calibrator_field, polang_calibrator_field], interp=['nearest','nearest','nearest'], solint='inf', parang=True)

#inspect cross-hand delay table
if doplots:
    plotcal(caltable=outputdir+split_ms_name+'_'+band+'_kcross.cal', xaxis='antenna', yaxis='delay', figfile=outputdir + 'kcross_'+band+'.png')

#solve for the Leakage terms
polcal(vis=split_ms_file, caltable=outputdir+split_ms_name+'_'+band+'_leakage.cal', field=leakage_calibrator_field, spw=bandspw, solint='inf', refant=refant, poltype='Df', gaintable=[outputdir+split_ms_name+'_'+band+'_delay.cal', outputdir+split_ms_name+'_'+band+'_bandpass.cal', outputdir+split_ms_name+'_'+band+'_gain.cal', outputdir+split_ms_name+'_'+band+'_kcross.cal'], gainfield=[flux_calibrator_field,flux_calibrator_field,leakage_calibrator_field, polang_calibrator_field], interp=['nearest','nearest', 'nearest','nearest'])

#inspect the leakage terms
if doplots:
    plotcal(caltable=outputdir+split_ms_name+'_'+band+'_leakage.cal', xaxis='freq', yaxis='amp', iteration='antenna', subplot=331, figfile=outputdir + 'leakage_'+band+'.png') 


###################################### Polarization Angle Calibration ##############################
#Calibrate PA, note that we **CANNOT add refant in this step
polcal(vis=split_ms_file, caltable=outputdir+split_ms_name+'_'+band+'_polx.cal', field=polang_calibrator_field, solint='inf', spw=bandspw, poltype='Xf', gaintable=[outputdir+split_ms_name+'_'+band+'_delay.cal', outputdir+split_ms_name+'_'+band+'_bandpass.cal', outputdir+split_ms_name+'_'+band+'_gain.cal', outputdir+split_ms_name+'_'+band+'_kcross.cal', outputdir+split_ms_name+'_'+band+'_leakage.cal'], gainfield=[flux_calibrator_field,flux_calibrator_field,polang_calibrator_field,polang_calibrator_field,leakage_calibrator_field], interp=['nearest','nearest','nearest','nearest','nearest'])

#inspect the absolute angle calibration solution
if doplots:
    plotcal(caltable=outputdir+split_ms_name+'_'+band+'_polx.cal', spw=bandspw, xaxis='freq', yaxis='phase', figfile=outputdir+'polx_'+band+'.png')


###################################### Absolute Flux Density Scale #################################
# boostrap flux density scale from flux_calibrator_field
fluxscale(vis=split_ms_file, fluxtable=outputdir+split_ms_name+'_'+band+'_fluxscaled.cal', caltable=outputdir+split_ms_name+'_'+band+'_gain.cal', reference=flux_calibrator_field)
if doplots:
    plotcal(caltable=outputdir+split_ms_name+'_'+band+'_gain.cal', xaxis='time', yaxis='amp', poln='R', figfile=outputdir+'flux_R'+band+'.png')
    plotcal(caltable=outputdir+split_ms_name+'_'+band+'_gain.cal', xaxis='time', yaxis='amp', poln='L', figfile=outputdir+'flux_L'+band+'.png')
