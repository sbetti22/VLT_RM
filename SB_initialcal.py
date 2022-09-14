#initial calibration for 12B-336

##############################Start calibration/pre-flagging step###############################
# Only run if it is the first time:
if firstrun == True:  
    #produce listobs in the outputdir directory
    DONElistobs(vis=ms_file,listfile=outputdir+raw_ms_name+'.ms.listobs')     
    #produce antennae plot in the outputdir directory
    DONEplotants(vis=ms_file,figfile=outputdir+raw_ms_name+'.plotants.png')
    clearstat()

#make antenna position correction if needed (may give error if no correction is needed)
    DONEgencal(vis=ms_file, caltable=outputdir+raw_ms_name+'.antpos', caltype='antpos')
#
################################# FLAG DATA #############################################
##Only run if it is the first time flag all data
if firstrun == True: 
    ## skip stuff that was done in importevla anyway
    flagdata(vis=ms_file,mode='unflag',flagbackup=False)
    #flag shadowed data
    flagdata(vis=ms_file,mode='shadow',flagbackup=False)
    #clip data with 0 vis
    flagdata(vis=ms_file,mode='clip',clipzeros=True,flagbackup=False)
    #remove the setup dummy scans
    flagdata(vis=ms_file,mode='manual', scan=setup_scan,flagbackup=T)
    #flag initial samples from start of each scan
    flagdata(vis=ms_file, mode='quack', quackinterval=avg_interval, quackmode='beg')
    # NOTE THAT HANNING SMOOTHING: if no outputvis is specified then the input vis will be overwritten
    hanningsmooth(vis=ms_file,datacolumn='data')

#from observing log, flag bad antennas
    flagdata(vis=ms_file,mode='manual', antenna=bad_antennae, spw=bandspw,flagbackup=True)

#Look at data to see obvious problems
     plotms(vis=ms_file, selectdata=True, field='0')
     plotms(vis=fil, field='', correlation='RR,LL', timerange='', antenna='ea01', spw='0:31', xaxis='time', yaxis='antenna2', plotrange=[-1,-1,0,26], coloraxis='field')

# rfi_free_channels='0:10~30,1:20~28,2:30~45,3:10~15,4:20~30,5~7:35~50,8:10~20,9:50~55,10~13:35~50,14:30~45,15:20~35'
#solve for initial phase using field only if we use the phase calibrator to get the bandpass... skip to the next step
# gaincal(vis=ms_file, field=gaincal_fields, caltable='initial_phase_'+band+'.cal', solint='60s', spw=rfi_free_channels, refant=refant, calmode='p')


#directly use flux_calibrator scan to get an initial bandpass, do not bin this one since we want to get rid of RFIs
    bandpass(vis=ms_file, caltable=outputdir + 'initial_bpass_'+band + '_spw' + bandspw +'.cal', field=flux_calibrator_field, spw=bandspw, solint='inf', refant=refant, solnorm=F)

# Examine initial bandpass
if doplots:
    plotcal(caltable=outputdir+'initial_bpass_'+band+'_spw'+bandspw+'.cal', xaxis='freq', yaxis='amp', spw=bandspw, iteration='antenna', antenna='0~8', subplot=331, figfile=outputdir + 'initial_bpass_'+band+'_p1.png', showgui=F)
    plotcal(caltable=outputdir+'initial_bpass_'+band+'_spw'+bandspw+'.cal', xaxis='freq', yaxis='amp', spw=bandspw, iteration='antenna', antenna='9~17', subplot=331, figfile=outputdir + 'initial_bpass_'+band+'_p2.png', showgui=F)
    plotcal(caltable=outputdir+'initial_bpass_'+band+'_spw'+bandspw+'.cal', xaxis='freq', yaxis='amp', spw=bandspw, iteration='antenna', antenna='18~26', subplot=331, figfile=outputdir + 'initial_bpass_'+band+'_p3.png', showgui=F)
    plotcal(caltable=outputdir+'initial_bpass_'+band+'_spw'+bandspw+'.cal', xaxis='freq', yaxis='amp', spw=bandspw, iteration='antenna', antenna='27~35', subplot=331, figfile=outputdir + 'initial_bpass_'+band+'_p4.png', showgui=F)

# now apply the antpos (if exists) and the initial bandpass calibration to the measurement set
# If we want to look at the data now, would need to look at the 'corrected data' column

if antpos:
    applycal(vis=ms_file, gaintable=[outputdir+raw_ms_name+'.antpos',
    outputdir+'initial_bpass_'+band+'_spw'+bandspw+'.cal'], spw=bandspw, calwt=False)
else:
    applycal(vis=ms_file, gaintable=[outputdir+'initial_bpass_'+band+'_spw'+ bandspw+'.cal'], spw=bandspw, calwt=False)

#flag bad antennae
flagdata(vis=ms_file, mode='manual', spw=flagspw)

#AUTOFLAG first then extend to all polarizations
#Just to examine first with action='calculate'

#flagdata(vis=ms_file,mode='rflag',correlation='rr,ll',datacolumn='corrected',spw=bandspw,timedevscale=2,freqdevscale=2,action='calculate',display='report')

flagdata(vis=ms_file, mode='rflag', correlation='rr,ll', datacolumn='corrected', spw=bandspw, timedevscale=2, freqdevscale=2, action='apply', display='report', flagbackup=True)

#First run of rflag can leave some leftover residuals due to individual large spikes corrupting mean and rms, so run again
bandpass(vis=ms_file, caltable=outputdir + 'second_bpass_'+band + '_spw' + bandspw +'.cal', field=flux_calibrator_field, spw=bandspw, solint='inf', refant=refant, solnorm=F)

applycal(vis=ms_file, gaintable=[outputdir+raw_ms_name+'.antpos', outputdir+'second_bpass_'+band+'_spw'+bandspw+'.cal'], spw=bandspw, calwt=False)

flagdata(vis=ms_file, mode='rflag', correlation='rr,ll', datacolumn='corrected', spw=bandspw, timedevscale=2, freqdevscale=2, action='apply', display='report', flagbackup=True)

flagdata(vis=ms_file, mode='extend', extendpols=True, action='apply', spw=bandspw, display='report', growtime=100, growfreq=100, flagbackup=True) 
# Inspect data now to see if additional flagging is needed

#summarize  how much data has been flagged
#tflagdata(vis=ms_file,mode='rflag',correlation='rr,ll',spw=bandspw,datacolumn='corrected',timedevscale=2,freqdevscale=2,action='calculate',display='both')

################################SPLIT##########################################
#split out data to make file smaller. Splitting must be done in the current directory because casa creates a temp file

split_ms_name = projnum + '.sb' + sbnum + '.autoflag.' + gain_solution_interval + '.' + str(binned_chan_width_MHz) + 'mhz'

cwd = os.getcwd()
os.chdir(datadir)
split(vis=raw_ms_name + '.ms', outputvis=split_ms_name + '.ms', datacolumn='data', width=str(binned_chan_width_num), timebin=gain_solution_interval)
os.chdir(cwd)
split_ms_file = datadir + split_ms_name + '.ms'

if firstrun == True:
    listobs(vis=split_ms_file,listfile=outputdir+split_ms_name+'.ms.listobs')

#reinitialize the corrected data column (to the data column) and the model_data column (total intensity=1,pol intensity=0)  now that the flagging procedure is done

cwd = os.getcwd()
os.chdir(datadir)
clearcal(vis=split_ms_name + '.ms')
os.chdir(cwd)
