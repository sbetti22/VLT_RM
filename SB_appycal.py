######################################  Apply Calbiration ######################

for mycalibrator in phasecal_target_list:
    applycal(vis=split_ms_file, field=mycalibrator, spw=bandspw, gaintable=[outputdir+split_ms_name+'_'+band+'_delay.cal', outputdir+split_ms_name+'_'+band+'_bandpass.cal', outputdir+split_ms_name+'_'+band+'_fluxscaled.cal', outputdir+split_ms_name+'_'+band+'_kcross.cal', outputdir+split_ms_name+'_'+band+'_leakage.cal', outputdir+split_ms_name+'_'+band+'_polx.cal'], gainfield=['','',mycalibrator,'','',''], interp=['','','nearest','','',''], parang=T, calwt=F)
    for mytarget in phasecal_target_list[mycalibrator]:
        applycal(vis=split_ms_file, field=mytarget, spw=bandspw, gaintable=[outputdir+split_ms_name+'_'+band+'_delay.cal', outputdir+split_ms_name+'_'+band+'_bandpass.cal', outputdir+split_ms_name+'_'+band+'_fluxscaled.cal', outputdir+split_ms_name+'_'+band+'_kcross.cal', outputdir+split_ms_name+'_'+band+'_leakage.cal', outputdir+split_ms_name+'_'+band+'_polx.cal'], gainfield=['','',mycalibrator,'','',''], interp=['linear','','','','',''], parang=T, calwt=F)
    split(vis=split_ms_file, spw=bandspw, outputvis=datadir + split_ms_name + '_' + mycalibrator+'_'+band+'.ms', datacolumn='corrected', field=phasecal_target_list[mycalibrator])

split(vis=split_ms_file, spw=bandspw, outputvis=datadir + split_ms_name + '_calibrators_' + band + '.ms', datacolumn='corrected', intent='CALIB*')
listobs(vis=datadir + split_ms_name + '_calibrators_' + band + '.ms', listfile=outputdir+split_ms_name + '_calibrators_' + band +'.ms.listobs')     
