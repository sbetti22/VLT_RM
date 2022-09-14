#run from 


casa
import sys

for my_calibrator in phasecal_target_list:
    #fill in your directory where the calibrators are
    fil = datadir + split_ms_name + '_' + my_calibrator+'_'+band+'.ms'
    #fil = '/students/sbetti/12B-336/12B-336.sb13594278.autoflag.30s.8mhz_calibrators_L.ms'
    for mytarget in phasecal_target_list[my_calibrator]:
 	if os.path.exists(datadir + 'IDL_Scripts/Fits_Images/' + mytarget + '.fits'):
        	print "NOTE: " + mytarget + '.fits' + ' exists!'
	else:
        	clean(vis=fil, mode='mfs', niter=3000, gain=0.1, threshold='0.25mJy', psfmode='clark', interactive=True, imsize=[128], cell=['0.4arcsec', '0.4arcsec'], stokes='IQU', weighting='natural', imagename=datadir+'CASA_Images/MFS/'+mytarget+'mfs', field=mytarget)
        	
            clean(vis=fil, mode='frequency', niter=300000, gain=0.1, threshold='0.25mJy', psfmode='clark', interactive=True, imsize=[128], cell=['0.4arcsec', '0.4arcsec'], stokes='IQU', weighting='natural', imagename=datadir+'CASA_Images/FREQ/'+mytarget+'freq', field=mytarget)
        	exportfits(imagename=datadir+'CASA_Images/FREQ/'+mytarget+'freq.image', fitsimage=datadir+ 'IDL_Scripts/Fits_Images/' + mytarget+'.fits', overwrite=True)
		contine = raw_input("enter 'y' to continue: ")
		if contine != 'y':
			sys.exit()
		else:
			continue  

tclean(vis='12B-336.sb.autoflag.30s.8mhz.ms',
imagename='SMITH_02316',
field='SMITH_02316',
spw='',
specmode='mfs',
deconvolver='hogbom',
nterms=1,
gridder='standard',
imsize=[128,128],
cell=['0.4arcsec'],
weighting='natural',
threshold='0mJy',
niter=5000,
interactive=True)

name = 'SMITH_02728'
clean(vis='12B-336.sb.autoflag.30s.8mhz.ms', 
      mode='frequency',
      niter=300000, 
      gain=0.1, 
      threshold='0.1mJy', 
      psfmode='clark', 
      interactive=True, 
      imsize=[128], 
      cell=['0.4arcsec', '0.4arcsec'], 
      stokes='IQU', 
      weighting='natural',
      imagename='CASA_images/' + name, 
      field=name)

exportfits(imagename='CASA_images/' + name + '.image', fitsimage='fits_images/' + name+'.fits', overwrite=True)



