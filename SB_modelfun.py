# Written by S.A.M on 1.11.13, same as setqu_of_model.py but it is now a function 
# setqu takes the central pixel values x0, y0 and the name of the model and inject 
# the model IQU flux to the central pixel in the image cube
#Effectively making a point source model for the polang calibrator
# USES Perley & Butler 2010 denisty scale

import numpy as np

def extrap(x, xp, yp):
    """np.interp function with linear extrapolation"""
    y = np.interp(x, xp, yp)
    y = np.where(x < xp[0], yp[0]+(x-xp[0])*(yp[0]-yp[1])/(xp[0]-xp[1]), y)
    y = np.where(x > xp[-1], yp[-1]+(x-xp[-1])*(yp[-1]-yp[-2])/(xp[-1]-xp[-2]), y)
    return y

# Setqu makes a model image cube (modelname) with zero values everywhere and
# the correct IQU flux at x0,y0
#This is effectively a point source model
# for polang_calibrator_source takes in the central pixel location 

def setqu(x0, y0,modelname,polang_calibrator_source):

    ia.open(modelname)
    
    f = open("polcal_table.txt", 'r')
    freq_mhz=[];pp_3C48=[];pa_3C48=[];pp_3C138=[];pa_3C138=[];pp_3C147=[];pa_3C147=[];pp_3C286=[]
    for line in f:
        col1,col2,col3,col4,col5,col6,col7,col8= line.split()
        freq_mhz.append(float(col1)) 
        pp_3C48.append(float(col2))
        pa_3C48.append(float(col3))
        pp_3C138.append(float(col4))
        pa_3C138.append(float(col5))
        pp_3C147.append(float(col6))
        pa_3C147.append(float(col7))
        pp_3C286.append(float(col8))


        freq_ghz = [x/1e3 for x in freq_mhz]


        spectraldict = {'3C48':[1.3197,-0.7253,-0.2023,0.0540],
	        '3C138':[1.0053,-0.4384,-0.1855,0.0511], 
		'3C147':[1.4428,-0.6300,-0.3142,0.1032], 
		'3C286':[1.2361,-0.4127,-0.1864,0.0294]} 

        dimension=imhead(imagename=modelname,mode="get",hdkey="shape")
        tot_chan = (dimension['value'][3])

# freq in GHz

        freq_info=imhead(imagename=modelname,mode="get",hdkey="crval4")
        freq_lower = (freq_info['value'])/1e9 


        delta_info=imhead(imagename=modelname,mode="get",hdkey="cdelt4")
        delta_freq = (delta_info['value'])/1e9

        chan_num = range(tot_chan)
        
        for i in chan_num:


            f_ghz = freq_lower+i*delta_freq
        
            si = pow(10,(spectraldict[polang_calibrator_source][0]+spectraldict[polang_calibrator_source][1]*np.log10(f_ghz)+spectraldict[polang_calibrator_source][2]*pow(np.log10(f_ghz),2)+spectraldict[polang_calibrator_source][3]*pow(np.log10(f_ghz),3)))

            if polang_calibrator_source == '3C286':
                my_pa=33.0 
                my_pp = extrap(f_ghz,freq_ghz,pp_3C286)
      
            if polang_calibrator_source == '3C138':
                my_pa = extrap(f_ghz,freq_ghz,pa_3C138)
                my_pp = extrap(f_ghz,freq_ghz,pp_3C138)

            if polang_calibrator_source == '3C48':
                my_pa = extrap(f_ghz,freq_ghz,pa_3C48)
                my_pp = extrap(f_ghz,freq_ghz,pp_3C48)

            if polang_calibrator_source == '3C147':
                my_pa = extrap(f_ghz,freq_ghz,pa_3C147)
                my_pp = extrap(f_ghz,freq_ghz,pp_3C147)


            sq = si*my_pp*0.01*np.cos(my_pa*2.0*np.pi/180.0)
            su = si*my_pp*0.01*np.sin(my_pa*2.0*np.pi/180.0)

        #now do polcal options=x
            ia.open(modelname)
            array_i=[si]
            array_q=[sq]
            array_u=[su]

            ia.putchunk(array_i,blc=[x0,y0,0,i])
            ia.putchunk(array_q,blc=[x0,y0,1,i])
            ia.putchunk(array_u,blc=[x0,y0,2,i])
        

    ia.close()

