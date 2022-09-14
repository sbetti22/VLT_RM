#Reduction script written on 1.7.13 by S.A.M for CASA 4.0
#based on calibration.py, which is specific to 12OCT07
#This script reduces data for all Stokes parameters (IQU) including as much automatically flagging .
#Require the following files to be in the same working directory: setqu_of_model_func.py, pol, polcal_table.txt
# Script is to be ran within the subdirectories for the different bands
# main MS file is in the datadir directory
####
#Modified by Alex Hill for Smith Cloud point source RMs 2013 Jan-Feb
#Modified by Sarah Betti for Smith Cloud point source RMs 2015 September

#import functions from setqu_of_model_func.py
cwd = os.getcwd()
os.chdir('scripts')
from setqu_of_model_func import *
os.chdir(cwd)

##################PICK the band to reduce###########################
#Set parameters specific to band
####################################################################
band = 'L' 
delnu = 1 # L band individual channel width is 1MHz (from listsobs output) 
bandspw='0~63'
bandspw_list=['0','1','2','3','4','5','6','7','8','9','10','11','12','13','14','15']
flagspw = '8,9'     # spws from which nothing can be recovered due to RFI
flagspw = ' '
# antenna ea18 causes rippling in the baseline in the RR correlation, so flag
bad_antennae='ea18' #bad antennae to flag from observing log
avg_interval = 1.0

antpos = False  # are there any antpos corrections to load?

phasecal_target_list={ #a dictionary matching calibrator and corresponding targets, refer to observing log for this, interpolation is nearest
     # calculated in IDL with create_phasecal_target_list.idl
     'J1939-1002': ('SMITH_00734', 'SMITH_00831', 'SMITH_00843', 'SMITH_00845', 'SMITH_00853', 'SMITH_00893', 'SMITH_00900', 'SMITH_00902', 'SMITH_00906', 'SMITH_00909', 'SMITH_00910', 'SMITH_00926', 'SMITH_00928', 'SMITH_00935', 'SMITH_00939', 'SMITH_00941', 'SMITH_00958', 'SMITH_00961', 'SMITH_00971', 'SMITH_00978', 'SMITH_00990', 'SMITH_00995', 'SMITH_01002', 'SMITH_01010', 'SMITH_01011', 'SMITH_01021', 'SMITH_01086', 'SMITH_01093', 'SMITH_01094', 'SMITH_01111', 'SMITH_01151', 'SMITH_01163', 'SMITH_01167', 'SMITH_01173', 'SMITH_01177', 'SMITH_01188', 'SMITH_01191', 'SMITH_01203', 'SMITH_01211', 'SMITH_01218', 'SMITH_01225', 'SMITH_01229', 'SMITH_01230', 'SMITH_01233', 'SMITH_01251', 'SMITH_01255', 'SMITH_01259', 'SMITH_01260', 'SMITH_01262', 'SMITH_01265', 'SMITH_01280', 'SMITH_01281', 'SMITH_01284', 'SMITH_01299', 'SMITH_01308', 'SMITH_01311', 'SMITH_01320', 'SMITH_01328', 'SMITH_01332', 'SMITH_01334', 'SMITH_01367', 'SMITH_01371', 'SMITH_01381', 'SMITH_01383', 'SMITH_01389', 'SMITH_01390', 'SMITH_01391', 'SMITH_01396', 'SMITH_01405', 'SMITH_01417', 'SMITH_01427', 'SMITH_01430', 'SMITH_01440', 'SMITH_01447', 'SMITH_01452', 'SMITH_01472', 'SMITH_01478', 'SMITH_01487', 'SMITH_01491', 'SMITH_01497', 'SMITH_01509', 'SMITH_01511', 'SMITH_01517', 'SMITH_01520', 'SMITH_01553', 'SMITH_01559', 'SMITH_01566', 'SMITH_01567', 'SMITH_01570', 'SMITH_01575', 'SMITH_01577', 'SMITH_01578', 'SMITH_01608', 'SMITH_01612', 'SMITH_01624', 'SMITH_01627', 'SMITH_01631', 'SMITH_01639', 'SMITH_01646', 'SMITH_01665', 'SMITH_01676', 'SMITH_01678', 'SMITH_01681', 'SMITH_01695', 'SMITH_01704', 'SMITH_01707', 'SMITH_01723', 'SMITH_01736', 'SMITH_01745', 'SMITH_01749', 'SMITH_01753', 'SMITH_01754', 'SMITH_01764', 'SMITH_01765', 'SMITH_01768', 'SMITH_01770', 'SMITH_01774', 'SMITH_01781', 'SMITH_01799', 'SMITH_01806', 'SMITH_01810', 'SMITH_01815', 'SMITH_01818', 'SMITH_01821', 'SMITH_01824', 'SMITH_01827', 'SMITH_01864', 'SMITH_01873', 'SMITH_01880', 'SMITH_01889', 'SMITH_01894', 'SMITH_01904', 'SMITH_01907', 'SMITH_01917', 'SMITH_01937', 'SMITH_01947', 'SMITH_01972', 'SMITH_01986', 'SMITH_01999', 'SMITH_02000', 'SMITH_02014', 'SMITH_02019', 'SMITH_02024', 'SMITH_02028', 'SMITH_02032', 'SMITH_02037', 'SMITH_02049', 'SMITH_02059', 'SMITH_02062', 'SMITH_02075', 'SMITH_02076', 'SMITH_02081', 'SMITH_02087', 'SMITH_02092', 'SMITH_02104', 'SMITH_02107', 'SMITH_02109', 'SMITH_02110', 'SMITH_02118', 'SMITH_02120', 'SMITH_02130', 'SMITH_02133', 'SMITH_02137', 'SMITH_02143', 'SMITH_02156', 'SMITH_02159', 'SMITH_02214', 'SMITH_02217', 'SMITH_02224', 'SMITH_02225', 'SMITH_02229', 'SMITH_02230', 'SMITH_02244', 'SMITH_02250', 'SMITH_02261', 'SMITH_02268', 'SMITH_02274', 'SMITH_02281', 'SMITH_02288', 'SMITH_02291', 'SMITH_02293', 'SMITH_02294', 'SMITH_02346', 'SMITH_02353', 'SMITH_02360', 'SMITH_02375', 'SMITH_02376', 'SMITH_02390', 'SMITH_02392', 'SMITH_02394', 'SMITH_02402', 'SMITH_02448', 'SMITH_02454', 'SMITH_02459', 'SMITH_02557', 'SMITH_02574', 'SMITH_02582', 'SMITH_02620', 'SMITH_02688', 'SMITH_02690', 'SMITH_02702', 'SMITH_02743', 'SMITH_02768', 'SMITH_02789', 'SMITH_02806', 'SMITH_02864', 'SMITH_02872', 'SMITH_03005', 'SMITH_03080', 'SMITH_03143'),      'J2025+0316': ('SMITH_02278', 'SMITH_02279', 'SMITH_02349', 'SMITH_02358', 'SMITH_02362', 'SMITH_02397', 'SMITH_02426', 'SMITH_02430', 'SMITH_02433', 'SMITH_02455', 'SMITH_02467', 'SMITH_02479', 'SMITH_02492', 'SMITH_02509', 'SMITH_02514', 'SMITH_02526', 'SMITH_02534', 'SMITH_02550', 'SMITH_02562', 'SMITH_02569', 'SMITH_02577', 'SMITH_02621', 'SMITH_02670', 'SMITH_02676', 'SMITH_02679', 'SMITH_02686', 'SMITH_02689', 'SMITH_02693', 'SMITH_02721', 'SMITH_02726', 'SMITH_02742', 'SMITH_02745', 'SMITH_02746', 'SMITH_02759', 'SMITH_02764', 'SMITH_02781', 'SMITH_02792', 'SMITH_02801', 'SMITH_02813', 'SMITH_02824', 'SMITH_02831', 'SMITH_02832', 'SMITH_02842', 'SMITH_02885', 'SMITH_02892', 'SMITH_02893', 'SMITH_02894', 'SMITH_02906', 'SMITH_02907', 'SMITH_02912', 'SMITH_02915', 'SMITH_02956', 'SMITH_02981', 'SMITH_02986', 'SMITH_02987', 'SMITH_02994', 'SMITH_03004', 'SMITH_03006', 'SMITH_03007', 'SMITH_03020', 'SMITH_03042', 'SMITH_03046', 'SMITH_03072', 'SMITH_03079', 'SMITH_03090', 'SMITH_03092', 'SMITH_03094', 'SMITH_03104', 'SMITH_03110', 'SMITH_03111', 'SMITH_03112', 'SMITH_03117', 'SMITH_03120', 'SMITH_03128', 'SMITH_03135', 'SMITH_03146', 'SMITH_03151', 'SMITH_03154', 'SMITH_03172', 'SMITH_03177', 'SMITH_03179', 'SMITH_03181', 'SMITH_03206', 'SMITH_03208', 'SMITH_03213', 'SMITH_03217', 'SMITH_03220', 'SMITH_03259', 'SMITH_03272', 'SMITH_03295', 'SMITH_03300', 'SMITH_03307', 'SMITH_03311', 'SMITH_03313', 'SMITH_03317', 'SMITH_03318', 'SMITH_03335', 'SMITH_03357', 'SMITH_03369', 'SMITH_03384', 'SMITH_03388', 'SMITH_03391', 'SMITH_03393', 'SMITH_03394', 'SMITH_03399', 'SMITH_03404', 'SMITH_03409', 'SMITH_03411', 'SMITH_03428', 'SMITH_03444', 'SMITH_03474', 'SMITH_03514', 'SMITH_03538', 'SMITH_03545', 'SMITH_03552', 'SMITH_03562', 'SMITH_03568', 'SMITH_03575', 'SMITH_03577', 'SMITH_03581', 'SMITH_03582', 'SMITH_03597', 'SMITH_03600', 'SMITH_03608', 'SMITH_03626', 'SMITH_03635', 'SMITH_03654', 'SMITH_03660', 'SMITH_03664', 'SMITH_03665', 'SMITH_03666', 'SMITH_03699', 'SMITH_03721', 'SMITH_03732', 'SMITH_03739', 'SMITH_03746', 'SMITH_03747', 'SMITH_03750', 'SMITH_03756', 'SMITH_03759', 'SMITH_03760', 'SMITH_03764', 'SMITH_03773', 'SMITH_03780', 'SMITH_03781', 'SMITH_03787', 'SMITH_03795', 'SMITH_03796', 'SMITH_03797', 'SMITH_03798', 'SMITH_03814', 'SMITH_03816', 'SMITH_03830', 'SMITH_03850', 'SMITH_03857', 'SMITH_03858', 'SMITH_03864', 'SMITH_03877', 'SMITH_03889', 'SMITH_03897', 'SMITH_03900', 'SMITH_03902', 'SMITH_03920', 'SMITH_03922', 'SMITH_03923', 'SMITH_03924', 'SMITH_03937', 'SMITH_03961', 'SMITH_03965', 'SMITH_03972', 'SMITH_03975', 'SMITH_04042', 'SMITH_04077', 'SMITH_04079', 'SMITH_04080', 'SMITH_04095', 'SMITH_04127', 'SMITH_04141', 'SMITH_04206', 'SMITH_04242', 'SMITH_04254', 'SMITH_04257', 'SMITH_04262', 'SMITH_04281', 'SMITH_04286', 'SMITH_04303', 'SMITH_04357', 'SMITH_04359', 'SMITH_04363', 'SMITH_04364', 'SMITH_04383', 'SMITH_04386', 'SMITH_04413', 'SMITH_04443', 'SMITH_04462', 'SMITH_04466', 'SMITH_04468', 'SMITH_04501', 'SMITH_04547', 'SMITH_04570', 'SMITH_04576', 'SMITH_04594', 'SMITH_04599', 'SMITH_04601', 'SMITH_04612', 'SMITH_04615', 'SMITH_04629', 'SMITH_04681', 'SMITH_04699', 'SMITH_04731', 'SMITH_04732', 'SMITH_04763', 'SMITH_04774', 'SMITH_04777', 'SMITH_04785', 'SMITH_04802', 'SMITH_04819', 'SMITH_04843', 'SMITH_04868', 'SMITH_04869', 'SMITH_04874', 'SMITH_04878', 'SMITH_04891', 'SMITH_04895', 'SMITH_04920', 'SMITH_04924', 'SMITH_04933', 'SMITH_04938', 'SMITH_04947', 'SMITH_04955', 'SMITH_05006', 'SMITH_05015', 'SMITH_05055', 'SMITH_05062', 'SMITH_05095', 'SMITH_05131', 'SMITH_05178', 'SMITH_05202', 'SMITH_05209', 'SMITH_05262', 'SMITH_05266', 'SMITH_05270', 'SMITH_05351', 'SMITH_05381', 'SMITH_05401'), 'J2134-0153': ('SMITH_05313', 'SMITH_05343', 'SMITH_05362', 'SMITH_05363', 'SMITH_05416', 'SMITH_05437', 'SMITH_05445', 'SMITH_05456', 'SMITH_05515', 'SMITH_05534', 'SMITH_05547', 'SMITH_05599', 'SMITH_05604', 'SMITH_05647', 'SMITH_05668', 'SMITH_05683', 'SMITH_05693', 'SMITH_05728', 'SMITH_05730', 'SMITH_05738', 'SMITH_05781', 'SMITH_05834', 'SMITH_05837', 'SMITH_05934', 'SMITH_06009', 'SMITH_06031', 'SMITH_06077', 'SMITH_06113', 'SMITH_06123', 'SMITH_06130', 'SMITH_06138', 'SMITH_06140', 'SMITH_06157', 'SMITH_06178', 'SMITH_06354', 'SMITH_06395', 'SMITH_06417', 'SMITH_06436', 'SMITH_06539', 'SMITH_06568', 'SMITH_06584', 'SMITH_06597', 'SMITH_06670', 'SMITH_06701', 'SMITH_06852', 'SMITH_06854', 'SMITH_06887', 'SMITH_06999', 'SMITH_07025', 'SMITH_07095', 'SMITH_07202', 'SMITH_07228', 'SMITH_07737'), }









#####################################################################
#Set paramters that should be kept the same for the same night of observation
#####################################################################


datadir = '/students/sbetti/12B-336/'

sbnum='13594278'
projnum='12B-336'
#raw_ms_name = '12B-336.sb13594278.eb13717631.56219.96819922454.ash2012-1-18'
raw_ms_name = '12B-336.sb13594278.eb13717631.56219.96819922454'


refant='ea16' #reference antenna near the center of the array

#first time running this script? If turned on, will perform hanning smoothing, otherwise will skip that step
#firstrun=True 
firstrun=False
doplots=False

#wish to use a polarization model image to calibrate polang? 
#True: make model iamge(for 3C286, 3C138 and 3C48)
#False: use the known flux of 3C286, PA ,PF and SETJY to set the IQU flux
#polcal_using_model=True
polcal_using_model=False 

setup_scan ='1' #dummy scan to remove 

without_edgechan='5~58' # each SPW have 64 channels, this selection excludes the edge channels

binned_chan_width_MHz = 8 # Desired final output channel width in units of  MHz, float type
binned_chan_width_num = binned_chan_width_MHz/delnu # Desired final output channel width in numbers of raw channel, integer type

#Define calibrators 

myflux_standard='Perley-Butler 2013' # The flux standard to be used (Make sure when comparing fluxes taken different epoch that they are on the same flux density scale)
flux_calibrator_field='1331+305=3C286' # field ID of the absolute flux calibrator at the specific band
polang_calibrator_field='1331+305=3C286' #field ID of the polarization angle calibrator at the specific band
leakage_calibrator_field='J2355+4950' #field ID of the leakage calibrator at the specific band
gaincal_fields=flux_calibrator_field + ',' + leakage_calibrator_field + ',J1939-1002,J2025+0316,J2134-0153,0137+331=3C48' # different band have different gaincal_fields and targets..

gain_solution_interval='30s'

flux_calibrator_source ='3C286'
polang_calibrator_source='3C286'

cal_model=flux_calibrator_source+'_'+band+'.im' # model used in SETJY

### END OF INPUTS ####

### filenames derived from inputs
outputdir = datadir + 'output-' + sbnum + '/'
ms_file = datadir + raw_ms_name + '.ms'

split_ms_name = projnum + '.sb' + sbnum + '.autoflag.' + gain_solution_interval + '.' + str(binned_chan_width_MHz) + 'mhz'

split_ms_file = datadir + split_ms_name + '.ms'

####### DO CALIBRATION #########

execfile('calibration_initial.py')

#execfile('calibration_calib.py')

#execfile('calibration_apply.py')
