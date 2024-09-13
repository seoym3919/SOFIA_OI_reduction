"""
 SOFIA data reduction for upGREAT 7 pixel observation
 written by Youngmin Seo

 This script generates the summary of fits files, which is needed to prepare for the data reduction 

 project.par contains the project name and the directory of data
 If there is no project.par, one should write the file
 for example,
 
 MONOCEROSR2
 /home/xxxx/WORK/SOFIA/GREAT/OC6N/20181204_F533/raw/r1/data/si/cycle6/20181204_NADINE_HFA_LFA/
  
"""

import numpy as np
import scipy as sp
import pandas as pd
import glob
from astropy.io import fits

try:
    with open('./project.par') as file:
        pass
except IOError as e:
    print("Unable to open project.par file") #Does not exist OR no read permissions
    print("Please write project.par file. Line 1 is the project name and Line 2 is the absolute path to the data file")
    print("for example") 
    print("MONOCEROSR2")
    print("/home/xxxx/WORK/SOFIA/GREAT/OC6N/20181204_F533/raw/r1/data/si/cycle6/20181204_NADINE_HFA_LFA/")

par = np.loadtxt('./project.par',dtype = np.str)
project_name = par[0]
dir = par[1]

dir_OTFT_L1b = './TA_OTF_L1b/'

file = glob.glob(dir_OTFT_L1b+'*.fits')



file_names = []
for i0 in range(len(file)):
	file_names.append(file[i0].split(dir_OTFT_L1b)[1])

scan_id = []
sample_id = []
band_id = []
beam_id = []
for i0 in range(len(file)):
	scan_id.append(file_names[i0].split('_')[0])
	sample_id.append(file_names[i0].split('_')[1])
	band_id.append(file_names[i0].split('_')[2])
	beam_id.append(file_names[i0].split('_')[3])


PLANID_id = []
TRACMODE_id = []
CHPSYM_id = []
CHOPPING_id = []
NODDING_id = []
MAPPING_id = []
INSTRUME_id = []
OBSRA_id = []
OBSDEC_id = []
INSTMODE_id = []
SOBSMODE_id = []
LLOADSN_id = []
NODN_id = []
LFRCASN_id = []
LSKYSN_id = []
OBJECT_id = []
COORDSYS_id = []
LAMBDA_id = []
BETA_id = []
LAMDEL_id = []
BETDEL_id = []
LINE_id = []
RESTFREQ_id = []
NUMCHAN_id = []
VELRES_id = []
REFPOSX_id = []
REFPOSY_id = []
REFRXDX_id = [] 
REFRXDY_id = [] 
RXDX_id = []
RXDY_id = []
THOT_id = []
TCOLD_id = []
TAMB_id = []
TEMP_OUT_id = []
EXPTIME_id = []
SCAN_id = []
OTFSLAM_id = []
OTFSBET_id = []
OTFVLAM_id = []
OTFVBET_id = []
OTFNDMP_id = []
POSANGLE_id = []
BEAMANG_id = []


for i0 in range(len(file)):
	hdu = fits.open(file[i0])

	PLANID_id.append(hdu[0].header['PLANID'])
	TRACMODE_id.append(hdu[0].header['TRACMODE'])
	MAPPING_id.append(hdu[0].header['MAPPING'])
	INSTRUME_id.append(hdu[0].header['INSTRUME'])
	OBSRA_id.append(hdu[0].header['OBSRA'])
	OBSDEC_id.append(hdu[0].header['OBSDEC'])
	INSTMODE_id.append(hdu[0].header['INSTMODE'])
	SOBSMODE_id.append(hdu[0].header['SOBSMODE'])
	LLOADSN_id.append(hdu[0].header['LLOADSN'])
	LFRCASN_id.append(hdu[0].header['LFRCASN'])
	OBJECT_id.append(hdu[0].header['OBJECT'])
	COORDSYS_id.append(hdu[0].header['COORDSYS'])
	LAMBDA_id.append(hdu[0].header['LAMBDA'])
	BETA_id.append(hdu[0].header['BETA'])
	LAMDEL_id.append(hdu[0].header['LAMDEL'])
	BETDEL_id.append(hdu[0].header['BETDEL'])
	LINE_id.append(hdu[0].header['LINE'])
	RESTFREQ_id.append(hdu[0].header['RESTFREQ'])
	NUMCHAN_id.append(hdu[0].header['NUMCHAN'])
	VELRES_id.append(hdu[0].header['VELRES'])
	REFPOSX_id.append(hdu[0].header['REFPOSX'])
	REFPOSY_id.append(hdu[0].header['REFPOSY'])
	REFRXDX_id.append(hdu[0].header['REFRXDX'])
	REFRXDY_id.append(hdu[0].header['REFRXDY'])
	RXDX_id.append(hdu[0].header['RXDX'])
	RXDY_id.append(hdu[0].header['RXDY'])
	EXPTIME_id.append(hdu[0].header['EXPTIME'])
	SCAN_id.append(hdu[0].header['SCAN'])
	OTFSLAM_id.append(hdu[0].header['OTFSLAM'])
	OTFSBET_id.append(hdu[0].header['OTFSBET'])
	OTFVLAM_id.append(hdu[0].header['OTFVLAM'])
	OTFVBET_id.append(hdu[0].header['OTFVBET'])
	OTFNDMP_id.append(hdu[0].header['OTFNDMP'])
	POSANGLE_id.append(hdu[0].header['POSANGLE'])
	BEAMANG_id.append(hdu[0].header['BEAMANG'])


column_names = ['SCAN','BAND','BEAM','PLANID','TRACMODE','MAPPING','INSTRUME','OBSRA','OBSDEC','INSTMODE','SOBSMODE','LLOADSN','LFRCASN','LAMBDA','BETA','LAMDEL','BETDEL','LINE','RESTFREQ','NUMCHAN','VELRES','REFPOSX','REFPOSY','REFRXDX','REFRXDY','RXDX','RXDY','OTFSLAM','OTFSBET','OTFVLAM','OTFVBET','OTFNDMP','POSANGLE', 'BEAMANG', 'EXPTIME','FILE NAME']

data = np.array([SCAN_id,band_id,beam_id,PLANID_id,TRACMODE_id,MAPPING_id,INSTRUME_id,OBSRA_id,OBSDEC_id,INSTMODE_id,SOBSMODE_id,LLOADSN_id,LFRCASN_id,LAMBDA_id,BETA_id,LAMDEL_id,BETDEL_id,LINE_id,RESTFREQ_id,NUMCHAN_id,VELRES_id,REFPOSX_id,REFPOSY_id,REFRXDX_id,REFRXDY_id,RXDX_id,RXDY_id,OTFSLAM_id,OTFSBET_id,OTFVLAM_id,OTFVBET_id,OTFNDMP_id,POSANGLE_id, BEAMANG_id, EXPTIME_id,file_names]).T

line_types = np.unique(np.array(LINE_id))
mask =  np.zeros([line_types.shape[0],data.shape[0]],dtype='bool')
for i0 in range(line_types.shape[0]):
	mask[i0,:] = (np.array(LINE_id) == line_types[i0])

writer = pd.ExcelWriter('./'+project_name+'_OTFT_summary.xlsx', engine='openpyxl')
df = pd.DataFrame(data, columns = column_names)
df.to_excel(writer, sheet_name='OTFT_summary')

writer.save()

