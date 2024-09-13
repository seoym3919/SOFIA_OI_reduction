"""
 SOFIA data reduction for upGREAT 7 pixel observation
 written by Youngmin Seo

 This script generates the summary of fits files, which is needed to prepare for the data reduction 

 project.par contains the project name and the directory of data
 If there is no project.par, one should write the file
 for example,
 
 MONOCEROSR2
 /home/xxxx/WORK/SOFIA/GREAT/OC6N/20181204_F533/raw/r1/data/si/cycle6/20181204_NADINE_HFA_LFA/
 OI_145_U
  
"""

from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
import numpy as np
import scipy as sp
import pandas as pd
import glob
import os
import shutil


try:
    with open('./project.par') as file:
        pass
except IOError as e:
    print("Error: Unable to open project.par file") #Does not exist OR no read permissions
    print("Please write project.par file. Line 1 is the project name, Line 2 is the absolute path to the data file, Line 3 is spectral line name" )
    print("for example" )
    print("MONOCEROSR2" )
    print("/home/xxxx/WORK/SOFIA/GREAT/OC6N/20181204_F533/raw/r1/data/si/cycle6/20181204_NADINE_HFA_LFA/" )
    print('OI_145_U' )

par = np.loadtxt('./project.par',dtype = np.str)
project_name = par[0]
dir = par[1]
line_name = par[2]

print('Start Ta conversion for '+project_name+'. The conversion will only for the total power observation. For OTF observation, please use TA_conversion_OI_OTF.py. ... \n')


os.system('rm -rf TA_OTF_L1b')
os.system('mkdir TA_OTF_L1b')

dir_out = './TA_OTF_L1b/'
dir_Tsys = './Tsys/'

CHECK_dir_Tsys = os.path.isdir(dir_Tsys)

if not CHECK_dir_Tsys:
	print('Error: Tsys is not calculated.')
	print('Please run Tsys_cal_OI.py first')
	exit()

try:
	with open('./'+project_name+'_summary.xlsx') as file:
		pass
except IOError as e:
	print("summary file does not exist") #Does not exist OR no read permissions
	print("running summary routine")
	print("summary generated \n")

print('Reading summary... \n')

summary = pd.read_excel('./'+project_name+'_summary.xlsx',sheet_name=line_name)

SCAN = np.array(summary['SCAN'])
INSTMODE = np.array(summary['INSTMODE'])
BEAM = np.array(summary['BEAM'])
SOBSMODE = np.array(summary['SOBSMODE'])
FNAME = np.array(summary['FILE NAME'])
LLOADSN = np.array(summary['LLOADSN'])
NUMCHAN = np.array(summary['NUMCHAN'])

print('These are the observation modes:')
print(np.unique(INSTMODE))
print('	')

if (INSTMODE == 'TP').any() == False:
	print('Error: Total power mode is in this data.')
	print('These are the observation modes:')
	print(np.unique(INSTMODE))
	print('For other observation mode please use different TA conversion scripts.')
	exit()

BEAM_names = np.unique(BEAM)
print('These are the beam to be reduced:')
print(BEAM_names)
print('	')

mask_ON = (SOBSMODE == 'ON')
mask_OFF = (SOBSMODE == 'OFF')


mask_OTFT = (INSTMODE == 'OTFT')
SCAN_OTFT = np.unique(SCAN[mask_OTFT])

mask_OTFT_OFF = (INSTMODE == 'OTFT-OFF')
SCAN_OTFT_OFF = np.unique(SCAN[mask_OTFT_OFF])


Tsys = np.zeros(np.int64(NUMCHAN[0])) 
for i0 in SCAN_OTFT:
	for ib in BEAM_names:
		print('Coverting Counts/s to Ta in SCAN = '+str(i0)+',  BEAM = '+ib+' ... \n')
		
		mask_i0 = (SCAN == i0)
		mask_ib = (BEAM == ib)
		mask_i0_ib = np.logical_and(mask_i0, mask_ib)

		mask_i0_ib_ON = np.logical_and(mask_i0_ib,mask_OTFT)
		mask_i0_ib_OFF = np.logical_and(mask_i0_ib,mask_OTFT_OFF)

		NUM_SAMPLE_ON = len(FNAME[mask_i0_ib_ON])
		NUM_SAMPLE_OFF = len(FNAME[mask_i0_ib_OFF])
		
		#print(NUM_SAMPLE_ON, NUM_SAMPLE_OFF)
		if (NUM_SAMPLE_ON == NUM_SAMPLE_OFF):
			for isam in range(0,NUM_SAMPLE_ON):
				#print('Processing '+FNAME[mask_i0_ib_ON][isam]+' & '+FNAME[mask_i0_ib_OFF][isam]+' ... \n')

				hdu_OFF = fits.open(dir+FNAME[mask_i0_ib_OFF][isam])
				int_OFF = np.array((hdu_OFF[0].data).reshape(hdu_OFF[0].header['NUMCHAN']))/hdu_OFF[0].header['EXPTIME']

				SCAN_Tsys =np.int64(hdu_OFF[0].header['LLOADSN'])
				file_Tsys = dir_Tsys+'Tsys_'+'{:06n}'.format(SCAN_Tsys)+'_'+ib
				CHECK_file_Tsys = os.path.isfile(file_Tsys)
				Tsys[:] = np.fromfile(file_Tsys)

				
				hdu_ON = fits.open(dir+FNAME[mask_i0_ib_ON][isam])
				int_ON = np.array((hdu_ON[0].data).reshape([hdu_ON[0].header['OTFNDMP'],hdu_ON[0].header['NUMCHAN']]))/hdu_ON[0].header['EXPTIME']
				int_TA_OTFT_L1b = np.zeros(int_ON.shape)
				for iOTF in range(0,hdu_ON[0].header['OTFNDMP']):
					int_TA_OTFT_L1b[iOTF,:] = Tsys * (int_ON[iOTF,:] - int_OFF[:])/int_OFF[:]
					#print((int_TA_OTFT_L1b[iOTF,:] == 0.).all())

		#		
				file_OTFT_L1b = dir_out+FNAME[mask_i0_ib_ON][isam]
				hdu_ON[0].header['PROCSTAT'] = 'LEVEL_1b'
				hdu_out = fits.PrimaryHDU(int_TA_OTFT_L1b, header = hdu_ON[0].header)
				hdu_out.writeto(file_OTFT_L1b)
		#
		#		#print('Ta-converted file saved to '+FNAME[mask_i0_ib_ON][isam]+' in '+dir_out+' ... \n')

print('Ta conversion done ... \n')

