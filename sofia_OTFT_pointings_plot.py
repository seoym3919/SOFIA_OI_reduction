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
import math
import matplotlib.pyplot as plt
import os
import shutil
from astropy.time import Time
from astropy.coordinates import ICRS, Galactic, FK4, FK5

def silentremove(filename):
# remove files without raising error
	try:
		os.remove(filename)
	except OSError as e: # this would be "except OSError, e:" before Python 2.6
		if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
			raise # re-raise exception if a different error occurred

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

os.system('rm -rf OTF_COORD')
os.system('mkdir OTF_COORD')

dir_out = './OTF_COORD/'


try:
	with open('./'+project_name+'_OTFT_summary.xlsx') as file:
		pass
except IOError as e:
	print("summary file does not exist") #Does not exist OR no read permissions
	print("running summary routine")
	print("summary generated \n")

print('Reading summary... \n')

summary_OTFT = pd.read_excel('./'+project_name+'_OTFT_summary.xlsx')

SCAN    = np.array(summary_OTFT['SCAN'])
BEAM    = np.array(summary_OTFT['BEAM'])
OBSRA   = np.array(summary_OTFT['OBSRA'])
OBSDEC  = np.array(summary_OTFT['OBSDEC'])
LAMBDA  = np.array(summary_OTFT['LAMBDA'])
BETA    = np.array(summary_OTFT['BETA'])
LAMDEL  = np.array(summary_OTFT['LAMDEL'])
BETDEL  = np.array(summary_OTFT['BETDEL'])
RXDX    = np.array(summary_OTFT['RXDX'])
RXDY    = np.array(summary_OTFT['RXDY'])
OTFSLAM = np.array(summary_OTFT['OTFSLAM'])
OTFSBET = np.array(summary_OTFT['OTFSBET'])
OTFVLAM = np.array(summary_OTFT['OTFVLAM'])
OTFVBET = np.array(summary_OTFT['OTFVBET'])
OTFNDMP = np.array(summary_OTFT['OTFNDMP'])

BEAM_names = np.unique(BEAM)

FIG_NAME = dir_out+'OTFT_COORD.png'
fig = plt.figure(figsize=(8,9))
BEAM_color = ['r','g','b','magenta','orange','k','cyan']
for ib in BEAM_names:
	print('Estimating pointings in BEAM = '+ib+' ... ')

	mask_BEAM = (BEAM == ib)
	NSAMPLE = len(SCAN[mask_BEAM])
	print('number of scans = '+str(NSAMPLE)+' in BEAM = '+ib+' ... ')

	TOTAL_NUM_POINTS = OTFNDMP[mask_BEAM].sum()
	print('number of OTF dumps = '+str(TOTAL_NUM_POINTS)+' in BEAM = '+ib+' ... ')
	
	NDMP_PER_SAMPLE = OTFNDMP[mask_BEAM]

	OBSRA_BEAM = OBSRA[mask_BEAM]
	LAMBDA_BEAM = LAMBDA[mask_BEAM]
	LAMDEL_BEAM = LAMDEL[mask_BEAM]
	OTFSLAM_BEAM = OTFSLAM[mask_BEAM] 
	OTFVLAM_BEAM = OTFVLAM[mask_BEAM]
	RXDX_BEAM = RXDX[mask_BEAM]


	OBSDEC_BEAM = OBSDEC[mask_BEAM]
	BETA_BEAM = BETA[mask_BEAM]
	BETDEL_BEAM = BETDEL[mask_BEAM]
	OTFSBET_BEAM = OTFSBET[mask_BEAM]
	OTFVBET_BEAM = OTFVBET[mask_BEAM]
	RXDY_BEAM = RXDY[mask_BEAM] 
		
		
	COORD_L = np.array([])
	COORD_B = np.array([])
	RA_coord = np.array([])
	DEC_coord = np.array([])
	for isam in range(0,NSAMPLE):
		OTFT_raster = np.arange(NDMP_PER_SAMPLE[isam])

		radec_skycoord  =SkyCoord(OBSRA_BEAM[isam],OBSDEC_BEAM[isam], frame = FK5, unit=(u.hourangle, u.deg)) 

		LAMBDA_off_start = LAMDEL_BEAM[isam] + RXDX_BEAM[isam]
		BETA_off_start = BETDEL_BEAM[isam] + RXDY_BEAM[isam]

		COORD_L_new = LAMBDA_BEAM[isam]+(LAMBDA_off_start +(OTFVLAM_BEAM[isam]*OTFT_raster))*1./3600.
		COORD_B_new = BETA_BEAM[isam]+(BETA_off_start +(OTFVBET_BEAM[isam]*OTFT_raster))*1./3600.
		COORD_L = np.append(COORD_L, COORD_L_new) 
		COORD_B = np.append(COORD_B, COORD_B_new) 

	
	FNAME_OUT1 = dir_out+'RA_coord_'+ib+'_'+'{:06n}'.format(len(COORD_L))
	RA_coord.tofile(FNAME_OUT1)

	FNAME_OUT2 = dir_out+'DEC_coord_B_'+ib+'_'+'{:06n}'.format(len(COORD_B))
	DEC_coord.tofile(FNAME_OUT2)
	print('Saved OTFT coordinate data to '+FNAME_OUT1+' & '+FNAME_OUT2+' ... ')
	print('The number in file name denotes number of elements')
	print('\n')
	
	plt.scatter(COORD_L,COORD_B, s = 9, c=BEAM_color[np.int64(ib.split('PX')[1])])

plt.xlabel(r'RA (Degree)')
plt.ylabel(r'Dec (Degree)')
plt.title('OTF coordinates')
fig.tight_layout()
fig.savefig(FIG_NAME)
plt.close('all')

