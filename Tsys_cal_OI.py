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

print('Start Tsys cal for '+project_name+' ... \n')


os.system('rm -rf Tsys')
os.system('mkdir Tsys')

dir_out = './Tsys/'

try:
	with open('./'+project_name+'_summary.xlsx') as file:
		pass
except IOError as e:
	print("summary file does not exist") #Does not exist OR no read permissions
	print("running summary routine")
	print("summary generated \n")

print('Reading summary... \n')

summary = pd.ExcelFile('./'+project_name+'_summary.xlsx')
tab_names = summary.sheet_names

if (np.array(tab_names) == np.array(line_name)).any() == False:
	print('Error: Target line is not observed in the data')
	print('observed lines')
	print(tab_names)
	print('Target line name = '+line_name)
	exit()

summary = pd.read_excel('./'+project_name+'_summary.xlsx',sheet_name=line_name)

SCAN = np.array(summary['SCAN'])
INSTMODE = np.array(summary['INSTMODE'])
BEAM = np.array(summary['BEAM'])
SOBSMODE = np.array(summary['SOBSMODE'])
FNAME = np.array(summary['FILE NAME'])


print('These are the observation modes:')
print(np.unique(INSTMODE))
print('	')



BEAM_names = np.unique(BEAM)
print('These are the beams:')
print(BEAM_names)
print('	')


mask_HOT = (SOBSMODE == 'HOT')
mask_COL = (SOBSMODE == 'COL')


SCAN_HOT = SCAN[mask_HOT]
SCAN_COL = SCAN[mask_COL]
SCAN_HOTCOL = np.intersect1d(SCAN_HOT,SCAN_COL)
print('These are the calibration scans to be reduced:')
print(SCAN_HOTCOL)
print('	')

for i0 in SCAN_HOTCOL:
	for ib in BEAM_names:
		print('Calculating Tsys in SCAN = '+str(i0)+' BEAM = '+ib+' ... \n')

		mask_i0 = (SCAN == i0)
		mask_ib = (BEAM == ib)
		mask_i0_ib = np.logical_and(mask_i0, mask_ib)
		
		mask_i0_ib_HOT = np.logical_and(mask_i0_ib,mask_HOT)
		mask_i0_ib_COL = np.logical_and(mask_i0_ib,mask_COL)
		
		
		hdu_HOT = fits.open(dir+FNAME[mask_i0_ib_HOT][0])
		int_HOT = np.array((hdu_HOT[0].data).reshape(hdu_HOT[0].header['NUMCHAN']),dtype='float')/hdu_HOT[0].header['EXPTIME']
		THOT = hdu_HOT[0].header['THOT']

		hdu_COL = fits.open(dir+FNAME[mask_i0_ib_COL][0])
		int_COL = np.array((hdu_COL[0].data).reshape(hdu_COL[0].header['NUMCHAN']),dtype='float')/hdu_COL[0].header['EXPTIME']
		TCOL = hdu_COL[0].header['TCOLD']
		
		yy = int_HOT/int_COL
		Tsys_temp = (THOT - yy * TCOL)/(yy - 1.)
		
		FNAME_OUT = 'Tsys_'+'{:06n}'.format(i0)+'_'+ib
		Tsys_temp.tofile(dir_out+FNAME_OUT)
		print('Saved Tsys array to '+FNAME_OUT+' ... ')
		print('\n')

	
print('Tsys calculation done ... \n')


