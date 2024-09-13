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
from scipy import constants
import pandas as pd
import glob
from astropy.io import fits
import matplotlib.pyplot as plt

def rebin(array,binsize):
	nx = len(array)
	nx_new = np.int64(nx/binsize)
	array_new = np.zeros(nx_new)
	for i0 in range(0,nx,binsize):
		i1 = np.int64(i0/binsize)
		array_new[i1] = np.mean(array[i0:i0+binsize])
	return array_new

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

dir_out = dir_OTFT_L1b 

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


NFILE = len(file)
BEAM_names = np.unique(beam_id)
for ib in BEAM_names:
	mask_BEAM = (np.array(beam_id) == ib)
	FILE_BEAM = np.array(file)[mask_BEAM]
	NFILE_BEAM = len(FILE_BEAM)

	print('Plotting spectra samples in BEAM = '+ib+' ... \n')
	
	FIG_NAME = dir_out+'Samples_OTFT_L1b_'+ib+'.png'
	fig = plt.figure(figsize=(8,NFILE_BEAM))
	for i0 in range(0, NFILE_BEAM):
		hdu = fits.open(FILE_BEAM[i0])
		OTFNDMP = np.float64(hdu[0].header['OTFNDMP'])
		idump = np.int64(OTFNDMP*0.5)
		spec = hdu[0].data
		
		NCHAN = np.float64(hdu[0].header['NUMCHAN'])
		REFCHAN = np.float64(hdu[0].header['STARTREF'])
		RESTFREQ = np.float64(hdu[0].header['RESTFREQ'])
		FREQOFFS = np.float64(hdu[0].header['FREQOFFS'])
		VELRES = np.float64(hdu[0].header['VELRES'])
		vel = (np.arange(NCHAN)-REFCHAN)*VELRES - FREQOFFS/RESTFREQ*constants.c/1.e3 +hdu[0].header['RVSYS']
		vel_new = rebin(vel,8)
		spec_new = rebin(spec[idump,:],8)
		
		ax = plt.subplot(NFILE_BEAM, 1, i0+1)
		plt.plot(vel_new, spec_new, drawstyle='steps')
		plt.xlim([-10.,30])
		plt.ylim([-10.,30])
	
	fig.tight_layout()
	fig.savefig(FIG_NAME)
	plt.close('all')




