"""
 SOFIA data reduction for upGREAT 7 pixel observation
 written by Youngmin Seo
 Created 08-01-2020
 Last updated 09-23-2021
 
 This script makes Level 2 OTF observation data by applying beam efficiency.  

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
import os
from astropy.io import fits
import matplotlib.pyplot as plt

# rebin module
def rebin(array,binsize):
	nx = len(array)
	nx_new = np.int64(nx/binsize)
	array_new = np.zeros(nx_new)
	for i0 in range(0,nx,binsize):
		i1 = np.int64(i0/binsize)
		array_new[i1] = np.mean(array[i0:i0+binsize])
	return array_new

# check whether or not there is project.par
# project.par is required file to setup the reduction
try:
    with open('./project.par') as file:
        pass
except IOError as e:
    print("Unable to open project.par file") #Does not exist OR no read permissions
    print("Please write project.par file. Line 1 is the project name and Line 2 is the absolute path to the data file")
    print("for example") 
    print("MONOCEROSR2")
    print("/home/xxxx/WORK/SOFIA/GREAT/OC6N/20181204_F533/raw/r1/data/si/cycle6/20181204_NADINE_HFA_LFA/")

# Load project parameters
par = np.loadtxt('./project.par',dtype = np.str)

project_name = par[0]
dir = par[1]
line_name = par[2]
beam2process = par[3].split(',')
basline_model_names = par[4].split(',')
nPCA_comp = np.int64(par[5])
nchan_new, nchan_bin, ichan_start = np.int64(par[6].split(','))
beam_eff = np.array(par[7].split(','),dtype = np.float)

# directory of Level 1c data
dir_OTFT_L1c = './TA_OTF_L1c/'


# Create directory for Level 2 data
os.system('rm -rf TA_OTF_L2')
os.system('mkdir TA_OTF_L2')
dir_out = './TA_OTF_L2/'

# Read Level 1c data
file = glob.glob(dir_OTFT_L1c+'*.fits')

# Get the file names
file_names = []
for i0 in range(len(file)):
	file_names.append(file[i0].split(dir_OTFT_L1c)[1])


# Get the scan number, sample number, band ID, and beam ID
scan_id = []
sample_id = []
band_id = []
beam_id = []
for i0 in range(len(file)):
	scan_id.append(file_names[i0].split('_')[0])
	sample_id.append(file_names[i0].split('_')[1])
	band_id.append(file_names[i0].split('_')[2])
	beam_id.append(file_names[i0].split('_')[3])

# total number of files
NFILE = len(file)

#
# Start of the beam efficiency
# The process will be done in beam by beam
#
spec_total = np.zeros(1024)
for ib in beam2process:

	# masking files with the beam ID of "ib"
	mask_BEAM = (np.array(beam_id) == ib)
	FILE_BEAM = np.array(file)[mask_BEAM]
	NFILE_BEAM = len(FILE_BEAM)
	FNAME_BEAM = np.array(file_names)[mask_BEAM]

	print('beam efficiency correction in BEAM = '+ib+' ... ')
#	
	# iteration through files with the beam ID of "ib"
	for i0 in range(0, NFILE_BEAM):
		hdu = fits.open(FILE_BEAM[i0])   # read in Level 1c data
		OTFNDMP = np.int64(hdu[0].header['OTFNDMP'])		# get number of OTF dumps within a single OTF file
		spec = hdu[0].data    # get OTF spectra
		NCHAN = np.int64(hdu[0].header['NUMCHAN']) # get number of channels

		mask_BEAM_name = (np.array(beam2process) == ib)  # get mask for beam
		beam_eff_single = np.array(beam_eff)[mask_BEAM_name] # get beam efficiency for Beam "ib"
		
		
		hdu[0].header['PROCSTAT'] = 'LEVEL_2'  # change the process state to L2
		file_OTFT_L2 = dir_out+FNAME_BEAM[i0]   # get L2 file name 
		hdu_out = fits.PrimaryHDU(spec/beam_eff_single, header = hdu[0].header)   # apply beam efficiency 
		hdu_out.writeto(file_OTFT_L2)		# save L2 spectra
		print('Saving beam efficiency corrected data. File name = '+ file_OTFT_L2)
	
		spec_total[:] = spec_total[:] + (spec/beam_eff_single).sum(axis=0) # accumulate spectra for plotting later
	print('\n')
#

# generate sample plot of OTF spectrum averaged over the entire region
FIG_NAME = './averaged_OTFT_total.png'
fig = plt.figure(figsize=(8,6))
vel = (np.arange(1024)-(hdu[0].header['CREFPIX']-1))*hdu[0].header['VELRES']+hdu[0].header['CREFVAL']

plt.plot(vel, spec_total/np.float64(NFILE_BEAM*OTFNDMP*6), drawstyle='steps')
plt.xlim([-10.,30])
plt.ylim([-1.,5])
plt.xlabel('LSR velocity [km/s]')
plt.ylabel('Main Beam Temperature [K]')

fig.tight_layout()
fig.savefig(FIG_NAME)
plt.close('all')

