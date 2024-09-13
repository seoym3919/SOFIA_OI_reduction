"""
 SOFIA data reduction for upGREAT 7 pixel observation
 written by Youngmin Seo

 This script makes Level 1c data by carrying out baseline correction using neural networks.  

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
import tensorflow as tf
from tensorflow import keras
from sklearn.decomposition import PCA
import pickle as pk
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


# Load project parameters
par = np.loadtxt('./project.par',dtype = np.str)

project_name = par[0]
dir = par[1]
line_name = par[2]
beam2process = par[3].split(',')
basline_model_names = par[4].split(',')
nPCA_comp = np.int64(par[5])
nchan_new, nchan_bin, ichan_start = np.int64(par[6].split(','))


# Load machine learning modules
model_autoencoder = keras.models.load_model(basline_model_names[0])
model_pca = pk.load(open(basline_model_names[1],'rb'))
model_eigen_correct = keras.models.load_model(basline_model_names[2])

# directory for Level 1b data
dir_OTFT_L1b = './TA_OTF_L1b/'

# Create directory for Level 1c data
os.system('rm -rf TA_OTF_L1c')
os.system('mkdir TA_OTF_L1c')
dir_out = './TA_OTF_L1c/'

# Read Level 1b data
file = glob.glob(dir_OTFT_L1b+'*.fits')

# Get the file names
file_names = []
for i0 in range(len(file)):
	file_names.append(file[i0].split(dir_OTFT_L1b)[1])

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
# Start of the baseline correction
# The process will be done in beam by beam
#
for ib in beam2process:
	# masking files with the beam ID of "ib"
	mask_BEAM = (np.array(beam_id) == ib)
	FILE_BEAM = np.array(file)[mask_BEAM]
	NFILE_BEAM = len(FILE_BEAM)
	FNAME_BEAM = np.array(file_names)[mask_BEAM]

	print('Correcting baselines of spectra in BEAM = '+ib+' ... \n')
#	
	# generate sample plots to check the quality of baseline corrections
	FIG_NAME = dir_out+'Baseline_correct_OTFT_L1c_'+ib+'.png'
	fig = plt.figure(figsize=(8,NFILE_BEAM))
	
	# iteration through files with the beam ID of "ib"
	for i0 in range(0, NFILE_BEAM):
		hdu = fits.open(FILE_BEAM[i0])
		OTFNDMP = np.int64(hdu[0].header['OTFNDMP'])
		spec = hdu[0].data
		NCHAN = np.int64(hdu[0].header['NUMCHAN'])
		
		# preparing for binning and baseline correction using neural networks  
		spec_new = np.zeros([OTFNDMP, np.int64(NCHAN/nchan_bin)])
		int_max = np.zeros(OTFNDMP)
		int_min = np.zeros(OTFNDMP)

		# iteration within the OTF dumps and normalize spectra
		for idmp in range(0,OTFNDMP):
			temp = rebin(spec[idmp,:],nchan_bin) # rebinning
			# The machine learning algoritms process normalized data whose values are from 0 to 1.
			int_max[idmp] = np.max(temp[ichan_start:ichan_start+nchan_new]) # maximum intensity
			int_min[idmp] = np.min(temp[ichan_start:ichan_start+nchan_new]) # minimum intensity
			spec_new[idmp,:] = (temp-int_min[idmp])/(int_max[idmp]-int_min[idmp]) # normalized spectra
		

		spectra_chunk = spec_new[:,ichan_start:ichan_start+nchan_new].reshape(OTFNDMP,nchan_new,1)
		spectra_decoded = model_autoencoder.predict(spectra_chunk)
		spectra_pca = (model_pca.transform(spectra_decoded.reshape([OTFNDMP,nchan_new]))).reshape(OTFNDMP, nPCA_comp)
		inputs_prediction = spectra_pca[0:OTFNDMP,:]
		baseline_decoded = (model_eigen_correct.predict(inputs_prediction)).reshape([OTFNDMP,nchan_new])

		spec_temp = spec_new.copy()
		spec_temp[:,ichan_start:ichan_start+nchan_new] = spec_new[:,ichan_start:ichan_start+nchan_new]-baseline_decoded[:,:]
		spec_temp[:,0:ichan_start] = 0.
		spec_temp[:,ichan_start+nchan_new:np.int64(NCHAN/nchan_bin)] = 0.

		spec_corrected = spec_temp.copy()
		for i1 in range(0,OTFNDMP):
			spec_corrected[i1,:] = spec_temp[i1,:]*(int_max[i1]-int_min[i1])
		
		REFCHAN = np.float64(hdu[0].header['STARTREF'])
		RESTFREQ = np.float64(hdu[0].header['RESTFREQ'])
		FREQOFFS = np.float64(hdu[0].header['FREQOFFS'])
		VELRES = np.float64(hdu[0].header['VELRES'])
		vel = (np.arange(NCHAN)-REFCHAN)*VELRES - FREQOFFS/RESTFREQ*constants.c/1.e3 +hdu[0].header['RVSYS']
		vel_new = rebin(vel,nchan_bin)[ichan_start:ichan_start+nchan_new]

		NCHAN_NEW_vel = np.int64(len(vel_new))
		VELRES_NEW = vel_new[np.int64(NCHAN_NEW_vel/2.)+1]-vel_new[np.int64(NCHAN_NEW_vel/2.)]
		header_new = (hdu[0].header).copy()


		# making new FITS header
		header_new['NAXIS1']  = NCHAN_NEW_vel
		header_new['PROCSTAT'] = 'LEVEL_1c'
		header_new['NUMCHAN'] = NCHAN_NEW_vel
		header_new['CREFPIX'] = np.int64(NCHAN_NEW_vel/2)
		header_new['STARTREF'] = np.int64(NCHAN_NEW_vel/2)
		header_new['VELRES'] = VELRES_NEW
		header_new.set('CREFVAL',value = vel_new[np.int64(NCHAN_NEW_vel/2)], after = 'CREFPIX')

		file_OTFT_L1c = dir_out+FNAME_BEAM[i0]
		hdu_out = fits.PrimaryHDU(spec_corrected[:,ichan_start:ichan_start+nchan_new], header = header_new)
		hdu_out.writeto(file_OTFT_L1c)		
		print('Saving baseline corrected data. File name = '+ file_OTFT_L1c)


		# plotting samples and baseline correctios for the quality check of baseline correction
		idump = np.int64(OTFNDMP*0.5)
		ax = plt.subplot(NFILE_BEAM, 1, i0+1)
		plt.plot(vel_new, spec_new[idump,ichan_start:ichan_start+nchan_new]*(int_max[idump]-int_min[idump])+int_min[idump], drawstyle='steps')
		plt.plot(vel_new, baseline_decoded[idump,:]*(int_max[idump]-int_min[idump])+int_min[idump], drawstyle='steps')
		plt.xlim([-30.,30])
		plt.ylim([-20.,40.])
	
	fig.tight_layout()
	fig.savefig(FIG_NAME)
	plt.close('all')

		
