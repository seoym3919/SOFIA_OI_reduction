from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import wcs
import numpy as np
import scipy as sp
import glob
import math
from astropy import constants as const
from grid_otf import grid_otf
import sys
from progressbar import ProgressBar

import math
import numpy
import sys
from scipy import interpolate
import scipy
import time
import os,errno

def silentremove(filename):
# remove files without raising error
    try:
        os.remove(filename)
    except OSError as e: # this would be "except OSError, e:" before Python 2.6
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred

#
# begine of main program
#	

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
line_name = par[2]
beam2process = par[3].split(',')
basline_model_names = par[4].split(',')
nPCA_comp = np.int64(par[5])
nchan_new = np.int64(par[6])
beam_eff = np.array(par[7].split(','),dtype = np.float)

dir_OTFT_L2 = './TA_OTF_L2/'

os.system('rm -rf TA_OTF_L3')
os.system('mkdir TA_OTF_L3')
dir_out = './TA_OTF_L3/'


file = glob.glob(dir_OTFT_L2+'*.fits')

file_names = []
for i0 in range(len(file)):
	file_names.append(file[i0].split(dir_OTFT_L2)[1])

scan_id = []
sample_id = []
band_id = []
beam_id = []
for i0 in range(len(file)):
	scan_id.append(file_names[i0].split('_')[0])
	sample_id.append(file_names[i0].split('_')[1])
	band_id.append(file_names[i0].split('_')[2])
	beam_id.append(file_names[i0].split('_')[3])


NFILE = np.int64(len(file))
hdu_temp = fits.open(file[0])
NDMP = np.int64(hdu_temp[0].header['OTFNDMP'])
NCHAN = np.int64(hdu_temp[0].header['NUMCHAN'])

vel = (np.arange(NCHAN)-(hdu_temp[0].header['CREFPIX']-1))*hdu_temp[0].header['VELRES']+hdu_temp[0].header['CREFVAL']


spec_array = np.zeros([NFILE*NDMP, NCHAN])
xpos = np.zeros([NFILE*NDMP])
ypos = np.zeros([NFILE*NDMP])
nchan_array = np.zeros([NFILE*NDMP])

for i0 in range(0,NFILE):
	hdu_temp = fits.open(file[i0])
	spec_array[i0*NDMP:(i0+1)*NDMP,:] = hdu_temp[0].data[:,:]

	OBSRA  = hdu_temp[0].header['OBSRA']
	OBSDEC  = hdu_temp[0].header['OBSDEC']
	LAMBDA  = hdu_temp[0].header['LAMBDA']
	BETA    = hdu_temp[0].header['BETA']
	LAMDEL  = hdu_temp[0].header['LAMDEL']
	BETDEL  = hdu_temp[0].header['BETDEL']
	RXDX    = hdu_temp[0].header['RXDX']
	RXDY    = hdu_temp[0].header['RXDY']
	OTFSLAM = hdu_temp[0].header['OTFSLAM']
	OTFSBET = hdu_temp[0].header['OTFSBET']
	OTFVLAM = hdu_temp[0].header['OTFVLAM']
	OTFVBET = hdu_temp[0].header['OTFVBET']

	OTFT_raster = np.arange(NDMP)
	COORD_L = OBSRA + LAMBDA+ LAMDEL + RXDX + OTFSLAM + OTFVLAM * OTFT_raster
	COORD_B = OBSDEC+ BETA  + BETDEL + RXDY + OTFSBET + OTFVBET * OTFT_raster

	xpos[i0*NDMP:(i0+1)*NDMP] = COORD_L[:]/3600.
	ypos[i0*NDMP:(i0+1)*NDMP] = COORD_B[:]/3600.


nchan_array[:] = NCHAN

header_sample = (hdu_temp[0].header).copy()


# dish size of SOFIA in cm
dish_diam = 250.
#
#
# wavelength of lines in cm
restfreq = header_sample['RESTFREQ']*1.e6
wavelength = const.c.cgs.value/restfreq
freq = restfreq
#
# beam size in array 
beam_fwhm = 1.22 * wavelength/dish_diam * np.rad2deg(1.)
#
#
# pixel size
pixPerBeam = 4.
pix_scale = beam_fwhm/pixPerBeam

#
# create header, wcs, and image size from header and given parameters
# image size
xRange = np.max(xpos)-np.min(xpos)
yRange = np.max(ypos)-np.min(ypos)
xsize = np.int64(math.ceil(xRange*1.1/pix_scale))+20
ysize = np.int64(math.ceil(yRange*1.1/pix_scale))+20

# set image center
refXsky = np.min(xpos)+0.5*xRange
refYsky = np.min(ypos)+0.5*yRange
refXpix = math.ceil(xsize*0.5)
refYpix = math.ceil(ysize*0.5)

# set coordinates and projection
xcoord = 'RA'
ycoord = 'DEC'
specSysDict = {'OBS':'TOPOCENT','GEO':'GEOCENTR','BAR':'BARYCENT','HEL':'HELIOCEN','GAL':'GALACTOC','LSD':'LSRD','LSR':'LSRK','LGR':'LOCALGRP','COB':'CMBDIPOL'}
coordType = [xcoord,ycoord]
radesys = ''
equinox = 0.
veldef = 'RADI'
specsys = specSysDict['LSR']


hdr = fits.Header()
# BASIC stuff, the WCS code needs this
hdr['SIMPLE'] = True
hdr['BITPIX'] = -64
hdr['NAXIS'] = 4
hdr['NAXIS1'] = xsize
hdr['NAXIS2'] = ysize
hdr['NAXIS3'] = NCHAN
hdr['NAXIS4'] = 1
hdr['BUNIT']   = 'K (Tmb)     '                           

ctypeDashes = '----'
###
xctype = coordType[0] + ctypeDashes[len(coordType[0]):]
yctype = coordType[1] + ctypeDashes[len(coordType[1]):]

# MAKE THE POSITION AXES##
proj = 'GLS'
hdr['CTYPE1'] = xctype + '-' + proj
hdr['CRVAL1'] = refXsky
hdr['CRPIX1'] = refXpix
hdr['CDELT1'] = -1.0*pix_scale
hdr['CROTA1'] = 0.
#
hdr['CTYPE2'] = yctype + '-' + proj
hdr['CRVAL2'] = refYsky
hdr['CRPIX2'] = refYpix
hdr['CDELT2'] = pix_scale
hdr['CROTA2'] = 0.
#
hdr['CTYPE3'] = 'VELO-LSR'
hdr['CUNIT3'] = 'km/s'
hdr['CRVAL3'] = vel[0]
hdr['CRPIX3'] = 1.0
hdr['CDELT3'] = (vel[1]-vel[0])
hdr['CROTA3'] = 0.
#
# STOKES axis - always I
hdr['CTYPE4'] = 'STOKES'
hdr['CRVAL4'] = 1.0
hdr['CRPIX4'] = 1.0
hdr['CDELT4'] = 1.0
hdr['CROTA4'] = 0.
#
#
#
hdr['RA']       = hdu_temp[0].header['OBSRA']                        
hdr['DEC']      = hdu_temp[0].header['OBSDEC']                            
hdr['EQUINOX']  = 0.2000000000000E+04
hdr['LINE']     = line_name                                      
hdr['OBJECT']   = project_name                                      
hdr['RESTFRQ'] = restfreq
hdr['VELO-LSR'] =  0.5000000000000E+04
hdr['VELREF']   =                  257
hdr['VELRES']   =     0.20000000298023                    
hdr['VELDEF']   = 'RADI-LSR'                              
hdr['PROCSTAT'] = 'LEVEL_4 '                              
hdr['SPECSYS'] = specsys

# create wcs object from STO2 header (non-trivial header)
wcsObj = wcs.WCS(hdr,relax=True)

#
# create spectral map 
cube, weight, beam_size = grid_otf(spec_array, xpos, ypos, wcsObj, NCHAN, xsize, ysize, pix_scale, beam_fwhm)



silentremove(dir_out+project_name+'_OTFT.fits')
hdu_cube_out = fits.PrimaryHDU(cube, header = hdr)
hdu_cube_out.writeto(dir_out+project_name+'_OTFT.fits')

