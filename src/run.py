#!/usr/bin/env python
# coding: utf-8

# This file executes the software the of the FWHM calculation.

#Denis Varise Bernardes.
#22/10/2020.

import FWHM
import os
from sys import exit

file_dir = os.path.dirname(os.getcwd()) + '\Sample images\\'
file_name = 'star.fits'
bias_name = 'bias.fits'

FWHM_obj = FWHM.fwhm(img_name = file_dir + '\\' +  file_name,                     
                     xy_star = (507,369),
                     sky_radius = 20,
                     bias_name = file_dir + '\\' +  bias_name,
                     ccd_gain = 3.36)

FWHM_obj.read_star_img()
FWHM_obj.get_max_count()
FWHM_obj.set_centroid()
fwhm, star_radius, x, y = FWHM_obj.calc_FWHM()

FWHM_obj.read_bias_img()
FWHM_obj.calc_dark_current()
FWHM_obj.read_exp_time()
FWHM_obj.read_em_gain()
FWHM_obj.calc_star_sky_flux()
snr, rn, sky_flux, star_flux, n_pixels, bias_level = FWHM_obj.calc_SNR()

print('\n')
print('FWHM: ', round(fwhm,2), 'pixels')
print('Star Radius:', round(star_radius), 'pixels')
print('Centroide: %i,%i'%(x,y))
print('SNR: ', round(snr,2))
print('Sky Flux: ', round(sky_flux,2), 'photons/s')
print('Star Flux:', round(star_flux,2), 'photons/s')
print('Star Pixels:', n_pixels)
