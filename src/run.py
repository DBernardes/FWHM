#!/usr/bin/env python
# coding: utf-8

# This file executes the software the of the FWHM calculation.

#Denis Varise Bernardes.
#22/10/2020.

import FWHM

file_dir = r'C:\Users\observer\Desktop\FWHM\Sample images'
file_name = 'star.fits'
bias_name = 'bias.fits'

FWHM_obj = FWHM.fwhm(img_name = file_dir + '\\' +  file_name,                     
                     xy_star = (507,369),
                     sky_radius = 20,
                     bias_name = file_dir + '\\' +  bias_name,
                     ccd_gain = 3.36,
                     em_gain = 0)

FWHM_obj.read_star_img()
FWHM_obj.get_max_count()
FWHM_obj.set_centroid()
fwhm, star_radius, x, y = FWHM_obj.calc_FWHM()

FWHM_obj.read_bias_img()
FWHM_obj.calc_dark_current()
FWHM_obj.read_exp_time()
FWHM_obj.read_em_gain()
FWHM_obj.calc_star_sky_flux()
snr, rn, sky_flux, star_flux, n_pixels = FWHM_obj.calc_SNR()

print('\n')
print('FWHM: ', fwhm)
print('Star Radius:', star_radius)
print('Centroide: %i,%i'%(x,y))
print('SNR: ', snr)
print('Sky Flux: ', sky_flux, 'photons/s')
print('Star Flux:', star_flux, 'photons/s')
print('Star Pixels:', n_pixels)
