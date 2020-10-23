#!/usr/bin/env python
# coding: utf-8

#ahuthor:Denis Varise Bernardes.
#26/11/2019.

'''
This library contains the class of the calculation of the star FWHM. Initially, it is needed to provide to the class the name
of the FITS image, the (x,y) corrdinates of the star, and the maximum sky radius. The star image will be read and the function
get_max_count() will seek for the pixel within the sky radius with the maximum counts value. The function set_centroid() is
optional and it resets the star centroid based on the maximum pixel value found in previous step. Then, the function
calc_FWHM() will calculate the star FWHM. The star radius length will be three times the FWHM length.
'''

from math import exp, sqrt
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
from sys import exit
import pandas as pd
from scipy.interpolate import UnivariateSpline

class fwhm:

    def __init__(self, img_name, xy_star, sky_radius, bias_name='', exp_time = 0, dark_current_noise = 0, ccd_gain = 0, em_gain = 0, ccd_serial = 9916, ccd_temp = -70):
        #FWHM calculation
        self.img_name = img_name
        self.x = xy_star[0]
        self.y = xy_star[1]
        self.sky_radius = sky_radius
        
        #SNR calculation
        self.dark_noise = dark_current_noise
        self.ccd_gain = ccd_gain
        self.em_gain = em_gain
        self.ccd_serial = ccd_serial
        self.ccd_temp = ccd_temp
        self.bias_name = bias_name
        self.exp_time = exp_time
        
        self.n_pixels_star = 0
        self.sky_flux = 0        
        self.read_noise = 0                    
        self.img_data = 0
        self.sky_flux = 0                
        self.star_flux = 0             
        self.bias_data = 0        
        self.max_star_flux = 0
        self.fwhm = 0
        self.SNR = 0


        
    def read_star_img(self):  
        #read the star image
        self.img_data = fits.getdata(self.img_name).astype(float)
        self.header = fits.getheader(self.img_name)
        img_shape = self.img_data.shape
        if img_shape[0] == 1: self.img_data = self.img_data[0]            
           
        

    def get_max_count(self):         
        #read image dimensions
        img_shape = self.img_data.shape
        #create a mask of 1s
        working_mask = np.ones(img_shape,bool)
        #create indices arrays
        ym, xm = np.indices(img_shape, dtype='float32')
        #create a array with the center in the (x,y) provided values
        r = np.sqrt((xm - self.x)**2 + (ym - self.y)**2)

        #create a mask, setting the value 1 for those pixels between the radius 2*sky_radius and 3*sky_radius
        mask = (r > 2* self.sky_radius) * (r < 3 * self.sky_radius) * working_mask       
        #calculate the madian of the pixels
        self.sky_flux = np.median(self.img_data[np.where(mask)])
        #subtract the bias image from the science image
        self.img_data = self.img_data - self.sky_flux        
        
        #create a mask of 1s for those pixels within the sky radius
        mask = (r < self.sky_radius) * working_mask
        #get all the pixels of the science image within the sky radius
        star_flux = self.img_data[np.where(mask)]
        #return the maximum value of the star
        self.max_star_flux = max(star_flux)        

    

    def set_centroid(self):        
        #return the coordenates of the pixel with the maximum value
        y_max, x_max = np.where(self.img_data==self.max_star_flux)
        #set the found coordenates as the new star coordenates
        self.x, self.y  = x_max[0], y_max[0]        
        
         

    def calc_FWHM(self):                        
        radius = int(self.sky_radius/2)
        x,y = self.x, self.y
        img_data = self.img_data[y-radius:y+radius,x-radius:x+radius]        
        #get the pixels values for the line in the x coordinate of the centroid
        light_profile = np.take(img_data, radius, axis=0)
        #calculate the half the maximum value found
        half_max = self.max_star_flux/2        
        #create a vector with the pixels number of the cetroid line of the star
        x = np.linspace(0, len(light_profile), len(light_profile))
        #interpolate the found profile subtracted by the half maximum
        spline = UnivariateSpline(x, light_profile-half_max, s=0)
        #calculate where the interpolation reaches the x axis
        r1, r2 = spline.roots()
        #The FWHM will be the distance between the two roots found
        self.fwhm = r2-r1
        #the star radius will be 3*fwhm found
        self.star_radius = 3*self.fwhm        
        x,y = self.x+1, self.y+1
        return self.fwhm, self.star_radius, x, y


    def read_bias_img(self):
        self.bias_data = fits.getdata(self.bias_name).astype(float)
        bias_shape = self.bias_data.shape        
        if bias_shape[0] == 1: self.bias_data = self.bias_data[0]
        self.read_noise = np.std(self.bias_data)*self.ccd_gain


    def calc_dark_current(self):     
        #equations of the dark current characterization of the SPARC4 CCDs
        #Access https://github.com/DBernardes/FWHM.git for more details
        if self.dark_noise == 0:
            ccd_temp = self.ccd_temp
            if self.ccd_serial == 9914:
                self.dark_noise = 24.66*exp(0.0015*ccd_temp**2+0.29*ccd_temp) 
            if self.ccd_serial == 9915:
                self.dark_noise = 35.26*exp(0.0019*ccd_temp**2+0.31*ccd_temp)
            if self.ccd_serial == 9916:
                self.dark_noise = 9.67*exp(0.0012*ccd_temp**2+0.25*ccd_temp)
            if self.ccd_serial == 9917:
                self.dark_noise = 5.92*exp(0.0005*ccd_temp**2+0.18*ccd_temp)

            

    def read_exp_time(self):        
        #try to find the exposure time in the image header
        if self.exp_time == 0:
            try:
                string_texp = self.header['exposure'].split(',')
                string_texp = string_texp[0] + '.' + string_texp[1]        
                self.exp_time = float(string_texp)
            except:
                string_texp = self.header['exposure']                 
                self.exp_time = float(string_texp)


    def read_em_gain(self):        
        #try to find the EM gain in the image header
        if self.em_gain == 0:
            self.em_gain = 1
            try:
                string_em_gain = self.header['EMGAIN'].split(',')
                string_em_gain = string_texp[0] + '.' + string_texp[1]        
                self.em_gain = float(string_texp)
            except:1
            try:
                string_em_gain = self.header['EMGAIN']                 
                self.em_gain = float(string_texp)
            except:1
                  

        

    def calc_star_sky_flux(self):
        #the image needs to be read again as its centroid has been previously adjusted
        self.img_data = fits.getdata(self.img_name).astype(float)
        header = fits.getheader(self.img_name)
        img_shape = self.img_data.shape
        if img_shape[0] == 1: self.img_data = self.img_data[0]

        img_shape = self.img_data.shape
        #creat mask with 1s
        working_mask = np.ones(img_shape,bool)
        #create two indice arrays
        ym, xm = np.indices(img_shape, dtype='float32')        
        #create a array with the center in the (x,y) provided values
        r = np.sqrt((xm - self.x)**2 + (ym - self.y)**2)        
        #create a mask seeting as 1 those pixels within 2 star radius and 3 star radius
        mask = (r > 2* self.star_radius) * (r < 3 * self.star_radius) * working_mask
        #calculate the median of the pixels for the sku flux
        self.sky_flux = np.median(self.img_data[np.where(mask)])

        #create a mask for those pixels within the star radius     
        mask = (r < self.star_radius) * working_mask
        #get the star flux and subtracts it from the sky flux
        self.star_flux = (self.img_data[np.where(mask)] - self.sky_flux).sum()
        #number of the star pixels
        self.n_pixels_star = len(self.img_data[np.where(mask)])
        


    def calc_SNR(self):        
        gain = self.ccd_gain
        star = self.star_flux*gain
        self.sky_flux -= np.median(self.bias_data)
        sky = self.sky_flux*gain
        n_pix = self.n_pixels_star        
        exp_time = self.exp_time
        dc = self.dark_noise*exp_time
        rn = self.read_noise        

        #if the CCD is a EMCCD, apply the em_gain
        em_gain = self.em_gain 
        star*= em_gain
        sky *= em_gain
        dc  *= em_gain

        #SNR equation
        aux = np.sqrt(star + n_pix * (sky + dc + rn**2))
        self.SNR = star/aux

        return self.SNR, rn, sky/exp_time, star/exp_time, n_pix
    

   
