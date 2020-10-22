#!/usr/bin/env python
# coding: utf-8

# Esta biblioteca contem a classe que calcula a relacao sinal-ruido de uma das imagens do janderson,
# assim como o fluxo de fotons da estrela e do ceu.
#Denis Varise Bernardes.
#26/11/2019.

from math import exp, sqrt
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np
from sys import exit
import Read_Noise_Calculation_Bib as RNC
import pandas as pd
from scipy.interpolate import UnivariateSpline

class PhotonFluxCalc:

    def __init__(self, img_name, bias_name, xy_star, sky_radius, ccd_serial = 9916):        
        self.exp_time = 0
        self.n_pixels_object = 1
        self.sky_flux = 0
        self.dark_noise = 0
        self.read_noise = 0
        #self.SNR = 0
        self.img_name = img_name
        self.x = xy_star[0]
        self.y = xy_star[1]
        self.sky_radius = sky_radius        
        self.img_data = 0
        self.sky_flux = 0                
        self.star_flux = 0        
        self.ganho = 0
        self.ccd_serial = ccd_serial        
        self.bias_data = 0
        self.bias_name = bias_name
        self.fwhm = 0
        self.max_star_flux = 0


    def read_bias_img(self):
        self.bias_data = fits.getdata(self.bias_name).astype(float)
        bias_shape = self.bias_data.shape        
        if bias_shape[0] == 1: self.bias_data = self.bias_data[0]        
        
        
    def read_img_get_info(self):
        #Esta funcao faz a leitura da imagem fornecida,
        # cria uma mascara com o tamanho da imagem
        # e dois arrays de indices com o tamanho da imagem
        # Tambem sao calculados os valores do ruido de RN e DC para o modo de operacao        
        self.img_data = fits.getdata(self.img_name).astype(float)
        header = fits.getheader(self.img_name)
        img_shape = self.img_data.shape
        if img_shape[0] == 1: self.img_data = self.img_data[0]       

        #------------------------------------------------------------------------------
        try:
            string_ccd_temp = header['temp'].split(',')
            string_ccd_temp = string_ccd_temp[0] + '.' + string_ccd_temp[1]
            ccd_temp = float(string_ccd_temp)
        except:
            string_ccd_temp = header['temp']            
            ccd_temp = float(string_ccd_temp)
        
        
        #equacao tirada do artigo de caract. dos CCDs
        if self.ccd_serial == 9914:
            self.dark_noise = 24.66*exp(0.0015*ccd_temp**2+0.29*ccd_temp) 
        if self.ccd_serial == 9915:
            self.dark_noise = 35.26*exp(0.0019*ccd_temp**2+0.31*ccd_temp)
        if self.ccd_serial == 9916:
            self.dark_noise = 9.67*exp(0.0012*ccd_temp**2+0.25*ccd_temp)
        if self.ccd_serial == 9917:
            self.dark_noise = 5.92*exp(0.0005*ccd_temp**2+0.18*ccd_temp)

        #------------------------------------------------------------------------------
        try:
            string_texp = header['exposure'].split(',')
            string_texp = string_texp[0] + '.' + string_texp[1]        
            self.exp_time = float(string_texp)
        except:
            string_texp = header['exposure']                 
            self.exp_time = float(string_texp)
        #------------------------------------------------------------------------------
                
        hss = int(1/(float(header['readtime'])*1e6)) #em MHz
        
        try:preamp = int(header['preamp'].split('x')[0])
        except:preamp = int(header['preamp'])
        
        binn = int(header['hbin'])
        em_mode = header['outptamp']        
        if em_mode == 'Conventional':em_mode=0
        else:em_mode=1
                
        RN = RNC.ReadNoiseCalc()        
        RN.write_operation_mode(em_mode, 2, hss, preamp, binn)        
        RN.calc_read_noise()        
        self.read_noise = float(RN.noise)        

        #------------------------------------------------------------------------------
        gain = 0
        if em_mode == 1:
            if hss == 30:
                if preamp == 1:
                    gain = 17.2
                if preamp == 2:
                    gain = 5.27
            if hss == 20:
                if preamp == 1:
                    gain = 16.4
                if preamp == 2:
                    gain = 4.39
            if hss == 10:
                if preamp == 1:
                    gain = 16.0
                if preamp == 2:
                    gain = 3.96
            if hss == 1:
                if preamp == 1:
                    gain = 15.9
                if preamp == 2:
                    gain = 3.88
        else:
            if hss == 1:
                if preamp == 1:
                    gain = 3.37
                if preamp == 2:
                    gain = 0.8
            if hss == 0.1:
                if preamp == 1:
                    gain = 3.35
                if preamp == 2:
                    gain = 0.8
        self.gain = gain        
        

    def get_max_count(self):         
        #lê as dimensões da imagem
        img_shape = self.img_data.shape
        #cria uma máscara de 1s
        working_mask = np.ones(img_shape,bool)
        #cria dois arrays de indices
        ym, xm = np.indices(img_shape, dtype='float32')
        #cria um array com um círculo com centro nas coordenadas (x,y) fornecidas
        r = np.sqrt((xm - self.x)**2 + (ym - self.y)**2)

        #cria uma máscara configurando como valor 1 os pixeis dentro dos dois cículos        
        mask = (r > 2* self.sky_radius) * (r < 3 * self.sky_radius) * working_mask       
        #calcula a mediana dos pixeis de fundo
        self.sky_flux = np.median(self.img_data[np.where(mask)])
        #subtrai uma imagem de bias da imagem de ciencia 
        self.img_data = self.img_data - self.sky_flux        
        
        #cria uma máscara com valor 1 para os pixeis dentro do raio do céu
        mask = (r < self.sky_radius) * working_mask
        #pega todos os píxeis dentro do raio estipulado
        star_flux = self.img_data[np.where(mask)]
        #retorna o valor máximo dos pixeis dentro do raio
        self.max_star_flux = max(star_flux)
        #print(self.max_star_flux, self.sky_flux),exit()       

    

    def set_centroid(self):        
        #retorna as coordenadas do pixel de valor máximo
        y_max, x_max = np.where(self.img_data==self.max_star_flux)
        #seta das coordenadas encontradas como o novo centróide da estrela
        self.x, self.y  = x_max[0], y_max[0]
        #print(self.x, self.y),exit()       
        
         

    def calc_FWHM(self):                        
        radius = int(self.sky_radius/2)
        x,y = self.x, self.y
        img_data = self.img_data[y-radius:y+radius,x-radius:x+radius]
        #obtem os valores dos pixeis para a linha da coordenada x do centroide
        light_profile = np.take(img_data, radius, axis=0)
        #meia_altura = metade do valor máximo encontrado
        half_max = self.max_star_flux/2
        #cria um vetor com o número de pixeis da linha do centróide da estrela
        x = np.linspace(0, len(light_profile), len(light_profile))
        #Faz uma interpolação do perfil obtido
        spline = UnivariateSpline(x, light_profile-half_max, s=0)
        #Calcula onde da interpolação toda o eixo X
        r1, r2 = spline.roots()
        #A FWHM será a distância entre as duas raízes encontradas
        self.fwhm = r2-r1
        # O raio da estrela será 3x o valor de fwhm encontrado
        self.star_radius = 3*self.fwhm
        #print(self.star_radius)
        

    def calc_star_sky_flux(self):
        self.img_data = fits.getdata(self.img_name).astype(float)
        header = fits.getheader(self.img_name)
        img_shape = self.img_data.shape
        if img_shape[0] == 1: self.img_data = self.img_data[0]

        img_shape = self.img_data.shape
        #cria uma máscara de 1s
        working_mask = np.ones(img_shape,bool)
        #cria dois arrays de indices
        ym, xm = np.indices(img_shape, dtype='float32')
        #cria um array de pixeis dentro de dois círculos concêntricos com raios 2r e 3r
        r = np.sqrt((xm - self.x)**2 + (ym - self.y)**2)
        #cria uma máscara configurando como valor 1 os pixeis dentro dos dois cículos        
        mask = (r > 2* self.star_radius) * (r < 3 * self.star_radius) * working_mask
        #calcula a mediana dos pixeis de fundo
        self.sky_flux = np.median(self.img_data[np.where(mask)])

        #cria uma máscara para os pixeis dentro do raio da estrela           
        mask = (r < self.star_radius) * working_mask
        #fluxo da estrela subtraído do fluxo de fundo
        self.star_flux = (self.img_data[np.where(mask)] - self.sky_flux).sum()
        #número de pixeis da estrela
        self.n_pixels_object = len(self.img_data[np.where(mask)])


    def calc_SNR(self):
        #Esta funcao calcula o SNR da estrela
        gain = self.gain
        star = self.star_flux
        self.sky_flux -= np.median(self.bias_data)
        sky = self.sky_flux
        n_pix = self.n_pixels_object
        dc = self.dark_noise
        t_exp = self.exp_time
        rn = self.read_noise        

        em_gain = 1 #testes da SNR com o emgain        
        self.star_flux*=em_gain
        self.sky_flux*=em_gain
        self.dark_noise*=em_gain        
        aux = np.sqrt(star*gain + n_pix * (sky*gain + dc*t_exp + rn**2))
        self.SNR = star*gain/aux            


    def get_results(self):
        fwhm = self.fwhm
        snr = self.SNR
        n_pixels = self.n_pixels_object
        t_exp = self.exp_time
        dc = self.dark_noise*t_exp
        rn = self.read_noise           
        std_bias = 276.608517296231 #ADU        
        sky_flux = self.gain * (self.sky_flux) / t_exp - dc        
        star_flux = self.gain * self.star_flux / t_exp        
        star_radius = self.star_radius
        x,y = self.x+1, self.y+1 #precisa somar 1 por causa da diferença pro DS9
        return snr, fwhm, dc, rn, sky_flux, star_flux, n_pixels, star_radius, x, y


file_dir = r'C:\Users\denis\Desktop\UNIFEI\Projeto_Mestrado\Dados_HATS24b'
file_name = 'hats-24_I_transito_001.fits'
bias_name = 'bias_final_002.fits'

PFC = PhotonFluxCalc(img_name = file_dir + '\\' +  file_name, bias_name = file_dir + '\\' +  bias_name, xy_star = (507,369), sky_radius=20, ccd_serial=9916)
PFC.read_bias_img()
PFC.read_img_get_info()
PFC.get_max_count()
PFC.set_centroid()
PFC.calc_FWHM() 
PFC.calc_star_sky_flux()
PFC.calc_SNR()
snr, fwhm, dc, rn, sky_flux, star_flux, n_pixels, star_radius, x, y = PFC.get_results()

print('\n')
print('SNR: ', snr)
print('Sky Flux: ', sky_flux, 'photons/s')
print('Star Flux:', star_flux, 'photons/s')
print('FWHM: ', fwhm)
print('Star Pixels:', n_pixels)
print('Star Radius:', star_radius)
print('Centroide: %i,%i'%(x,y))
