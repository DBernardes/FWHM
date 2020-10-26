# Introduction
This repo presents the library to calculate the Full Width at Half Maximum (FWHM) of a star image. The software seeks for the pixel within the sky radius with the maximum count value. The star centroid is reset for the (x,y) coordinates of the maximum pixel found. With this information, the star FWHM is calculated through interpolation of the light profile, subtracted by half of the maximum pixel value. The star radius is given by three times the length of the FWHM found. Also, the software provides the option to calculate the signal-to-noise ratio (SNR) of the star image. This README presents a description of the developed software, and a simple example to run it. The figure below presents an image with the star radius obtained through the software.

<p align="center">
  <img src="https://github.com/DBernardes/FWHM/blob/main/star_image.png" />
</p>


# Description

This is the software for the calculation of the FWHM of a star image. Initially, it is needed to provide to the software the image name, the (x,y) coordinates of the star, and the maximum sky radius. Then, the software seeks for the pixel with the maximum count value within a circle centered in the provided coordinates and with a radius equal to the maximum star radius. The software provides the option to reset the centroid of the star for the maximum pixel. Based on this information, the software calculates the FHWM of the star. To accomplish that, it is done an interpolation for the light profile subtracted by half of the maximum pixel value along the x-direction of the star centroid. The FWHM is given by the distance between the point where the interpolation touches the x-axis. Also, the star radius equals three times the FWHM found.

Beyond the FWHM, the software allows us to calculate the SNR of the star. For this purpose, it is needed to provide the name of a bias image, the CCD gain (in e-/ADU), the exposure time, the dark current noise, and the EM gain (for EMCCDs). A bias image is an image acquired with the same CCD operation mode used to acquire the star image, with no light incidence, and for an exposure time equal to zero. The exposure time, dark current noise, and the EM gain are optional parameters. The software seeks the exposure time and the EM gain values in the image header. The dark current noise is estimated based on a model presented by Bernardes et al. ([link](https://arxiv.org/abs/1806.02191)) for the [SPARC4](https://www.spiedigitallibrary.org/proceedings/Download?fullDOI=10.1117/12.924976) EMCCD cameras. This model is based on the CCD temperature. For this reason, the software has the option to provide the CCD temperature used to acquire the star image. 

Then, the software uses the information obtained so far so calculate the sky and star flux in photons. The sky flux is given by the median of those pixels within the region between two circles with 2 times and 3 times the star radius found in the previous step. The star flux is given by the sum of those pixels within a circle with the star radius, subtracted by the sky flux. Therefore, the calculation of the SNR is given by

<p align="center">
  <img src="https://latex.codecogs.com/svg.latex?SNR&space;=&space;\frac{S}{\sqrt{S&space;&plus;&space;n_p&space;\times&space;[S_{sky}&space;&plus;&space;S_{dc}&space;&plus;&space;(\sigma&space;\times&space;G&space;/&space;G_{em})^2]}}" title="SNR = \frac{S}{\sqrt{S + n_p \times [S_{sky} + S_{dc} + (\sigma \times G / G_{em})^2]}}" />
</p>

where S is the star flux in photons, n<sub>p</sub> is the pixels number of the star, S<sub>sky</sub> is the sky flux in photon, S<sub>dc</sub> is the dark current noise in electrons, &sigma; is the image read noise in analogical-to-digital unit (ADU), G is the CCD gain in e-/ADU, and G<sub>em</sub> is the CCD Electrion Multiplying gain. 

 
# Running the software

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

## Prerequisites
There are some packages that need to be installed before running the software.They are:

* [math](https://docs.python.org/3/library/math.html)
* [astropy](https://www.astropy.org/)
* [numpy](https://numpy.org/)
* [scipy](https://www.scipy.org/)

To install these packages it is suggested to use the pip command as follows
```
pip install <package_name>
```

## Installing
Clone this repo using ```https://github.com/DBernardes/FWHM.git```

## Running the tests
1. Open the file run.py. In this file, it is needed to provide to the ```file_dir``` variable the path of the sample images in your current machine. 
2. Execute the run.py file. The calculation results will be printed on the screen.   


## Authors and Contact

* **Denis Bernardes**: 

email: denis.bernardes099@gmail.com 

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/DBernardes/FWHM/blob/main/LICENSE) file for details
