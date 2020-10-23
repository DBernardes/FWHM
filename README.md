# Introduction
This repo presents the library to calculate the FWHM of a star image. For that, it is needed to provide to the software the image name, the star (x,y) coordinates in pixels, and the maximum sky radius. The sofware will seek for the pixel within the sky radius with the maximum count value. The star centroid will be reset for the (x,y) coordinates of the maximum pixel found. With this information, the star FWHM will be calculated through the light profile. The star radius will be the length of the FWHM found. This README presents a description of the developed software, and a simple example of its working. Figure below presents an image with the star radius obtained through the software.

<p align="center">
  <img src="https://github.com/DBernardes/FWHM/blob/main/star_image.png" />
</p>


## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites
There are some packages that need to be installed before running the software. The first one is the Software Development Kit (SDK) developed by Andor Technology to control the CCDs. The second one is the GFITSIO package, used to save the data acquired by the camera in FITS format. 

![Software Development Kit (SDK)](https://andor.oxinst.com/products/software-development-kit/)

![GFITSIO](https://github.com/USNavalResearchLaboratory/GFITSIO)


### Installing
Clone this repo using ```https://github.com/DBernardes/FWHM.git```

## Running the tests
1. Before running the software, you need EMCCD to be connected to your PC.
2. Open the project SPARC4_AC.lvproj.
3. Run the VI SPARC4_GUI.vi.
4. Wait until the camera starts.
5. Set the night directory where the acquired images should be saved.
6. Press the Acquire button to start an acquisition. This would allow you to obtain a FITS files in your directory with the data acquired by the camera.

## Authors and Contact

* **Denis Bernardes**: 

email: denis.bernardes099@gmail.com 

## License

This project is licensed under the MIT License - see the [LICENSE.md](https://github.com/DBernardes/FWHM/blob/main/LICENSE) file for details
