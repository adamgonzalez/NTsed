# NTsed
A spectral energy distribution (SED) model for active galactic nuclei (AGNs) that includes emission from the accretion disc (UV/optical) and an energetically linked hot corona (X-ray).

The disc emission is computed following Novikov & Thorne (1973), specifically the formulation provided in Page & Thorne (1974). The hot corona emission is computed using the thermal Comptonization model of Zdziarski, Johnson & Magdziarz (1996) and Zycki, Done & Smith (1999) called nthComp, the luminosity of which is linked to the accretion flow following Kubota & Done (2018).

## Installation
This is an XSPEC model, and therefore requires a working installation of HEASOFT (tested successfully on versions >=6.30). After either downloading the `srccode.cxx` and `model.dat` files to a model directory or cloning the repository, the model can be installed by opening XSPEC and executing `initpackage ntsed model.dat .`. A successful installation will not produce any errors when loading the model using `lmod ntsed .`. To avoid loading the model manually each time you want to use it, you can add `load /path/to/model/libntsed.so` to your `~/.xspec/xspec.rc` file.

## Model parameters
The model has the following physically relevant variable parameters:
1. `logMBH`
