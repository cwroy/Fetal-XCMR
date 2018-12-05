Fetal XCMR v1.0 2018-12-04
============================================================
This package contains scripts designed to simulate fetal cardiac MRI data.

Example usage of the main script "Fetal_XCMR_Main.m" is demonstrated by two demos for Cartesian and Radial data acquisitions respectively contained in "Fetal_XCMR_Demo.m".

Dependencies
============================================================
This code was tested in Matlab R2018a

All code is self contained and subfunctions are located in the Utilities folder with the exception of the non-uniform fourier transform which is required for radial data generation and reconstruction. It can easily be implemented using the iGRASP code from Li Feng and Ricardo Otazo.

Feng L, Grimm R, Tobias Block K, Chandarana H, Kim S, Xu J, Axel L, Sodickson DK, Otazo R. Golden-angle radial sparse parallel MRI: Combination of compressed sensing, parallel imaging, and golden-angle radial sampling for fast and flexible dynamic volumetric MRI. Magn Reson Med. 2013 Oct 18. doi: 10.1002/mrm.24980.

Which can be downloaded here:
https://cai2r.net/resources/software/grasp-matlab-code

Generating new source images using XCAT
============================================================
Source images are included representing dynamic fetal cardiac and maternal respiratory (50 frames each).

XCAT is a separate software, which allows you to create additional volumes which could be used by Fetal XCMR
The XCAT software is available here: http://www.hopkinsradiology.org/DMIP/Research/xcat

Distributing Fetal XCMR
============================================================
Fetal XCMR source code is distributed as Matlab files to support its use for research purposes. Anyone is allowed to use the original and/or modified Fetal XCMR code for non-commercial purposes. If you publish research using Fetal XCMR, please cite the Fetal XCMR publication currently submitted to JCMR.

Contact
============================================================
For any questions, comments and contributions, please contact
Chris Roy 
fetal.xcmr@gmail.com
https://github.com/cwroy/Fetal-XCMR/
(c) Christopher W. Roy 2018


