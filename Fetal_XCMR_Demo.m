clear;clc;close all
DemoStartTime=tic;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This demo generates fetal XCMR data using either Cartesian or radial
% sampling. It takes about 20-30 minutes to run both demos.
%
% Simulation and reconstruction of radial data requires the iGRASP code from Li Feng and Ricardo Otazo.
% Feng L, Grimm R, Tobias Block K, Chandarana H, Kim S, Xu J, Axel L, Sodickson DK, Otazo R. Golden-angle radial sparse parallel MRI: Combination of compressed sensing, parallel imaging, and golden-angle radial sampling for fast and flexible dynamic volumetric MRI. Magn Reson Med. 2013 Oct 18. doi: 10.1002/mrm.24980.
% It can be downloaded here:
% https://cai2r.net/resources/software/grasp-matlab-code
%
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo 1 Cartesian Data Simualtion and CINE reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
addpath('Utilities')
[ACQ,PHYSIO] = Fetal_XCMR_Cartesian_Demo_Parameters;
[KSPACE,TIME,NOISE,ACQ,PHYSIO,CART_GROUNDTRUTH] = Fetal_XCMR_Main(ACQ,PHYSIO);

% Reconstruct ground truth and Cartesian CINEs
CART_CINE = CART_CINE_Reconstruction(KSPACE+NOISE, TIME, cumsum([0,PHYSIO.RRInterval]), 15);
CART_CINE0 = GROUNDTRUTH_CINE_Reconstruction(CART_GROUNDTRUTH,ACQ,PHYSIO, 15);


% Display Results
if ~isempty(CART_CINE)
    figure(1)
    imCropY=101:228;    
    imCropX=108:235;
    subplot(2,1,1)
    imagesc(abs([CART_CINE0(imCropY,imCropX,1),CART_CINE0(imCropY,imCropX,8),imresize(squeeze(CART_CINE0(imCropY,176,:)),[length(imCropY),30])]));axis image;colormap gray
    title('Ground Truth')
    subplot(2,1,2)
    imagesc(abs([CART_CINE(imCropY,imCropX,1),CART_CINE(imCropY,imCropX,8),imresize(squeeze(CART_CINE(imCropY,176,:)),[length(imCropY),30])]));axis image;colormap gray
    title('Cartesian Reconstruction')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Demo 2 Cartesian Data Simualtion and CINE reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[ACQ,PHYSIO] = Fetal_XCMR_Radial_Demo_Parameters;
[KSPACE,TIME,NOISE,ACQ,PHYSIO,RAD_GROUNDTRUTH] = Fetal_XCMR_Main(ACQ,PHYSIO);

% Reconstruct ground truth and radial CINEs
% Note that only a gridding reconstruction is provided in this package but
% the data can readily be used with available compressed sensing
% reconstructions to supress streaking artifact i.e. the iGRASP code mentioned above.
RAD_CINE = RADIAL_CINE_Reconstruction(KSPACE,ACQ,TIME,PHYSIO, 15);
RAD_CINE0 = GROUNDTRUTH_CINE_Reconstruction(RAD_GROUNDTRUTH,ACQ,PHYSIO, 15);

% Display Results
if ~isempty(RAD_CINE)
    figure(2)
    imCropY=101:228;    
    imCropX=108:235;
    subplot(2,1,1)
    imagesc(abs([RAD_CINE0(imCropY,imCropX,1),RAD_CINE0(imCropY,imCropX,8),imresize(squeeze(RAD_CINE0(imCropY,176,:)),[length(imCropY),30])]));axis image;colormap gray
    title('Ground Truth')
    subplot(2,1,2)
    imagesc(abs([RAD_CINE(imCropY,imCropX,1),RAD_CINE(imCropY,imCropX,8),imresize(squeeze(RAD_CINE(imCropY,176,:)),[length(imCropY),30])]));axis image;colormap gray
    title('Radial Gridding Reconstruction')
end
display(['Demo completed in: ',num2str(round(toc(DemoStartTime))),'s'])