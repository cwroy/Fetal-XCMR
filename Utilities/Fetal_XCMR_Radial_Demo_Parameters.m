function [ACQ,PHYSIO] = Fetal_XCMR_Radial_Demo_Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script should be edited according to the desired values within
% the range of user selected parameters. Alternatively, these values can be
% changed in the main XCMR script directly.
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ACQ.GroundTruthFlag = 1; %
ACQ.Trajectory           = 'Rad'; % ['CART' or 'RAD']
ACQ.SliceOrientation     = 'SAx';  % ['TRA' or 'SAG' or 'COR' or 'SAx']
ACQ.FlipAngle            = 70;     % [Use reasonable scanner values i.e. 70] degrees
ACQ.TR                   = 4.06;   % [Use reasonable scanner values i.e. 3-5] ms
ACQ.TE                   = 2;      % [Use reasonable scanner values i.e. 1-3] ms

ACQ.SliceThickness       = 4;      % [Use reasonable scanner values i.e. 3-8] mm
ACQ.nSlices              = 1;      % [Choose according to slice thickness and heart dimensions which are aprox. 50mm x 57mm x 46mm)]
ACQ.SpatialResolution    = 1.5;      % [Use reasonable scanner values i.e. 1-2] mm

ACQ.nCoils               = 4;      % [Use reasonable (even numbered) scanner values: 4-32]
ACQ.SNR                  = 20;     % [Use reasonable scanner values]

% For radial this determines the number of spokes
ACQ.nMeasurements        = 3000;%3000;%3000;    % [Use reasonable scanner values 500-5000]
% or For Cartesian this determines the number of frames which is then multiplied by the number of lines to get the equivalent number of measurements. This is ignored for radial
ACQ.nFrames              = 20;     % [Use reasonable scanner values 1-30] 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User Selected Physiological Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Displacement of the fetus in 3 directions based on maternal respiration
% Can choose any value but the code will restrict fetal movement within the region of the maternal abdomen
PHYSIO.RespiratoryMotionAmplitude = [1,1,1];
% Displacement of the fetus in 3 directions based on random gross fetal movement
% Can choose any value but the code will restrict fetal movement within the region of the maternal abdomen
PHYSIO.FetalMotionAmplitude = [0,0,0];
% Fetal cardiac Motion (Yes = 1, No = 0)
PHYSIO.CardiacMotionFlag = 1;
% Maternal respiratory Motion (Yes = 1, No = 0 (i.e. breath-hold) )
PHYSIO.RespiratoryMotionFlag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% End of User selected Parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
