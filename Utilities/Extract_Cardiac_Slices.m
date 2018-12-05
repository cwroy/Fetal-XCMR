function Slices = Extract_Cardiac_Slices(I)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to identify the region containing only the heart to speed up processing
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,zmin]=find((squeeze(sum(I,2)))>0,1,'first');
[~,zmax]=find((squeeze(sum(I,2)))>0,1,'last');
Slices=zmin:zmax;
end