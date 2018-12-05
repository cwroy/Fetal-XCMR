function CINE = CART_CINE_Reconstruction(KSpace, Times, RWaveTimes, nFrames)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates CINE images of nFrames from Cartesian KSpace data and
% related physiological signals
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize K-Space
nRows = size(KSpace,2);
nCols = size(KSpace,1);
nCoils = size(KSpace,4);
rKSpace = zeros(nCols, nRows, nFrames, nCoils, 'single'); 
[LeftIndices,RightIndices,LeftWeights,RightWeights] = CardiacPhase_Interpolation(Times,RWaveTimes,nFrames);
LeftWeights(isnan(LeftWeights))=0;RightWeights(isnan(RightWeights))=0;
for iFrame = 1:nFrames
    for iRow = 1:nRows
        rKSpace(:,iRow,iFrame,:) = KSpace(:,iRow,LeftIndices(iFrame,iRow),:).*LeftWeights(iFrame,iRow)+KSpace(:,iRow,RightIndices(iFrame,iRow),:).*RightWeights(iFrame,iRow);
    end
end
CINE=sqrt(sum(ifft2c(rKSpace).^2,4));

end

function [LeftIndices,RightIndices,LeftWeights,RightWeights] = CardiacPhase_Interpolation(Times,RWaveTimes,nFrames)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calculates weights and indices for linear nearest neighboor interpolation
%
% nRows: Number of Rows in Original Data Array
% nCPs: Maximum number of cardiac phases in Original Data Array
% Times: Timing Data Matrix (nRows, nCardiacPhases)
% RWaveTimes: RWaveTimes Vector
% nFrames: Number of desired reconstructed cardiac phases
%
%
% LeftIndices & RightIndices: Indices for nearest neighboors data corresponding to each desired cardiac phase (nFrames, nRows)
% LeftWeight & RightWeights: Weightings for nearest neighboor interpolation corresponding to each desired cardiac phase (nFrames, nRows)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nRows=size(Times,1);nCP=size(Times,2);
 
% Lay out desired cardiac phases
ReconPhases = linspace(0,1,nFrames+1);ReconPhases = ReconPhases(1:end-1);
 
% Calculate original cardiac phases
CP=Calculate_Phases(Times,RWaveTimes)';
 
% The minimum and maximum original cardiac phases are used to interpolate endpoint CPs (wrap)
[maxCP,maxIndex]=max(CP);
[minCP,minIndex]=min(CP);
wrapIndices=[maxIndex;repmat((1:nCP)',[1,nRows]);minIndex];
 
% Sort cardiac phases to facilitate nearest neighboor calculations bellow
[sortCP,sortIndices]=sort([maxCP-1;CP;minCP+1],1);
% sub2ind is used to map indices
isort=sub2ind([nCP+2,nRows,1,1],sortIndices,repmat(1:nRows,[nCP+2,1]));
  
% Find nearest neighboors
nearestIndices=zeros(nFrames,size(CP,2));
for loop=1:nFrames
    tempCP=sortCP-ReconPhases(loop);
        tempCP(tempCP>0)=nan;
    [~,i]=min(abs(tempCP));
    nearestIndices(loop,:)=i;  
end
 
% sub2ind is used to map indices
imin=sub2ind([nCP+2,nRows,1,1],nearestIndices,repmat(1:nRows,[nFrames,1]));
 
% map nearrest neighboors indices to original indices
LeftIndices=wrapIndices(isort(imin));
RightIndices=wrapIndices(isort(imin+1));
 
% Calculate nearest neighboors distances
LeftDistances=abs(repmat(ReconPhases',[1 nRows])-sortCP(imin(:,:,1)));
RightDistances=abs(sortCP(imin(:,:,1)+1)-repmat(ReconPhases',[1 nRows]));
 
% Check to make sure Left and Right Distances are not both zero
 
i=find((LeftDistances+RightDistances)==0);
if ~isempty(i)
temp=abs(sortCP(imin(:,:,1)+2)-repmat(ReconPhases',[1 nRows]));
RightDistances(i)=temp(i);
temp=wrapIndices(isort(imin+2));
RightIndices(i)=temp(i);
end
 
% Calculate nearest neighboors interpolation weightings
LeftWeights=1-LeftDistances./(LeftDistances+RightDistances);
RightWeights=1-RightDistances./(LeftDistances+RightDistances);
end