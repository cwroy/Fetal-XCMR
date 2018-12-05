function CP=Calculate_Phases(Times,PhysioSignal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script calculates the cardiac or respiratory phase at a given time point given a
% simulated heart rate or respiratory rate(PhysioSignal)
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make sure the PhysioSignal spans the length of desired time points
while PhysioSignal(end)<max(Times(:))
    PhysioSignal(end+1)=PhysioSignal(end)+median(diff(PhysioSignal));%ScanLength+10;
end


% Times=single(Times);
if size(PhysioSignal,1)<size(PhysioSignal,2)
PhysioSignal=PhysioSignal';
end

if 10*size(Times,1)<size(Times,2)
Times=Times';
end


Before_Index=zeros(size(Times,1),size(Times,2),length(PhysioSignal)-1,'double');
After_Index=zeros(size(Times,1),size(Times,2),length(PhysioSignal)-1,'double');
% Determine which time points occur before and after each RWaveTime
% (indices)
for loop=1:(length(PhysioSignal)-1)
Before_Index(:,:,loop)=Times<=PhysioSignal(loop+1);
After_Index(:,:,loop)=Times>PhysioSignal(loop+1);
end
Before_Index=abs(sum(Before_Index,3)-length(PhysioSignal));
After_Index=sum(After_Index,3)+1;
%Keeps track of indices where data was not collected
Before_Index(isnan(Times))=1;
After_Index(isnan(Times))=1;
% For each time point, determine which RWaves come before and after
Last_RWaves=PhysioSignal(After_Index);
Next_RWaves=PhysioSignal(Before_Index+1);
% Calculate Cardiac Phase for each time point
CP=(Times-Last_RWaves)./(Next_RWaves-Last_RWaves);
%Keeps track of indices where data was not collected
CP(isnan(Times))=nan;
CP(CP==1)=0;
end