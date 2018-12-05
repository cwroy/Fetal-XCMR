function Motion=Interp_Motion(MeasuredPhases,InterpolatedPhases)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates motion at simulated time points by
% linearly interpolating the measured positions from the XCAT volumes.
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Measured_Phases=linspace(0,1,size(MeasuredPhases,1)+1);
MeasuredPhases=cat(1,MeasuredPhases(end,:),MeasuredPhases,MeasuredPhases(1,:));
Measured_Phases=[-Measured_Phases(2),Measured_Phases];

Motion=zeros(length(InterpolatedPhases),3);
for iPhase=1:length(InterpolatedPhases)
i=find(Measured_Phases==InterpolatedPhases(iPhase));
if ~isempty(i)
    Motion(iPhase,:)=(MeasuredPhases(i,:));
else
    i=find(Measured_Phases<InterpolatedPhases(iPhase),1,'last');
    LeftWeight=InterpolatedPhases(iPhase)-Measured_Phases(i);
    RightWeight=Measured_Phases(i+1)-InterpolatedPhases(iPhase);
    TotalWeight=LeftWeight+RightWeight;
    Motion(iPhase,:)=(MeasuredPhases(i,:)).*(1-LeftWeight./TotalWeight)+(MeasuredPhases(i+1,:)).*(1-RightWeight./TotalWeight);
end
end
end
