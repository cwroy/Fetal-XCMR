function Interp_KSpace=Interp_Phantom(Phantom,InterpolatedPhase)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates a phantom volume at simulated time points by
% linearly interpolating the XCAT volumes.
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Measured_Phases=linspace(0,1,size(Phantom,4)-1);
Measured_Phases=[-Measured_Phases(2),Measured_Phases];

i=find(Measured_Phases==InterpolatedPhase);
if ~isempty(i)
    Interp_KSpace=(Phantom(:,:,:,i));
else
    i=find(Measured_Phases<InterpolatedPhase,1,'last');
    LD=InterpolatedPhase-Measured_Phases(i);
    RD=Measured_Phases(i+1)-InterpolatedPhase;
    TD=LD+RD;
    Interp_KSpace=(Phantom(:,:,:,i)).*(1-LD./TD)+(Phantom(:,:,:,i+1)).*(1-RD./TD);
end
end
