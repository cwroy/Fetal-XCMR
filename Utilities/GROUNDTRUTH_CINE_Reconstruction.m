function IMG = GROUNDTRUTH_CINE_Reconstruction(GROUNDTRUTH,ACQ,PHYSIO, nFrames)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates a CINE image with nFrames from the GROUNDTRUTH images and
% simulated cardiac intervals
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try    
    CP=Calculate_Phases(ACQ.TIME,cumsum([0,PHYSIO.RRInterval]));
    rCP=linspace(0,1,nFrames+1);
    IMG=zeros([ACQ.ImSize,nFrames]);
    for iFrame=1:nFrames
        i=CP>=rCP(iFrame)&CP<rCP(iFrame+1);
        IMG(:,:,iFrame)=mean(GROUNDTRUTH(:,:,i),3);
    end
catch
    disp('Error...')
    IMG=[];
end
end
