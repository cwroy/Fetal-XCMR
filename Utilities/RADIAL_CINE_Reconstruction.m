function IMG = RADIAL_CINE_Reconstruction(KSpace,ACQ,Times,PHYSIO, nFrames)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates CINE images of nFrames by sorting and gridding Radial KSpace
%
% Reconstruction of radial data requires the iGRASP code
% from Li Feng and Ricardo Otazo.
% Feng L, Grimm R, Tobias Block K, Chandarana H, Kim S, Xu J, Axel L, Sodickson DK, Otazo R. Golden-angle radial sparse parallel MRI: Combination of compressed sensing, parallel imaging, and golden-angle radial sampling for fast and flexible dynamic volumetric MRI. Magn Reson Med. 2013 Oct 18. doi: 10.1002/mrm.24980.
% It can be downloaded here:
% https://cai2r.net/resources/software/grasp-matlab-code
%
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(1)
    if exist('MCNUFFT.m','file')==0
        inputdlg('Simulation and reconstruction of radial data requires the iGRASP code. Add to path or downloaded from: https://cai2r.net/resources/software/grasp-matlab-code','Error');
        IMG=[];
        break
    end
    
    try
        % Calculate cardiac phases
        CP=Calculate_Phases(Times,cumsum([0,PHYSIO.RRInterval]));
        rCP=linspace(0,1,nFrames+1);
        w = Radial_Weighting(size(KSpace,1));
        IMG=zeros([ACQ.ImSize,nFrames]);
        for iFrame=1:nFrames
            % Bin the data according to cardiac phase
            i=CP>=rCP(iFrame)&CP<rCP(iFrame+1);
            D = Sampling_Density(ACQ.kx(:,i)+1i*ACQ.ky(:,i));
            % Create non-uniform Fourier transform operator
            FT = MCNUFFT(ACQ.kx(:,i)+1i*ACQ.ky(:,i),(bsxfun(@times,D.^2,repmat(w,[1,sum(i)]).^2)),Resize_Volume(ACQ.Coils,ACQ.ImSize));
            % Perform gridding reconstruction
            IMG(:,:,iFrame)=(FT'*double(KSpace(:,i,:)))/mean(D(:));
        end
    catch
        disp('Error...')
        IMG=[];
    end
    break
end
end
function D = Sampling_Density(k)
% Calculate the sampling density according to the distance between sorted
% radial spokes
Angles=(squeeze(atand(real(k(1,:))./imag(k(1,:)))));
[sAngles,i]=sort(Angles);
D=zeros(1,length(Angles));
D(1,i(1))=(sAngles(2)-(sAngles(end)-180))/2;
D(1,i(end))=((180+sAngles(1))-sAngles(end-1))/2;
for iAngle=2:length(sAngles)-1
    D(1,i(iAngle))=(sAngles(iAngle+1)-sAngles(iAngle-1))/2;
end

end