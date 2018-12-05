function NOISE = Iterative_Noise(KSpace,ACQ)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iteratively adjust noise to match desired SNR level
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(ACQ.Trajectory,'CART')
    IMG=abs(sqrt(sum(ifft2c(KSpace).^2,4)));
    S=XCAT_to_MR(13,ACQ.FlipAngle,ACQ.TR,ACQ.TE,'Maternal');
    i=find((abs(IMG-S)/S)<0.01);
    Noise=col(mean(abs(IMG(i))))/ACQ.SNR;
    NOISE=Noise.*randn(size(KSpace)).*exp(1i*2*pi*(rand(size(KSpace))));
    IMG=abs(sqrt(sum(ifft2c(KSpace+NOISE).^2,4)));
    MeasuredSNR=mean(IMG(i))/std(IMG(i));
    iter=0;
    while abs((ACQ.SNR-MeasuredSNR)/ACQ.SNR)>0.05
        iter=iter+1;
        if MeasuredSNR>ACQ.SNR
            Noise=Noise*1.01;
            NOISE=Noise.*randn(size(KSpace)).*exp(1i*2*pi*(rand(size(KSpace))));
            IMG=abs(sqrt(sum(ifft2c(KSpace+NOISE).^2,4)));
            MeasuredSNR=mean(IMG(i))/std(IMG(i));
        end
        if MeasuredSNR<ACQ.SNR
            Noise=Noise*0.99;
            NOISE=Noise.*randn(size(KSpace)).*exp(1i*2*pi*(rand(size(KSpace))));
            IMG=abs(sqrt(sum(ifft2c(KSpace+NOISE).^2,4)));
            MeasuredSNR=mean(IMG(i))/std(IMG(i));
        end
        if iter>1000
            break;
        end
    end
    
elseif strcmpi(ACQ.Trajectory,'RAD')
    w = Radial_Weighting(size(KSpace,1));
    FT = MCNUFFT(ACQ.kx+1i*ACQ.ky,(w).^2,Resize_Volume(ACQ.Coils,ACQ.ImSize));
    NoiseROIy=0.5*ACQ.ImSize(1)-5:0.5*ACQ.ImSize(1)+5;
    NoiseROIx=0.5*ACQ.ImSize(1)-5:0.5*ACQ.ImSize(1)+5;
    IMG=FT'*KSpace;
    S=mean(col(abs(IMG(NoiseROIy,NoiseROIx))));
    Noise=double(XCAT_to_MR(13,ACQ.FlipAngle,ACQ.TR,ACQ.TE,'Maternal')/ACQ.SNR);
    NOISE=Noise.*randn(size(KSpace)).*exp(1i*2*pi*(rand(size(KSpace))));
    IMG = FT'*(KSpace+NOISE);
    MeasuredSNR=S/std(col(abs(IMG(NoiseROIy,NoiseROIx))));
    iter=0;
    while abs((ACQ.SNR-MeasuredSNR)/ACQ.SNR)>0.05
        iter=iter+1;
        if MeasuredSNR>ACQ.SNR
            Noise=Noise*1.01;
            NOISE=Noise.*randn(size(KSpace)).*exp(1i*2*pi*(rand(size(KSpace))));
            IMG = FT'*(KSpace+NOISE);
            MeasuredSNR=S/std(col(abs(IMG(NoiseROIy,NoiseROIx))));
        end
        if MeasuredSNR<ACQ.SNR
            Noise=Noise*0.99;
            NOISE=Noise.*randn(size(KSpace)).*exp(1i*2*pi*(rand(size(KSpace))));
            IMG = FT'*(KSpace+NOISE);
            MeasuredSNR=S/std(col(abs(IMG(NoiseROIy,NoiseROIx))));
        end
        if iter>1000
            break;
        end
    end
end

end