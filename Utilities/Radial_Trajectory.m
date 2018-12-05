function [kx,ky] = Radial_Trajectory(Echolength, nSpokes, gatype)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script is used to generate a radial trajectory for nSpokes of length
% Echolength. The trajectory is based on the golden where gatype:
%(1 - standard GA, 7 - tiny golden angle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

k0 = [zeros(1,Echolength); -floor(Echolength/2):ceil(Echolength/2-1)];
dPhi=pi/(gatype-1+(1+sqrt(5))/2);
kx=zeros(Echolength,nSpokes);
ky=zeros(Echolength,nSpokes);
for i=1:nSpokes
    rot_angle = (i-1)*dPhi;    
    R = [cos(rot_angle), -sin(rot_angle);
        sin(rot_angle),  cos(rot_angle)];
    ktmp       = (R*k0).';
    kx(:,i) = ktmp(:,1);
    ky(:,i) = ktmp(:,2);
end
kx = kx./Echolength;
ky = ky./Echolength;
end



