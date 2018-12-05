function CSM = Simulate_Coils(imSize,nCoils)
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates (nCoils) coil sensitivity maps (CSM) for an image of size (imSize)
%
%
% This code is based off of code from Santiago Aja-Fernandez
% www.lpi.tel.uva.es/~santi
 
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 

vX=repmat((1:imSize(1))',[1,imSize(2)]);
vY=repmat(1:imSize(2),[imSize(1),1]);
vX=pi/2.*vX./max(vX(:));
vY=pi/2.*vY./max(vY(:));
 
CSM=zeros([imSize,nCoils]);
Theta=0:(2*pi/nCoils):(2*pi-(2*pi/nCoils));
Theta=Theta(1:end-1);
for ii=1:nCoils/2
  if (Theta(ii)<=pi/2)
      N1=vX.*cos(Theta(ii))+vY.*sin(Theta(ii));
      N1=N1./max(N1(:)).*pi/2; 
      CSM(:,:,ii)=cos(N1);
      CSM(:,:,nCoils./2+ii)=sin(N1);
  else %>pi/2
      N1=(vX).*abs(cos(Theta(ii)))+(pi/2-vY).*sin(Theta(ii));
      N1=N1./max(N1(:)).*pi/2; 
      CSM(:,:,ii)=sin(N1);
      CSM(:,:,nCoils./2+ii)=cos(N1);       
  end
end
 
CSM=CSM./sqrt(nCoils/2);