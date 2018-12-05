function Motion = Generate_Motion(N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script generates bounded random motion
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Motion=zeros(N,1);
Motion(1,1) = rand(1,1)-0.5;
for i=2:N
   Motion(i,1) = Motion(i-1,1) + (rand(1,1)-0.5);
end
Motion=Motion-mean(Motion);
Motion=Motion./max(abs(Motion(:)));

end