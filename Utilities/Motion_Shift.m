function NewVolume = Motion_Shift(PhantomVolume,Displacement)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Shift PhantomVolume by desired displacement. Can be fractional
%
% Based off of fraccircshift
% https://www.mathworks.com/matlabcentral/fileexchange/45950-fraccircshift
% 
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int = floor(Displacement);     %integer portions of Displacement
fra = Displacement - int;      %fractional portions of Displacement
dim = numel(Displacement);
NewVolume = PhantomVolume;
for n = 1:numel(Displacement)  %The dimensions are treated one after another.
    intn = int(n);
    fran = fra(n);
    Displacement1 = zeros(dim,1);
    Displacement1(n) = intn;
    Displacement2 = zeros(dim,1);
    Displacement2(n) = intn+1;
    %Linear intepolation:
    NewVolume = (1-fran)*circshift(NewVolume,Displacement1) + fran*circshift(NewVolume,Displacement2);
end