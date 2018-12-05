function [y,x,z]=Shrink_Volume(Images)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This script shrinks a multi-dimensional image to minimize outer regions
% of zero signal. 
%
% Christopher W. Roy 2018-12-04
% fetal.xcmr@gmail.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Images=sum(sum(Images,4),5);
if size(Images,3)>1
st = regionprops(Images~=0, 'BoundingBox' );
Coord=[st(1).BoundingBox(1),st(1).BoundingBox(2),st(1).BoundingBox(3),st(1).BoundingBox(4),st(1).BoundingBox(5),st(1).BoundingBox(6)];
Coord(1:3)=floor(Coord(1:3));
Coord(4:6)=ceil(Coord(4:6));
y=Coord(2):Coord(2)+Coord(5);
x=Coord(1):Coord(1)+Coord(4);
z=Coord(3):Coord(3)+Coord(6);
x(x==0)=1;
y(y==0)=1;
z(z==0)=1;
else

st = regionprops(Images~=0, 'BoundingBox' );
Coord=[st(1).BoundingBox(1),st(1).BoundingBox(2),st(1).BoundingBox(3),st(1).BoundingBox(4)];
Coord(1:2)=floor(Coord(1:2));
Coord(3:4)=ceil(Coord(3:4));
y=Coord(2):Coord(2)+Coord(4);
x=Coord(1):Coord(1)+Coord(3);
x(x==0)=1;
y(y==0)=1;
z=1;
end
x=unique(x);
y=unique(y);
z=unique(z);
if rem(length(x),2),x=x(1:end-1);end
if rem(length(y),2),y=y(1:end-1);end
end