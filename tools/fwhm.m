%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A sinple function for determining the full width at half maximum
% Input: 
% x: the array of x-values
% f: the array of function values
% Output:
% fwhm: FWHM of f in the scale given by x
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% old code
% % function FWHM = fwhm(x,f)
% % if nargin < 2
% %     f = x;
% %     x=[1:length(f)];
% % end
% % f1 = f-0.5*max(f);
% % ind = find(f1(1:end-1).*f1(2:end) <=0);
% % FWHM = x(ind(2))-x(ind(1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[Fullw]=fwhm(data)
% This function determines  full width at half 
% maximum of a peak if the input data  has two columns:
% Column 1 = x
% Column 2 = y
%Coded by Ebo Ewusi-Annan
%University of Florida 
%August 2012
x = data(:,1);
y= data(:,2);
maxy = max(y); 
f = find(y==maxy); 
cp = x(f);% ignore Matlabs suggestion to fix!!!
y1= y./maxy;
ydatawr(:,1) = y1;
ydatawr(:,2) = x;
newFit1=find(x>= cp);
newFit2=find(x < cp);
ydatawr2 = ydatawr(min(newFit1):max(newFit1),:);
ydatawr3 = ydatawr(min(newFit2):max(newFit2),:);
sp1 = spline(ydatawr2(:,1),ydatawr2(:,2),0.5);
sp2 = spline(ydatawr3(:,1),ydatawr3(:,2),0.5);
Fullw = sp1-sp2;
 end