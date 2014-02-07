function [px py] = doSelectEllipsePoint(ccc, orient ,k)    
%%%%%% smaples the number of points within the ellipse currently set to 1
%%%%%% but can also be increased standart values are out commented
%phi = linspace(0,2*pi,50);
%%%%%% Orient set to 0 or 180 will mark the long axis -+ 90 the short axis
%%%%%% Input parameters is the centroid and orientation ouptut from
%%%%%% regionprops. The parameter k selects the region identified by
%%%%%% regionprops.
%orient =+70

phi = linspace(0,2*pi,1);
cosphi = cos(phi);
sinphi = sin(phi);
    
xbar = ccc(k).Centroid(1);
ybar = ccc(k).Centroid(2);

    a = ccc(k).MajorAxisLength/2;
    b = ccc(k).MinorAxisLength/2;

    theta = pi*((ccc(k).Orientation-orient)/180);
    
    %theta = pi*(ccc(k).Orientation-orient/180);
    
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)];

    xy = [a*cosphi; b*sinphi];
    xy = R*xy;

    px = xy(1,:) + xbar;
    py = xy(2,:) + ybar;
    
%    imshow(BW2)
 %   hold on
 %  plot(px,py,'ro')
 %  plot(ccc(k).Centroid(1),ccc(k).Centroid(2),'xr')