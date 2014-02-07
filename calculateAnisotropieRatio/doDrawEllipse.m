function [px py] = doDrawEllipse(ccc, orientationStore ,k,BW2)    
%%%%%% This function draws an ellipse which can be rotated. It expects
%%%%%% output from regionprops. k selects the desired region.

phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);
    
xbar = ccc(k).Centroid(1);
ybar = ccc(k).Centroid(2);

    a = ccc(k).MajorAxisLength/2;
    b = ccc(k).MinorAxisLength/2;

    %theta = pi*((ccc(k).Orientation-orient)/180);
    
    theta = pi*(orientationStore/180);
    
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)];

    xy = [a*cosphi; b*sinphi];
    xy = R*xy;

    px = xy(1,:) + xbar;
    py = xy(2,:) + ybar;
    
    imshow(BW2)
    hold on
   plot(px,py,'ro')
   plot(ccc(k).Centroid(1),ccc(k).Centroid(2),'xr')