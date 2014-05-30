function [flank1 flank2 flankMid1 flankMid2] = doSplitWatershedInTwo(t3Store, greenStack)

lauf = 1;

L2 = t3Store(:,:,lauf);
  cc = regionprops(L2,'all');
  B = bwboundaries(L2,8); 
  boundary = B{1};

  bStore{lauf} = B{1};
  
k=1
phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);
    
xbar = cc(k).Centroid(1);
    ybar = cc(k).Centroid(2);

    a = cc(k).MajorAxisLength/2;
    b = cc(k).MinorAxisLength/2;

    theta = pi*(cc(k).Orientation/180);
    
    R = [ cos(theta)   sin(theta)
         -sin(theta)   cos(theta)];

    xy = [a*cosphi; b*sinphi];
    xy = R*xy;

    x = xy(1,:) + xbar;
    y = xy(2,:) + ybar;

    imshow(L2,[])
    hold on
    plot(x,y,'r','LineWidth',2);


   % tan((cc.Orientation+90)/360*2*pi)*250
   
   %%%%% left
   
          
       
   
    orientationSet  =90;
   
   
   
    [xCoord,yCoord,subRegionIntersect,selectedRoi] = splitFurrowsCleavageFurrow(L2,x,y,boundary, orientationSet, 300,greenStack(:,:,lauf),cc)    
    
   
    cNewMid1(lauf) = xCoord;
    rNewMid1(lauf) = yCoord;
    
      flank1 = selectedRoi;
      flankMid1 = subRegionIntersect;
      
    
   %%%%% right
    orientationSet  =-90;
   
   
   
    [xCoord,yCoord,subRegionIntersect,selectedRoi] = splitFurrowsCleavageFurrow(L2,x,y,boundary, orientationSet, 300,greenStack(:,:,lauf),cc)    
    
   % [xCoord, yCoord, regionIntersect,selectedRoi] = detectAxis_v1(L2, cc,x,y,boundary,400)
    cNewMid2(lauf) = xCoord;
    rNewMid2(lauf) = yCoord;
    
    
    
    flank2 = selectedRoi;
    flankMid2 = subRegionIntersect;
    
    
      
end