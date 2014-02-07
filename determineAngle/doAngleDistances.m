
%%%%%%%%%% Determine the angel of the image by fitting an ellipse around
%%%%%%%%%% the image and measure the chromatin chromatin, as well as the
%%%%%%%%%% cortex diamter.


function [distChrom distContractileRing angle] = doAngleDistances(imgMid,voxelX_mumMid,tifFilename)

test  = imgMid(:,:,2);
figure(1)
imshow(test(:,:,:),[]);
%marker = testNorm;
testNorm = normalizedImage(test);
marker = false(size(testNorm));


    for lauf = 1:2

         fh = figure(1);
         title('Mark the Background')
        frapBkg = roipoly(testNorm(:,:,1));
        close(fh);
   
        marker(frapBkg) = 255;
                
    end

%testNorm = normalizedImage(redStack(:,:,15));
%%%%% seeded wathershed to segment the image. This step is done to
%%%%% determine the orientation of the cell

    h = fspecial('gaussian', 10, 10) 
    redStackGauss = imfilter(test,h);
    imshow(redStackGauss(:,:,1),[])
     testNorm = normalizedImage(redStackGauss);
     mp = imimposemin(testNorm,marker);
     L2 = watershed(mp);
                %mpReg = imregionalmin(mp);
     t3 = false(size(testNorm));
     t3(L2 == 0 ) = 255;
     imshow(t3,[])          
     L2 = t3;
     cc = regionprops(L2,'all');
     B = bwboundaries(L2,8); 
     boundary = B{1};
     
     %%%%%%%%%%%
     %%%% This part will calculate and plot the ellipse, for visual
     %%%% inspection
     %%%%
     
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

    h= figure(1)
        imshow(L2,[])
          hold on
         plot(x,y,'r','LineWidth',2);
    print(h,'-dpdf', [ tifFilename,'_Ellipse.pdf']);%tifCurvetifFilename);
    close all
 
    
    
    %%%%%% This section smoothes the images with a gaussian filter and
    %%%%%% applies a canny edge detection to get the outline of the
    %%%%%% chromatin. Then an ellipse is fitted around the outline and the
    %%%%%% center of mass determined.
    orientation = cc(k).Orientation;
    
    %mip = max(imgMid, [], 3);
    
   
    imshow(imgMid(:,:,1),[])
    hold on
    plot(center(1),center(2),'x')
    
    G = fspecial('gaussian',[5 5],2);
  
    imgFilter = imfilter(imgMid(:,:,1),G,'same');
  
    BW2 = edge(imgFilter,'canny',0.4,10);
    imshow(BW2)
    
   
    
    %%%%%%% Dilation and erosion is required to close gaps in the outline
    
   se = strel('disk',6);
   
  % BW3 = imdilate(BW2,se);
  
   % BW4 = imerode(BW3,se,6);
    BW4 = BW2;
    ccc = regionprops(BW4,'all');
   
   
    
    
    BW2(:,:) = 0;
   
    BW2(ccc(1).PixelIdxList) = 255;
    imshow(BW2)
   
    
    %%%%% Compute the weighted center of mass
    
     s = regionprops(BW4, imgMid(:,:,1), {'Centroid','WeightedCentroid'});
   
    
    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% This section determines the contractile ring diameter
    %%%%%%%%%
    
    
        imshow(imgMid(:,:,2),[])
   
        hold on
        [xEq,yEq] = ginput(2)
        close all
        
      
        dist1 = pdist2([xEq(1),yEq(1)],[xEq(2),yEq(2)]);
                   
       
        distContractileRing = dist1*voxelX_mumMid;

    %%%%%%%%%%%%%%%
    %%%%% This section tests whether a metaphase plate or divided
    %%%%% chromatin was identified.

    dividedChromatin = length(ccc);
    
    if dividedChromatin > 1
         
            chromatin1 = s(1).WeightedCentroid;
         
            chromatin2 = s(2).WeightedCentroid;
    
    
   
  %  imshow(BW4)
  %  hold on
  %  plot(s(1).Centroid(1),s(1).Centroid(2),'xr')
  %  plot(s(2).Centroid(1),s(2).Centroid(2),'xr')
  %  plot(chromatin1(1),chromatin1(2),'xb')
  %   plot(chromatin2(1),chromatin2(2),'xb')
  %  line([chromatin1(1) chromatin2(1)],[chromatin1(2) chromatin2(2)])
     
   % close all
    

   % Calculate the angle between the points and convert it from rad to the
   % degree system
   
    dx = chromatin2(1) - chromatin1(1);
    dy = chromatin2(2) - chromatin1(2);
    
    angle = (abs(atan2(dy,dx)))*180/pi
    
    
                   
    
    BW_test = BW2;
    BW_test(:,:) = 0;
   imshow(BW_test)
hold on

                   
         dist2 = pdist2([ chromatin1(1),chromatin1(2)],[ chromatin2(1), chromatin2(2)]);
                   
         distChrom = dist2*voxelX_mumMid;
  
         h=figure(1)  
         imshow(imgMid(:,:,1),[])
         hold on
    
             plot(chromatin1(1),chromatin1(2),'xb')
             plot(chromatin2(1),chromatin2(2),'xb')
                 for subrun = 1:2
       
                   shapeTmp = ccc(subrun).PixelList ;
                   shapeX = shapeTmp(:,1);
                   shapeY = shapeTmp(:,2);
    
                    plot(shapeX, shapeY,'xr') 
                 end

     
            print(h,'-dpdf', [ tifFilename,'_Chromatin.pdf']);%tifCurvetifFilename);
             close all
    
    else
        
        
        chromatin1 = s(1).WeightedCentroid;
                   
        dist2 =0;
      
                  angle = 0; 
      
         distChrom = dist2*voxelX_mumMid;
  
         h=figure(1)  
         imshow(imgMid(:,:,1),[])
         hold on
    
            plot(chromatin1(1),chromatin1(2),'xb')
   
          subrun = 1
          shapeTmp = ccc(subrun).PixelList ;
          shapeX = shapeTmp(:,1);
          shapeY = shapeTmp(:,2);
    
         plot(shapeX, shapeY,'xr') 
     
     
          print(h,'-dpdf', [ tifFilename,'_Chromatin.pdf']);%tifCurvetifFilename);
         close all
     end
    
end
