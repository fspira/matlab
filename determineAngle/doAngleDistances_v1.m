
%%%%%%%%%% Determine the angel of the image by fitting an ellipse around
%%%%%%%%%% the image and measure the chromatin chromatin, as well as the
%%%%%%%%%% cortex diamter.


function [distChrom distContractileRing angle] = doAngleDistances(imgMid,voxelX_mumMid,tifFilename)

    %mip = max(imgMid, [], 3);
    
   
   
  [cChromatin1 rChromatin1 cChromatin2 rChromatin2] = doChromatinChromatinDistance(imgMid(:,:,1))
    
   
  %%%%% Test whether input cell is a metaphase cell
  if cChromatin1 == cChromatin2
  
    G1 = fspecial('gaussian',[75 75],15);
    imgFilter = imfilter(imgMid(:,:,1),G1,'same');
    BW2 = edge(imgFilter,'canny',0.7,10);
    BW2 = bwmorph(BW2,'bridge',8);
    imshow(BW2)
    cc = regionprops(BW2,'orientation');
    
    angle = cc(1).Orientation;
    
    
  else
          dx = cChromatin2(1) - cChromatin1(1);
          dy = rChromatin2(1) - rChromatin1(1);

          angle = (abs(atan2(dy,dx)))*180/pi

      
  end
    
    %%%%% Compute the weighted center of mass
    
   % s = regionprops(BW4, imgMid(:,:,1), {'Centroid','WeightedCentroid'});
   
    
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
    %%%%% This section computes the chromatin chromatin distance
    %%%%%%%%%%%%%%%

                   
         dist2 = pdist2([ cChromatin1(1),rChromatin1(1)],[ cChromatin2(1), rChromatin2(1)]);
                   
         distChrom = dist2*voxelX_mumMid;
  
         h=figure(1)  
         imshow(imgMid(:,:,1),[])
         hold on
    
             plot(cChromatin1(1),rChromatin1(1),'xb')
             plot(cChromatin2(1),rChromatin2(1),'xb')
             %    for subrun = 1:2
       
       %            shapeTmp = ccc(subrun).PixelList ;
      %             shapeX = shapeTmp(:,1);
        %           shapeY = shapeTmp(:,2);
    
         %           plot(shapeX, shapeY,'xr') 
         %        end

     
          %  print(h,'-dpdf', [ tifFilename,'_Chromatin.pdf']);%tifCurvetifFilename);
             close all
    
  
        
    
   
     end
    

