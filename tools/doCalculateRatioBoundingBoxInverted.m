function [AvgRatio,S1Avg,S2Avg] = doCalculateRatioBoundingBoxInverted(ratio, BWBox, BWBoxPixel, S1, S2)

%BWBox = BWHor;
%BWBoxPixel = BW_Box;

    %%%%%%% Calculate D and ratio for horizontal box
   % roiMaskD(:,:,lauf) = D_Anisotropie(:,:,lauf).*BWHor;
   % averageRatioD(lauf) = (sum(sum(roiMaskD(:,:,lauf)))/BW_Box);
   
   lauf = 1
   
    roiMask(:,:,lauf) = ratio(:,:,lauf).*BWBox;
    averageRatio(lauf) = (sum(sum(roiMask(:,:,lauf)))/BWBoxPixel);
    
    
    
    %%%% This part averages intensities within the box and subsequently
    %%%% calculates the ratio
    
  % roiMask(:,:,lauf) = ratio(:,:,lauf).*BWHor;
    
   S1Mask(:,:,lauf) = S1(:,:,lauf) .*BWBox;
   S2Mask(:,:,lauf) = S2(:,:,lauf) .*BWBox;
   S1Avg(lauf) = (sum(sum(S1Mask(:,:,lauf)))/BWBoxPixel);
   S2Avg(lauf) = (sum(sum(S2Mask(:,:,lauf)))/BWBoxPixel);
   
   
   AvgRatio(lauf) = S2Avg(lauf)/S1Avg(lauf);
     
  

