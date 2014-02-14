function [first_css second_css] = doCalculateCenter(B_FirstBoundary,B_SecondBoundary,BW2,img)

    BW_first = BW2;
    BW_first(:,:) = 0;
    
    BW_second = BW2;
    BW_second(:,:) = 0;
    
    for subLauf = 1:length(B_FirstBoundary)
    
     BW_first(B_FirstBoundary(subLauf,1),B_FirstBoundary(subLauf,2)) = 255;
    end
    
     for subLauf = 1:length(B_SecondBoundary)
    
   
     BW_second(B_SecondBoundary(subLauf,1),B_SecondBoundary(subLauf,2)) = 255;
    end
    
    imshow(BW_first)
    firstBoundary = imfill(BW_first,'holes');
    secondBoundary =  imfill(BW_second,'holes');
    
    imshow(secondBoundary)
    
    firstCenterImg = uint16(firstBoundary) .* uint16(img);
    secondCenterImg = uint16(secondBoundary) .* uint16(img);
    
    first_css = regionprops(BW_first, firstCenterImg, {'WeightedCentroid'});
    second_css = regionprops(BW_second, secondCenterImg, {'WeightedCentroid'});
    
   
   