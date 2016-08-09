function [ redMax1_new ] = detectCenterInterpolated_lower(redStackMid,flank1Tmp,rMark1,cMark1)
redStackMid(:,:) = 0;
redStackMid_center = redStackMid;
rMark1 = rMark1(1)
cMark1 = cMark1(1)
%flank1Tmp = flank2Tmp
        for subrun = 1:length(flank1Tmp)

            redStackMid(round(flank1Tmp(2,subrun)),round(flank1Tmp(1,subrun))) = 25000;

        end
    

    [mm nn pp] = size(redStackMid)
   vLine = rMark1 +15
   %hLine=1
    for subrun = 1: 45
        
        redStackMid_center(vLine, cMark1) = 25000;
        redStackMid_center(vLine+1, cMark1+1) = 25000;
        vLine = vLine -1

    end

    mergeFileCenter = redStackMid_center + redStackMid;

    imshow(mergeFileCenter(:,:,1),[])
    
    max_I = max(max(mergeFileCenter(:,:,1)))
   
    [xCoord yCoord] = find(mergeFileCenter == max_I)
    
    redMax1x_Tmp = find(round(flank1Tmp(2,:)) == xCoord(1))  
    redMax1y_Tmp = find(round(flank1Tmp(1,:)) == yCoord(1))
    
   
     redMax1_newTmp = intersect(redMax1x_Tmp, redMax1y_Tmp)
    
        
    redMax1_new = redMax1_newTmp(1)

end

