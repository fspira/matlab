function imgMergeSaveVertical = doDrawROI16Bit(imgMerge,BWVer,analysisFrame)

BW2 = edge(BWVer,'canny',0.1,0.1);
BW3 = bwboundaries(BW2);
BW3Work = BW3{1};

BWTest = BW2;
BWTest(:,:) =0;


imgMergeSaveVertical = imgMerge(:,:,:,analysisFrame);


for lauf = 1:length(BW3Work)
    
    imgMergeSaveVertical(BW3Work(lauf,1),BW3Work(lauf,2),:,:) = 65000;
    
end