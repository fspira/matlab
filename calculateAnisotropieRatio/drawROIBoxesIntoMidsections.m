%%%% Plot the ROI used for analysis and save the file


distvec = linspace(0,1,256)';
red = zeros(256,3);
red(:,1) = distvec;

green = zeros(256,3);
green(:,2) = distvec;


for lauf =1 :length(xFurrow1)

    greenStackRGB(:,:,:,lauf) = ind2rgb(S1(:,:,lauf),green);
    redStackRGB(:,:,:,lauf) = ind2rgb(S2(:,:,lauf),red);
             
    imgMergeTmp(:,:,:,lauf) = greenStackRGB(:,:,:,lauf) + redStackRGB(:,:,:,lauf);
    imgMerge(:,:,:,lauf) = normalizedImage3D(imgMerge(:,:,:,lauf));
    
    
end


for lauf = 1:length(xFurrow1)

    x = xFurrow1{lauf}
    y = yFurrow1{lauf}
    
    S1_In = imgMerge(:,:,:,lauf);    
    
 [BW_Box_Out BWHor_Out BWVer,imgMergeSaveVertical]=doDrawBoxAtCoordinates_Hor_1_1(x,y,S1_In,voxelX_mum);
 
    x = xFurrow2{lauf}
    y = yFurrow2{lauf}
 
 [BW_Box_Out BWHor_Out BWVer,imgMergeSaveVertical]=doDrawBoxAtCoordinates_Hor_1_1(x,y,imgMergeSaveVertical,voxelX_mum);
 
    x = xPole1{lauf}
    y = yPole{lauf}
 
 [BW_Box_Out BWHor_Out BWVer,imgMergeSaveVertical]=doDrawBoxAtCoordinates_Ver_1_1(x,y,imgMergeSaveVertical,voxelX_mum);
 
 
    x = xPole2{lauf}
    y = yPole2{lauf}
 
 [BW_Box_Out BWHor_Out BWVer,imgMergeSaveVertical]=doDrawBoxAtCoordinates_Ver_1_1(x,y,imgMergeSaveVertical,voxelX_mum);
 
    saveImg(:,:,:,lauf) = imgMergeSaveVertical;
 
end

tiffwrite_RBG(saveImg,[tifFilename,'_ROI_1_1.tif']);


clear all