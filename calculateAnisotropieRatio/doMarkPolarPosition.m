function [ x_Out y_Out] = doMarkPolarPosition(imgMergeIn,analysisFrame)

imgMerge = normalizedImage(imgMergeIn);

imshow(imgMerge(:,:,:,analysisFrame),'InitialMagnification',1000)


%colormap(jet)
[xx,yy] = ginput(2)
close all

    for subrun = 1:2

    %%%% This command draws an hourglass
    %yTest = [y-10 y-10 y+10 y+10]
    %xTest = [x+10 x-10 x+10 x-10];

    %%%%%%%% setup of the ROI - 1.5µm x 4µm

    x = xx(subrun);
    y = yy(subrun);

    x_Out{subrun} = x;
    y_Out{subrun} = y;
    
    end
end
