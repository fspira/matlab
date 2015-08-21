function [BW_Box_Out BWHor_Out BWVer_Out x_Out y_Out] = doDrawAnalysisBox12x1_ver(imgMerge,analysisFrame,ratio,voxelX_mum)

imshow(imgMerge(:,:,:,analysisFrame),'InitialMagnification',1000)


%colormap(jet)
[xx,yy] = ginput(4)
close all

    for subrun = 1:4

    %%%% This command draws an hourglass
    %yTest = [y-10 y-10 y+10 y+10]
    %xTest = [x+10 x-10 x+10 x-10];

    %%%%%%%% setup of the ROI - 1.5µm x 4µm

    x = xx(subrun);
    y = yy(subrun);

    shortAxis = round((1 /voxelX_mum)/2);
    longAxis = round((1.2 / voxelX_mum)/2);


    %%%%%% Edges of the horizontal box
    yTest =    [y-shortAxis y+shortAxis y+shortAxis y-shortAxis y-shortAxis]'
    xTest =    [x-longAxis x-longAxis x+longAxis x+longAxis x-longAxis]'

    yTestVer = [y-longAxis y-longAxis y+longAxis y+longAxis y-longAxis]'
    xTestVer = [x-shortAxis x+shortAxis x+shortAxis x-shortAxis x-shortAxis]'


    %%%%%% Total number of pixels within the box
    BW_Box = (yTest(2)-yTest(1)) * (xTest(3) - xTest(2));

    %%%%%% Draw the horizontal box
    BWHor = roipoly(ratio(:,:,analysisFrame),xTest, yTest);

    %%%%% Draw the vertical box
    BWVer = roipoly(ratio(:,:,analysisFrame),xTestVer, yTestVer);

    BW_Box_Out{subrun} = BW_Box;
    BWHor_Out{subrun} =BWHor;
    BWVer_Out{subrun} = BWVer;
    x_Out{subrun} = x;
    y_Out{subrun} = y;
    
    end
end
