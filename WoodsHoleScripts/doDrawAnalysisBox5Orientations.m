function [BW_Box BWI0 BWI45 BWI90 BWI135 BWI135B x_Out y_Out...
    I0Avg I45Avg I90Avg I135Avg I0BleachAvg] = doDrawAnalysisBox5Orientations(imgMerge,I0, I45, I90, I135, I0Bleach ,voxelX_mum)

imgMerge = normalizedImage(imgMerge);

imshow(imgMerge(:,:),'InitialMagnification',1000)


%colormap(jet)
[xx,yy] = ginput(4)
close all

    for subrun = 1:4

    %%%% This command draws an hourglass
    %yTest = [y-10 y-10 y+10 y+10]
    %xTest = [x+10 x-10 x+10 x-10];

    %%%%%%%% setup of the ROI - 1.5µm x 4µm

    x = xx(subrun)
    y = yy(subrun);

    shortAxis = round((0.8 /voxelX_mum)/2);
    longAxis = round((0.8 / voxelX_mum)/2);


    %%%%%% Edges of the horizontal box
    yTest =    [y-shortAxis y+shortAxis y+shortAxis y-shortAxis y-shortAxis]'
    xTest =    [x-longAxis x-longAxis x+longAxis x+longAxis x-longAxis]'

   


    %%%%%% Total number of pixels within the box
    BW_Box{subrun} = (yTest(2)-yTest(1)) * (xTest(3) - xTest(2));
    BWBoxPixel = BW_Box{1};
    %%%%%% Draw the horizontal box
    BWI0{subrun} = roipoly(I0,xTest, yTest);
    BWI45{subrun} = roipoly(I45,xTest, yTest);
    BWI90{subrun} = roipoly(I90,xTest, yTest);
    BWI135{subrun} = roipoly(I135,xTest, yTest);
    BWI135B{subrun} = roipoly(I0Bleach,xTest, yTest);

    
    
    x_Out{subrun} = x;
    
    y_Out{subrun} = y;

    
    averageMask = uint16(BWI0{subrun}) .* uint16(I0);
    I0Avg(subrun) = (sum(sum(averageMask))/BWBoxPixel);
    
    averageMask = uint16(BWI45{subrun}) .* uint16(I45);
    I45Avg(subrun) = (sum(sum(averageMask))/BWBoxPixel);
    
    averageMask = uint16(BWI90{subrun}) .* uint16(I90);
    I90Avg(subrun) = (sum(sum(averageMask))/BWBoxPixel);
    
    averageMask = uint16(BWI135{subrun}) .* uint16(I135);
    I135Avg(subrun) = (sum(sum(averageMask))/BWBoxPixel);
    
    averageMask = uint16(BWI135B{subrun}) .* uint16(I0Bleach);
    I0BleachAvg(subrun) = (sum(sum(averageMask))/BWBoxPixel);
    
   
   
    
    end
end
