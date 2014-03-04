function [BW_Box BWHor BWVer x y] = doDrawAnalysisBox(imgMerge,analysisFrame,ratio,voxelX_mum)

imshow(imgMerge(:,:,:,analysisFrame),'InitialMagnification',1000)


%colormap(jet)
[x,y] = ginput(1)
close all
%%%% This command draws an hourglass
%yTest = [y-10 y-10 y+10 y+10]
%xTest = [x+10 x-10 x+10 x-10];

%%%%%%%% setup of the ROI - 1.5µm x 4µm

shortAxis = round((0.5 /voxelX_mum)/2);
longAxis = round((1 / voxelX_mum)/2);


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
end
