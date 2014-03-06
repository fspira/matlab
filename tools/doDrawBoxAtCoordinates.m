function [BW_Box_Out BWHor_Out BWVer]=doDrawBoxAtCoordinates(x,y,S1)
    
    yTest =    [y-shortAxis y+shortAxis y+shortAxis y-shortAxis y-shortAxis]'
    xTest =    [x-longAxis x-longAxis x+longAxis x+longAxis x-longAxis]'

    yTestVer = [y-longAxis y-longAxis y+longAxis y+longAxis y-longAxis]'
    xTestVer = [x-shortAxis x+shortAxis x+shortAxis x-shortAxis x-shortAxis]'


    %%%%%% Total number of pixels within the box
    BW_Box = (yTest(2)-yTest(1)) * (xTest(3) - xTest(2));

    %%%%%% Draw the horizontal box
    BWHor = roipoly(S1(:,:,analysisFrame),xTest, yTest);

    %%%%% Draw the vertical box
    BWVer = roipoly(S1(:,:,analysisFrame),xTestVer, yTestVer);

    BW_Box_Out = BW_Box;
    BWHor_Out =BWHor;
    BWVer_Out = BWVer;
    
    
    
    
    
end