function  slidingWindowIntensities = getIntensityUnderSlidingWindow(redLinescan,sWSet,furrowCenter)

yRed1 = redLinescan;
%sWSet = 4
redMax = furrowCenter




%sWSet = 40;
sW = sWSet;
sW1 = sW; %%% in pixels
sW2 = sW;

lauf = 1;

    redValuesSliding_1tmp = redLinescan;
    
        %%%%% This section ensures that the sliding windows fits the
        %%%%% selection. If not the selection is cropped to fit the window
    
    %    redMax1(dist1 == 1) = 2;
        
        if redMax(lauf)+ (sW1/2) > length(redValuesSliding_1tmp)
            
            sW1 = round(length(redValuesSliding_1tmp)-(redMax(lauf)))-1
            if  mod(sW1,2) ==1
                sW1 = sW1 -1;
            end
            
        elseif redMax(lauf)-(sW1/2) < 1
            sW1 = round(redMax(lauf)/2)-1
            if  mod(sW1,2) ==1
                sW1 = sW1 -1;
            end
        end
        
    slidingWindowIntensities = redValuesSliding_1tmp(round(redMax-(sW1/2)):round(redMax+(sW1/2)));
    
    
