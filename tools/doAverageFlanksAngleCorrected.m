function [greenOut midOut] = doAverageFlanksAngleCorrected(redMax1,redMax2,yGreen1,yGreen2)
%%%%% Not adjusted still orginal file
for subrun = 1:length(yGreen1)

    maxPos1 = redMax1(subrun);
    maxPos2 = redMax2(subrun);
    
    greenSub1 = yGreen1{subrun};
    greenSub2 = yGreen2{subrun};
    
    
    lengthAB(1) = length(greenSub1) - maxPos1
    lengthAB(2) = length(greenSub2) - maxPos2
    
    maxPosVec(1) = maxPos1;
    maxPosVec(2) = maxPos2;
    
    %%%% Average curves bigger than the middle peak
    index = 1;
    for subsubRun = 0:(min(lengthAB))
       
        greenMergeRight(index) = (greenSub1(maxPos1+subsubRun)+greenSub2(maxPos2+subsubRun))/2;
     
        index = index+1;
        
    end
        
    
    %%%% Average curves smaller than the middle peak
    index = 1;
    for subsubRun = 1:(min(maxPosVec))-1
       
        greenMergeLeft((min(maxPosVec))-subsubRun) = (greenSub1(maxPos1-subsubRun)+greenSub2(maxPos2-subsubRun))/2;
      
        
        index = index+1;
        
    end    
    
    
    greenOut{subrun} = cat(2,greenMergeLeft, greenMergeRight)';
    
   
    
    midOut(subrun) = min(maxPosVec)
         
    clear greenMergeRight redMergeRight greenMergeLeft redMergeLeft;
    
end