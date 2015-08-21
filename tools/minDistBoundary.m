function out = minDistBoundary(in)


inTmp = in;

for iSub = 1: length(in)
    
    if inTmp(iSub) <5
        inTmp(iSub) = 5;
        
        
    end
    
    out = inTmp;




end