function out = doEliminateNegativeValues(in)

inTmp = in;

for iSub = 1: length(in)
    
    if inTmp(iSub) ==0
    inTmp(iSub) = 1;
    elseif inTmp(iSub) < 0
        inTmp(iSub) = sqrt((inTmp(iSub)^2));
    end
    
    out = inTmp;
    
end