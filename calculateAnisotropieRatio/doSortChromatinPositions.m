
%%%%%% Test whether the segmented region was flipped and therefore whether
%%%%%% the position weighted centroid was positioned correctly

function [cChromatin1 rChromatin1 cChromatin2 rChromatin2] = doSortChromatinPositions(cChromatin1, rChromatin1, cChromatin2, rChromatin2)

for lauf = 2:length(cChromatin1)
    
  
    
        distSort1 = pdist2([cChromatin1(lauf-1),rChromatin1(lauf-1)],[cChromatin1(lauf),rChromatin1(lauf)]);
        distSort2 = pdist2([cChromatin2(lauf-1),rChromatin2(lauf-1)],[cChromatin2(lauf),rChromatin2(lauf)]);
        
        distCrossSort1 = pdist2([cChromatin1(lauf-1),rChromatin1(lauf-1)],[cChromatin2(lauf),rChromatin2(lauf)]);
        distCrossSort2 = pdist2([cChromatin2(lauf-1),rChromatin2(lauf-1)],[cChromatin1(lauf),rChromatin1(lauf)]);
        
    if distSort1 > distCrossSort1 || distSort2 > distCrossSort2
        
        cChromatin1(lauf) =  cChromatin2(lauf);
        rChromatin1(lauf) =  rChromatin2(lauf);
        
        cChromatin2(lauf) =  cChromatin1(lauf);
        rChromatin2(lauf) =  rChromatin1(lauf);
        
        
    else
        
        cChromatin1(lauf) =  cChromatin1(lauf);
        rChromatin1(lauf) =  rChromatin1(lauf);
        
        cChromatin2(lauf) =  cChromatin2(lauf);
        rChromatin2(lauf) =  rChromatin2(lauf);
        
        
    end
  
    
end