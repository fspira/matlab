
%%%%%% Test whether the segmented region was flipped and therefore whether
%%%%%% the position weighted centroid was positioned correctly

function [cChromatin1FuncOut rChromatin1FuncOut cChromatin2FuncOut rChromatin2FuncOut] = doSortChromatinPositions(cChromatin1Func, rChromatin1Func, cChromatin2Func, rChromatin2Func)

%%%%% for debuggin variables
%cChromatin1Func = cChromatin1
%rChromatin1Func = rChromatin1
%cChromatin2Func = cChromatin2
%rChromatin2Func = rChromatin2

cChromatin1FuncOut(1) = cChromatin1Func(1)
rChromatin1FuncOut(1) = rChromatin1Func(1)
cChromatin2FuncOut(1) = cChromatin2Func(1)
rChromatin2FuncOut(1) = rChromatin2Func(1)


for lauf = 2:length(cChromatin1Func)
    
  
    
        distSort1 = pdist2([cChromatin1Func(lauf-1),rChromatin1Func(lauf-1)],[cChromatin1Func(lauf),rChromatin1Func(lauf)]);
        distSort2 = pdist2([cChromatin2Func(lauf-1),rChromatin2Func(lauf-1)],[cChromatin2Func(lauf),rChromatin2Func(lauf)]);
        
        distCrossSort1 = pdist2([cChromatin1Func(lauf-1),rChromatin1Func(lauf-1)],[cChromatin2Func(lauf),rChromatin2Func(lauf)]);
        distCrossSort2 = pdist2([cChromatin2Func(lauf-1),rChromatin2Func(lauf-1)],[cChromatin1Func(lauf),rChromatin1Func(lauf)]);
        
    if distSort1 - distCrossSort1 == 0
        
         cChromatin1FuncOut(lauf) =  cChromatin1Func(lauf)
        rChromatin1FuncOut(lauf) =  rChromatin1Func(lauf)
        
        cChromatin2FuncOut(lauf) =  cChromatin2Func(lauf)
        rChromatin2FuncOut(lauf) =  rChromatin2Func(lauf)
        
    elseif  distCrossSort1 > distSort1 || distCrossSort2 > distSort2
        
        cChromatin1FuncOut(lauf) =  cChromatin2Func(lauf);
        rChromatin1FuncOut(lauf) =  rChromatin2Func(lauf);
        
        cChromatin2FuncOut(lauf) =  cChromatin1Func(lauf);
        rChromatin2FuncOut(lauf) =  rChromatin1Func(lauf);
        
       % loop = 1
    else
        
        cChromatin1FuncOut(lauf) =  cChromatin1Func(lauf)
        rChromatin1FuncOut(lauf) =  rChromatin1Func(lauf)
        
        cChromatin2FuncOut(lauf) =  cChromatin2Func(lauf)
        rChromatin2FuncOut(lauf) =  rChromatin2Func(lauf)
        
      %  loop = 2
    end
end
  
    
end