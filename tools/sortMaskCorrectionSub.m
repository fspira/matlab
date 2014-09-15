function [minIndex, subFunc] = sortMaskCorrectionSub(d1Shape,dist1,dist1Min, minIndex,idxSort)


            if length(find(minIndex  == idxSort))  >= 2  ;
                
               
                            
                d1Shape(minIndex) = NaN;
                dist1Min = min(d1Shape);
                minIndex = find(dist1  == dist1Min);
            
                subFunc =1;
            
            else
                subFunc =0;
            end
            
end
    