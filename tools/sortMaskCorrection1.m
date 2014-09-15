function [cNew rNew minIndexStore] = sortMaskCorrection1(c,r,statfix,minIndexOld, idxSort)

%minIndexOld = 36;
      %statfix = 37;
                    
            for subrun = 1:length(r)
                    
            
                dist1(subrun) = pdist2([c(statfix),r(statfix)],[c(subrun),r(subrun)]);
              
            
            end
                     
                     
                        d1Shape = dist1;
                
                        d1Shape(dist1 == 0) = NaN;
                
                        dist1Min = min(d1Shape);
    
                        minIndex = find(dist1  == dist1Min);
                        
                        if length(minIndex) > 1
                            
                            
                            minIndex = minIndex(length(minIndex))
                            
                        elseif minIndex == minIndexOld
                            
                            
                            
                            d1Shape(minIndexOld) = NaN;
                           
                             dist1Min = min(d1Shape);
    
                             minIndex = find(dist1  == dist1Min);
                          
                        end
                        
                      %  subFunc = 1;
                        
                      %   while subFunc == 1
                        
                      %      [minIndex, subFunc] = sortMaskCorrectionSub(d1Shape,dist1,dist1Min, minIndex(1),idxSort)
                      %      subFunc
                      %   end

    
                       cNew = c(minIndex(1));
                       rNew = r(minIndex(1));
                       minIndexStore = minIndex(1); 
                       
                   
                       
                       
                       %subrun
                       %statfix
                       
       