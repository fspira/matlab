function [cNew rNew minIndexStore] = sortMaskCorrection(c,r)
inkrLauf =1;
        for statfix = 2:length(r)
                    
            for subrun = 1:length(r)
                    
            
                dist1(subrun) = pdist2([c(statfix),r(statfix)],[c(subrun),r(subrun)]);
              
            
            end
                     
                     
                        d1Shape = dist1;
                
                        d1Shape(dist1 == 0) = NaN;
                
                        dist1Min = min(d1Shape);
    
                        minIndex = find(dist1  == dist1Min);
    
                       cNew(inkrLauf) = c(minIndex(1));
                       rNew(inkrLauf) = r(minIndex(1));
                       minIndexStore(inkrLauf) = minIndex(1); 
                       
                       minIndex
                       
                       inkrLauf = inkrLauf +1;
                       
                       
                       %subrun
                       %statfix
                       
       	end