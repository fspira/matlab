function [cSort1, rSort1, idxSort]  =   RefineRoi(r,c,lauf)
              
           
             
                 
              
                    idxSort =1;
                    idxNew = 1;
                    idxOld  =1 ;
                    
                    for subLauf = 1: length(r)
                    
                        [cNew1 rNew1 idxNew] = sortMaskCorrection1(c,r,idxNew(1),idxOld, idxSort);
                    
                        cSort1(subLauf) = cNew1;
                        rSort1(subLauf) = rNew1;
                        idxSort(subLauf) = idxNew
                        if subLauf == 1
                            idxOld = idxSort(subLauf)
                        else
                            idxOld = idxSort(subLauf-1)
                        end
                        
                    
                    end
               
                   