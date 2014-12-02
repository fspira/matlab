
   

%%%%%% Test whether the segmented region was flipped and therefore whether
%%%%%% the position weighted centroid was positioned correctly

function [chromatin1Mid chromatin2Mid] = doSortChromatinPositions(chromatin1Mid, chromatin2Mid, AnalysisEnd)

%%%%% debugging variables
%cChromatin1Func = cChromatin1
%rChromatin1Func = rChromatin1
%cChromatin2Func = cChromatin2
%rChromatin2Func = rChromatin2
   chromatin1TmpUpdate =[];
   chromatin2TmpUpdate =[];
    
  
    
    chromatin1FirstCenter = chromatin1Mid(1,1:2);
    chromatin2FirstCenter = chromatin2Mid(1,1:2);
    
    for lauf = 1:AnalysisEnd

        if lauf == AnalysisEnd

                     chromatin1Dist = pdist2([  chromatin1FirstCenter],[chromatin1Mid(lauf,1:2)])
                     chromatin2Dist = pdist2([  chromatin2FirstCenter],[chromatin2Mid(lauf, 1:2)])

                    if chromatin1Dist >chromatin2Dist

                       
                        chromatin1MidTmp = chromatin1Mid(lauf)
                        chromatin2MidTmp = chromatin2Mid(lauf)

                        chromatin1Mid(lauf) = chromatin2MidTmp;
                        chromatin2Mid(lauf) = chromatin1MidTmp;

                    else

                    end  


            else

           chromatin1Dist = pdist2([  chromatin1FirstCenter],[chromatin1Mid(lauf,1:2)])
           chromatin2Dist = pdist2([  chromatin2FirstCenter],[chromatin2Mid(lauf,1:2)])

            if chromatin1Dist  > chromatin2Dist

                     

                        chromatin1MidTmp = chromatin1Mid(lauf)
                        chromatin2MidTmp = chromatin2Mid(lauf)

                        chromatin1Mid(lauf) = chromatin2MidTmp;
                        chromatin2Mid(lauf) = chromatin1MidTmp;

                    
            else

            end
            
        end
        
    end
    
