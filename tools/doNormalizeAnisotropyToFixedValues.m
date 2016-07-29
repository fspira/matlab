function [anisotropyNormOut] = doNormalizeAnisotropyToFixedValues(time, anisotropy, metaphaseEquatorFixedCells)
[m n p] = size(anisotropy)


for lauf =1:n
   
    timeNorm = find(time(:,lauf) == 0)
  
   
    anisoTmp = anisotropy(:,lauf)
        
    anisoTmp(isnan(anisoTmp(:,1)),:)=[];

    timeTmp = time(:,lauf);
    timeTmp(isnan(timeTmp(:,1)),:)=[];
    
    timeTmp = find(timeTmp(:,1) == 0)
    
     if timeTmp >= 3

      normValue = mean(anisoTmp(timeNorm-3:timeNorm)) - metaphaseEquatorFixedCells

        elseif  timeTmp >= 2
           
      normValue = mean(anisoTmp(timeNorm-1:timeNorm)) - metaphaseEquatorFixedCells


        elseif  timeTmp >= 1

        normValue = mean(anisoTmp(timeNorm:timeNorm)) - metaphaseEquatorFixedCells


     end
        
    
    
    
    anisotropyNormOut(:,lauf) = anisotropy(:,lauf) - normValue;
    
    
end