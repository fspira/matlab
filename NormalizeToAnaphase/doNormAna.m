
function [normGreenStore2] = doNormAna(greenStore2,anaOnset)

%greenStore2 = meanLifeactMRLC

[m n p] = size(greenStore2);





    greenEvalTmp = greenStore2;
    
    greenStore2Tmp = greenEvalTmp(isfinite(greenEvalTmp(:, 1)), :)
    

    normGreenStore2  =  greenStore2Tmp ./ greenStore2Tmp(anaOnset); 
   
   
    
    

    
    
  

