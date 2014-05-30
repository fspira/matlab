function [red_furrow_mean_ana]  = doNormaliseToAnaphase(red_furrow_mean,time)
   

%red_furrow_mean= redFurrow;
%%%%%%%%% normalize to Anaphase Onset
[m n p] = size(time);

    for lauf = 1:n
        zeroPos = find(time(:,lauf)==0)
        
        if zeroPos >= 2
        
          
            red_furrow_mean_ana(:,lauf) =  red_furrow_mean(:,lauf) ./ mean(red_furrow_mean(zeroPos,lauf));
          
            
        else
            
            red_furrow_mean_ana(:,lauf) =  red_furrow_mean(:,lauf) ./ mean(red_furrow_mean(zeroPos-2:zeroPos,lauf));
           
        end
        
      lauf
    end