function [greenVariablesNorm]=doNormalizeDataToAnaphase(time, greenVariables)

[m n p] = size(greenVariables)

    for lauf = 1:n

        timeTmp = time(:,lauf)
        timeTmp(isnan(timeTmp(:,1)),:)=[]

         zeroTimeTmp = find(timeTmp(:,1)==0);
                


         zeroTime = find(time(:,lauf)==0);

         if zeroTimeTmp >= 3

             greenVariablesNorm(:,lauf) = greenVariables(:,lauf) ./ mean(greenVariables(zeroTime-2:zeroTime,lauf))

         elseif zeroTimeTmp == 2
             greenVariablesNorm(:,lauf) = greenVariables(:,lauf) ./ mean(greenVariables(zeroTime-1:zeroTime,lauf))

         elseif zeroTimeTmp == 1
             greenVariablesNorm(:,lauf) = greenVariables(:,lauf) ./ mean(greenVariables(zeroTime:zeroTime,lauf))


         end

    end
end


%greenTest = greenVariables

 % greenIdx = find(greenTest<0);
 % greenTest(greenIdx) = 0