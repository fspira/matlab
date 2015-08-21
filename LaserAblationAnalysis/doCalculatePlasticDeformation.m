
function [matOut thalfOut] = doCalculatePlasticDeformation(timeRatio,matIn) 

%matIn = lateAna_parallel_MinorAxisManual
%matIn = earlyAna_parallel_MinorAxisManual(1:9,16)

[m n p] = size(matIn)
matInZero_Speed=zeros(m-length(m-1),n)
zeroArray = zeros(1,n)

  
matIn_zero = cat(1,zeroArray, matIn)

  
  
timeRatioZero = cat(1,0,timeRatio)

matInZero_SpeedFirst = (matIn_zero(1+1,1:end) - matIn_zero(1,1:end)) / timeRatio(1)

    

for lauf = 1:m
   
   matInZero_Speed(lauf,1:end) = (matIn_zero(lauf+1,1:end) - matIn_zero(lauf,1:end))/0.78

end

  matInZero_Speed = cat(1, matInZero_SpeedFirst,matInZero_Speed(2:m,:) )


for lauf = 1:n
    [fitresult gof] = fitExp2964(timeRatioZero(1:length(matInZero_Speed)), matInZero_Speed(:,lauf))
    
    thalbeSingle_mat_Speed = log(0.5) / (fitresult.b)
   thalbeStore(lauf) = thalbeSingle_mat_Speed
             
                
    confidenceIntervalsSingle_mat_Speed = log(0.5) ./ (confint(fitresult))
end


matOut =   matInZero_Speed;
thalfOut = thalbeStore

