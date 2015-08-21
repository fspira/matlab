function [stdGreenLine,greenMean] = calculateStd(greenStore)

[m n p] = size(greenStore);

for lauf = 1:m
    greenStoreLine(:,lauf) = greenStore(lauf,1:n);
    greenTmp = greenStoreLine(:,lauf)
    greenTmp = greenTmp(isfinite(greenTmp(:, 1)), :)

    stdGreenLine(lauf) = std(greenTmp);
    greenMean(lauf) = median(greenTmp);

end

end