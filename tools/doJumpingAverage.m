function [matOut] = doJumpingAverage(matIn, windowSize,matEnd)

 % filter length
 %matEnd = 246;
 %testData1 = testData;
%matIn = testData;
 
 matLength = length(matIn)

 %testData(isfinite(testData(:, 1)), :)
 
matDimension = matLength/windowSize

 
xx = reshape(matIn(1:matEnd)',windowSize,[]);
yy = sum(xx,1)./size(xx,1);
y = reshape(repmat(yy, size(xx,1),1),1,[]);
y = y'

matOut = y
