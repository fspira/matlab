function [density flagOut predicted_time] = doMapTime(timeToPredict)

load('/Users/spira/Desktop/Rpe1AnisotropiePhalloidin/analysisWorkflow/modelK9.mat');
load('/Users/spira/Desktop/Rpe1AnisotropiePhalloidin/analysisWorkflow/x1.mat');
load('/Users/spira/Desktop/Rpe1AnisotropiePhalloidin/analysisWorkflow/x2.mat');

predicted_time = polyvaln(m,timeToPredict);
[bandwidth,density,X,Y]=kde2d([x1 x2], 256, [0 0], [30 30]);

ff = @(xxx) xxx/30*255+1;
density(ff(timeToPredict(1)),ff(timeToPredict(2))) > 0.004;

df = @(xx1,yy1) density(ff(xx1),ff(yy1)) > 0.004;
flagOut = df(timeToPredict(1),timeToPredict(2));






end