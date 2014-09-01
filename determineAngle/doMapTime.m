function [density flagOut predicted_time] = doMapTime(timeToPredict)

load('/Users/spira/Desktop/programme/staging/modelK4New.mat');
load('/Users/spira/Desktop/programme/staging/x1New.mat');
load('/Users/spira/Desktop/programme/staging/x2New.mat');
%%% timeToPredict = stagingParameters
predicted_time = polyvaln(m,timeToPredict);
[bandwidth,density,X,Y]=kde2d([x1 x2], 256, [0 0], [30 30]);


timeToPredict = uint16(timeToPredict)

ff = @(xxx) xxx/30*255+1;
density(ff(timeToPredict(1)),ff(timeToPredict(2))) > 0.0001;

df = @(xx1,yy1) density(ff(xx1),ff(yy1)) > 0.0001;
flagOut = df(timeToPredict(1),timeToPredict(2));






end