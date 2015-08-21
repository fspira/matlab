   
function [fitresult, gof] = fitExp2815(timeRatio,mean_Curve) 
%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( timeRatio, mean_Curve );

% Set up fittype and options.
ft = fittype( 'exp1' );
opts = fitoptions( ft );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.Lower = [-Inf -Inf];
opts.Robust = 'LAR';
opts.StartPoint = [max(mean_Curve) min(mean_Curve)];
opts.Upper = [Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );
thalbeSingle = log(0.5) / (fitresult.b)
confidenceIntervalsSingle = log(0.5) ./ (confint(fitresult))

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'mean_meta_VerticalSpeed vs. timeRatio', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel( 'timeRatio' );
ylabel( 'mean_meta_VerticalSpeed' );
grid on