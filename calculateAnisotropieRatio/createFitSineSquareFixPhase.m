function [fitresult, gof] = createFitSineSquareFixPhase(xComponent, yComponent)
%CREATEFIT(XCOMPONENT,YCOMPONENT)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : xComponent
%      Y Output: yComponent
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 10-Jun-2016 13:06:49


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( xComponent, yComponent );

% Set up fittype and options.
ft = fittype( 'a0 + a1*sin(x*pi/180)^2', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Display = 'Off';
opts.Lower = [-Inf -Inf];
opts.Robust = 'Bisquare';
opts.StartPoint = [0.171186687811562 0.706046088019609];
opts.Upper = [Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'yComponent vs. xComponent', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel( 'xComponent' );
ylabel( 'yComponent' );
%grid on


