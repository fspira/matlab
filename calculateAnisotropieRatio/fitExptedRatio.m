function [fitresult, gof] =fitExptedRatio(test1, test)
%CREATEFIT(TEST1,TEST)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : test1
%      Y Output: test
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 30-Oct-2014 14:33:38


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( test1, test );

% Set up fittype and options.
ft = fittype( 'poly6' );
opts = fitoptions( ft );
opts.Lower = [-Inf -Inf -Inf -Inf -Inf -Inf -Inf];
opts.Robust = 'LAR';
opts.Upper = [Inf Inf Inf Inf Inf Inf Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'test vs. test1', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel( 'Distance [�m]' );
ylabel( 'Normalized intensity [A.U.}' );
grid on


