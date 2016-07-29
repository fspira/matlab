function [fitresult, gof] = createFitPolyFitAngleIntensity_sine(alphaDeg, flankCrop)

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( alphaDeg, flankCrop );

ft = fittype( 'a0 + a1*sin(x*w+ p)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( ft );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';
opts.Lower = [-Inf -Inf -Inf 0];
opts.Robust = 'Bisquare';
opts.Lower = [-Inf -Inf -Inf 0.034];
opts.StartPoint = [0.959492426392903 0.655740699156587 1.636 0.0349];
opts.Upper = [Inf Inf Inf 0.035];


% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData,'.k'); hold on
%set(h, 'MarkerFaceColor', 'c')


%legend( h, 'Normalized intensity', 'Fit', 'Location', 'NorthEast' );
% Label axes
xlabel( 'Angle (deg)' );
ylabel( 'Anisotropy (a.u.)' );
%grid on


