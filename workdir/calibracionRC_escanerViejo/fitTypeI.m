function [fitresult, gof] = fitTypeI(dosesGy, aa, aa_errors)
%CREATEFIT(DOSESGY,AA)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : dosesGy
%      Y Output: aa
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 06-Aug-2020 12:22:06


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( dosesGy, aa );

% Set up fittype and options.
ft = fittype( 'alpha+beta./(x-gamma)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Robust = 'LAR';
opts.StartPoint = [0.04 0.24 -0.8];
% Fit model to data.
[fitresult] = fit( xData, yData, ft, opts);

opts.StartPoint = coeffvalues(fitresult);
opts.Weights = aa_errors.^(-2);
[fitresult, gof] = fit( xData, yData, ft, opts);
% % Plot fit with data.

% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% hold on
% errorbar(xData, yData, aa_errors, 'b.');
% legend( h, 'aa vs. dosesGy', 'untitled fit 1', 'Location', 'NorthEast', 'Interpreter', 'none' );
% % Label axes
% xlabel( 'dosesGy', 'Interpreter', 'none' );
% ylabel( 'aa', 'Interpreter', 'none' );
% grid on


