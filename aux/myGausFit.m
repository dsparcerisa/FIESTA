function [fitresult, gof] = myGausFit(x, sumX, wx)
% myGausFit(x, sumX, wx)
%  Create a fit to a Gaussian + background

if nargin==2
    wx = ones(size(x));
end

%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( x, sumX );

% Set up fittype and options.
ft = fittype( 'a0 + a*exp(-((x-b)/c)^2)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Trust-Region';
opts.Weights = wx;
opts.Display = 'Off';
opts.Lower = [0 0 -100 0];
opts.Robust = 'LAR';
opts.StartPoint = [0.3619 0 0 0.7449];
opts.Upper = [Inf 0 100 100];

% Fit model to data after fixing baseline
[fitresult, gof] = fit( xData, yData, ft, opts );


