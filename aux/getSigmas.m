function [sigmaX, sigmaY, deltasigmaX, deltasigmaY, xCenter, yCenter, FWHMx, FWHMy] = getSigmas(x, y, dose, deltaDose)
% [sigmaX, sigmaY, deltasigmaX, deltasigmaY, xCenter, yCenter] = getSigmas(x,y,dose,deltaDose*)
% Fit to a single gaussian in X and Y and calculate center position, sigmaX
% and sigmaY
% Inputs: X, Y positions in cm
% Dose matrix
% deltaDose matrix*

% 0. Default weights if no errors are assigned
if nargin == 3
    deltaDose = 0.01*mean(dose(:))*ones(size(dose));
end

% 0. Transpose matrices if necessary
if size(x,1) < size(x,2)
    x = x';
end

if size(y,1) < size(y,2)
    y = y';
end

% 1. Filler for nan values
dose(isnan(dose)) = 0;
deltaDose(isnan(deltaDose)) = 0.3;

% 2. Create baselines as the mean of the N lowest values
NlowVal = min([5,floor(numel(x)/10),floor(numel(y)/10)]);
sumX = mean(dose, 1, 'omitnan');
sortedSumX = sort(sumX);
baseLineX = mean(sortedSumX(1:NlowVal), 'omitnan');
sumX = sumX - baseLineX;
sumY = mean(dose, 2, 'omitnan');
sortedSumY = sort(sumY);
baseLineY = mean(sortedSumY(1:NlowVal), 'omitnan');
sumY = sumY - baseLineY;

% 3. Determine weights of each point in X and Y and masks to avoid errors
errSumX = sqrt(mean(deltaDose.^2, 1, 'omitnan'));
errSumY = sqrt(mean(deltaDose.^2, 2, 'omitnan'));
wx = errSumY.^(-2);
wy = errSumX.^(-2);
xmask = ~isnan(wx) & ~isinf(wx) & ~isnan(sumY);
ymask = ~isnan(wy) & ~isinf(wy) & ~isnan(sumX);

% 4. Fits in X and Y and output information
FX = myGausFit(x(xmask), sumY(xmask), wx(xmask));
xCenter = FX.b;
sigmaX = FX.c / sqrt(2);
FXerrors = confint(FX);
deltasigmaX = max(abs(FXerrors(:,4)-FX.c)) / sqrt(2);

FY = myGausFit(y(ymask), sumX(ymask), wy(ymask));
yCenter = FY.b;
sigmaY = FY.c / sqrt(2);
FYerrors = confint(FY);
deltasigmaY = max(abs(FYerrors(:,4)-FY.c)) / sqrt(2);

end

