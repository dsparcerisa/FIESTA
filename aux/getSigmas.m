function [sigmaX, sigmaY, deltasigmaX, deltasigmaY, xCenter, yCenter] = getSigmas(x,y,dose,deltaDose,offCenterSigmas)

if size(x,1) < size(x,2)
    x = x';
end

if size(y,1) < size(y,2)
    y = y';
end

if nargin <5 
    offCenterSigmas = 3;
end
    
dose(isnan(dose)) = 0;
deltaDose(isnan(deltaDose)) = 0.3;

% Baseline as 0.1% of max value
%sortedDose = sort(dose(:));
% baseLineDose = mean(sortedDose(1:10), 'omitnan')
% dose(dose<0) = 0;
% baseLineDose = mean(dose(dose<0.01*max(dose(:))))
% dose = dose - baseLineDose;

% figure(1);
sumX = mean(dose, 1, 'omitnan');
sortedSumX = sort(sumX);
baseLineX = mean(sortedSumX(1:5), 'omitnan');
%baseLineX = max(min(sumX),0.5*(sumX(1) + sumX(end)));
%baseLineX = min(sumX);
%baseLineX = 0;
sumX = sumX - baseLineX;
errSumX = sqrt(mean(deltaDose.^2, 1, 'omitnan'));
sumY = mean(dose, 2, 'omitnan');
sortedSumY = sort(sumY);
baseLineY = mean(sortedSumY(1:5), 'omitnan');
%baseLineY = max(min(sumY),0.5*(sumY(1) + sumY(end)));
%baseLineY = min(sumY);
%baseLineY = 0;
sumY = sumY - baseLineY;
errSumY = sqrt(mean(deltaDose.^2, 2, 'omitnan'));

wx = errSumY.^(-2);
xmask = ~isnan(wx) & ~isinf(wx) & ~isnan(sumY);

FX = myGausFit(x(xmask), sumY(xmask), wx(xmask));
%FX = fit(x(xmask), sumY(xmask), 'gauss2', 'Weight', wx(xmask), 'Lower', [0 -10 0 -2 -10 100], 'Upper', [500 10 10 2 10 10000])
% FXerrors = confint(FX);
xCenter = FX.b;
sigmaX = FX.c / sqrt(2);
% deltasigmaX = max(abs(FXerrors(:,3)-FX.c1)) / sqrt(2);

subplot(1,2,1); hold off
% errorbar(x(xmask), sumY(xmask), errSumY(xmask), 'b.'); hold on
errorbar(x, sumY, errSumY, 'b.'); hold on
plot(x,FX(x),'r-');
title('FX');

wy = errSumX.^(-2);
ymask = ~isnan(wy) & ~isinf(wy) & ~isnan(sumX);
%FY = fit(y(ymask), sumX(ymask)', 'gauss2', 'Weight', wy(ymask)', 'Lower', [0 -10 0 -2 -10 100], 'Upper', [500 10 10 2 10 10000])
FY = myGausFit(y(ymask), sumX(ymask), wy(ymask));
%FYerrors = confint(FY);
yCenter = FY.b;
sigmaY = FY.c / sqrt(2);
%deltasigmaY = max(abs(FYerrors(:,3)-FY.c1)) / sqrt(2);

subplot(1,2,2); hold off
%errorbar(y(ymask), sumX(ymask), errSumX(ymask), 'b.'); hold on
errorbar(y, sumX, errSumX,'b.'); hold on
plot(y,FY(y),'r-');
title('FY');
% pause

% Erase points outside N sigmas
distances = (x-xCenter).^2 + (y'-yCenter).^2;
outMask = distances > (offCenterSigmas * max(sigmaX, sigmaY));
% figure(2)
% subplot(1,3,1);
% imagesc(dose);
% subplot(1,3,2)
% imagesc(outMask)
% subplot(1,3,3);
% dose(outMask) = 0;
% deltaDose(outMask) = 0.3; % TODO: mirar un valor razonable
% imagesc(dose);

% FIND LINE
% xLine = interp2(y,x,dose,yCenter*ones(size(x)),x);
% yLine = interp2(y,x,dose,xCenter*ones(size(y)),y);
% % Repeat fit
% figure(1);
sumX = mean(dose, 1, 'omitnan');
sortedSumX = sort(sumX);
baseLineX = mean(sortedSumX(1:5), 'omitnan');
%baseLineX = max(min(sumX),0.5*(sumX(1) + sumX(end)));
%baseLineX = min(sumX);
%baseLineX = 0;
sumX = sumX - baseLineX;
errSumX = sqrt(mean(deltaDose.^2, 1, 'omitnan'));
sumY = mean(dose, 2, 'omitnan');
sortedSumY = sort(sumY);
baseLineY = mean(sortedSumY(1:5), 'omitnan');
%baseLineY = max(min(sumY),0.5*(sumY(1) + sumY(end)));
%baseLineY = min(sumY);
%baseLineY = 0;
sumY = sumY - baseLineY;
errSumY = sqrt(mean(deltaDose.^2, 2, 'omitnan'));

wx = (errSumY).^(-2);
xmask = ~isnan(wx) & ~isinf(wx) & ~isnan(sumY);

%FX = fit(x(xmask), sumY(xmask), 'gauss1', 'Weight', wx(xmask), 'Lower', [0 -10 0], 'Upper', [500 10 10])
% miFun = @(x) a1*exp(-((x-b1)/c1)^2) + a2
%FX = fit(x(xmask), sumY(xmask), 'gauss2', 'Weight', wx(xmask), 'Lower', [0 -10 0 -2 0 10000], 'Upper', [500 10 10 2 0 10000])
FX = myGausFit(x(xmask), sumY(xmask), wx(xmask));
% FX2 = myGausFit(x, xLine)
FXerrors = confint(FX);
xCenter = FX.b;
sigmaX = FX.c / sqrt(2);
deltasigmaX = max(abs(FXerrors(:,4)-FX.c)) / sqrt(2);

subplot(1,2,1); hold on
%errorbar(x(xmask), sumY(xmask), errSumY(xmask), 'g.'); hold on
errorbar(x, sumY, errSumY, 'g.'); hold on
plot(x,FX(x),'m-');

wy = (errSumX).^(-2);
ymask = ~isnan(wy) & ~isinf(wy) & ~isnan(sumX);
%FY = fit(y(ymask), sumX(ymask)', 'gauss1', 'Weight', wy(ymask)', 'Lower', [0 -10 0], 'Upper', [500 10 10])
%FY = fit(y(ymask), sumX(ymask)', 'gauss2', 'Weight', wy(ymask)', 'Lower', [0 -10 0 -2 0 10000], 'Upper', [500 10 10 2 0 10000])
FY = myGausFit(y(ymask), sumX(ymask), wy(ymask));
%FY2 = myGausFit(y, yLine)
FYerrors = confint(FY);
yCenter = FY.b;
sigmaY = FY.c / sqrt(2);
deltasigmaY = max(abs(FYerrors(:,4)-FY.c)) / sqrt(2);

subplot(1,2,2); hold on
%errorbar(y(ymask), sumX(ymask), errSumX(ymask), 'g'); hold on
errorbar(y, sumX, errSumX, 'g.'); hold on
plot(y,FY(y),'m-');

end

