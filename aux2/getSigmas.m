function [sigmaX, sigmaY, deltasigmaX, deltasigmaY, xCenter, yCenter] = getSigmas(x,y,dose,deltaDose)

if size(x,1) < size(x,2)
    x = x';
end

if size(y,1) < size(y,2)
    y = y';
end

sumX = sum(dose, 1);
sumX = sumX - min(sumX);
errSumX = rssq(deltaDose, 1);
sumY = sum(dose, 2);
sumY = sumY - min(sumY);
errSumY = rssq(deltaDose, 2);

wx = errSumY.^(-2);
xmask = ~isnan(wx) & ~isinf(wx);
FX = fit(x(xmask), sumY(xmask), 'gauss1', 'Weight', wx(xmask));
FXerrors = confint(FX);
xCenter = FX.b1;
sigmaX = FX.c1 / sqrt(2);
deltasigmaX = max(abs(FXerrors(:,3)-FX.c1)) / sqrt(2);

wy = errSumX.^(-2);
ymask = ~isnan(wy) & ~isinf(wy);
FY = fit(y(ymask), sumX(ymask)', 'gauss1', 'Weight', wy(ymask)');
FYerrors = confint(FY);
yCenter = FY.b1;
sigmaY = FY.c1 / sqrt(2);
deltasigmaY = max(abs(FYerrors(:,3)-FY.c1)) / sqrt(2);

end

