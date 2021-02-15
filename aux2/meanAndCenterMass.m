function [mask, xcentre, ycentre, sigmaX, sigmaY, meanvalue, maxX, maxY, deltasigmaX, deltasigmaY, stdDose] = meanAndCenterMass(img,radius_pixels)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Xvalues = 1:size(img, 2);
Yvalues = 1:size(img, 1);
[x, y] = meshgrid(Xvalues, Yvalues);
%weightedx = x .* img;
%weightedy = y .* img;
%xcentre = sum(weightedx(:)) / sum(img(:));
%ycentre = sum(weightedy(:)) / sum(img(:));
sumX = sum(img, 1, 'omitnan');
sumX = - (sumX - max(sumX));
sumY = sum(img, 2, 'omitnan');
sumY = - (sumY - max(sumY));
F1 = fit(Xvalues', sumX', 'gauss1');
F1errors = confint(F1);
xcentre = F1.b1;
sigmaX = F1.c1 / sqrt(2);
deltasigmaX = max(abs(F1errors(:,3)-F1.c1)) / sqrt(2);
F2 = fit(Yvalues', sumY, 'gauss1');
F2errors = confint(F2);
ycentre = F2.b1;
sigmaY = F2.c1 / sqrt(2);
deltasigmaY = max(abs(F2errors(:,3)-F2.c1)) / sqrt(2);
distance2 = (x-xcentre).^2 + (y-ycentre).^2;
mask = distance2 <= (radius_pixels^2);
img = -img;
meanvalue = mean2(img(mask==1));
stdDose = std(img(mask==1));
maxX = F1.a1;
maxY = F2.a1;

end

