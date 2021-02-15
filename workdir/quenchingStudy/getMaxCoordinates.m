function [Xcenter, Ycenter] = getMaxCoordinates(I)
figure(77);
[croppedI, rect] = imcrop(I);
Xvalues = 1:size(croppedI, 2);
Yvalues = 1:size(croppedI, 1);
croppedI = sum(croppedI,3);
sumX = sum(croppedI, 1);
sumX = - (sumX - max(sumX));
sumY = sum(croppedI, 2);
sumY = - (sumY - max(sumY));
F1 = fit(Xvalues', sumX', 'gauss1');
Xcenter = F1.b1 + rect(1);

F2 = fit(Yvalues', sumY, 'gauss1');
Ycenter = F2.b1 + rect(2);
close(77);

end

