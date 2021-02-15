dxy = 0.1;
gaussXY = @(Xval, Yval, Xc, Yc, sigmaX, sigmaY) exp(-((Xval-Xc)./(sqrt(2)*sigmaX)).^2) * exp(-((Yval-Yc)./(sqrt(2)*sigmaY)).^2);
Xval = (-3:dxy:3)';
Yval = Xval';
A = gaussXY(Xval, Yval, -1.2, -1.2, 0.78, 0.89);
A = A + gaussXY(Xval, Yval, -1.2, 0, 0.78, 0.89);
A = A + gaussXY(Xval, Yval, -1.2, 1.2, 0.78, 0.89);
A = A + gaussXY(Xval, Yval, 0, -1.2, 0.78, 0.89);
A = A + gaussXY(Xval, Yval, 0, 0, 0.78, 0.89);
A = A + gaussXY(Xval, Yval, 0, 1.2, 0.78, 0.89);
A = A + gaussXY(Xval, Yval, 1.2, -1.2, 0.78, 0.89);
A = A + gaussXY(Xval, Yval, 1.2, 0, 0.78, 0.89);
A = A + gaussXY(Xval, Yval, 1.2, 1.2, 0.78, 0.89);

rVal = 2; % mm
inMask = (Xval.^2 + Yval.^2) <= rVal;

subplot(2,1,1);
imagesc(Xval, Yval, A')
xlabel('X'); ylabel('Y');

subplot(2,1,2);
imagesc(Xval, Yval, inMask')
xlabel('X'); ylabel('Y');

meanDose = mean(A(inMask==1));
pctStdDose = 100 * std(A(inMask==1)) / meanDose 