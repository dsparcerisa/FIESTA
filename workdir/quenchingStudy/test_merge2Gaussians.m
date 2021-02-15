NXY = 80;

% Función para crear gaussiana
createG = @(center, sigma, Xval, Yval) exp(-((Xval-center(1)).^2 + (Yval-center(2)).^2)./(2*sigma.^2))
Xval = 1:NXY;
Yval = Xval';

% Gaussian 1
center1 = [34.343 45.112];
sigma1 = 20;
S1 = createG(center1, sigma1, Xval, Yval);

% Gaussian 2
center2 = [44.343 37.112];
sigma2 = 23;
S2 = createG(center2, sigma2, Xval, Yval);



% Merge 2 gaussians
mG = merge2Gaussians(S1, S2, center1, center2, sigma1, sigma2);
% Make plot
subplot(2,2,1);
imagesc(S1); hold on;plot(center1(1), center1(2), 'bo');
colorbar
subplot(2,2,2);
imagesc(S2); hold on;plot(center2(1), center2(2), 'bo');
colorbar
subplot(2,1,2);
imagesc(Xval-40.5, Yval-40.5, mG); hold on;plot(0, 0, 'bo');
colorbar
