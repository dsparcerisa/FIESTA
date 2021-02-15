function [mG, XvalG, YvalG] = merge2Gaussians(S1, S2, center1, center2, sigma1, sigma2, sigFactor)

if nargin<7 || ~exist('sigFactor')
    sigFactor = 2;
end

% Centrar en X e Y e interpolar
size1X = size(S1,1);
size1Y = size(S1,2);
size2X = size(S2,1);
size2Y = size(S2,2);

Xval1 = (1:size1X) - center1(1);
Yval1 = (1:size1Y)' - center1(2);
Xval2 = (1:size2X) - center2(1);
Yval2 = (1:size2Y)' - center2(2);

sigmaG = max(sigma1, sigma2);
maxR = sigFactor*sigmaG;

sizeGX = max(size1X,size2X);
sizeGY = max(size1Y,size2Y);

XvalG = -((sizeGX-1)/2):((sizeGX-1)/2);
YvalG = (-((sizeGY-1)/2):((sizeGY-1)/2))';

S1i = interp2(Xval1,Yval1,S1,XvalG,YvalG,'cubic');
S2i = interp2(Xval2,Yval2,S2,XvalG,YvalG,'cubic');

% Delete extra-radius points
distance = sqrt(XvalG.^2 + YvalG.^2);
outerMask = distance > maxR;
S1i(outerMask) = 0;
S2i(outerMask) = 0;

% Merge using MAX function
mG = max(S1i, S2i);

% Make plot
figure(7)
subplot(2,2,1); 
imagesc(S1i); hold on;plot(0, 0, 'bo');
colorbar
subplot(2,2,2);
imagesc(S2i); hold on;plot(0, 0, 'bo');
colorbar
subplot(2,2,3);
hold off
plot(sum(S1i,'omitnan'),'r');
hold on
plot(sum(S2i,'omitnan'),'g');
plot(sum(mG,'omitnan'),'b:');
subplot(2,2,4);
hold off
plot(sum(S1i,2,'omitnan'),'r');
hold on
plot(sum(S2i,2,'omitnan'),'g');
plot(sum(mG,2,'omitnan'),'b:');

end

