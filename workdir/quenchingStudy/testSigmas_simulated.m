clear all

% Problem matrix
dxy = 0.01; % 0.1 mm
sizeXY = 2; % cm
sigmaXY = 0.2; % 2 mm
GP = createGaussProfile(dxy, dxy, sizeXY, sizeXY, sigmaXY, sigmaXY);
noiseLevel = 1/100;
noiseStd = noiseLevel * max(GP.data(:));
noiseMat = noiseStd * randn(size(GP.data));
noisyGP = CartesianGrid2D(GP);
noisyGP.data = GP.data + noiseMat;
% Crop it to uncenter it
%noisyGP.crop([-0.5 0.5 -0.5 0.5])
noisyGP.plotSlice

%% Find plot
figure(3)
noisyGP.plotLinear('X',0)
hold on
noisyGP.plotLinear('Y',0,'r')


%% Use methods to recover the sigma
%[~, xcentre, ycentre, sigmaX1, sigmaY1, ~, ~, ~, deltaSigmaX1, deltaSigmaY1, ~] = meanAndCenterMass(-noisyGP.data,5);

%% GetSigmas
[sigmaX2, sigmaY2, deltaSigmaX2, deltaSigmaY2, xCenter, yCenter] = getSigmas(noisyGP.getAxisValues('X'),noisyGP.getAxisValues('Y'),noisyGP.data,noiseStd*ones(size(noisyGP.data)),4);
100*[sigmaX2 deltaSigmaX2; sigmaY2 deltaSigmaY2]
