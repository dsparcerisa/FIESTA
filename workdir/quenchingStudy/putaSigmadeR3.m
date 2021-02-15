close all; clear all

%% Load
filmType = 'EBT2';
E0 = 4; % MeV
radName = 'Rad5b';

calFile = ['/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo/fitCoef_FQS_' filmType '.mat'];
load(calFile);
addpath '/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo'



basePath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_07_08 - Estudio quenching CMAM/scan-dani-6sept';
fullPath = [basePath filesep filmType filesep sprintf('%i',E0) 'MeV'];
imgPath1 = [fullPath filesep radName '.tif'];
imgPath2 = [fullPath filesep radName ' (2).tif'];
imgPath3 = [fullPath filesep radName ' (3).tif'];



[pixCM, maxBits] = getImgMetaInfo(imgPath1);

I1 = imread(imgPath1);
I2 = imread(imgPath2);
I3 = imread(imgPath3);

% Merge all three scans of the image
I = (1/3) * (double(I1) + double(I2) + double(I3));
I = uint16(I);
warning('Select last point of highest dose')
I = imcrop(I);
deltas = -0.5:0.001:0.5;
d1 = getDoseMicke(double(I), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);

%% Probar con la suma
sumX = mean(d1.data,2);
baseLineX = min(sumX);
sumX = sumX - baseLineX;
X = d1.getAxisValues('X');
subplot(1,2,1);
plot(X, sumX, 'b.')
hold on
title('X');
    [FX, gofx] = fit(X', sumX, 'gauss2', 'Lower', [0 -Inf 0 0 -Inf 0])
    %[FX, gofx] = fit(X', sumX, 'gauss1', 'Lower', [0 -Inf 0])

    FXerrors = confint(FX);

    if FX.a1 > FX.a2
        SX = 10*FX.c1/sqrt(2);
        dSX = 10*max(abs(FXerrors(:,3)-FX.c1)) / sqrt(2);
    else
        SX = 10*FX.c2/sqrt(2);
        dSX = 10*max(abs(FXerrors(:,6)-FX.c2)) / sqrt(2);
    end
    plot(FX)

    sumY = mean(d1.data);
    baseLineY = min(sumY);
    sumY = sumY - baseLineY;
    Y = d1.getAxisValues('Y');
    subplot(1,2,2);
    plot(Y, sumY, 'b.')
    hold on
    title('Y');
    [FY, gofy] = fit(Y', sumY', 'gauss2', 'Lower', [0 -Inf 0 0 -Inf 0])
    FYerrors = confint(FY);

    if FY.a1 > FY.a2
        SY = 10*FY.c1/sqrt(2);
        dSY = 10*max(abs(FYerrors(:,3)-FY.c1)) / sqrt(2);
    else
        SY = 10*FY.c2/sqrt(2);
        dSY = 10*max(abs(FYerrors(:,6)-FY.c2)) / sqrt(2);
    end
    plot(FY)

fprintf('\tSigma\tdSigma\tR2\n');
fprintf('Xlo\t%3.4f\t%3.4f\t%3.4f\n', SX, dSX, gofx.rsquare);
fprintf('Ylo\t%3.4f\t%3.4f\t%3.4f\n', SY, dSY, gofy.rsquare);