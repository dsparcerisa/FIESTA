clear all; close all;

filmType = 'EBT2';
%% Load calibration file
addpath '/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo'
calFile = ['/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo/fitCoef_FQS_' filmType '.mat'];
load(calFile);
deltas = -0.5:0.001:0.5;


basePath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_12_01 Samuel y Ana Espinosa';
filmName = 'samu2';
fullPath = [basePath filesep filmName];
imgPath1 = [fullPath '_1.tif'];
imgPath2 = [fullPath '_2.tif'];
imgPath3 = [fullPath '_3.tif'];

[pixCM, maxBits] = getImgMetaInfo(imgPath1);

I1 = imread(imgPath1);
I2 = imread(imgPath2);
I3 = imread(imgPath3);

% Merge all three scans of the image
I = (1/3) * (double(I1) + double(I2) + double(I3));
I = uint16(I);

%% Apply general cropping
warning('Apply general cropping');

I = imcrop(I);
close all
imshow(I);

%% Load
d1 = getDoseMicke(double(I), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);

%% Apply RE correction
E0 = 10;
% Load LET values
LETpath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/25-TOPAS/LETScorer2/results';
LETFileName = sprintf('LETvalues_%s_%iMeV',filmType,E0);
fullLETPath = [LETpath filesep LETFileName];
load(fullLETPath);
Zval = 5.3 + (0:10);
LET_16cm = interp1(Zval, meanLET, 33.0, 'linear', 'extrap')

% Load Quenching function
load('/Users/dani/Documents/FIESTA/workdir/quenchingStudy/finalResults_EBT2.mat', 'F2');

d1_corrected = CartesianGrid2D(d1);
d1_corrected.data = d1.data ./ F2(LET_16cm);

%% Save
save('samu2.mat','d1','d1_corrected');