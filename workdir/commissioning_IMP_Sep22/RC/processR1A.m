clear all
addpath '/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo'
addpath '/Users/dani/Documents/FIESTA/workdir/quenchingStudy'
basePath = '/Users/dani/Documents/FIESTA/workdir/commissioning_IMP_Sep22/RC';
fullPath = [basePath filesep 'todasTC_conCristal.tif'];

I  = imread(fullPath);
I = uint16(I);

[pixCM, maxBits] = getImgMetaInfo(fullPath);

%% Apply general cropping
warning('Apply general cropping');

I = imcrop(I);
close all
imshow(I);

%% Base processing data
deltas = -0.005:0.001:0.005;
load('/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo/fitCoef_HF_EBT3unl.mat')
radius_max_cm = 0.03;
radius_max_pixels = radius_max_cm*pixCM;
radius_pocillo_cm = 0.2;
radius_pocillo_pixels = radius_pocillo_cm*pixCM;

%% Primera parte: estabilidad del halo de 10s entre CDEFG de R1

%% Number of items

Npoints = 5; 
allI = {};

close all;

% Process all elements in allI
meanDoses = nan(Npoints,1);
stdDoses = nan(Npoints,1);

for i=1:Npoints
    warning('Select subimage %i / %i', i, Npoints);
    figure(1);
    allI{i} = imcrop(I);
    close(1)
    
    %figure(2)
    [dose, varMat, dr, dg, db] = getDoseMicke(double(allI{i}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
    
    % Corrección de valores negativos con canal azul
    negValuesMask = dose.data<0 | dose.data>30;
    dose.data(negValuesMask) = db(negValuesMask);
       
    meanDoses(i) = mean(dose.data, 'all');
    stdDoses(i) = std(dose.data(:));
    fprintf('Mean dose (all): %3.3f +- %3.3f\n', meanDoses(i), stdDoses(i));
    
end

