% Ejemplo de buscar ángulo con datos simulados.
clear all

%% Load Rad16
addpath '/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo'
basePath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_07_08 - Estudio quenching CMAM/scan-dani-6sept';
fullPath = [basePath filesep 'EBT3unl' filesep '10MeV'];
radName = 'Rad16';
imgPath1 = [fullPath filesep radName '.tif'];
imgPath2 = [fullPath filesep radName ' (2).tif'];
imgPath3 = [fullPath filesep radName ' (3).tif'];

I1 = imread(imgPath1);
I2 = imread(imgPath2);
I3 = imread(imgPath3);

% Merge all three scans of the image

I = (1/3) * (double(I1) + double(I2) + double(I3));
I = uint16(I);
warning('Select comet area');
I = imcrop(I);
close all

[pixCM, maxBits] = getImgMetaInfo(imgPath1);
deltas = -0.5:0.001:0.5;
load('/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo/fitCoef_HF_EBT3unl.mat')

[dose, varMat] = getDoseMicke(double(I), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);

%% Función para crear gaussiana
createG = @(center, sigma, Xval, Yval) exp(-((Xval-center(1)).^2 + (Yval-center(2)).^2)./(2*sigma.^2))
Xval_mm = -10:0.1:10;
Yval_mm = Xval_mm';

% Probaremos con Rad16 (EBT3unl a 10MeV)

deviation_deg = 16;

sigmas = [0.3727
0.5067
0.6791
0.815
1.0629];

sigmas = sigmas * 0.9

% Gaussian 1
Zvalues = (0:2:8) - 1;
S = {};
for i=1:numel(Zvalues)
    S{i} = createG([Zvalues(i)*tand(deviation_deg) 0], sigmas(i), Xval_mm, Yval_mm);
    subplot(3,2,i);
    imagesc(Xval_mm, Yval_mm, S{i});
    grid on
    hold on
    plot(0,Zvalues(i)*tand(deviation_deg), 'bo');
    if i==1
        suma = S{i};
    else
        suma = suma + S{i};
    end
    axis([-4 4 -4 4]);
    caxis([0 3]);
    colorbar
end
subplot(3,2,6)
imagesc(Xval_mm, Yval_mm, suma)
axis([-4 4 -4 4]);
caxis([0 3]);

grid on
colorbar

%% Compare the two
subplot(2,1,1);
dose.plotSlice
subplot(2,1,2);
imagesc(Xval_mm, Yval_mm, suma)

%% Prueba de regionprops
testProps = regionprops(dose.data,'Centroid','Orientation');
centroids = cat(1,testProps.Centroid);
imagesc(dose.data)
hold on
for i=1:size(centroids,1)
    plot(centroids(i,1),centroids(i,2),'bo')
    % pause
end
% --> El último valor de centroids da el punto central
% --> con un fit podemos quizá sacar el angulo phi y comparar con
% orientations (media de los dos primeros).
F = fit(centroids(:,1), centroids(:,2), 'poly1');
plot(F);
theAngle = atand(-F.p1)
allOrientations = [testProps.Orientation];
theAngle2 = mean(allOrientations(1:3))

%% Complete image rotation (does it preserve distances?)
rotDose = imrotate(dose.data, -theAngle, 'bicubic', 'crop');
prop2 = regionprops(rotDose,'Centroid');
centroids2 = cat(1,prop2.Centroid);
theCentroid2 = centroids2(end-1,:)
rotDose_Xval = ((1:size(rotDose,2)) - theCentroid2(1))/pixCM*10
rotDose_Yval = ((1:size(rotDose,1)) - theCentroid2(2))/pixCM*10
imagesc(rotDose_Xval, rotDose_Yval, rotDose);
grid on

%% Plot both and compare
close all
subplot(2,1,1);
imagesc(rotDose_Xval, rotDose_Yval, rotDose);
axis([-2 5 -3 3])
caxis([ 0 10]);
colorbar
title('Measured dose');
subplot(2,1,2);
fac = 0.9*max(rotDose(:)) / max(suma(:));
imagesc(Xval_mm, Yval_mm, fac*suma)
colorbar
axis([-2 5 -3 3])
caxis([ 0 10]);
title('Estimated');