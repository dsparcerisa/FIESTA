clear all; close all
addpath '/Users/dani/FIESTA/workdir/calibracionRC_nuevo'
addpath '/Users/dani/FIESTA/workdir/quenchingStudy'
basePath = '/Users/dani/OneDrive - Universidad Complutense de Madrid (UCM)/UCM/01-Proyectos/40-TFG Laura/Escaneos/22_12_22';
fullPath = [basePath filesep '1Escaner2.tif'];
I  = imread(fullPath);
I = uint16(I);

[pixCM, maxBits] = getImgMetaInfo(fullPath);

%% Apply general cropping
warning('Apply general cropping');

I = imcrop(I);
close all
imshow(I);

%% Base processing data
deltas = -0.02:0.001:0.02;
load('/Users/dani/FIESTA/workdir/calibracionRC_nuevo/fitCoef_FQS_EBT3.mat')

[dose, varMat, dr, dg, db] = getDoseMicke(double(I), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);

% Corrección de valores n>10 con canal rojo
%corrValuesMask = dose.data>1;
%dose.data(corrValuesMask) = dr(corrValuesMask);

figure(1)
subplot(2,2,1);
dose.plotSlice;
title('Dose');
colorbar
subplot(2,2,2);
varMat.plotSlice;
colorbar
title('Variations');
subplot(2,2,3);
plot(dose.getAxisValues('Y'),mean(dose.data),'b')
grid on
title('Sum in X');
subplot(2,2,4);
plot(dose.getAxisValues('X'),mean(dose.data, 2),'b')
grid on
title('Sum in Y');


figure(2)
subplot(2,2,1);
imagesc(dose.getAxisValues('X'), dose.getAxisValues('Y'), dr')
caxis([-1 12]);
title('Red Channel');
colorbar
subplot(2,2,2);
imagesc(dose.getAxisValues('X'), dose.getAxisValues('Y'), dg')
caxis([-1 12]);
title('Green Channel');
colorbar
subplot(2,2,3);
imagesc(dose.getAxisValues('X'), dose.getAxisValues('Y'), db')
caxis([-1 12]);
title('Blue Channel');
colorbar
subplot(2,2,4);
histogram(dose.data(:),'EdgeColor',[0 0 0]);
hold on
histogram(dr(:),'EdgeColor',[1 0 0]);
histogram(dg(:),'EdgeColor',[0 1 0]);
histogram(db(:),'EdgeColor',[0 0 1]);
histogram(dose.data(:),'EdgeColor',[0 0 0]);
title('Histogram');

mean(dose.data(:))
std(dose.data(:))