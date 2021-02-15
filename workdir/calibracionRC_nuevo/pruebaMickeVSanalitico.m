%% Abrir todos los archivos
clear all; close all

% Choose dataset
filmIndex = 3; % EBT2,2=EBT3,3=EBT3unl
dataIndex = 1; % 1=HF fotones, 2=FQS protones

%% Load data
if dataIndex==1
    load('basicData_HF.mat');
elseif dataIndex==2
    load('basicData_FQS.mat');
else
    error('No data available.');
end

%% Cargar coeficientes necesarios
if dataIndex==1
    load(sprintf('fitCoef_HF_%s.mat',filmOrder{filmIndex}));
elseif dataIndex==2
    load(sprintf('fitCoef_FQS_%s.mat',filmOrder{filmIndex}));
else
    error('Cannot load fit data.');
end


%% Estudio sobre las imágenes completas (iterativa)
autoDoseValues = nan(Nsamples,1);
autoDoseErrValues = nan(Nsamples,1);
autoDoseValues2 = nan(Nsamples,1);
autoDoseErrValues2 = nan(Nsamples,1);
fittedDoseValues = nan(Nsamples,1);
fittedDoseErrValues = nan(Nsamples,1);

for j=1:Nsamples    
    I1a = imageSubsets{filmIndex,j,1};
    I1b = fliplr(I1a);
    I1c = flipud(I1a);
    I1d = fliplr(I1c);
    
    I2a = imageSubsets{filmIndex,j,2};
    I2b = fliplr(I2a);
    I2c = flipud(I2a);
    I2d = fliplr(I2c);
    
    I3a = imageSubsets{filmIndex,j,3};
    I3b = fliplr(I3a);
    I3c = flipud(I3a);
    I3d = fliplr(I3c);
    
    I = (double(I1a) + double(I1b) + double(I1c) + double(I1c) ...
      + double(I2a) + double(I2b) + double(I2c) + double(I2c) ...
      + double(I3a) + double(I3b) + double(I3c) + double(I3d)) / 12;
    
    pxmax = double(intmax(class(I1a)));
    deltasMicke = -0.5:0.001:0.5;
    [dose, varMat] = getDoseMicke(I, CoefR1, CoefG1, CoefB1, pixelsXcm, deltasMicke, pxmax);
    autoDoseValues(j) = mean2(dose.data);
    autoDoseErrValues(j) = std2(dose.data);
    
    [dose2, varMat2] = getDoseMicke_mult(I, CoefR1, CoefG1, CoefB1, pixelsXcm, deltasMicke, pxmax);
    autoDoseValues2(j) = mean2(dose2.data);
    autoDoseErrValues2(j) = std2(dose2.data);    
    
    % Calcular dosis directamente con el fit
    R = I(:,:,1) / pxmax;
    G = I(:,:,2) / pxmax;
    B = I(:,:,3) / pxmax;    
    
    stackedDoses = cat(3, getChannelDoseT1(CoefR1, R), getChannelDoseT1(CoefG1, G), getChannelDoseT1(CoefB1, B));
    fittedDose = mean(stackedDoses, 3);
    fittedDoseValues(j) = mean2(fittedDose);
    fittedDoseErrValues(j) = std2(fittedDose);
    
% 	figure(10+j);

%     subplot(1,2,1);
%     dose.plotSlice
%     title(j);
%     caxis([0 round(max(dosesGy))]);
%     colorbar
%     subplot(1,2,2);
%     varMat.plotSlice  
%     caxis([min(deltasMicke) max(deltasMicke)]);
%     colorbar
%    fprintf('%3.3f +- %3.3f\n',autoDoseValues(j),autoDoseErrValues(j));
    
end

%% PLOTTEAR
subplot(3,1,1)
errorbar(dosesGy, autoDoseValues, autoDoseErrValues, 'b.'); hold on
errorbar(dosesGy, fittedDoseValues, fittedDoseErrValues, 'r.');
errorbar(dosesGy, autoDoseValues2, autoDoseErrValues2, 'g.'); hold on
plot(dosesGy, dosesGy, 'k:');
legend({'Micke aditivo','Average of 3 components', 'Micke multiplicativo'},'Location','SouthEast');
xlabel('Nominal doses (Gy)');
ylabel('Measured doses (Gy)');
if dataIndex==1
    title(sprintf('Photons (%s)',filmOrder{filmIndex}));
elseif dataIndex==2
    title(sprintf('Protons (%s)',filmOrder{filmIndex}));
end
grid on

subplot(3,1,2)
errorbar(dosesGy, (autoDoseValues-dosesGy'), autoDoseErrValues, 'b.'); hold on
errorbar(dosesGy, (fittedDoseValues-dosesGy'), fittedDoseErrValues, 'r.');
errorbar(dosesGy, (autoDoseValues2-dosesGy'), autoDoseErrValues2, 'g.'); hold on

axis([0 12 -0.5 0.5]);
xlabel('Nominal doses (Gy)');
ylabel('Absolute Error (Gy)');
grid on

subplot(3,1,3)
errorbar(dosesGy, 100*(autoDoseValues-dosesGy')./(dosesGy' + eps), 100*autoDoseErrValues./(dosesGy' + eps), 'b.'); hold on
errorbar(dosesGy, 100*(fittedDoseValues-dosesGy')./(dosesGy' + eps), 100*fittedDoseErrValues./(dosesGy' + eps), 'r.');
errorbar(dosesGy, 100*(autoDoseValues2-dosesGy')./(dosesGy' + eps), 100*autoDoseErrValues2./(dosesGy' + eps), 'g.'); hold on

axis([0 12 -10 10]);
xlabel('Nominal doses (Gy)');
ylabel('Relative Error (%)');
grid on