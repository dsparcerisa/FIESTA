
%% Abrir todos los archivos
clear all

%% Primero lo hacemos con la calibración de FQS
basePath = '/Users/dani/Google Drive/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/21_04_16 Calibracion en medicina'
Nsamples = 12;
dosesGy = [0 0.5 1 1.5 2 3 4 6 8 10 12 15]
filenames = {'0','0_5','1', '1_5', '2', '3', '4', '6', '8', '10', '12', '15'};

% O bien, para la de HF
% basePath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_08_11 - Calibracion en Fuenlabrada'
% Nsamples = 11;
% MU = [0 52 105 157 209 419 629 838 1048 1257 2095];
% dosesGy = MU ./ 104.8;
    
% Valores generales
scansPerSample = 1;
filmsPerSample = 2;
filmOrder = {'EBT3', 'EBT3unl'};

%% Cropping rectangle: strip of 1x3 cm
cropHt_cm = 3;
cropWd_cm = 1;

exampleFilePath = [basePath filesep 'EBT3' filesep '0.tif'];
[pixelsXcm, maxBits] = getImgMetaInfo(exampleFilePath);

cropHt_px = pixelsXcm * cropHt_cm;
cropWd_px = pixelsXcm * cropWd_cm;

%% Load all 
imageSubsets = {}; % {film, doseSample, scan}
for i=1:filmsPerSample
    filmFolder = [basePath filesep filmOrder{i}];
    for j=1:Nsamples
        fileName = filenames{j};
        filePath = [filmFolder filesep fileName '.tif'];
        I = imread(filePath);
        imageSubsets{i, j} = cropCentralStrip(I, cropWd_px, cropHt_px);
    end
end

%% Plot all images
for i=1:filmsPerSample
   for j=1:Nsamples
           subplot(filmsPerSample, Nsamples, j+(i-1)*Nsamples )
           imshow(imageSubsets{i, j});
           title(sprintf('%i %i', i, j));
   end
end

% save('basicData_HF.mat');
save('basicData_Medicina_SCANNED.mat');

%% Study 1. Variability between scans
meanValues = nan(Nsamples, filmsPerSample, 3);
stdValues = nan(Nsamples, filmsPerSample, 3);
for i=1:filmsPerSample
    for j=1:Nsamples
        for l=1:3
            meanValues(i,j,l) = mean2(imageSubsets{i,j}(:,:,l));
            stdValues(i,j,l) = std2(imageSubsets{i,j}(:,:,l));
        end
    end
end



%% Plot results
for i=1:filmsPerSample
    subplot(2,1,i);
    pixelValuesR = meanValues(i,:,1) / maxBits;
    pixelValuesG = meanValues(i,:,2) / maxBits;
    pixelValuesB = meanValues(i,:,3) / maxBits;
    
    errValuesR = stdValues(i,:,1) / maxBits;
    errValuesG = stdValues(i,:,2) / maxBits;
    errValuesB = stdValues(i,:,3) / maxBits;
    
    errorbar(dosesGy, pixelValuesR, errValuesR, 'r.');
    hold on
    errorbar(dosesGy, pixelValuesG, errValuesG, 'g.');
    errorbar(dosesGy, pixelValuesB, errValuesB, 'b.');
    title(filmOrder{i})
    grid on
    ylim([0 0.7]);
    xlim([0 round(max(dosesGy))]);
    xlabel('Dose [Gy]');
    ylabel('Relative pixel value');
end

%% Create single value table
finalPixelValues = nan(filmsPerSample, Nsamples, 3);
finalPixelErrValues = nan(filmsPerSample, Nsamples, 3);
for i=1:filmsPerSample
    for j=1:Nsamples
        for l=1:3
            finalPixelValues(i,j,l) = mean(meanValues(i,j,l)) / maxBits;
            finalPixelErrValues(i,j,l) = rssq(stdValues(i,j,l)) / scansPerSample / maxBits;
        end
    end
end

%% Plot results
plotPoints = 0:0.1:max(dosesGy);
allCoefs = {};
for i=1:filmsPerSample
    subplot(1,filmsPerSample,i);
        
    errorbar(dosesGy, finalPixelValues(i,:,1), finalPixelErrValues(i,:,1), 'r.');
    hold on
    errorbar(dosesGy, finalPixelValues(i,:,2), finalPixelErrValues(i,:,2), 'g.');
    errorbar(dosesGy, finalPixelValues(i,:,3), finalPixelErrValues(i,:,3), 'b.');

    % Find calibration
    FR = fit(dosesGy', finalPixelValues(i,:,1)', 'rat11');
    CoefR = [FR.p2 FR.p1 FR.q1];
    FG = fit(dosesGy', finalPixelValues(i,:,2)', 'rat11');
    CoefG = [FG.p2 FG.p1 FG.q1];
    FB = fit(dosesGy', finalPixelValues(i,:,3)', 'rat11'); 
    CoefB = [FB.p2 FB.p1 FB.q1];
    
    allCoefs{i, 1} = CoefR;
    allCoefs{i, 2} = CoefG;
    allCoefs{i, 3} = CoefB;

    plot(plotPoints, FR(plotPoints), 'r-');
    plot(plotPoints, FG(plotPoints), 'g-');
    plot(plotPoints, FB(plotPoints), 'b-');
    
    title(filmOrder{i})
    grid on
    ylim([0 0.7]);
    xlim([0 round(max(dosesGy))]);
    xlabel('Dose [Gy]');
    ylabel('Relative pixel value');
    
    %saveFileName = ['CoefFitFQS_' filmOrder{i} '.mat'];  
    %saveFileName = ['CoefFitHF_' filmOrder{i} '.mat'];  
    saveFileName = ['CoefFitMedicina_' filmOrder{i} '.mat']; 

    save(saveFileName, 'CoefR', 'CoefG', 'CoefB');
end

%% Test self-consistency (without using delta Method)
autoDoses = nan(filmsPerSample, Nsamples, scansPerSample);
autoDoseErrors = nan(filmsPerSample, Nsamples, scansPerSample);
tic
for i=1:filmsPerSample
    for j=1:Nsamples
        
        % PV(x) = (p1*D + p2) / (D + q1)
        % D(PV) = (q1 * PV - p2) ./ (p1 - PV);
        % COEF= [p2 p1 q1];
        
        p2R = allCoefs{i, 1}(1);
        p1R = allCoefs{i, 1}(2);
        q1R = allCoefs{i, 1}(3);
        
        p2G = allCoefs{i, 2}(1);
        p1G = allCoefs{i, 2}(2);
        q1G = allCoefs{i, 2}(3);
        
        p2B = allCoefs{i, 3}(1);
        p1B = allCoefs{i, 3}(2);
        q1B = allCoefs{i, 3}(3);
        
        doseR = @(PV) (q1R.*PV - p2R) ./ (p1R - PV);
        doseG = @(PV) (q1G.*PV - p2G) ./ (p1G - PV);
        doseB = @(PV) (q1B.*PV - p2B) ./ (p1B - PV);
        
        autoDoses(i,j) = doseR(meanValues(i,j,1) / maxBits);
        %autoDoseErrors(i,j,k) = std2(dose.data);
        plot(dosesGy, autoDoses(i, :), 'r.'); hold on
        
        autoDoses(i,j) = doseG(meanValues(i,j,2) / maxBits);
        %autoDoseErrors(i,j,k) = std2(dose.data);
        plot(dosesGy, autoDoses(i, :), 'g.'); hold on
        
        autoDoses(i,j) = doseB(meanValues(i,j,3) / maxBits);
        %autoDoseErrors(i,j,k) = std2(dose.data);
        plot(dosesGy, autoDoses(i, :), 'b.'); hold on
        
        grid on
    end
end
toc





