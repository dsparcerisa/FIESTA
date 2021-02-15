%% Abrir todos los archivos
clear all
basePath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_07_23 - Calibracion en FQS'
Nsamples = 10;
scansPerSample = 3;
filmsPerSample = 3;
filmOrder = {'EBT2', 'EBT3', 'EBT3-unl'};

%% Find all cropping rectangles
croppingRectangles = {};
for i = 1:Nsamples
    
    fileName = sprintf('Sample%i-scan1.tif', i-1);
    filePath = [basePath filesep fileName];
    
    for j = 1:filmsPerSample
        msg = sprintf('Draw rectangle for %s film', filmOrder{j});        
        croppingRectangles{i, j} = getCroppingRect(filePath, msg);
    end
end

save('croppingRectangles.mat');

%% Autocrop
imageSubsets = {}; % {sample, film, scan}
for i=1:Nsamples
   for j=1:filmsPerSample
       for k=1:scansPerSample
           fileName = sprintf('Sample%i-scan%i.tif', i-1, k);
           filePath = [basePath filesep fileName];
           I = imread(filePath);
           imageSubsets{i, j, k} = imcrop(I, croppingRectangles{i,j});
       end
   end
end

%% Plot all images
for i=1:Nsamples
   for j=1:filmsPerSample
       for k=1:scansPerSample
           subplot(filmsPerSample*scansPerSample, Nsamples, sub2ind([10 3 3],i,j,k))
           imshow(imageSubsets{i, j, k});
           title(sprintf('%i %i %i', i, j, k));
       end
   end
end

%% Study 1. Variability between scans
meanValues = nan(Nsamples, filmsPerSample, scansPerSample, 3);
stdValues = nan(Nsamples, filmsPerSample, scansPerSample, 3);
for i=1:Nsamples
   for j=1:filmsPerSample
       for k=1:scansPerSample
            for l=1:3
                meanValues(i,j,k,l) = mean2(imageSubsets{i,j,k}(:,:,l));
                stdValues(i,j,k,l) = std2(imageSubsets{i,j,k}(:,:,l));
            end
       end
   end
end

%% Plot each channel for each type of film:
MU = [0 443 887 1330 1774 3548 5322 7095 8869 10643];
% 1681 MU = 1.8953 Gy
dosesGy = MU*1.8953/1681;
[pixelsXcm, maxBits] = getImgMetaInfo(filePath);

%% Plot results
for j=1:filmsPerSample
    subplot(1,3,j);
    for k=1:scansPerSample
        pixelValuesR = meanValues(:,j,k,1) / maxBits;
        pixelValuesG = meanValues(:,j,k,2) / maxBits;
        pixelValuesB = meanValues(:,j,k,3) / maxBits;

        errValuesR = stdValues(:,j,k,1) / maxBits;
        errValuesG = stdValues(:,j,k,2) / maxBits;
        errValuesB = stdValues(:,j,k,3) / maxBits;        
            
        errorbar(dosesGy, pixelValuesR, errValuesR, 'r.');
        hold on
        errorbar(dosesGy, pixelValuesG, errValuesG, 'g.');
        errorbar(dosesGy, pixelValuesB, errValuesB, 'b.');
    end
    title(filmOrder{j})
    grid on
    ylim([0 0.7]);
    xlim([0 12]);
    xlabel('Dose [Gy]');
    ylabel('Relative pixel value');
end

%% Create single value table
finalPixelValues = nan(Nsamples, filmsPerSample, 3);
finalPixelErrValues = nan(Nsamples, filmsPerSample, 3);
for i=1:Nsamples
    for j=1:filmsPerSample
        for l=1:3
            finalPixelValues(i,j,l) = mean(meanValues(i,j,:,l)) / maxBits;
            finalPixelErrValues(i,j,l) = rssq(stdValues(i,j,:,l)) / scansPerSample / maxBits;
        end
    end
end

%% Plot results
plotPoints = 0:0.1:12;
allCoefs = {};
for j=1:filmsPerSample
    subplot(1,filmsPerSample,j);
        
    errorbar(dosesGy, finalPixelValues(:,j,1), finalPixelErrValues(:,j,1), 'r.');
    hold on
    errorbar(dosesGy, finalPixelValues(:,j,2), finalPixelErrValues(:,j,2), 'g.');
    errorbar(dosesGy, finalPixelValues(:,j,3), finalPixelErrValues(:,j,3), 'b.');

    % Find calibration
    FR = fit(dosesGy', finalPixelValues(:,j,1), 'rat11');
    CoefR = [FR.p2 FR.p1 FR.q1];
    FG = fit(dosesGy', finalPixelValues(:,j,2), 'rat11');
    CoefG = [FG.p2 FG.p1 FG.q1];
    FB = fit(dosesGy', finalPixelValues(:,j,3), 'rat11'); %% Probar con otra función para el azul!
    CoefB = [FB.p2 FB.p1 FB.q1];
    
    allCoefs{j, 1} = CoefR;
    allCoefs{j, 2} = CoefG;
    allCoefs{j, 3} = CoefB;

    plot(plotPoints, FR(plotPoints), 'r-');
    plot(plotPoints, FG(plotPoints), 'g-');
    plot(plotPoints, FB(plotPoints), 'b-');
    
    title(filmOrder{j})
    grid on
    ylim([0 0.7]);
    xlim([0 12]);
    xlabel('Dose [Gy]');
    ylabel('Relative pixel value');
    
    saveFileName = ['CoefFitFQS_' filmOrder{j} '.mat'];    
    save(saveFileName, 'CoefR', 'CoefG', 'CoefB');
end

%% Test self-consistency (without using delta Method)
autoDoses = nan(Nsamples, filmsPerSample, scansPerSample);
autoDoseErrors = nan(Nsamples, filmsPerSample, scansPerSample);
tic
for i=1:Nsamples
    for j=1:filmsPerSample
        
        % PV(x) = (p1*D + p2) / (D + q1)
        % D(PV) = (q1 * PV - p2) ./ (p1 - PV);
        % COEF= [p2 p1 q1];
        
        p2R = allCoefs{j, 1}(1);
        p1R = allCoefs{j, 1}(2);
        q1R = allCoefs{j, 1}(3);
        
        p2G = allCoefs{j, 2}(1);
        p1G = allCoefs{j, 2}(2);
        q1G = allCoefs{j, 2}(3);
        
        p2B = allCoefs{j, 3}(1);
        p1B = allCoefs{j, 3}(2);
        q1B = allCoefs{j, 3}(3);
        
        doseR = @(PV) (q1R.*PV - p2R) ./ (p1R - PV);
        doseG = @(PV) (q1G.*PV - p2G) ./ (p1G - PV);
        doseB = @(PV) (q1B.*PV - p2B) ./ (p1B - PV);

        for k=1:scansPerSample
            autoDoses(i,j,k) = doseR(meanValues(i,j,k,1) / maxBits);
            %autoDoseErrors(i,j,k) = std2(dose.data);
            plot(dosesGy, autoDoses(:, j, k), 'r.'); hold on
            
            autoDoses(i,j,k) = doseG(meanValues(i,j,k,2) / maxBits);
            %autoDoseErrors(i,j,k) = std2(dose.data);
            plot(dosesGy, autoDoses(:, j, k), 'g.'); hold on
            
            autoDoses(i,j,k) = doseB(meanValues(i,j,k,3) / maxBits);
            %autoDoseErrors(i,j,k) = std2(dose.data);
            plot(dosesGy, autoDoses(:, j, k), 'b.'); hold on
            
        end
    end
    i
end
toc


%% Test self-consistency
autoDoses = nan(Nsamples, filmsPerSample, scansPerSample);
autoDoseErrors = nan(Nsamples, filmsPerSample, scansPerSample);
deltas = 0.8:0.005:1.2;
tic
parfor i=1:Nsamples
    for j=1:filmsPerSample
        for k=1:scansPerSample
            dose = getDoseFromRCv2(imageSubsets{i,j,k}, allCoefs{j, 1}, allCoefs{j, 2}, pixelsXcm, deltas);
            autoDoses(i,j,k) = mean2(dose.data);
            autoDoseErrors(i,j,k) = std2(dose.data);
        end
    end
    i
end
toc

%% Plot result of self-consistency test
figure(10);    
for k=1:scansPerSample
    errorbar(dosesGy, autoDoses(:,1,k), autoDoseErrors(:,1,k), 'bo'); hold on
    errorbar(dosesGy, autoDoses(:,2,k), autoDoseErrors(:,2,k), 'ro'); hold on
    errorbar(dosesGy, autoDoses(:,3,k), autoDoseErrors(:,3,k), 'go'); hold on
end
grid on
ylim([0 12]);
xlim([0 12]);
xlabel('Dose [Gy]');
ylabel('Measured dose [Gy]');
legend(filmOrder);
title('Self-consistency test');





