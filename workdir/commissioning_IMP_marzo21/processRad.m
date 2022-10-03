clear all;  close all

addpath '/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo'
addpath '/Users/dani/Documents/FIESTA/workdir/quenchingStudy'
basePath = '/Users/dani/Library/CloudStorage/OneDrive-UniversidadComplutensedeMadrid(UCM)/UCM/01-Proyectos/37-Commissioning linea IMP - Nueva carpeta/3. Junio 21/Radiocromicas';
fullPath = [basePath filesep 'R123.tif'];
I  = imread(fullPath);
I = uint16(I);

[pixCM, maxBits] = getImgMetaInfo(fullPath);

%% Apply general cropping
warning('Apply general cropping');

I = imcrop(I);
close all
imshow(I);

%% Base processing data
deltas = -0.000001:0.001:0.000001;
load('/Users/dani/Documents/FIESTA/workdir/calibracionRC_nuevo/fitCoef_FQS_EBT3.mat')
radius_max_cm = 0.03;
radius_max_pixels = radius_max_cm*pixCM;
radius_pocillo_cm = 0.2;
radius_pocillo_pixels = radius_pocillo_cm*pixCM;

%% Number of items

Npoints = input('Select number of subimages: ')
allI = {};

close all;

% Process all elements in allI
sigmaX_mm = nan(Npoints,1);
sigmaY_mm = nan(Npoints,1);
deltasigmaX_mm = nan(Npoints,1);
deltasigmaY_mm = nan(Npoints,1);
meanDoses_max = nan(Npoints,1);
stdDoses_max = nan(Npoints,1);
meanDoses_4mm = nan(Npoints,1);
stdDoses_4mm = nan(Npoints,1);
xcenters = nan(Npoints,1);
ycenters = nan(Npoints,1);
FWHMx_mm = nan(Npoints,1);
FWHMy_mm = nan(Npoints,1);

for i=1:Npoints
    warning('Select subimage %i / %i', i, Npoints);
    figure(1);
    allI{i} = imcrop(I);
    close(1)
    
    figure(2)
    [dose, varMat, dr, dg, db] = getDoseMicke(double(allI{i}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
    
    % Corrección de valores negativos con canal azul
    negValuesMask = dose.data<0 | dose.data>30;
    dose.data(negValuesMask) = db(negValuesMask);
    
    try
        [~, ~, ~, ~, ~, meanvalueMax, ~, ~, ~, ~, stdDoseMax] = meanAndCenterMass(-dose.data,radius_max_pixels);       
        [~, ~, ~, ~, ~, meanvalue, ~, ~, ~, ~, stdDose] = meanAndCenterMass(-dose.data,radius_pocillo_pixels);
        [sigmaX, sigmaY, deltasigmaX, deltasigmaY, xCenter, yCenter, FWHMx, FWHMy] = getSigmas2(10*dose.getAxisValues('X'),10*dose.getAxisValues('Y'),dose.data,0.05*dose.data+0.01);  
    catch
        warning('Could not process meanAndCenterMass, using nans');
        mask = nan;
        xCenter = nan;
        yCenter  = nan;
        sigmaX  = nan;
        sigmaY = nan;
        meanvalue = nan;
        stdDose = nan;
        meanvalueMax = nan;
        stdDoseMax = nan;        
        deltasigmaX = nan;
        deltasigmaY = nan;
        FWHMx = nan;
        FWHMy = nan;

    end    
    
    sigmaX_mm(i) = sigmaX;
    sigmaY_mm(i) = sigmaY;
    deltasigmaX_mm(i) = deltasigmaX;
    deltasigmaY_mm(i) = deltasigmaY;
    meanDoses_max(i) = meanvalueMax;
    stdDoses_max(i) = stdDoseMax;    
    meanDoses_4mm(i) = meanvalue;
    stdDoses_4mm(i) = stdDose;
    xcenters(i) = xCenter; %xcentre;
    ycenters(i) = yCenter; %ycentre;    
    FWHMx_mm(i) = FWHMx;
    FWHMy_mm(i) = FWHMy;    
    
    fprintf('SigmaX (mm): %3.3f +- %3.3f\n', sigmaX, deltasigmaX);
    fprintf('SigmaY (mm): %3.3f +- %3.3f\n', sigmaY, deltasigmaY);
    fprintf('FWHMx (mm): %3.3f (%3.2f sigmaX)\n', FWHMx, FWHMx/sigmaX);
    fprintf('FWHMy (mm): %3.3f (%3.2f sigmaY)\n', FWHMy, FWHMy/sigmaY);
    fprintf('Mean dose (top): %3.3f +- %3.3f\n', meanvalueMax, stdDoseMax);
    fprintf('Mean dose (r2mm): %3.3f +- %3.3f\n\n', meanvalue, stdDose);
    
    figure(3);
    imshow(allI{i});
    title(i)
    
    figure(4);
    subplot(1,2,1);
    dose.plotSlice;
    hold on
    %midVal = meanvalueMax*0.5;
    %contour(dose.data,[0 midVal],'EdgeColor','k')
    title('Dose');
    colorbar
    subplot(1,2,2);
    varMat.plotSlice;
    title('Multichannel correction');
    colorbar
    pause
end

finalSigmas = 0.5*(sigmaX_mm + sigmaY_mm);


%% Surf plot
figure(6);
subplot(1,2,1);
surf(dose.getAxisValues('X'),dose.getAxisValues('Y'),dose.data')
xlabel('X');
ylabel('Y');
subplot(1,2,2);
imagesc(dose.getAxisValues('X'),dose.getAxisValues('Y'),dose.data')
colorbar
xlabel('X');
ylabel('Y');
set(gca,'ColorScale','log')

%% Get sigmas
YNum= dose.getAxisValues('Y');
YVal=sum(dose.data);
figure(7)
plot(YNum,YVal)