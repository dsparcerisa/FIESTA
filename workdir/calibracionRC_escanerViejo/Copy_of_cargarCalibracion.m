% 22 enero 2020
clear all

% Info for processing
pixelsXCM = 157.48;

% File paths
fileCal1 = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/20-FLASH/radiocromicas/scans/cal1.tif';
fileCal2 = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/20-FLASH/radiocromicas/scans/cal2.tif';
N = 8;
N0 = 3; % añadir 3 puntos más con dosis 0, uno en cada imagen

%% 0. Leer y croppear manualmente
cropRects = {};
allI = {};
I1=imread(fileCal1);
% Leer las del file1
for i=1:N
    figure(1);
    [I, rect] = imcrop(I1);
    close(1);
    cropRects{i} = rect;
    allI{i} = I;
end

I2=imread(fileCal2);
% Leer las del file1
for i=1:N
    figure(1);
    [I, rect] = imcrop(I2);
    close(1);
    cropRects{N+i} = rect;
    allI{N+i} = I;
end

for i=1:N0
    figure(1);
    [I, rect] = imcrop(I1);
    cropRects{2*N + i} = rect;
    allI{2*N+i} = I;    
    close(1);    
end

for i=1:N0
    figure(1);
    [I, rect] = imcrop(I2);
    cropRects{2*N + N0 + i} = rect;
    allI{2*N+ N0 + i} = I;    
    close(1);    
end

figure(30);
for i=1:numel(allI); subplot(6,4,i); imshow(allI{i}); end

save('Step0.mat');

%% 1. Hallar sigmaX y sigmaY de todos
clear all
load('Step0.mat');

sigmasX = nan(2*N, 1);
sigmasY = nan(2*N, 1);

% Hacerlo con el modelo inicial:
for i=1:2*N
    I = allI{i};
    IR = double(I(:,:,1));
    IG = double(I(:,:,2));    
    IB = double(I(:,:,3));
    IS = IR + IG + IB;
    [mask, xcentre, ycentre, sigmaX, sigmaY, meanvalue, mX, mY] = meanAndCenterMass(IR, 1);
    sigmasX(i) = sigmaX;
    sigmasY(i) = sigmaY;  
end

newSigmasX = nan(2*N, 1);
newSigmasY = nan(2*N, 1);
oldCal = load('CoefFitCMAM2.mat');

% Hacerlo con las dosis
for i=1:2*N
    I = allI{i};
    doseRGB = getCalDose(I, oldCal.CoefR, oldCal.CoefG, oldCal.CoefB);    
    [mask, xcentre, ycentre, sigmaX, sigmaY, meanvalue, mX, mY] = meanAndCenterMass(-doseRGB, 1);
    newSigmasX(i) = sigmaX;
    newSigmasY(i) = sigmaY;  
end


sigmas_mm = 0.5*(sigmasX+sigmasY) / pixelsXCM * 10;
newSigmas_mm = 0.5*(newSigmasX+newSigmasY) / pixelsXCM * 10;

meanSigma = mean(sigmas_mm);
stdSigma = std(sigmas_mm);

meanNewSigma = mean(newSigmas_mm);
stdNewSigma = std(newSigmas_mm);

save('Step1.mat');

%% 2. Calcular la energía depositada por disparo (asumiendo 1nA en MUESTRA)

clear all
load('Step1.mat');

sigma = meanNewSigma; % mm
radiusInMM = 0.25;
radiusInVoxels = radiusInMM * pixelsXCM / 10;

fractionNprot = gaussEffect(sigma, radiusInMM, 2e8);

t_shot_ms = 34;
I_nA = 1;
Q_pC = t_shot_ms * I_nA;
rhoW_gcm3 = 1;

e = 1.602e-19; % Carga elemental C
Nprot = Q_pC * 1e-12 / e; % Número total de protones por nA
Nprot_in = Nprot * fractionNprot;

Sw = 164.7; % en MeV/cm, a una distancia de 7.7 cm, por tanto con E = 1.90 MeV
MeV2J = e*1e6;

thickness_um = 30;
thickness_cm = thickness_um*1e-4;
Edep_cyl_J = Nprot_in * Sw * thickness_cm * MeV2J;

Vcyl_cm3 = thickness_cm * pi *((radiusInMM*0.1)^2);
mcyl_kg = Vcyl_cm3 * rhoW_gcm3 * 1e-3;

Dcyl_Gy = Edep_cyl_J / mcyl_kg;

% La capa activa son 30 µm. Asumamos agua.
save('Step2.mat');

%% 3. Calcular la dosis en un determinado radio en cada uno de los puntos
clear all
load('Step2.mat');
factorImuestra = 0.81;
disparos = 1:16;
I_pA = nan(16,1);
I_pA_shots = factorImuestra*linspace(48,45,36);
I_pA(1) = I_pA_shots(1);
I_pA(2) = mean(I_pA_shots(2:3));
I_pA(3) = mean(I_pA_shots(4:6));
I_pA(4) = mean(I_pA_shots(7:10));
I_pA(5) = mean(I_pA_shots(11:15));
I_pA(6) = mean(I_pA_shots(16:21));
I_pA(7) = mean(I_pA_shots(22:28));
I_pA(8) = mean(I_pA_shots(29:36));

I_pA_shots2 = factorImuestra*linspace(43,42,100);
I_pA(9) = mean(I_pA_shots2(1:9));
I_pA(10) = mean(I_pA_shots2(10:19));
I_pA(11) = mean(I_pA_shots2(20:30));
I_pA(12) = mean(I_pA_shots2(31:42));
I_pA(13) = mean(I_pA_shots2(43:55));
I_pA(14) = mean(I_pA_shots2(56:69));
I_pA(15) = mean(I_pA_shots2(70:84));
I_pA(16) = mean(I_pA_shots2(85:100));

I_nA = I_pA' / 1000;
D_final = Dcyl_Gy*disparos.*I_nA;

D_final = [D_final zeros(1,6)];
save('Step3.mat');

%% 4. Calcular el fit para los valores de Dosis en: maxX, maxY, meanValue para R, G, B, RG, RGB
% Obtener todos los datos
clear all
load('Step3.mat');
R = nan(numel(D_final),1);
G = nan(numel(D_final),1);
B = nan(numel(D_final),1);

pxmax = double(intmax(class(allI{1}))); % maximum pixel value 2^16

for i=1:numel(D_final)
    I = allI{i};  
    IR = double(I(:,:,1));
    IG = double(I(:,:,2));    
    IB = double(I(:,:,3));
    
    if D_final(i)==0
        R(i) = mean2(IR) / pxmax;
        G(i) = mean2(IG) / pxmax;
        B(i) = mean2(IB) / pxmax;
    else
        mask = meanAndCenterMass(IG, radiusInVoxels);
        R(i)= mean2(IR(mask==1)) / pxmax;
        G(i)= mean2(IG(mask==1)) / pxmax;
        B(i)= mean2(IB(mask==1)) / pxmax;
    end
end

% % Plot results
% % mean Values
% figure(1);
% plot(D_final, results(:,1), 'r');
% hold on
% plot(D_final, results(:,4), 'g');
% plot(D_final, results(:,7), 'm');
% plot(D_final, results(:,10), 'k');
% xlabel('Dose (Gy)');
% ylabel('Mean values');
% title('Mean values');
% 
% % maxX
% figure(2);
% plot(D_final, results(:,2), 'r');
% hold on
% plot(D_final, results(:,5), 'g');
% plot(D_final, results(:,8), 'm');
% plot(D_final, results(:,11), 'k');
% xlabel('Dose (Gy)');
% ylabel('Max values (X)');
% title('Max values (X)');
% 
% % maxY
% figure(3);
% plot(D_final, results(:,3), 'r');
% hold on
% plot(D_final, results(:,6), 'g');
% plot(D_final, results(:,9), 'm');
% plot(D_final, results(:,12), 'k');
% xlabel('Dose (Gy)');
% ylabel('Max values (Y)');
% title('Max values (Y)');
% 
% % maxXY
% figure(4);
% plot(D_final, results(:,3)+results(:,2), 'r');
% hold on
% plot(D_final, results(:,6)+results(:,5), 'g');
% plot(D_final, results(:,9)+results(:,8), 'm');
% plot(D_final, results(:,12)+results(:,11), 'k');
% xlabel('Dose (Gy)');
% ylabel('Max values (X+Y)');
% title('Max values (X+Y)');

save('Step4.mat');

%% 5. Fits for meanValue
clear all
load('Step4.mat');

%aa = D_final(12);
%D_final(12) = D_final(11);
%D_final(11) = aa;

%aa = D_final(3);
%D_final(3) = D_final(4);
%D_final(4) = aa;

[FR, gofR] = fit(D_final', R, 'rat11');
CoefR = [FR.p2 FR.p1 FR.q1];
[FG, gofG] = fit(D_final', G, 'rat11');
CoefG = [FG.p2 FG.p1 FG.q1];
[FB, gofB] = fit(D_final', B, 'rat11');
CoefB = [FB.p2 FB.p1 FB.q1];

% [ gofR.rsquare gofG.rsquare gofRG.rsquare gofRGB.rsquare ];
% 
% % Fits for maxX and maxY
% normValue = 10000; % At least the maximum value since R must be [0-1]
% R2 = results(:,2) / normValue;
% G2 = results(:,5) / normValue;
% RG2 = results(:,8) / normValue;
% RGB2 = results(:,11) / normValue;
% 
% [FR2, gofR2] = fit(D_final', R2, 'rat11');
% [FG2, gofG2] = fit(D_final', G2, 'rat11');
% [FRG2, gofRG2] = fit(D_final', RG2, 'rat11');
% [FRGB2, gofRGB2] = fit(D_final', RGB2, 'rat11');
% [ gofR2.rsquare gofG2.rsquare gofRG2.rsquare gofRGB2.rsquare ];
% 
% R3 = results(:,3) / normValue;
% G3 = results(:,6) / normValue;
% RG3 = results(:,9) / normValue;
% RGB3 = results(:,12) / normValue;
% 
% [FR3, gofR3] = fit(D_final', R3, 'rat11');
% [FG3, gofG3] = fit(D_final', G3, 'rat11');
% [FRG3, gofRG3] = fit(D_final', RG3, 'rat11');
% [FRGB3, gofRGB3] = fit(D_final', RGB3, 'rat11');
% [ gofR3.rsquare gofG3.rsquare gofRG3.rsquare gofRGB3.rsquare ];
% 
% R23 = 0.5*(R2 + R3);
% G23 = 0.5*(G2 + G3);
% RG23 = 0.5*(RG2 + RG3);
% RGB23 = 0.5*(RGB2 + RGB3);
% [FR23, gofR23] = fit(D_final', R23, 'rat11');
% [FG23, gofG23] = fit(D_final', G23, 'rat11');
% [FRG23, gofRG23] = fit(D_final', RG23, 'rat11');
% [FRGB23, gofRGB23] = fit(D_final', RGB23, 'rat11');
% [ gofR23.rsquare gofG23.rsquare gofRG23.rsquare gofRGB23.rsquare ];
% 
% % --> Best FITS for Green channel: both for the sum 

save('CoefFitCMAM3.mat', 'CoefR', 'CoefG', 'CoefB');
%% 6. Calculate doses using basic code
% Choice of parameters for meanValue
% p1 = FG.p1;
% p2 = FG.p2;
% q1 = FG.q1;
% 
% D_mv = (q1 * G - p2) ./ (p1 - G);
% 
% Choice of parameters for maxXY
% p1b = FR.p1;
% p2b = FR.p2;
% q1b = FR.q1;
% D_mv2 = (q1b * R - p2b) ./ (p1b - R);
% 
% hold off
% plot(D_final, D_mv, 'go');
% hold on
% plot(D_final, D_mv2, 'ro');
% plot(D_final, D_final, 'k');
% 
% sum((D_mv-D_final').^2)
% sum((D_mv2-D_final').^2)

%% 7. Calculate dose using minimization
doseRG = {};
doseRGB = {};
doseR = {};
doseG = {};
deltas = {};
for ii=1:numel(allI)
    Image = allI{ii};
    ImageODR = log10(pxmax./double(Image(:,:,1)));
    ImageODG = log10(pxmax./double(Image(:,:,2)));
    ImageODB = log10(pxmax./double(Image(:,:,3)));
    
    % Perturbed dose
    delta = 0.8:0.001:1.2;
    dev_min = 1.e30*ones(size(Image,1),size(Image,2));
    delta0 = nan(size(Image,1),size(Image,2));    
    
    for k = 1:numel(delta)
        for i = 1:size(Image,1)
            for j = 1:size(Image,2)
                DR = (CoefR(3)*10^(-ImageODR(i,j)*delta(k))-CoefR(1))/(CoefR(2)-10^(-ImageODR(i,j)*delta(k)));
                DG = (CoefG(3)*10^(-ImageODG(i,j)*delta(k))-CoefG(1))/(CoefG(2)-10^(-ImageODG(i,j)*delta(k)));
                DB = (CoefB(3)*10^(-ImageODB(i,j)*delta(k))-CoefB(1))/(CoefB(2)-10^(-ImageODB(i,j)*delta(k)));
                
                dev = (DR-DG)^2+(DR-DB)^2+(DB-DG)^2;
                
                if dev < dev_min(i,j)
                    dev_min(i,j) = dev;
                    delta0(i,j) = delta(k);
                end
            end
        end
    end
    
    DR = (CoefR(3)*10.^(-ImageODR.*delta0)-CoefR(1))./(CoefR(2)-10.^(-ImageODR.*delta0));
    DG = (CoefG(3)*10.^(-ImageODG.*delta0)-CoefG(1))./(CoefG(2)-10.^(-ImageODG.*delta0));
    DB = (CoefB(3)*10.^(-ImageODB.*delta0)-CoefB(1))./(CoefB(2)-10.^(-ImageODB.*delta0));
    dev = (DR-DG).^2+(DR-DB).^2+(DB-DG).^2;
    
    doseR{ii} = DR;
    doseG{ii} = DG;
    doseRG{ii} = (DR+DG)./2;
    doseRGB{ii} = (DR+DG+DB)./3;
    deltas{ii} = delta0;    
end

% Plots
figure(30);
for i=1:numel(allI); subplot(6,4,i); imagesc(doseRG{i}); caxis([0 11]); colorbar; end

figure(31);
for i=1:numel(allI); subplot(6,4,i); imagesc(deltas{i}); caxis([0.8 1.2]); colorbar; end
save('Step5.mat');

%% 8. Calculate mean dose in area R=radiusInVoxels for each image (or total mean dose for empties)
clear all
load('Step5.mat');
dValRG = nan(1, numel(allI));
dValRGB = nan(1, numel(allI));
dValR = nan(1, numel(allI));
dValG = nan(1, numel(allI));

for i = 1:(2*N)

    %[~, ~, ~, ~, ~, dValRG(i), ~, ~] = meanAndCenterMass( -imgaussfilt(doseRG{i},10), radiusInVoxels);
 
    [~, ~, ~, ~, ~, dValRGB(i), ~, ~] = meanAndCenterMass( -imgaussfilt(doseRGB{i},10), radiusInVoxels);

    %[~, ~, ~, ~, ~, dValR(i), ~, ~] = meanAndCenterMass( -imgaussfilt(doseR{i},10), radiusInVoxels);
      
    %[~, ~, ~, ~, ~, dValG(i), ~, ~] = meanAndCenterMass( -imgaussfilt(doseG{i},10), radiusInVoxels);    

end
for i = (2*N+1):numel(doseRG)   
%        dValRG(i) = mean2(doseRG{i});
        dValRGB(i) = mean2(doseRGB{i});
%        dValR(i) = mean2(doseR{i});        
 %       dValG(i) = mean2(doseG{i});                
end
%% Plot
sum((D_final-dValRGB).^2)
hold off
%plot(D_final, dValR, 'ro');
hold on
%plot(D_final, dValG, 'go');
%plot(D_final, dValRG, 'mo');
plot(D_final, dValRGB, 'bo');
plot(D_final, D_final, 'k');
axis([0 11 0 11]);
xlabel('Reference dose [Gy]');
ylabel('Calculated dose [Gy]');
legend('Measured dose', 'Calibration dose', 'Location', 'Northwest');
set(gca, 'FontSize', 14)