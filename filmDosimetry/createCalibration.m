%% create Calibration file - generic code
clear all
FIESTAConfig

% Info for processing
% pixelsXCM = 157.48; % 400 ppp
% pixelsXCM = 118.11; % 300 ppp

% File paths
fileCal1 = fullfile(Fcfg.filmPath,'cal0006.tif');
fileCal2 = fullfile(Fcfg.filmPath,'cal0005.tif');
fileCal3 = fullfile(Fcfg.filmPath,'cal0002.tif');

N = 4; % Number of areas per image
N0 = 2; % añadir 2 puntos más con dosis 0

[pixelsXcm, maxInt] = getImgMetaInfo(fileCal1);

%% 0. Read and crop files
filePaths = {fileCal1, fileCal2, fileCal3};
nCrops = [N+N0, N+N0, N+N0];
allI = loadNcropFiles(filePaths, nCrops);

%% 1. Specify doses
D_final = [
    0.9910
    2.0162
    3.0072
    3.9558
    0
    0
    5.0088
    6.0136
    7.0046
    7.9302
    0
    0
    9.0215
   10.0114
   11.0029
   11.8858
   0
   0]; 

% Temporary fix (positions changed): must see logs before
aa = D_final(15);
D_final(15) = D_final(16);
D_final(16) = aa;

%% 2. Calculate RGB relative values

R = nan(numel(D_final),1);
G = nan(numel(D_final),1);
B = nan(numel(D_final),1);

pxmax = double(maxInt); % maximum pixel value

for i=1:numel(allI)
    I = allI{i};
    IR = double(I(:,:,1));
    IG = double(I(:,:,2));
    IB = double(I(:,:,3));
    
    R(i) = mean2(IR) / pxmax;
    G(i) = mean2(IG) / pxmax;
    B(i) = mean2(IB) / pxmax;
    
end

%% 3. Fit and save calibration
[FR, gofR] = fit(D_final, R, 'rat11');
CoefR = [FR.p2 FR.p1 FR.q1];
[FG, gofG] = fit(D_final, G, 'rat11');
CoefG = [FG.p2 FG.p1 FG.q1];
[FB, gofB] = fit(D_final, B, 'rat11');
CoefB = [FB.p2 FB.p1 FB.q1];

plot(FR, 'r-'); hold on
plot(D_final,R,'ro'); 
plot(FG, 'g-');
plot(D_final,G,'go'); 
plot(FB, 'b-');
plot(D_final,B,'bo'); 
title('Fit');
xlabel('Dose (Gy)');
ylabel('RGB value');
grid on

save('CoefFitCMAM_HG.mat', 'CoefR', 'CoefG', 'CoefB');

%% 4. Test calibration
D_retest = nan(numel(allI), 1);
std_retest = nan(numel(allI), 1);
allD = {};
deltas = 0.99:0.001:1.01;
for i=1:numel(allI)
    dose = getDoseFromRC(allI{i},CoefR, CoefG, CoefB, pixelsXcm, deltas);
    allD{i} = dose.data;
    D_retest(i) = mean(dose.data(:));
    std_retest(i) = std(dose.data(:));
end

% Temporary fix (positions changed): must see logs before
aa = D_final(15);
D_final(15) = D_final(16);
D_final(16) = aa;

figure(1)
plotNimages(allD);
figure(2)
errorbar(D_final, D_retest, std_retest, 'b.');
hold on
plot(D_final, D_final, 'ro');
legend('Measured dose', 'Nominal dose', 'Location', 'Northwest');
xlabel('Nominal dose (Gy)');
ylabel('Measured dose (Gy)');
title('Validation of calibration curve');
grid on
set(gca, 'FontSize', 14)