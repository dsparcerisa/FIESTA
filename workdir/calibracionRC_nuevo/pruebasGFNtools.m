%% Abrir todos los archivos
clear all; close all

% Choose dataset
filmIndex = 2; % EBT2,2=EBT3,3=EBT3unl
dataIndex = 2; % 1=HF fotones, 2=FQS protones

%% Load data
if dataIndex==1
    load('basicData_HF.mat');
elseif dataIndex==2
    load('basicData_FQS.mat');
else
    error('No data available.');
end

%% Create single value table
meanValues = nan(filmsPerSample, Nsamples, scansPerSample, 3);
stdValues = nan(filmsPerSample, Nsamples, scansPerSample, 3);
for i=1:filmsPerSample
   for j=1:Nsamples
       for k=1:scansPerSample
            for l=1:3
                meanValues(i,j,k,l) = mean2(imageSubsets{i,j,k}(:,:,l));
                stdValues(i,j,k,l) = std2(imageSubsets{i,j,k}(:,:,l));
            end
       end
   end
end

finalPixelValues = nan(filmsPerSample, Nsamples, 3);
finalPixelErrValues = nan(filmsPerSample, Nsamples, 3);
for i=1:filmsPerSample
    for j=1:Nsamples
        for l=1:3
            finalPixelValues(i,j,l) = mean(meanValues(i,j,:,l)) / maxBits;
            finalPixelErrValues(i,j,l) = rssq(stdValues(i,j,:,l)) / scansPerSample / maxBits;
        end
    end
end

%% Estudiar SOLAMENTE uno de los tipos
Rmask = 1:10;
Bmask = 1:10;
Gmask = 1:10;

R = finalPixelValues(filmIndex,:,1);
G = finalPixelValues(filmIndex,:,2);
B = finalPixelValues(filmIndex,:,3);
dR = finalPixelErrValues(filmIndex,:,1);
dG = finalPixelErrValues(filmIndex,:,2);
dB = finalPixelErrValues(filmIndex,:,3);
dosePoints = 0:0.1:(round(max(dosesGy)));

%% Cargar coeficientes necesarios
if dataIndex==1
    load(sprintf('fitCoef_HF_%s.mat',filmOrder{filmIndex}));
elseif dataIndex==2
    load(sprintf('fitCoef_FQS_%s.mat',filmOrder{filmIndex}));
else
    error('Cannot load fit data.');
end

%% Hacer la batería de pruebas

subplot(3,1,1);
% FIT 1
% En rojo
[D_1_R, dD_1_R] = getChannelDoseT1_wErrors(CoefR1, R, dCoefR1, dR)

% En verde
[D_1_G, dD_1_G] = getChannelDoseT1_wErrors(CoefG1, G, dCoefG1, dG)

% En azul
[D_1_B, dD_1_B] = getChannelDoseT1_wErrors(CoefB1, B, dCoefB1, dB)

errorbar(dosesGy,D_1_R,dD_1_R,'r.'); hold on
errorbar(dosesGy,D_1_G,dD_1_G,'g.');
errorbar(dosesGy,D_1_B,dD_1_B,'b.');

D_1 = (D_1_R ./ dD_1_R.^2 +  D_1_G ./ dD_1_G.^2 +  D_1_B ./ dD_1_B.^2) ...
    ./ (1 ./ dD_1_R.^2 +  1 ./ dD_1_G.^2 +  1 ./ dD_1_B.^2);
dD_1 = (1 ./ dD_1_R.^2 +  1 ./ dD_1_G.^2 +  1 ./ dD_1_B.^2).^(-0.5);

errorbar(dosesGy,D_1,dD_1,'k.'); hold on
plot(dosePoints,dosePoints,'y:');
grid on
axis([0 round(max(dosesGy)) 0 round(max(dosesGy))]);
residue_1 = rssq(D_1-dosesGy)
residue_1b = rssq(D_1(1:end-1)-dosesGy(1:end-1))
%%
% FIT 2
% En rojo
D_2_R = getChannelDoseT2(CoefR2, R);
dD_2_R = max(getChannelDoseT2(CoefR2, R-dR) - D_2_R, D_2_R - getChannelDoseT2(CoefR2, R+dR));
% En verde
D_2_G = getChannelDoseT2(CoefG2, G);
dD_2_G = max(getChannelDoseT2(CoefG2, G-dG) - D_2_G, D_2_G - getChannelDoseT2(CoefG2, G+dG));
% En azul
D_2_B = getChannelDoseT2(CoefB2, B);
dD_2_B = max(getChannelDoseT2(CoefB2, B-dB) - D_2_B, D_2_B - getChannelDoseT2(CoefB2, B+dB));
% errorbar(dosesGy,D_2_R,dD_2_R,'r.'); hold on
% errorbar(dosesGy,D_2_G,dD_2_G,'g.');
% errorbar(dosesGy,D_2_B,dD_2_B,'b.');
D_2 = (D_2_R ./ dD_2_R.^2 +  D_2_G ./ dD_2_G.^2 +  D_2_B ./ dD_2_B.^2) ...
    ./ (1 ./ dD_2_R.^2 +  1 ./ dD_2_G.^2 +  1 ./ dD_2_B.^2);
dD_2 = (1 ./ dD_2_R.^2 +  1 ./ dD_2_G.^2 +  1 ./ dD_2_B.^2).^(-0.5);
errorbar(dosesGy,D_2,dD_2,'g.'); 
residue_2 = rssq(D_2-dosesGy)
residue_2b = rssq(D_2(1:end-1)-dosesGy(1:end-1))
% FIT 3
% En rojo
D_3_R = getChannelDoseT3(CoefR3, R, R(1));
dD_3_R = max(getChannelDoseT3(CoefR3, R-dR, R(1)) - D_3_R, D_3_R - getChannelDoseT3(CoefR3, R+dR, R(1)));
% En verde
D_3_G = getChannelDoseT3(CoefG3, G, G(1));
dD_3_G = max(getChannelDoseT3(CoefG3, G-dG, G(1)) - D_3_G, D_3_G - getChannelDoseT3(CoefG3, G+dG, G(1)));
% En azul
D_3_B = getChannelDoseT3(CoefB3, B, B(1));
dD_3_B = max(getChannelDoseT3(CoefB3, B-dB, B(1)) - D_3_B, D_3_B - getChannelDoseT3(CoefB3, B+dB, B(1)));
% errorbar(dosesGy,D_3_R,dD_3_R,'r.'); hold on
% errorbar(dosesGy,D_3_G,dD_3_G,'g.');
% errorbar(dosesGy,D_3_B,dD_3_B,'b.');
D_3 = (D_3_R ./ dD_3_R.^2 +  D_3_G ./ dD_3_G.^2 +  D_3_B ./ dD_3_B.^2) ...
    ./ (1 ./ dD_3_R.^2 +  1 ./ dD_3_G.^2 +  1 ./ dD_3_B.^2);
dD_3 = (1 ./ dD_3_R.^2 +  1 ./ dD_3_G.^2 +  1 ./ dD_3_B.^2).^(-0.5);
errorbar(dosesGy,D_3,dD_3,'b.'); 
residue_3 = rssq(D_3-dosesGy)
residue_3b = rssq(D_3(1:end-1)-dosesGy(1:end-1))

subplot(3,1,2);
errorbar(dosesGy, D_1-dosesGy, dD_1, 'r.'); hold on
errorbar(dosesGy, D_2-dosesGy, dD_2, 'g.');
errorbar(dosesGy, D_3-dosesGy, dD_3, 'b.');
legend('Fit I','Fit II', 'Fit III');
subplot(3,1,1);
title('Photons')
xlabel('Nominal dose (Gy)');
ylabel('Dose averaging three channels (Gy)');
if dataIndex==1
    title(sprintf('Photons - %s', filmOrder{filmIndex}));
else
    title(sprintf('Protons - %s', filmOrder{filmIndex}));
end
subplot(3,1,2);
ylabel('Absolute error (Gy)');
xlabel('Nominal dose (Gy)');
legend('Fit I','Fit II', 'Fit III');
ylim([-0.5 0.5]);
grid on

subplot(3,1,3);
errorbar(dosesGy, 100*(D_1-dosesGy)./dosesGy, 100*dD_1./dosesGy, 'r.'); hold on
errorbar(dosesGy, 100*(D_2-dosesGy)./dosesGy, 100*dD_2./dosesGy, 'g.');
errorbar(dosesGy, 100*(D_3-dosesGy)./dosesGy, 100*dD_3./dosesGy, 'b.');
ylabel('Relative error (%)');
xlabel('Nominal dose (Gy)');
legend('Fit I','Fit II', 'Fit III');
ylim([-10 10]);
grid on