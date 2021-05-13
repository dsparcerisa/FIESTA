%% Abrir todos los archivos
clear all; close all

% Choose dataset
filmIndex = 2; % EBT2,2=EBT3,3=EBT3unl
dataIndex = 3; % 1=HF fotones, 2=FQS protones, 3=Medicina

for filmIndex=1:3
subplot(1,3,filmIndex)
%% Load data
if dataIndex==1
    load('basicData_HF.mat', 'filmsPerSample', 'Nsamples', 'scansPerSample', 'imageSubsets', 'maxBits', 'dosesGy', 'filmOrder');
elseif dataIndex==2
    load('basicData_FQS.mat', 'filmsPerSample', 'Nsamples', 'scansPerSample', 'imageSubsets', 'maxBits', 'dosesGy', 'filmOrder');
elseif dataIndex==3
    load('basicData_Medicina_SCANNED.mat', 'filmsPerSample', 'Nsamples', 'scansPerSample', 'imageSubsets', 'maxBits', 'dosesGy', 'filmOrder');
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

%% Fit type I
figure(1)

[fR1, gofR1] = fitTypeI(dosesGy(Rmask), R(Rmask), dR(Rmask));
plot(dosePoints, fR1(dosePoints), 'r-'); hold on
[fG1, gofG1] = fitTypeI(dosesGy, G, dG);
plot(dosePoints, fG1(dosePoints), 'g-');
[fB1, gofB1] = fitTypeI(dosesGy(Bmask), B(Bmask), dB(Bmask));
plot(dosePoints, fB1(dosePoints), 'b-');
title(sprintf('%s', filmOrder{filmIndex}));
errorbar(dosesGy, R, dR, 'r.')
hold on
errorbar(dosesGy, G, dG, 'g.')
errorbar(dosesGy, B, dB, 'b.')
grid on
xlabel('Doses (Gy)')
if (filmIndex==1)
    ylabel('Relative PV');
end
xlim([0 round(max(dosesGy))]);

CoefR1 = [fR1.alpha fR1.beta fR1.gamma];
CoefG1 = [fG1.alpha fG1.beta fG1.gamma];
CoefB1 = [fB1.alpha fB1.beta fB1.gamma];

ypos = max([max(R), max(G), max(B)]) * 0.8;
%text(8,ypos,sprintf('R^2 = %3.4f',gofR1.rsquare),'Color','red')
%text(8,ypos*0.9,sprintf('R^2 = %3.4f',gofG1.rsquare),'Color','green')
%text(8,ypos*0.8,sprintf('R^2 = %3.4f',gofB1.rsquare),'Color','blue')

axis([0 20 0 1]);
end

dataIndex = 2;
hold on

for filmIndex=1:3
subplot(1,3,filmIndex)
%% Load data
if dataIndex==1
    load('basicData_HF.mat', 'filmsPerSample', 'Nsamples', 'scansPerSample', 'imageSubsets', 'maxBits', 'dosesGy', 'filmOrder');
elseif dataIndex==2
    load('basicData_FQS.mat', 'filmsPerSample', 'Nsamples', 'scansPerSample', 'imageSubsets', 'maxBits', 'dosesGy', 'filmOrder');
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

%% Fit type I
figure(1)

[fR1, gofR1] = fitTypeI(dosesGy(Rmask), R(Rmask), dR(Rmask));
plot(dosePoints, fR1(dosePoints), 'r:', 'LineWidth', 2);
hold on
[fG1, gofG1] = fitTypeI(dosesGy, G, dG);
plot(dosePoints, fG1(dosePoints), 'g:','LineWidth', 2);
[fB1, gofB1] = fitTypeI(dosesGy(Bmask), B(Bmask), dB(Bmask));
plot(dosePoints, fB1(dosePoints), 'b:','LineWidth', 2);
errorbar(dosesGy, R, dR, 'k.')
hold on
errorbar(dosesGy, G, dG, 'k.')
errorbar(dosesGy, B, dB, 'k.')
title(sprintf('%s', filmOrder{filmIndex}));
grid on
xlabel('Doses (Gy)')
xlim([0 round(max(dosesGy))]);

CoefR1 = [fR1.alpha fR1.beta fR1.gamma];
CoefG1 = [fG1.alpha fG1.beta fG1.gamma];
CoefB1 = [fB1.alpha fB1.beta fB1.gamma];

ypos = max([max(R), max(G), max(B)]) * 0.8;
%text(8,ypos,sprintf('R^2 = %3.4f',gofR1.rsquare),'Color','red')
%text(8,ypos*0.9,sprintf('R^2 = %3.4f',gofG1.rsquare),'Color','green')
%text(8,ypos*0.8,sprintf('R^2 = %3.4f',gofB1.rsquare),'Color','blue')

axis([0 12 0 0.8]);
set(gca,'FontSize',12)
end
subplot(1,3,3)
title('EBT3 unlaminated');
set(gcf, 'Position',  [100, 100, 1200, 400])

% 
% %% Fit type II
% figure(2)
% 
% ODR = -log10(R);
% dODR = (dR./R./log(10));
% ODG = -log10(G);
% dODG = (dG./G./log(10));
% ODB = -log10(B);
% dODB = (dB./B./log(10));
% 
% %subplot(1,2,1);
% % errorbar(dosesGy, ODR, dODR, 'r.')
% % hold on
% % errorbar(dosesGy, ODG, dODG, 'g.')
% % errorbar(dosesGy, ODB, dODB, 'b.')
% % xlabel('Doses (Gy)')
% % ylabel('OD');
% % grid on
% % [fR2, gofR2] = fitTypeII(dosesGy(Rmask), ODR(Rmask), dODR(Rmask));
% % plot(dosePoints, fR2(dosePoints), 'r-');
% % [fG2, gofG2] = fitTypeII(dosesGy, ODG, dODG);
% % plot(dosePoints, fG2(dosePoints), 'g-');
% % [fB2, gofB2] = fitTypeII(dosesGy(Bmask), ODB(Bmask), dODB(Bmask));
% % plot(dosePoints, fB2(dosePoints), 'b-');[]
% % 
% % text(6,0.4,sprintf('R^2 = %3.4f',gofR2.rsquare),'Color','red')
% % text(6,0.45,sprintf('R^2 = %3.4f',gofG2.rsquare),'Color','green')
% % text(6,0.5,sprintf('R^2 = %3.4f',gofB2.rsquare),'Color','blue')
% % 
% % xlim([0 round(max(dosesGy))]);
% % title(sprintf('Type II - Poly (%s)', filmOrder{filmIndex}));
% % 
% % subplot(1,2,2);
% errorbar(dosesGy, ODR, dODR, 'r.')
% hold on
% errorbar(dosesGy, ODG, dODG, 'g.')
% errorbar(dosesGy, ODB, dODB, 'b.')
% xlabel('Doses (Gy)')
% ylabel('OD');
% grid on
% [fR2r, gofR2r] = fitTypeIIb(dosesGy(Rmask), ODR(Rmask), dODR(Rmask));
% plot(dosePoints, fR2r(dosePoints), 'r-');
% [fG2r, gofG2r] = fitTypeIIb(dosesGy, ODG, dODG);
% plot(dosePoints, fG2r(dosePoints), 'g-');
% [fB2r, gofB2r] = fitTypeIIb(dosesGy(Bmask), ODB(Bmask), dODB(Bmask));
% plot(dosePoints, fB2r(dosePoints), 'b-');
% 
% text(6,0.3,sprintf('R^2 = %3.4f',gofR2r.rsquare),'Color','red')
% text(6,0.4,sprintf('R^2 = %3.4f',gofG2r.rsquare),'Color','green')
% text(6,0.5,sprintf('R^2 = %3.4f',gofB2r.rsquare),'Color','blue')
% 
% xlim([0 round(max(dosesGy))]);
% title(sprintf('Type II - Rational (%s)', filmOrder{filmIndex}));
% 
% CoefR2 = [fR2r.p1 fR2r.p2 fR2r.q1];
% CoefG2 = [fG2r.p1 fG2r.p2 fG2r.q1];
% CoefB2 = [fB2r.p1 fB2r.p2 fB2r.q1];
% 
% %% Type 3
% NODR = log10(R(1) ./ R);
% dNODR = (1./R./log(10)).*dR;
% NODG = log10(G(1) ./ G);
% dNODG = (1./G./log(10)).*dG;
% NODB = log10(B(1) ./ B);
% dNODB = (1./B./log(10)).*dB;
% 
% figure(3)
% errorbar(dosesGy, NODR, dNODR, 'r.')
% hold on
% errorbar(dosesGy, NODG, dNODG, 'g.')
% errorbar(dosesGy, NODB, dNODB, 'b.')
% xlabel('Doses (Gy)')
% ylabel('NOD');
% grid on
% 
% [fR3, gofR3] = fitTypeIII(dosesGy(Rmask), NODR(Rmask), dNODR(Rmask));
% plot(dosePoints, fR3(dosePoints), 'r-');
% [fG3, gofG3] = fitTypeIII(dosesGy, NODG, dNODG);
% plot(dosePoints, fG3(dosePoints), 'g-');
% [fB3, gofB3] = fitTypeIII(dosesGy(Bmask), NODB(Bmask), dNODB(Bmask));
% plot(dosePoints, fB3(dosePoints), 'b-');
% title(sprintf('Type III - Potential (%s)', filmOrder{filmIndex}));
% text(1,0.5,sprintf('R^2 = %3.4f',gofR3.rsquare),'Color','red')
% text(1,0.45,sprintf('R^2 = %3.4f',gofG3.rsquare),'Color','green')
% text(1,0.4,sprintf('R^2 = %3.4f',gofB3.rsquare),'Color','blue')
% ylim([0 1]);
% 
% CoefR3 = [fR3.alpha fR3.beta fR3.gamma];
% CoefG3 = [fG3.alpha fG3.beta fG3.gamma];
% CoefB3 = [fB3.alpha fB3.beta fB3.gamma];
% 
% %% Save data
% 
% if dataIndex==1
%     save(sprintf('fitCoef_HF_%s.mat',filmOrder{filmIndex}),'CoefR1','CoefG1','CoefB1','CoefR2','CoefG2','CoefB2', 'CoefR3','CoefG3','CoefB3');
% elseif dataIndex==2
%     save(sprintf('fitCoef_FQS_%s.mat',filmOrder{filmIndex}),'CoefR1','CoefG1','CoefB1','CoefR2','CoefG2','CoefB2', 'CoefR3','CoefG3','CoefB3');
% else
%     error('Cannot save data.');
% end
