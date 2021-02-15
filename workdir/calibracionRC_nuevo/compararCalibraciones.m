MU_p = [0 443 887 1330 1774 3548 5322 7095 8869 10643];
% 1681 MU = 1.8953 Gy (valores de FQS)
dosesGy_p = MU_p*1.8953/1681;

MU_phot = [0 52 105 157 209 419 629 838 1048 1257 2095];
% 104.8 MU = 1 Gy (valores de HF)
dosesGy_phot  = MU_phot ./ 104.8;

plotPoints = 0:0.1:12;

% 1. Plottear los fits de todas ellas
% ans(x) = (p1*x + p2) / (x + q1)
applyRatFit = @(d, Coef) (Coef(2).*d + Coef(1)) ./ (d + Coef(3));

% EBT2
subplot(1,3,1);
load('fitCoef_FQS_EBT2.mat')
plot(plotPoints, applyRatFit(plotPoints, CoefR1), 'r-'); hold on
plot(plotPoints, applyRatFit(plotPoints, CoefG1), 'g-'); 
plot(plotPoints, applyRatFit(plotPoints, CoefB1), 'b-'); 
load('fitCoef_HF_EBT2.mat')
plot(plotPoints, applyRatFit(plotPoints, CoefR1), 'r:'); hold on
plot(plotPoints, applyRatFit(plotPoints, CoefG1), 'g:'); 
plot(plotPoints, applyRatFit(plotPoints, CoefB1), 'b:'); 
xlim([0 12]); grid on
ylim([0 0.8]);
ylabel('Pixel Value');
title('EBT2');
xlabel('Dose (Gy)');
% EBT3
subplot(1,3,2);
load('fitCoef_FQS_EBT3.mat')
plot(plotPoints, applyRatFit(plotPoints, CoefR), 'r-'); hold on
plot(plotPoints, applyRatFit(plotPoints, CoefG), 'g-'); 
plot(plotPoints, applyRatFit(plotPoints, CoefB), 'b-'); 
load('fitCoef_HF_EBT3.mat')
plot(plotPoints, applyRatFit(plotPoints, CoefR), 'r:'); hold on
plot(plotPoints, applyRatFit(plotPoints, CoefG), 'g:'); 
plot(plotPoints, applyRatFit(plotPoints, CoefB), 'b:'); 
xlim([0 12]); grid on
title('EBT3');
ylim([0 0.8]);
xlabel('Dose (Gy)');

% EBT3unl
subplot(1,3,3);
load('fitCoef_FQS_EBT3unl.mat')
plot([13 13],[0 0],'k-'); hold on
plot([13 13],[0 0],'k:'); hold on

plot(plotPoints, applyRatFit(plotPoints, CoefR), 'r-'); hold on
plot(plotPoints, applyRatFit(plotPoints, CoefG), 'g-'); 
plot(plotPoints, applyRatFit(plotPoints, CoefB), 'b-'); 
load('fitCoef_HF_EBT3unl.mat')
plot(plotPoints, applyRatFit(plotPoints, CoefR), 'r:'); hold on
plot(plotPoints, applyRatFit(plotPoints, CoefG), 'g:'); 
plot(plotPoints, applyRatFit(plotPoints, CoefB), 'b:'); 
title('EBT3unl');
xlim([0 12]); grid on
xlabel('Dose (Gy)');
ylim([0 0.8]);
legend('Calibration with protons', 'Calibration with photons');