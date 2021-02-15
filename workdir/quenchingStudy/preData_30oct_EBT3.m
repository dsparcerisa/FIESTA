clear all
load('preData_30oct_EBT3.mat');
E = myData(:,1);
LET = myData(:,2);
ratio = myData(:,3);
dratio = myData(:,4);
doses = myData(:,5);
vM = doses<12;

%% Fit to rational (Energy)
figure(1);
errorbar(E(vM),ratio(vM),dratio(vM),'.');
grid on
wgt = dratio.^(-2);
F = fit(E(vM), ratio(vM), 'rat22', 'Weights', wgt(vM))
hold on;
Evals = 0:0.1:10;
plot(Evals, F(Evals), 'R-');
xlabel('Energy (MeV)');
ylabel('RE');

%% Fit to rational (LET)
figure(2);
hold off
errorbar(LET(vM),ratio(vM),dratio(vM),'.'); hold on
%errorbar(LET(vM==0),ratio(vM==0),dratio(vM==0),'g.');

grid on
wgt = dratio.^(-2);
ft = fittype( 'D0 / (a*x + D0)', 'independent', 'x', 'dependent', 'y' );
ft2 = fittype( '1 - a*x^b', 'independent', 'x', 'dependent', 'y' );

F = fit(LET(vM), ratio(vM), ft, 'Weights', wgt(vM));% 'Lower', [45 0], 'Upper', [45 1])
F2 = fit(LET(vM), ratio(vM), ft2, 'Weights', wgt(vM), 'Lower', [0 0], 'Upper', [1 5])

hold on;
Evals = 10:1:400;
%plot(Evals, F(Evals), 'm-');
plot(Evals, F2(Evals), 'r-');
xlabel('LET (MeV/cm)');
ylabel('RE');
ylim([0 2]);
set(gca,'xscale','log')
axis([40 300 0 1.40]);
title('EBT3');