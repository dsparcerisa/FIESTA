clear all; close all

%% Energy at film surface

a = 2.08; % cmMev^-p
p = 1.75;
b = 0.265; % MeV*1.73
q = -0.73;

Es = @(z,E0) ((E0 - b*E0.^q).^p - z./a).^(1/p)

%% Average LETd inside film

LETa_3 = 4.08E05;
LETb_3 = 2.88;
LETc_3 = 22.5;
LETd_3 = 0.142;

LETa_2 = 8.66E03;
LETb_2 = 2.11;
LETc_2 = 19.3;
LETd_2 = 0.136;

LETa_u = 94.8;
LETb_u = 1.56;
LETc_u = 17.1;
LETd_u = 0.122;

LET_Es_3 = @(Es) LETa_3 .* exp(-LETb_3.*Es) + LETc_3 .* exp(-LETd_3.*Es)
LET_Es_2 = @(Es) LETa_2 .* exp(-LETb_2.*Es) + LETc_2 .* exp(-LETd_2.*Es)
LET_Es_u = @(Es) LETa_u .* exp(-LETb_u.*Es) + LETc_u .* exp(-LETd_u.*Es)

LET_E_3 = @(z,E0) LET_Es_3(Es(z,E0))
LET_E_2 = @(z,E0) LET_Es_2(Es(z,E0))
LET_E_u = @(z,E0) LET_Es_u(Es(z,E0))

%% Relative efficiency
a_RE = 0.012;
b_RE = 1.01;
RE_LET = @(L) 1 -  a_RE .* L .^ b_RE;
RE_E_3 = @(z,E0) RE_LET(LET_E_3(z,E0));
RE_E_2 = @(z,E0) RE_LET(LET_E_2(z,E0));
RE_E_u = @(z,E0) RE_LET(LET_E_u(z,E0));


%% Load data
distances = 5.3 + (0:1:10);
distPlot = 4:0.1:16;
radius_cm = 0.03;

EBT2_5_radName = 'Rad8';
EBT2_10_radName = 'Rad18';
[LET2_5, dLET2_5, RE2_5, dRE2_5] = processAnyRad(EBT2_5_radName, radius_cm);
[LET2_10, dLET2_10, RE2_10, dRE2_10] = processAnyRad(EBT2_10_radName, radius_cm);

EBT3_5_radName = 'Rad7';
EBT3_10_radName = 'Rad17';
[LET3_5, dLET3_5, RE3_5, dRE3_5] = processAnyRad(EBT3_5_radName, radius_cm);
[LET3_10, dLET3_10, RE3_10, dRE3_10] = processAnyRad(EBT3_10_radName, radius_cm);

EBT3unl_5_radName = 'Rad6';
EBT3unl_5b_radName = 'Rad9';
EBT3unl_10_radName = 'Rad16';
[LET3unl_5, dLET3unl_5, RE3unl_5, dRE3unl_5] = processAnyRad(EBT3unl_5_radName, radius_cm);
[LET3unl_5b, dLET3unl_5b, RE3unl_5b, dRE3unl_5b] = processAnyRad(EBT3unl_5b_radName, radius_cm);
[LET3unl_10, dLET3unl_10, RE3unl_10, dRE3unl_10] = processAnyRad(EBT3unl_10_radName, radius_cm);

%% Plot
close all
set(gcf, 'Position',  [100, 100, 900, 300])
ha = tight_subplot(1,3,0.02,[0.15 0.08],[0.06 0.02])
axes(ha(1))
errorbar(distances, RE3_5, dRE3_5, 'bo')
hold on
errorbar(distances, RE3_10, dRE3_10, 'rx')
plot(distPlot, RE_E_3(distPlot,5),'b-', 'LineWidth', 2);
plot(distPlot, RE_E_3(distPlot,10),'r-', 'LineWidth', 2);
legend({'5 MeV','10 MeV', 'Model 5 MeV', 'Model 10 MeV'},'Location','SouthWest');
title('EBT3');
xlabel('Air gap (cm)');
grid on
ylim([0 1.2])
ylabel('Relative efficiency');
set(gca,'FontSize',14)


axes(ha(2))
errorbar(distances, RE2_5, dRE2_5, 'bo')
hold on
errorbar(distances, RE2_10, dRE2_10, 'rx')
plot(distPlot, RE_E_2(distPlot,5),'b-', 'LineWidth', 2);
plot(distPlot, RE_E_2(distPlot,10),'r-', 'LineWidth', 2);
%legend({'5 MeV','10 MeV', 'Model 5 MeV', 'Model 10 MeV'},'Location','SouthWest');
title('EBT2');
xlabel('Air gap (cm)');
grid on
ylim([0 1.2])
set(gca,'FontSize',14)
yticklabels([]);

axes(ha(3))
%errorbar(distances, RE3unl_5, dRE3unl_5, 'b.')
errorbar(distances, RE3unl_5b, dRE3unl_5b, 'bo')
hold on
errorbar(distances, RE3unl_10, dRE3unl_10, 'rx')
plot(distPlot, RE_E_u(distPlot,5),'b-', 'LineWidth', 2);
plot(distPlot, RE_E_u(distPlot,10),'r-', 'LineWidth', 2);
%legend({'5 MeV','10 MeV', 'Model 5 MeV', 'Model 10 MeV'},'Location','SouthWest');
title('EBT3 unlaminated');
xlabel('Air gap (cm)');
grid on
ylim([0 1.2])
set(gca,'FontSize',14)
yticklabels([]);

%% Cálculo de incertidubres de ejemplo
L1 = LET_E_2(15,10)
L2 = LET_E_3(5,4)
RE1 = RE_LET(L1)
RE2 = RE_LET(L2)
da_RE = 0.002;
db_RE = 0.04;

dRE1 = sqrt( ( L1^b_RE * da_RE )^2 + ((1-RE1)*log(L1)*db_RE)^2  )
dRE2 = sqrt( ( L2^b_RE * da_RE )^2 + ((1-RE2)*log(L2)*db_RE)^2  )
