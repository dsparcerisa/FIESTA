clear all; close all;

%% Figura que muestre el LET promedio y su barra de error en función de la energía en cada punto.

Zpos = 5.3 + (0:10);
E_EBT3 = [];
LET_EBT3 = [];
dLET_EBT3 = [];
E_EBT2 = [];
LET_EBT2 = [];
dLET_EBT2 = [];
E_EBT3unl = [];
LET_EBT3unl = [];
dLET_EBT3unl = [];
LETw_EBT3 = [];
dLETw_EBT3 = [];
LETw_EBT2 = [];
dLETw_EBT2 = [];
LETw_EBT3unl = [];
dLETw_EBT3unl = [];

filmType = 'EBT3';

for E0=[4 5 6 8 10]
%E0 = 4;

%% Para calcular la eficiencia relativa mejor juntar todos los valores de todas las medidas
kaptonThickness_um = 8;
kaptonThickness_cm = 1e-4*kaptonThickness_um;
E0kap_vector = energyStoppingPowerKapton(E0, [0 kaptonThickness_cm]);
E0kap = E0kap_vector(2);
% 
airVecPos = 0:0.01:20;
[Eair_vec, ~, Sw_vec] = energyStoppingPower(E0kap, airVecPos);
validMask = ~isnan(Eair_vec);
Epos_preFilm = interp1(airVecPos(validMask), Eair_vec(validMask), Zpos');

% Load LET values
LETpath = '/Users/dani/Google Drive/UCM/01-Proyectos/25-TOPAS/LETScorer2/results';
LETFileName = sprintf('LETvalues_%s_%iMeV',filmType,E0);
fullLETPath = [LETpath filesep LETFileName];
load(fullLETPath);
Sw_preFilm = meanLET';
dSw_preFilm = errLET';
LETw = meanLETW';
dLETw = errLETW';
%% Make plot
if E0==4 
    Epos_preFilm = Epos_preFilm(1:3);
    Sw_preFilm = Sw_preFilm(1:3);
    dSw_preFilm = dSw_preFilm(1:3);
    LETw = LETw(1:3);
    dLETw = dLETw(1:3);
end

E_EBT3 = [E_EBT3; Epos_preFilm];
LET_EBT3 = [LET_EBT3; Sw_preFilm];
dLET_EBT3 = [dLET_EBT3; dSw_preFilm];
LETw_EBT3 = [LETw_EBT3; LETw]; 
dLETw_EBT3 = [LETw_EBT3; dLETw];

%errorbar(Epos_preFilm, Sw_preFilm, dSw_preFilm, 'b.')

hold on
end



filmType = 'EBT2';
for E0=[4 5 6 8 10]
%E0 = 4;

%% Para calcular la eficiencia relativa mejor juntar todos los valores de todas las medidas
kaptonThickness_um = 8;
kaptonThickness_cm = 1e-4*kaptonThickness_um;
E0kap_vector = energyStoppingPowerKapton(E0, [0 kaptonThickness_cm]);
E0kap = E0kap_vector(2);
% 
airVecPos = 0:0.01:20;
[Eair_vec, ~, Sw_vec] = energyStoppingPower(E0kap, airVecPos);
validMask = ~isnan(Eair_vec);
Epos_preFilm = interp1(airVecPos(validMask), Eair_vec(validMask), Zpos');

% Load LET values
LETpath = '/Users/dani/Google Drive/UCM/01-Proyectos/25-TOPAS/LETScorer2/results';
LETFileName = sprintf('LETvalues_%s_%iMeV',filmType,E0);
fullLETPath = [LETpath filesep LETFileName];
load(fullLETPath);
Sw_preFilm = meanLET';
dSw_preFilm = errLET';
LETw = meanLETW';
dLETw = errLETW';

%% Make plot
if E0==4 
    Epos_preFilm(9:10) = [];
    Sw_preFilm(9:10) = [];
    dSw_preFilm(9:10) = [];
    LETw(9:10) = [];
    dLETw(9:10) = [];
end

E_EBT2 = [E_EBT2; Epos_preFilm];
LET_EBT2 = [LET_EBT2; Sw_preFilm];
dLET_EBT2 = [dLET_EBT2; dSw_preFilm];
LETw_EBT2 = [LETw_EBT2; LETw]; 
dLETw_EBT2 = [LETw_EBT2; dLETw];
hold on
end

filmType = 'EBT3unl';
for E0=[3 4 5 6 8 10]
%E0 = 4;

%% Para calcular la eficiencia relativa mejor juntar todos los valores de todas las medidas
kaptonThickness_um = 8;
kaptonThickness_cm = 1e-4*kaptonThickness_um;
E0kap_vector = energyStoppingPowerKapton(E0, [0 kaptonThickness_cm]);
E0kap = E0kap_vector(2);
% 
airVecPos = 0:0.01:20;
[Eair_vec, ~, Sw_vec] = energyStoppingPower(E0kap, airVecPos);
validMask = ~isnan(Eair_vec);
Epos_preFilm = interp1(airVecPos(validMask), Eair_vec(validMask), Zpos');

% Load LET values
LETpath = '/Users/dani/Google Drive/UCM/01-Proyectos/25-TOPAS/LETScorer2/results';
LETFileName = sprintf('LETvalues_%s_%iMeV',filmType,E0);
fullLETPath = [LETpath filesep LETFileName];
load(fullLETPath);
Sw_preFilm = meanLET';
dSw_preFilm = errLET';
LETw = meanLETW';
dLETw = errLETW';


%% Make plot
if E0==3 
    Epos_preFilm(9:10) = [];
    Sw_preFilm(9:10) = [];
    dSw_preFilm(9:10) = [];
    LETw(9:10) = [];
    dLETw(9:10) = [];    
end

E_EBT3unl = [E_EBT3unl; Epos_preFilm];
LET_EBT3unl = [LET_EBT3unl; Sw_preFilm];
dLET_EBT3unl = [dLET_EBT3unl; dSw_preFilm];
LETw_EBT3unl = [LETw_EBT3unl; LETw]; 
dLETw_EBT3unl = [LETw_EBT3unl; dLETw];

end

%% Make plot
subplot(1,2,1);
errorbar(E_EBT3, LET_EBT3, dLET_EBT3, '.', 'Color', [0, 0.4470, 0.7410]);
hold on
errorbar(E_EBT2, LET_EBT2, dLET_EBT2, '.', 'Color', [0.8500, 0.3250, 0.0980]); 
errorbar(E_EBT3unl, LET_EBT3unl, dLET_EBT3unl, '.', 'Color', [0.4660 0.6740 0.1880])

%errorbar(E_EBT3, LETw_EBT3, dLET_EBT3, '.', 'Color', 0.5*[0, 0.4470, 0.7410]);
%errorbar(E_EBT2, LETw_EBT2, dLET_EBT2, '.', 'Color', 0.5*[0.8500, 0.3250, 0.0980]); 
%errorbar(E_EBT3unl, LETw_EBT3unl, dLET_EBT3unl, '.', 'Color', 0.5*[0.4660 0.6740 0.1880])

Eforplot=0:0.1:10;

F_ebt3 = fit(E_EBT3, LET_EBT3, 'exp2', 'Exclude', isnan(E_EBT3)| isnan(LET_EBT3))%, 'Lower', [0 -10 1 -10], 'Upper', [1e6 0 1 0])
plot(Eforplot, F_ebt3(Eforplot), '-', 'Color', [0, 0.4470, 0.7410]);

F_ebt2 = fit(E_EBT2, LET_EBT2, 'exp2', 'Exclude', isnan(E_EBT2)| isnan(LET_EBT2))%, 'Lower', [0 -10 1 -10], 'Upper', [1e6 0 1 0])
plot(Eforplot, F_ebt2(Eforplot), '-', 'Color', [0.8500, 0.3250, 0.0980]);

F_ebt3unl = fit(E_EBT3unl, LET_EBT3unl, 'exp2', 'Exclude', isnan(E_EBT3unl) | isnan(LET_EBT3unl))% 'Lower', [0 -10 1 -10], 'Upper', [1e6 0 1 0])
plot(Eforplot, F_ebt3unl(Eforplot), '-', 'Color', [0.4660 0.6740 0.1880]);
grid on
xlabel('E on film surface (MeV)')
set(gca, 'FontSize', 14);
ylabel('LET in active layer (keV/µm)')
legend('EBT3','EBT2', 'EBT3 unlaminated')
ylim([0 70]);

% Plot LET in medium VS LET in water
subplot(1,2,2);
ratios3 = LET_EBT3./LETw_EBT3;
[~, o3] = sort(LET_EBT3);
plot(LET_EBT3(o3),ratios3(o3), '-', 'Color', [0, 0.4470, 0.7410], 'LineWidth', 2);
hold on
ratios2 = LET_EBT2./LETw_EBT2;
[~, o2] = sort(LET_EBT2);
plot(LET_EBT2(o2),ratios2(o2), '-', 'Color', [0.8500, 0.3250, 0.0980], 'LineWidth', 2); 
ratios3u = LET_EBT3unl./LETw_EBT3unl;
[~, o3u] = sort(LET_EBT3unl);
plot(LET_EBT3unl(o3u),ratios3u(o3u), '-', 'Color', 0.5*[0.4660 0.6740 0.1880], 'LineWidth', 2);
grid on
xlabel('LET in medium (keV/µm)')
ylabel('LET in medium / LET in water')
set(gca, 'FontSize', 14);
xlim([0 70])
legend('EBT3','EBT2', 'EBT3 unlaminated')

%% Confidence interval
pctValue = 68.2;
confI3 = confint(F_ebt3,pctValue/100)
dI3 = 0.5*(confI3(2,:) - confI3(1,:))
confI2 = confint(F_ebt2,pctValue/100)
dI2 = 0.5*(confI2(2,:) - confI2(1,:))
confI3u = confint(F_ebt3unl,pctValue/100)
dI3u = 0.5*(confI3u(2,:) - confI3u(1,:))