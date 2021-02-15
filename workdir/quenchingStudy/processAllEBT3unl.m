clear all; close all
radNames = {'Rad1','Rad2','Rad3','Rad4','Rad6','Rad9','Rad11','Rad13','Rad16'};
allSw_preFilm = [];
alldSw_preFilm = [];
allratios = [];
alldratios = [];

radius_cm = 0.03;

for i = 1:numel(radNames)
    radName = radNames{i};
    [Sw_preFilm, dSw_preFilm, ratios, dratios] = processAnyRad(radName, radius_cm);
    allSw_preFilm = [allSw_preFilm, Sw_preFilm];
    alldSw_preFilm = [alldSw_preFilm, dSw_preFilm];
    allratios = [allratios, ratios];
    alldratios = [alldratios, dratios];    
end

allSw_preFilm = allSw_preFilm / 10;
alldSw_preFilm = alldSw_preFilm / 10;

%% Include uncertainties for fit
N_ITER = 100;
ft2 = fittype( '1 - a*x^b', 'independent', 'x', 'dependent', 'y' );

a_i = nan(N_ITER,1);
b_i = nan(N_ITER,1);

for i=1:N_ITER
    allSw_preFilm_i = allSw_preFilm + alldSw_preFilm.*randn(1,numel(allSw_preFilm));
    allratios_i = allratios + alldratios.*randn(1,numel(allSw_preFilm));
    F_i = fit(allSw_preFilm_i', allratios_i', ft2, 'Lower', [0 0])
    %plot(allSw_preFilm_i, allratios_i, '.'); hold on
    a_i(i) = F_i.a;
    b_i(i) = F_i.b;
end

%% Plots
errorbar(allSw_preFilm, allratios, alldratios, alldratios, alldSw_preFilm, alldSw_preFilm, 'k.');
set(gca, 'XScale', 'log');
xlabel('LET (keV/um)');
ylabel('RE');
ylim([0 1.2]);
xlim([1 100]);
grid on

% hold on
% wgt = alldratios.^(-2);
% wgt(isnan(wgt))=0;
% allratios(isnan(allratios))=0;
% allSw_preFilm(isnan(allSw_preFilm))=0;
% %ft = fittype( 'D0 / (a*x + D0)', 'independent', 'x', 'dependent', 'y' );
% ft2 = fittype( '1 - a*x^b', 'independent', 'x', 'dependent', 'y' );
% 
% %F = fit(allSw_preFilm', allratios', ft, 'Weights', wgt', 'Lower', [20 0], 'Upper', [80 1])
%F2 = fit(allSw_preFilm', allratios', ft2, 'Lower', [0 0])
% a = F2.a
% b = F2.b

a = mean(a_i);
da = std(a_i);
b = mean(b_i);
db = std(b_i);

%F2.a = a; F2.b = b;

% confI = confint(F2,0.682);
% da = 0.5*(confI(2,1) - confI(1,1))
% db = 0.5*(confI(2,2) - confI(1,2))

% Generar líneas de confianza
N_ITER = 10000;
ai = a + da*randn(N_ITER,1);
bi = b + db*randn(N_ITER,1);
Svals = 1:1:80;
REi = (1 - ai.*Svals.^bi);

hold on
%plot(Svals, F(Svals), 'm-');
plot(Svals, 1 - a.*Svals.^b, 'r-');
hold on
plot(Svals, prctile(REi, 68.2), 'r:');
%plot(Svals, prctile(REi, 50), 'k-');
plot(Svals, prctile(REi, 100-68.2), 'r:');

xlabel('LET (keV/µm)');
ylabel('RE');
set(gca,'xscale','log')
axis([5 80 0 1.20]);
title('EBT3 unlaminated');

%% Save
save('finalResults_EBT3unl.mat');