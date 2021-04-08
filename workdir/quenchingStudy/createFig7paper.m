clear all; close all;

set(gcf, 'Position',  [100, 100, 800, 400])
ha = tight_subplot(1,2,0.03,[0.12 0.05],[0.06 0.02])
axes(ha(2))
%% Load all
EBT2 = load('finalResults_EBT2.mat');
EBT3 = load('finalResults_EBT3.mat');
EBT3unl = load('finalResults_EBT3unl.mat');

%% Plot all three together
errorbar(EBT3.allSw_preFilm, EBT3.allratios, EBT3.alldratios, EBT3.alldratios, EBT3.alldSw_preFilm, EBT3.alldSw_preFilm, '.', 'Color', [0, 0.4470, 0.7410]);
hold on
errorbar(EBT2.allSw_preFilm, EBT2.allratios, EBT2.alldratios, EBT2.alldratios, EBT2.alldSw_preFilm, EBT2.alldSw_preFilm, '.', 'Color', [0.8500, 0.3250, 0.0980]); 
errorbar(EBT3unl.allSw_preFilm, EBT3unl.allratios, EBT3unl.alldratios, EBT3unl.alldratios, EBT3unl.alldSw_preFilm, EBT3unl.alldSw_preFilm, '.', 'Color', [0.4660 0.6740 0.1880])
plot([0],[0],'o', 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', [0.7 0.7 0.7]);
%plot(EBT2.Svals, EBT2.F2(EBT2.Svals), 'k-');
%plot(EBT2.Svals, EBT3.F2(EBT2.Svals), 'r-');
%plot(EBT2.Svals, EBT3unl.F2(EBT2.Svals), 'b-');

%% Fit to a single function 1-ax^b
% figure(2)
X = [EBT3.allSw_preFilm'; EBT2.allSw_preFilm'; EBT3unl.allSw_preFilm'];
dX = [EBT3.alldSw_preFilm'; EBT2.alldSw_preFilm'; EBT3unl.alldSw_preFilm'];
Y = [EBT3.allratios'; EBT2.allratios'; EBT3unl.allratios'];
dY = [EBT3.alldratios'; EBT2.alldratios'; EBT3unl.alldratios'];
% errorbar(X,Y,dY,dY,dX,dX, '.k')

wgt = dY.^(-2);
ft2 = fittype( '1 - a*x^b', 'independent', 'x', 'dependent', 'y' );
F2 = fit(X, Y, ft2, 'Weights', wgt', 'Lower', [0 0])
a = F2.a
b = F2.b
pctValue = 68.2;
confI = confint(F2,pctValue/100);
da = 0.5*(confI(2,1) - confI(1,1))
db = 0.5*(confI(2,2) - confI(1,2))

% Generar líneas de confianza
N_ITER = 10000;
ai = a + da*randn(N_ITER,1);
bi = b + db*randn(N_ITER,1);
Svals = 1:1:80;
REi = (1 - ai.*Svals.^bi);
hold on
plot(Svals, F2(Svals), 'b-');
hold on
Ydata_plus = prctile(REi, pctValue);
Ydata_minus = prctile(REi, 100-pctValue);
x_plot =[Svals, fliplr(Svals)];
y_plot=[Ydata_minus, fliplr(Ydata_plus)];


%% Include literature fit
LF = load('literaturefit.mat');
LFai = LF.a2 + LF.da2*randn(N_ITER,1);
LFbi = LF.b2 + LF.db2*randn(N_ITER,1);
LFSvals = 1:1:100;
LFREi = (1 - LFai.*LFSvals.^LFbi);
hold on
plot(LFSvals, 1-LF.a2.*LFSvals.^LF.b2, 'k-');
hold on
LFYdata_plus = prctile(LFREi, pctValue);
LFYdata_minus = prctile(LFREi, 100-pctValue);
LFx_plot =[LFSvals, fliplr(LFSvals)];
LFy_plot=[LFYdata_minus, fliplr(LFYdata_plus)];
fill(LFx_plot, LFy_plot,1,'facecolor', [0.7 0.7 0.7], 'edgecolor', 'black', 'edgealpha', 1, 'facealpha', 0.4, 'LineStyle',':');
fill(x_plot, y_plot,1,'facecolor', [0.6 0.6 1], 'edgecolor', 'blue', 'edgealpha', 1, 'facealpha', 0.4, 'LineStyle',':');

%% Hacer fit combinado

X = [EBT3.allSw_preFilm'; EBT2.allSw_preFilm'; EBT3unl.allSw_preFilm'];
Y = [EBT3.allratios'; EBT2.allratios'; EBT3unl.allratios'];

combX = [X; vertcat(LF.X{1:8}); LF.X{10}];
combY = [Y; vertcat(LF.Y{1:8}); LF.Y{10}];
%FC = fit(combX, combY, ft2, 'Lower', [0 1], 'Upper',[1 1]) LINEAR FIT
FC = fit(combX, combY, ft2)

aC = FC.a
bC = FC.b
confIC = confint(FC,pctValue/100);
daC = 0.5*(confIC(2,1) - confIC(1,1))
dbC = 0.5*(confIC(2,2) - confIC(1,2))
% 
% % Generar líneas de confianza
N_ITER = 10000;
aiC = aC + daC*randn(N_ITER,1);
biC = bC + dbC*randn(N_ITER,1);
Svals = 1:1:100;
REiC = (1 - aiC.*Svals.^biC);
hold on
% plot(Svals, FC(Svals), 'g-');
hold on
Ydata_plusC = prctile(REiC, pctValue);
Ydata_minusC = prctile(REiC, 100-pctValue);
x_plotC =[Svals, fliplr(Svals)];
y_plotC=[Ydata_minusC, fliplr(Ydata_plusC)];
% fill(x_plotC, y_plotC,1,'facecolor', [0.8 1 0.8], 'edgecolor', 'green', 'edgealpha', 0.6, 'facealpha', 0.2, 'LineStyle',':');

%% Plot literature data in grey
for i=1:10 
   plot(LF.X{i}, LF.Y{i},  LF.markers{i}, 'Color', [0.7 0.7 0.7], 'MarkerFaceColor', [0.7 0.7 0.7]); hold on
end

errorbar(EBT3.allSw_preFilm, EBT3.allratios, EBT3.alldratios, EBT3.alldratios, EBT3.alldSw_preFilm, EBT3.alldSw_preFilm, '.', 'Color', [0, 0.4470, 0.7410], 'MarkerSize', 10);
errorbar(EBT2.allSw_preFilm, EBT2.allratios, EBT2.alldratios, EBT2.alldratios, EBT2.alldSw_preFilm, EBT2.alldSw_preFilm, '.', 'Color', [0.8500, 0.3250, 0.0980], 'MarkerSize', 10); 
errorbar(EBT3unl.allSw_preFilm, EBT3unl.allratios, EBT3unl.alldratios, EBT3unl.alldratios, EBT3unl.alldSw_preFilm, EBT3unl.alldSw_preFilm, '.', 'Color', [0.4660 0.6740 0.1880], 'MarkerSize', 10);
% plot(Svals, FC(Svals), 'b-');
plot(LFSvals, 1-LF.a2.*LFSvals.^LF.b2, 'k-');
plot(Svals, F2(Svals), 'b-');

%% Retoques finales
legend({'EBT3','EBT2','EBT3unlaminated', 'Literature values', 'Data fit', 'Literature fit excluding Grilj (2018) data'}, 'Location', 'SouthWest')
set(gca, 'XScale', 'log');
xlabel('LET (keV/µm)');
ylim([0 1.4]);
xlim([1 100]);
grid on
title('Measured values');
set(gca, 'FontSize', 14);
xticks([1 10 100]);
xticklabels({'1','10','100'});
yticklabels([])
axes(ha(1))
compareWithLiterature
xticks([1 10 100]);
ylim([0 1.4]);
yticks([0 0.2 0.4 0.6 0.8 1 1.2 1.4]);
xticklabels({'1','10','100'});
ylabel('RE');
