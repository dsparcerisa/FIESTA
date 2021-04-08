clear all; close all;


%% Load all
EBT2 = load('finalResults_EBT2.mat');
EBT3 = load('finalResults_EBT3.mat');
EBT3unl = load('finalResults_EBT3unl.mat');

color3 = [0, 0.4470, 0.7410];
color2 = [0.8500, 0.3250, 0.0980]; 
color3u = [0.4660 0.6740 0.1880];
%% Plot EBT3
pctValue = 68.2;

set(gcf, 'Position',  [100, 100, 800, 800])
ha = tight_subplot(2,2,0.04,[0.07 0.05],[0.06 0.02])
axes(ha(1))

axes(ha(1))
errorbar(EBT3.allSw_preFilm, EBT3.allratios, EBT3.alldratios, EBT3.alldratios, EBT3.alldSw_preFilm, EBT3.alldSw_preFilm, 'k.', 'MarkerSize', 10);

hold on
Y3data_plus = prctile(EBT3.REi, pctValue);
Y3data_minus = prctile(EBT3.REi, 100-pctValue);
x3_plot =[EBT3.Svals, fliplr(EBT3.Svals)];
y3_plot=[Y3data_minus, fliplr(Y3data_plus)];
fill(x3_plot, y3_plot,1,'facecolor', color3, 'edgecolor', color3, 'edgealpha', 1, 'facealpha', 0.1, 'LineStyle',':');
plot(EBT3.Svals, 1 - EBT3.a.*EBT3.Svals.^EBT3.b, '-', 'Color', color3, 'LineWidth', 2);

set(gca, 'XScale', 'log');
ylabel('RE');
ylim([0 1.4]);
xlim([5 80]);
xticks([5 10 20 40 80])
xticklabels([]);
grid on
title('EBT3');
set(gca, 'FontSize', 14);

axes(ha(2))
errorbar(EBT2.allSw_preFilm, EBT2.allratios, EBT2.alldratios, EBT2.alldratios, EBT2.alldSw_preFilm, EBT2.alldSw_preFilm, 'k.', 'MarkerSize', 10);

hold on
Y2data_plus = prctile(EBT2.REi, pctValue);
Y2data_minus = prctile(EBT2.REi, 100-pctValue);
x2_plot =[EBT2.Svals, fliplr(EBT2.Svals)];
y2_plot=[Y2data_minus, fliplr(Y2data_plus)];
fill(x2_plot, y2_plot,1,'facecolor', color2, 'edgecolor', color2, 'edgealpha', 1, 'facealpha', 0.1, 'LineStyle',':');
plot(EBT2.Svals, 1 - EBT2.a.*EBT2.Svals.^EBT2.b, '-', 'Color', color2, 'LineWidth', 2);

set(gca, 'XScale', 'log');
ylim([0 1.4]);
xlim([5 80]);
xticks([5 10 20 40 80])
xticklabels([]);
yticklabels([]);
grid on
title('EBT2');
set(gca, 'FontSize', 14);

axes(ha(3))
errorbar(EBT3unl.allSw_preFilm, EBT3unl.allratios, EBT3unl.alldratios, EBT3unl.alldratios, EBT3unl.alldSw_preFilm, EBT3unl.alldSw_preFilm, 'k.', 'MarkerSize', 10);

hold on
Y3udata_plus = prctile(EBT3unl.REi, pctValue);
Y3udata_minus = prctile(EBT3unl.REi, 100-pctValue);
x3u_plot =[EBT3unl.Svals, fliplr(EBT3unl.Svals)];
y3u_plot=[Y3udata_minus, fliplr(Y3udata_plus)];
fill(x3u_plot, y3u_plot,1,'facecolor', color3u, 'edgecolor', color3u, 'edgealpha', 1, 'facealpha', 0.1, 'LineStyle',':');
plot(EBT3unl.Svals, 1 - EBT3unl.a.*EBT3unl.Svals.^EBT3unl.b, '-', 'Color', color3u, 'LineWidth', 2);

set(gca, 'XScale', 'log');
xlabel('LET (keV/µm)');
ylabel('RE');
ylim([0 1.4]);
xlim([5 80]);
xticks([5 10 20 40 80])
title('EBT3 unlaminated');
grid on
set(gca, 'FontSize', 14);

axes(ha(4))
X = [EBT3.allSw_preFilm'; EBT2.allSw_preFilm'; EBT3unl.allSw_preFilm'];
Y = [EBT3.allratios'; EBT2.allratios'; EBT3unl.allratios'];

plot(X,Y,'.','Color',[0.7 0.7 0.7], 'MarkerSize', 10);
hold on
fill(x3_plot, y3_plot,1,'facecolor', color3, 'edgecolor', color3, 'edgealpha', 1, 'facealpha', 0.1, 'LineStyle',':');
fill(x2_plot, y2_plot,1,'facecolor', color2, 'edgecolor', color2, 'edgealpha', 1, 'facealpha', 0.1, 'LineStyle',':');
fill(x3u_plot, y3u_plot,1,'facecolor', color3u, 'edgecolor', color3u, 'edgealpha', 1, 'facealpha', 0.1, 'LineStyle',':');
plot(EBT3.Svals, 1 - EBT3.a.*EBT3.Svals.^EBT3.b, '-', 'Color', color3, 'LineWidth', 2);
plot(EBT2.Svals, 1 - EBT2.a.*EBT2.Svals.^EBT2.b, '-', 'Color', color2, 'LineWidth', 2);
plot(EBT3unl.Svals, 1 - EBT3unl.a.*EBT3unl.Svals.^EBT3unl.b, '-', 'Color', color3u, 'LineWidth', 2);


set(gca, 'XScale', 'log');
xlabel('LET (keV/µm)');
ylim([0 1.4]);
xlim([5 80]);
xticks([5 10 20 40 80])
xticklabels('auto')
yticklabels([]);
grid on
title('Fit comparison');
set(gca, 'FontSize', 14);

%% Retoques finales
% legend({'EBT3','EBT2','EBT3unlaminated', 'Literature values', 'Data fit', 'Literature fit excluding Grilj (2018) data'}, 'Location', 'SouthWest')
% set(gca, 'XScale', 'log');
% xlabel('LET (keV/µm)');
% ylabel('RE');
% ylim([0 1.2]);
% xlim([1 100]);
% grid on
% title('Measured values');
% 
% axes(ha(2))
% %compareWithLiterature
% %yticklabels([])
% %ylabel('');

%% ANnotations

axes(ha(4))
E2_X = [0.6 0.6];
E2_Y = [0.26 0.33];
annotation('textarrow',E2_X,E2_Y,'String','EBT2', 'FontSize', 14)

E3_X = [0.85 0.92];
E3_Y = [0.12 0.12];
annotation('textarrow',E3_X,E3_Y,'String','EBT3', 'FontSize', 14)

E3u_X = [0.87 0.87];
E3u_Y = [0.33 0.24];
annotation('textarrow',E3u_X,E3u_Y,'String','EBT3 unlaminated', 'FontSize', 14)