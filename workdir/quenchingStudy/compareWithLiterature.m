%clear all; close all;
E3color = [0, 0.4470, 0.7410];
E2color = [0.8500, 0.3250, 0.0980];
E3ucolor = [0.4660 0.6740 0.1880];
E1color = [0.4940 0.1840 0.5560];

clear X Y

% Dataset 1
X{1} = [2.00235004
    2.369110327
2.916202626
3.283817965
4.574264294
5.493026715
6.862609615
8.573672252
11.76694937
15.5999123];
Y{1} = [0.999404018
0.979293486
0.999404018
0.939487072
0.919814226
0.939798059
0.909770478
0.869411197
0.819445854
0.67950143];
tag{1} = 'Kirby et al, 2009'
color{1} = E1color;

% Dataset 2
X{2} = [18.07
16.08
12.62
10.89];
Y{2} = [0.68
0.75
0.72
0.78];
tag{2} = 'Martisíková & Jäkel, 2010'
color{2} = E1color;

% Dataset 3
X{3} = [1.270
1.395
1.516
1.620
1.723
1.831
1.939
2.068
2.183
2.259
2.392
2.470
2.602
2.823
3.048
3.321
3.542
3.707
3.975
4.747
5.326
8.351
8.944
17.186
18.251];
Y{3} = [0.999
0.981
0.992
0.989
0.967
0.962
0.964
0.984
0.952
0.944
0.903
0.989
0.935
0.928
0.908
0.910
0.893
0.882
0.887
0.850
0.858
0.786
0.792
0.733
0.728];
tag{3} = 'Angellier et al, 2011'
color{3} = E2color;

% Dataset 4
X{4} = [11.0082
5.5025
3.5771
2.9155
2.3968
2.3091
2.1262
2.0157
1.9307
1.824
1.6022
1.3816
1.1558
1.0555
0.951
0.7291
0.6325
0.5832
0.5621];
Y{4} = [0.835
0.977
0.953
0.988
0.998
0.999
0.986
0.986
0.956
0.998
0.984
0.969
1.002
0.999
1.013
1.010
1.010
1.022
0.981];
tag{4} = 'Reinhardt et al, 2012 (EBT2)'
color{4} = E2color;

% Dataset 5
X{5} = [12.5034
5.6113
3.0111
2.4784
2.3975
2.235
2.0331
1.9179
1.7272
1.5168
1.2867
1.1746
1.0821
0.9478
0.8077
0.8058
0.7017];
Y{5} = [0.811
0.927
0.961
0.972
0.962
0.957
0.971
0.969
0.999
0.986
0.983
0.999
1.003
1.010
1.000
1.002
0.986];
tag{5} = 'Reinhardt et al, 2012 (EBT3)'
color{5} = E3color;

% Dataset 6
X{6} = [1.17
1.29
1.44
1.64
1.80
1.96
2.21
2.50
2.91
3.61
4.30
5.01
5.96
7.46
9.44
11.94
14.64];
Y{6} = [1.00
1.00
1.01
1.01
1.01
0.99
0.98
0.97
0.96
0.94
0.93
0.91
0.90
0.86
0.82
0.76
0.69];
tag{6} = 'Carnicer et al, 2013'
color{6} = E3color;

% Dataset 7
X{7} = [5.05
5.37
5.70
6.06
6.44
6.85
7.28
7.74
8.23
8.76
9.33
9.96
10.68
11.56
12.74
14.52
17.58
23.37
35.06];
Y{7} = [0.98
0.97
0.97
0.97
0.97
0.96
0.96
0.96
0.95
0.94
0.93
0.92
0.91
0.88
0.85
0.81
0.76
0.68
0.59];
tag{7} = 'Fiorini et al, 2014'
color{7} = E3color;

% Dataset 8
X{8} = [24.70
14.00
13.80
10.90];
Y{8} = [0.71
0.84
0.82
0.80];
tag{8} = 'Castriconi et al, 2017'
color{8} = E3color;

% % Dataset 9
X{9} = [12
16
20
22
28
40
78
96];
Y{9} = [0.904191617
0.859582543
0.819168174
0.796133568
0.743842365
0.67916042
0.545783133
0.478858351];
tag{9} = 'Grilj et al, 2018'
color{9} = E3ucolor;

% Dataset 10
% OLD VALUES
% X{10} = [1.30
% 1.50
% 2.03
% 3.05
% 3.48
% 4.49
% 5.39
% 7.45
% 1.17
% 1.88
% 2.19
% 2.50
% 3.23
% 4.23
% 5.44
% 6.56
% 7.21
% 7.90
% 9.27
% 0.86
% 1.16
% 1.61
% 2.33
% 3.08
% 4.07
% 5.01
% 6.28
% 7.23
% 8.35];
X{10} = [1.30
1.50
2.03
3.05
3.48
4.49
5.87
8.55
1.17
1.88
2.19
2.50
3.23
4.23
5.93
7.38
8.23
9.16
11.18
0.86
1.16
1.61
2.33
3.08
4.07
5.01
7.01
8.26
9.77];

Y{10} = [0.98
0.97
0.96
0.97
0.94
0.94
0.87
0.78
0.99
0.96
0.96
0.95
0.94
0.92
0.90
0.88
0.84
0.81
0.80
1.00
0.99
0.98
0.96
0.94
0.93
0.89
0.84
0.83
0.82];
tag{10} = 'Anderson et al, 2019'
color{10} = E3color;

markers = {'+','o','*','x','s','d','^','v','<','p'};

%% Plot
for i=1:10 
   plot(X{i}, Y{i},  markers{i}, 'Color', color{i}, 'MarkerFaceColor', color{i}); hold on
end

set(gca, 'XScale', 'log');
xlabel('LET (keV/µm)');
ylabel('RE');
ylim([0 1.2]);
xlim([1 100]);
grid on
set(gca, 'FontSize', 14);
 
%% Fit to a single function 1-ax^b
x = vertcat(X{:});
y = vertcat(Y{:});
ft2 = fittype( '1 - a*x^b', 'independent', 'x', 'dependent', 'y' );
F2 = fit(x, y, ft2, 'Lower', [0 0])
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
Svals = 1:1:100;
REi = (1 - ai.*Svals.^bi);
hold on
plot(Svals, F2(Svals), 'r-');
hold on
Ydata_plus = prctile(REi, pctValue);
Ydata_minus = prctile(REi, 100-pctValue);
x_plot =[Svals, fliplr(Svals)];
y_plot=[Ydata_minus, fliplr(Ydata_plus)];

%% Fit excluding Grilj
%% Fit to a single function 1-ax^b
x2 = vertcat(X{1:8},X{10});
y2 = vertcat(Y{1:8},Y{10});
F22 = fit(x2, y2, ft2, 'Lower', [0 0])
a2 = F22.a
b2 = F22.b
pctValue = 68.2;
confI2 = confint(F22,pctValue/100);
da2 = 0.5*(confI2(2,1) - confI2(1,1))
db2 = 0.5*(confI2(2,2) - confI2(1,2))

% Generar líneas de confianza
N_ITER = 10000;
ai2 = a2 + da2*randn(N_ITER,1);
bi2 = b2 + db2*randn(N_ITER,1);
Svals = 1:1:100;
REi2 = (1 - ai2.*Svals.^bi2);
hold on
plot(Svals, F22(Svals), 'k-');
hold on
Ydata_plus2 = prctile(REi2, pctValue);
Ydata_minus2 = prctile(REi2, 100-pctValue);
x_plot2 =[Svals, fliplr(Svals)];
y_plot2=[Ydata_minus2, fliplr(Ydata_plus2)];

%% Fit ONLY EBT2
x3 = vertcat(X{3:4});
y3 = vertcat(Y{3:4});
F23 = fit(x3, y3, ft2, 'Lower', [0 0])
a3 = F23.a
b3 = F23.b
pctValue = 68.2;
confI3 = confint(F23,pctValue/100);
da3 = 0.5*(confI3(2,1) - confI3(1,1))
db3 = 0.5*(confI3(2,2) - confI3(1,2))

% Generar líneas de confianza
N_ITER = 10000;
ai3 = a3 + da3*randn(N_ITER,1);
bi3 = b3 + db3*randn(N_ITER,1);
Svals = 1:1:100;
REi3 = (1 - ai3.*Svals.^bi3);
hold on
plot(Svals, F23(Svals), 'r-');
hold on
Ydata_plus3 = prctile(REi3, pctValue);
Ydata_minus3 = prctile(REi3, 100-pctValue);
x_plot3 =[Svals, fliplr(Svals)];
y_plot3=[Ydata_minus3, fliplr(Ydata_plus3)];

%% Fit ONLY EBT3
x4 = vertcat(X{5:8}, X{10});
y4 = vertcat(Y{5:8}, Y{10});
F24 = fit(x4, y4, ft2, 'Lower', [0 0])
a4 = F24.a
b4 = F24.b
pctValue = 68.2;
confI4 = confint(F24,pctValue/100);
da4 = 0.5*(confI4(2,1) - confI4(1,1))
db4 = 0.5*(confI4(2,2) - confI4(1,2))

% Generar líneas de confianza
N_ITER = 10000;
ai4 = a4 + da4*randn(N_ITER,1);
bi4 = b4 + db4*randn(N_ITER,1);
Svals = 1:1:100;
REi4 = (1 - ai4.*Svals.^bi4);
hold on
plot(Svals, F24(Svals), 'b-');
hold on
Ydata_plus4 = prctile(REi4, pctValue);
Ydata_minus4 = prctile(REi4, 100-pctValue);
x_plot4 =[Svals, fliplr(Svals)];
y_plot4=[Ydata_minus3, fliplr(Ydata_plus4)];

%% Plot confidence intervals
fill(x_plot, y_plot,1,'facecolor', [1 0.8 0.8], 'edgecolor', 'red', 'edgealpha', 0.6, 'facealpha', 0.2, 'LineStyle',':');
fill(x_plot2, y_plot2,1,'facecolor', [0.7 0.7 0.7], 'edgecolor', 'black', 'edgealpha', 0.6, 'facealpha', 0.2, 'LineStyle',':');
fill(x_plot3, y_plot3,1,'facecolor', [0.7 0 0], 'edgecolor', 'black', 'edgealpha', 0.6, 'facealpha', 0.2, 'LineStyle',':');
fill(x_plot4, y_plot4,1,'facecolor', [0.0 0 0.8], 'edgecolor', 'black', 'edgealpha', 0.6, 'facealpha', 0.2, 'LineStyle',':');

% Leyenda
tag{11} = 'Fit';
tag{12} = 'Fit excluding Grilj (2018) data';
tag{13} = 'Fit of EBT2 data';

legend(tag, 'Location', 'SouthWest')
title('Detail of values from the literature');
save('literaturefit.mat','a','da','b','db','X','Y', 'a2','da2','b2','db2','x2','y2','markers');