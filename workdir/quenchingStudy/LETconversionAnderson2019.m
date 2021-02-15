clear all; close all

LETa_3 = 4.08E05;
LETb_3 = 2.88;
LETc_3 = 22.5;
LETd_3 = 0.142;

LETa_u = 94.8;
LETb_u = 1.56;
LETc_u = 17.1;
LETd_u = 0.122;

LET_Es_3 = @(Es) LETa_3 .* exp(-LETb_3.*Es) + LETc_3 .* exp(-LETd_3.*Es)
LET_Es_u = @(Es) LETa_u .* exp(-LETb_u.*Es) + LETc_u .* exp(-LETd_u.*Es)

E = 0:0.1:10;
LET3 = LET_Es_3(E);
LET3u = LET_Es_u(E);

plot(LET3u,LET3)

LETproblem = [5.39
7.45
1.17
1.88
2.19
2.50
3.23
4.23
5.44
6.56
7.21
7.90
9.27
0.86
1.16
1.61
2.33
3.08
4.07
5.01
6.28
7.23
8.35];
LETproblem = LETproblem(LETproblem>5.3);
LETsolution = interp1(LET3u, LET3, LETproblem)