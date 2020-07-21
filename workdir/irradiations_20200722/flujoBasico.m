%% Comenzar
clear all; close all
FIESTAConfig_win
startShutter

%% Realizar irradiacion
numCiclos = 3;
t_open = 1; % segundos
t_wait = 3; % segundos

irradiatePlanOnlyShutter(Fcfg, COMshutter, numCiclos, t_open, t_wait)