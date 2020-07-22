%% Comenzar
clear all; close all
FIESTAConfig_win
startShutter

%% Realizar irradiacion
numCiclos = 1;
t_open = 10; % segundos
t_wait = 3; % segundos

irradiatePlanOnlyShutter(Fcfg, COMshutter, numCiclos, t_open, t_wait)

%% Abrir shutter en manual (NO LOGGEA)
t_manual = 0.5;
configureShutter(COMshutter, 't', t_manual);
shutter(COMshutter, 'n', 1);

%% Cálculo de la energía
separacionAire_cm = 8.1;
E0 = 8.22;

dz_cm = 0.001;
Zval = 0:dz_cm:separacionAire_cm;
energyA = energyStoppingPower(E0, Zval);
finalE = energyA(end)