%% Comenzar
clear all; close all
FIESTAConfig_win
startShutter

%% Realizar irradiacion
numCiclos = 1;
t_open = 0.03 % segundos
t_wait = 30; % segundos

irradiatePlanOnlyShutter(Fcfg, COMshutter, numCiclos, t_open, t_wait)

%% Abrir shutter en manual (NO LOGGEA)
t_manual = 200;
configureShutter(COMshutter, 't', t_manual);
shutter(COMshutter, 'n', 1);

%% Cálculo de la energía
separacionAire_cm = 9.2;
E0 = 8.755;

dz_cm = 0.001;
Zval = 0:dz_cm:separacionAire_cm;
energyA = energyStoppingPower(E0, Zval);
finalE = energyA(end)