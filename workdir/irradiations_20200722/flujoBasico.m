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

E0 = 6.75

%. 1. Caída en el kapton
kaptonWindowThickness_cm = 8e-4; % 8 um
E_temp = energyStoppingPowerKapton(E0, [0:1e-4:kaptonWindowThickness_cm]);
E0_afterK = E_temp(end)

% 2. Caída en el aire
separacionAire_cm = 4.37;

dz_cm = 0.001;
Zval = 0:dz_cm:separacionAire_cm;
energyA = energyStoppingPower(E0_afterK, Zval);
finalE = energyA(end)