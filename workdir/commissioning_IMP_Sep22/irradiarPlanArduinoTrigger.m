% Utilizar instrumentos ya conectados o abrirlos a mano con openInstruments
instrreset
clear all; close all
FIESTAConfig_win_sep22
startStage
startTrigger
%% Abrir la GUI
stageControlStart(COMstage)
%% stageControlStart(COMStage);
stageLimits = [-10 0.1 -9 0.1 -0.1 9.1]; % Básicos si está centrada
% La fibra en -7.5 -4.8 9
% La FC en -7.5 -4.3 6.5
%% Alinear en posición 0 y medir la distancia (ANOTAR!)
beamExit2RCdistance = input('Measure distance between beam exit and RC (cm): ');

%% Alinear manualmente con el aspa y la lámina radioluminiscente;
% disp('Press any key when beam is on');
% Configure_shutter(COMshutter,'t',10);
% Shutter(COMshutter,'n',10);

%% When beam is in position, continue
% readStatus(COMStage, 'N'); % Set zero position here

% Pocillo 0,0 desde el ASPA:
% X_2_poc00 = [9.55 6.3 0]; 

%% Cargar el plan de calibración de RC;
%planName = 'plan0.txt';
%planName = 'plan1.txt';
%planName = 'plan2.txt';
planName = 'planGFN.txt';
planPath = fullfile(Fcfg.planPath, planName);
thePlan = readPlan(planPath);
%thePlan.Y = thePlan.Y - 4.5;
%thePlan.Z = thePlan.Z + 4.5;
%% Irradiar plan
tic
shift = [-0.5 -4.5 4.5] 
irradiatePlan_Arduino(Fcfg, COMstage, COMtrigger, thePlan, shift, stageLimits, beamExit2RCdistance)
toc