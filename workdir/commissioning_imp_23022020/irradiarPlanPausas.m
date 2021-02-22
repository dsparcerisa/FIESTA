%% Comenzar
clear all; close all
FIESTAConfig_win_23feb2021
instrreset

% Establecer la velocidad de movimiento de la stage (0-6000)
% ANTES de llamar a startStage. Si no la encuentra, utilizará el valor
% máximo de 6000.
varSpeed = 6000; 
startStage


%% Abrir el controlador de la stage en modo manual
stageControlStart(COMstage)

%% Leer plan
planName = 'planW.txt'; % Cambiar aqui
planPath = [Fcfg.planPath filesep planName];
thePlan = readPlan(planPath);

%% Irradiar plan
maxIterations = 50; % NUMERO DE REPETICIONES (SOLO SE CAMBIA AQUI, NO ESTA EN EL PLAN)
vector2startPoint = [0 0 0]; % para mover todo el plan un vector constante
stageLimits = [-11 11 -11 11 -11 11]; % corregir antes de irradiar para evitar fallos de movimiento
irradiatePlanOnlyStage(Fcfg, COMstage, thePlan, vector2startPoint, stageLimits, maxIterations)

