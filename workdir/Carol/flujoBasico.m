%% Comenzar
clear all; close all
FIESTAConfig_win
instrreset
startStage

%% Abrir el controlador de la stage en modo manual
stageControlStart(COMstage)

%% Leer plan
planName = 'pruebaErrores.txt'; % Cambiar aqui
planPath = [Fcfg.planPath filesep planName];
thePlan = readPlan(planPath);

%% Irradiar plan
maxIterations = 1; % NUMERO DE REPETICIONES
vector2startPoint = [0 0 0]; % para mover todo el plan un vector constante
stageLimits = [-11 11 -11 11 -11 11]; % corregir antes de irradiar para evitar fallos de movimiento
irradiatePlanOnlyStage(Fcfg, COMstage, thePlan, vector2startPoint, stageLimits, maxIterations)
