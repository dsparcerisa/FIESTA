% Add path to consider s
% stageControlPath = '../stageControlGUI';
% addpath(stageControlPath);

% Start instruments
startInstruments

%% planPath
% Read plan
plan = readPlan('plan_COM_Conv.txt');
plan = readPlan('plan_COM_Flash.txt');

%% Data necessary for irradiation
% PPfactor = 5e-4;
% 
% Ifactor = 0.81;
% 
% I_FC1_nA = 2; 
% I_muestra_nA = I_FC1_nA * Ifactor * PPfactor;
% 
% plan.t_s = plan.Q ./ I_muestra_nA / 1000;

%% Center stage
%autoCenter(s2);   
%[~, ~, ~, Zpos] = monitorStatus(s2);

%% Align in Z position
% align manually and then set here (EN EL ASPA) y lo más cerca posible
%Zpos = 0;
%stageControlStart;

%% Calculo del vector global
% Aspa al Pocillo 0 0
%deltaX = +95.5 mm
%deltaY = +63 mm
X2poc00 = [9.55 6.3 0];

% Pocillo 0 0 al centro de la placa:
poc00toplateCtr = [-0.899*5.5 -0.899*3.5 0];

globalVector = X2poc00 + poc00toplateCtr
%% Irradiate plan for 
tic
irradiatePlan(COMstage, COMshutter, plan, globalVector)
toc