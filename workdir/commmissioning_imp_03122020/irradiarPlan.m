% Utilizar instrumentos ya conectados o abrirlos a mano con openInstruments

%% Abrir la GUI
% stageControlStart(COMStage);
% X_2_poc00 = [9.55 6.3 0]; 

stageLimits = [-10.3 0.3 -9.93 0.7 -7.5 -1]; % DEFINE!!


%% Alinear en posición 0 y medir la distancia (ANOTAR!)

% PLATES
if (modo == 1)
    %beamExit2WellDistance = input('Measure distance between beam exit and well edge (cm): ');
    beamExit2WellDistance = 0.6 + 2.7 + 0.2; % 2mm - distancia 6 feb
    beamExit2WellBottomDistance = beamExit2WellDistance + 1.1;
else
% CUBETAS
    % distanciaAFondoCubeta = input('Distancia de la salida del haz al FONDO de la cubeta (cm): ');
    distanciaAFondoCubeta = 0.6 + 2.7 + 0.5 + 0.8;
end

ConvPlanPath = 'plans/plan_CONV_feb6_V2_placa.txt';

ConvPlan = readPlan(ConvPlanPath);

fprintf('INTENDED INTENSITY (FC1): %5.3f nA\n', ConvPlan.I);

%% Irradiar plan
tic
if (modo==1) % Plate
    shiftPocilloaFondo = [0 0 0.3];
    irradiatePlan(COMStage, COMShutter, ConvPlan, shiftPocilloaFondo, stageLimits, beamExit2WellBottomDistance) % placa 
else
    shiftCubetas = [0 -0.3 0];
    irradiatePlan(COMStage, COMShutter, ConvPlan, shiftCubetas, stageLimits, distanciaAFondoCubeta) % cubeta
end
toc