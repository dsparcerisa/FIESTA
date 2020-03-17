% Utilizar instrumentos ya conectados o abrirlos a mano con openInstruments

%% Abrir la GUI
% stageControlStart(COMStage);
% X_2_poc00 = [9.55 6.3 0]; 

stageLimits = [-10.3 0.3 -9.93 0.7 -7.5 -1]; % DEFINE!!

%% Alinear en posición 0 y medir la distancia (ANOTAR!)
distNarizAGelatina_cm = 12;

GelPath = 'plans/plan_gelatina_feb17_prueba.txt'; 

GelPlan = readPlan(GelPath);

fprintf('INTENDED INTENSITY (FC1): %5.3f nA\n', GelPlan.I);

%% Irradiar plan
shiftGelatina = [0 0 0];
irradiatePlan(COMStage, COMShutter, GelPlan, shiftGelatina, stageLimits, distNarizAGelatina_cm) % gelatina
