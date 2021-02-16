% Utilizar instrumentos ya conectados o abrirlos a mano con openInstruments
clearvars -except COMshutter COMstage stageLimits Fcfg
%% Abrir la GUI
% stageControlStart(COMStage);
stageLimits = [-5.2 5.2 -5.2 5.2 -10.5 0]; % Definidos 5 Feb 11:30!

beamExit2WellBottomDistance = 4.4;

%% Alinear manualmente con el aspa y la lámina radioluminiscente;
% disp('Press any key when beam is on');
% Configure_shutter(COMShutter,'t',10);
% Shutter(COMShutter,'n',10);

%% When beam is in position, continue
% readStatus(COMStage, 'N'); % Set zero position here

% Pocillo 0 0 al centro de la placa:
poc00_2_PlateCtr = [-0.899*5.5 -0.899*3.5 0];

% ASPA al centro de la placa:
% X_2_PlateCtr = X_2_poc00 + poc00_2_PlateCtr;

%% COMMISSIONING DE PLAN FLASH

%% Cargar el plan de medida de haces:
commissioningPlanFlashPath = fullfile(Fcfg.planPath, 'plan_COM_Conv_feb6.txt');
%commissioningPlanFlashPath = fullfile(Fcfg.planPath, 'medioPlan.txt');
commissioningPlanFlash = readPlan(commissioningPlanFlashPath);

%% Irradiar plan
tic
irradiatePlan(Fcfg, COMstage, COMshutter, commissioningPlanFlash, [0 1.5 0], stageLimits, 0)
toc