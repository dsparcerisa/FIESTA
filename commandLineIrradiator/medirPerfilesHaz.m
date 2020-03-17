% Utilizar instrumentos ya conectados o abrirlos a mano con openInstruments
clearvars -except COMShutter COMStage stageLimits
%% Abrir la GUI
% stageControlStart(COMStage);
stageLimits = [-10.3 0.7 -9.9 0.7 -8 0]; % Definidos 5 Feb 11:30!

%% Alinear en posición 0 y medir la distancia (ANOTAR!)
%beamExit2WellDistance = input('Measure distance between beam exit and well edge (cm): ');
%beamExit2WellBottomDistance = beamExit2WellDistance + 1.1;

beamExit2WellBottomDistance = 4.4;

%% Alinear manualmente con el aspa y la lámina radioluminiscente;
% disp('Press any key when beam is on');
% Configure_shutter(COMShutter,'t',10);
% Shutter(COMShutter,'n',10);

%% When beam is in position, continue
% readStatus(COMStage, 'N'); % Set zero position here

% Pocillo 0,0 desde el ASPA:
% X_2_poc00 = [9.55 6.3 0]; 

% Pocillo 0 0 al centro de la placa:
poc00_2_PlateCtr = [-0.899*5.5 -0.899*3.5 0];

% ASPA al centro de la placa:
% X_2_PlateCtr = X_2_poc00 + poc00_2_PlateCtr;

%% COMMISSIONING DE PLAN FLASH

%% Cargar el plan de medida de haces:
commissioningPlanFlashPath = 'plans/plan_COM_Flash_feb6.txt';
commissioningPlanFlash = readPlan(commissioningPlanFlashPath);

%% Simular plan y crear radiocrómica de prueba
dxy = 0.01;
sizeX = 30;
sizeY = 30;
doseCanvas = createEmptyCG2D(dxy, sizeX, sizeY);
dz = 0.01;
targetTh = 0.001;
targetSPR = 1;
sigmaXY = 0.1;
N0 = createGaussProfile(dxy, dxy, sizeX, sizeY, sigmaXY, sigmaXY);
factorImuestra = 0.76; % Estimación para 3 MeV
dosePlanFlash = getDoseFromPlan(doseCanvas, commissioningPlanFlash, dz, targetTh, targetSPR, N0, factorImuestra, beamExit2WellBottomDistance);
dosePlanFlash.crop([-5 5 -4 4]);
%% Plot
subplot(2,1,1);
dosePlanFlash.plotSlice
set(gca,'YDir','reverse','XDir','reverse');
colorbar
title('Dose distribution FLASH (Gy)');
subplot(2,1,2);
imshow(simulateRC(flip(dosePlanFlash.data',2)));
title('Expected RC film');
%% Irradiar plan
tic
irradiatePlan(COMStage, COMShutter, commissioningPlanFlash, poc00_2_PlateCtr, stageLimits, beamExit2WellBottomDistance)
toc

%% COMMISSIONING DE PLAN CONV
%% Cargar el plan de medida de haces:
commissioningPlanConvPath = 'plans/plan_COM_Conv_feb6.txt';
commissioningPlanConv = readPlan(commissioningPlanConvPath);

dxy = 0.01;
sizeX = 30;
sizeY = 30;
doseCanvas = createEmptyCG2D(dxy, sizeX, sizeY);
dz = 0.01;
targetTh = 0.001;
targetSPR = 1;
sigmaXY = 0.01; % Agujero de 100 um
N0 = createGaussProfile(dxy, dxy, sizeX, sizeY, sigmaXY, sigmaXY);
factorImuestra = 0.76; % Estimación para 3 MeV
dosePlanConv = getDoseFromPlan(doseCanvas, commissioningPlanConv, dz, targetTh, targetSPR, N0, factorImuestra, beamExit2WellBottomDistance);
dosePlanConv.crop([-5 5 -4 4]);

%% Plot
subplot(2,1,1);
dosePlanConv.plotSlice
set(gca,'YDir','reverse','XDir','reverse');
colorbar
title('Dose distribution CONV (Gy)');
subplot(2,1,2);
imshow(simulateRC(flip(dosePlanConv.data',2)));
title('Expected RC film');

%% Perform plan irradiation (x2!)
tic
irradiatePlan(COMStage, COMShutter, commissioningPlanConv, poc00_2_PlateCtr, stageLimits, beamExit2WellBottomDistance)
toc