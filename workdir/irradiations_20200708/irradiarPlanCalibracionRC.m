% Utilizar instrumentos ya conectados o abrirlos a mano con openInstruments

%% Abrir la GUI
%% stageControlStart(COMStage);

stageLimits = [-10.3 0.3 -9.9 0.7 -8 0]; % Definidos 5 Feb 11:30!

%% Alinear en posición 0 y medir la distancia (ANOTAR!)
beamExit2RCdistance = input('Measure distance between beam exit and RC(cm): ');

%% Alinear manualmente con el aspa y la lámina radioluminiscente;
% disp('Press any key when beam is on');
% Configure_shutter(COMShutter,'t',10);
% Shutter(COMShutter,'n',10);

%% When beam is in position, continue
% readStatus(COMStage, 'N'); % Set zero position here

% Pocillo 0,0 desde el ASPA:
% X_2_poc00 = [9.55 6.3 0]; 

%% Cargar el plan de calibración de RC;
RCCalibrationPlanPath = 'plans/plan_calibrationRC.txt';
%RCCalibrationPlanPath = 'plans/plan_calibrationRCx10.txt';
%RCCalibrationPlanPath = 'plans/plan_calibrationRC_media.txt';
RCCalibrationPlan = readPlan(RCCalibrationPlanPath);

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
dosePlanRC = getDoseFromPlan(doseCanvas, RCCalibrationPlan, dz, targetTh, targetSPR, N0, factorImuestra, beamExit2RCdistance);
dosePlanRC.crop([-11 1 -8 1]);
%% Plot
subplot(2,1,1);
dosePlanRC.plotSlice
set(gca,'YDir','reverse','XDir','reverse');
colorbar
title('Dose distribution (Gy)');
subplot(2,1,2);
imshow(simulateRC(flip(dosePlanRC.data',2)));
title('Expected RC film');

%% Irradiar plan
tic
irradiatePlan(COMStage, COMShutter, RCCalibrationPlan, [0 0 0], stageLimits, beamExit2RCdistance)
toc