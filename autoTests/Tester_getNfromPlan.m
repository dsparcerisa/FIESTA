clear all

dxy = 0.01;
sizeX = 10;
sizeY = 10;
sigmaX = 0.01;
sigmaY = 0.01;

% limits = [-1 1 -1 1];
pC2pNumber = 1e7/1.602176634;
N0 = createGaussProfile(dxy, dxy, sizeX, sizeY, sigmaX, sigmaY);
plan = readPlan('fullplan.txt');

% Utilizar shift y crop de CG2D
doseCanvas = createEmptyCG2D(dxy, sizeX, sizeY);


% Target info
z = 8;
E0 = 3;
dz = 0.001;
targetTh = 0.001;
targetSPR = 1;

% Test code for getDoseFromPlan

dose = getDoseFromPlan(doseCanvas, plan, E0, z, dz, targetTh, targetSPR, N0);
dose.plotSlice
colorbar

