clear all

%% MW86 plates
wellDiam = 6.35; % mm
wellRadius_cm = 0.1 * wellDiam/2; % mm
well2wellDist_cm = 0.899; %cm

%% Study 3x3 canvas
dxy = 0.01;
sizeX = 5;
sizeY = 5;
sigmaX = 0.1; % Initially say it has a 1x1mm sigma.
sigmaY = 0.1;

N0 = createGaussProfile(dxy, dxy, sizeX, sizeY, sigmaX, sigmaY);
z = 11;
E0 = 3;
%N0 = createSquarePinhole(dxy, sizeX, sizeY, sigmaX);

%% Create well positions
doseCanvas = createEmptyCG2D(dxy, sizeX, sizeY);
well0 = getWell(doseCanvas, wellRadius_cm, [0 0]);

shifts = well2wellDist_cm*[-1 0 1];
[xs, ys] = meshgrid(shifts, shifts);
XYshifts = [xs(:) ys(:)];

wells = {};
allWells = doseCanvas.copy;
for i=1:9
    wells{i} = well0.copy;
    wells{i}.shift(XYshifts(i, :));
    wells{i}.data = double(wells{i}.data);
    allWells = allWells + wells{i};
end
%allWells.crop([doseCanvas.minX doseCanvas.maxX doseCanvas.minY doseCanvas.maxY]);

%% limits = [-1 1 -1 1];
% pC2pNumber = 1e7/1.602176634;
dz = 0.001;
targetTh = 0.001;
targetSPR = 1;

load(fullfile('TOPAS','polyFlu.mat'), 'polyFlu');
sigma_XY = getSigma(polyFlu, E0, z);

% Create new plan here
plan.E = E0;
plan.Z = z;
plan.I = 1;
plan.codFiltro = 'PP100';
plan.name = 'autoplan1';

% Utilizar shift y crop de CG2D
% Hacer esto a mano;
deltaXY = sqrt(sigma_XY^2+sigmaX^2)*1.5; 
% N = 3;
deltaShifts = [-deltaXY 0 deltaXY];
[deltaX, deltaY] = meshgrid(deltaShifts, deltaShifts);
plan.X = deltaX(:);
plan.Y = deltaY(:);

%plan.X = 0; plan.Y = 0;
%plan.Q = 1;
plan.Q = ones(size(plan.X));
plan.numSpots = numel(plan.X);

% Test code for getDoseFromPlan
dose = getDoseFromPlan(doseCanvas, plan, dz, targetTh, targetSPR, N0);
dose.plotSlice
colorbar
hold on
contour(allWells.getAxisValues('X'),allWells.getAxisValues('Y'),allWells.data)
caxis([0 max(dose.data(:))]);
axis([-1 1 -1 1]);

%% Get stats
dose4 = getStats(wells{4}, dose);
close all
dose5 = getStats(wells{5}, dose);
mean(dose5)
