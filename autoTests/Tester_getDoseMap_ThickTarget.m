% function doseMap = getDoseMap_ThickTarget(E0, z, Nprot, targetTh, materialList, N0, sigmaPoly)



% Create material list
Material = {'Air'; 'Copper'; 'Air'; 'Water'};
Thickness_cm = [4; 0.0020; 8; 0.0030];
materialTable = table(Material, Thickness_cm);

sigmaPoly_ct = [0 0 4]; % Valor en mm
E0 = 8;
Nprot = 6.25e6 * 400;

dxy = 0.01;
sizeX = 30;
sizeY = 30;
sigmaX = 0.001; % REVISAR
sigmaY = 0.001; % REVISAR
N0 = createGaussProfile(dxy, dxy, sizeX, sizeY, sigmaX, sigmaY);

doseMap = getDoseMap_ThickTarget(E0, Nprot, materialTable, N0, sigmaPoly_ct);
doseMap.plotSlice
colorbar