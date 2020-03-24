clear all

minE = 1;
maxE = 10;

%% Voxel dimensions
% # Results for scorer Fluence
% # R in 100 bins of 0.01 cm
% # Phi in 1 bin  of 360 deg
% # Z in 2000 bins of 0.02 cm
% # Fluence ( /mm2 ) : Sum   
% # Edep ( MeV ) : Sum

NR = 100;
NZ = 2000;
dR = 0.01;
dZ = 0.02;
RValues = dR * (1:NR) - dR/2;
ZValues = dZ * (1:NZ) - dZ/2;

sigFlu = {};
sigEdep = {};

for i = minE:maxE
    [sigFlu{i}, sigEdep{i}] = loadSingleEnergy(i);
end

save('allSigmas.mat');