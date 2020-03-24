clear all

%% Voxel dimensions
% # Results for scorer Fluence
% # R in 100 bins of 0.01 cm
% # Phi in 1 bin  of 360 deg
% # Z in 2000 bins of 0.02 cm
% # Fluence ( /mm2 ) : Sum   
% # Edep ( MeV ) : Sum

NR = 300;
NZ = 2000;
dR = 0.01;
dZ = 0.02;
RValues = dR * (1:NR) - dR/2;
ZValues = dZ * (1:NZ) - dZ/2;

[sigFlu16, sigEdep16, dose16] = loadSingleEnergyAl('16um');
[sigFlu16g, sigEdep16g, dose16g] = loadSingleEnergyAl('AuAlFoil');
[sigFlu32, sigEdep32, dose32] = loadSingleEnergyAl('32um');
[sigFlu0, sigEdep0, dose0] = loadSingleEnergy(3);
%% Plot
figure(11); hold off
plot(ZValues, sigFlu0, 'k'); hold on;
plot(ZValues, sigFlu16, 'r'); 
plot(ZValues, sigFlu16g, 'b'); 
plot(ZValues, sigFlu32,'m')
legend('0', '16', '16+gold', '32');
axis([0 max(ZValues) 0 0.5]);

figure(12); hold off;
plot(ZValues, sigEdep0, 'k'); hold on;
plot(ZValues, sigEdep16, 'r'); 
plot(ZValues, sigEdep16g, 'b'); 
plot(ZValues, sigEdep32, 'm')
legend('0', '16', '16+gold', '32');
axis([0 max(ZValues) 0 0.5]);

figure(13); hold off;
plot(ZValues, dose16, 'k'); hold on
plot(ZValues, dose16g, 'r'); hold on
plot(ZValues, dose32, 'b'); hold on
plot(ZValues, dose0, 'm'); hold on
axis([0 max(ZValues) 0 1.2*max(dose0)]);


save('allSigmasAl.mat');