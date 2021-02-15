clear all; close all;

%% Create variables
E0=1:0.1:10;
airVecPos = 0:0.1:20;
finalE = nan(numel(E0),numel(airVecPos));

%% Calculate finalE
for iE = 1:numel(E0)
    
kaptonThickness_um = 8;
kaptonThickness_cm = 1e-4*kaptonThickness_um;
E0kap_vector = energyStoppingPowerKapton(E0(iE), [0 kaptonThickness_cm]);
E0kap = E0kap_vector(2);
Eair_vec = energyStoppingPower(E0kap, airVecPos);
finalE(iE,:) = Eair_vec;
end

%% Plot finalE
figure(1)
subplot(1,2,1);
imagesc(airVecPos,E0,finalE)
set(gca, 'YDir','normal')
subplot(1,2,2);
plot(airVecPos, finalE')
ylim([0 10]);

%% Calculate for real E
airVecPos2 = 0:0.1:120;
Eair_vec2 = max(0,energyStoppingPower(10, airVecPos2));
%plot(airVecPos2, Eair_vec2)
vM = Eair_vec2>0;
Z = airVecPos2(vM);
E = Eair_vec2(vM);
Rmax = Z(end);
R = Rmax - Z;
R(1) = eps;
F1 = fit(E', R', 'Power1', 'Robust', 'LAR')
a = F1.a;
p = F1.b;
Efun = @(E0,z) real((E0.^p - z./a).^(1/p))

%% Calculate fData
fData = Efun(E0', airVecPos); 

%% Plot2
figure(2)
subplot(1,2,1);
imagesc(airVecPos,E0,fData)
set(gca, 'YDir','normal')
subplot(1,2,2);
plot(airVecPos, fData')
ylim([0 10]);

%% Comparison plot
plot(airVecPos, finalE'); 
hold on
plot(airVecPos, fData', ':');

%% Analysis
figure(3)
vM = finalE(:)>0.5;
histogram((fData(vM) - finalE(vM)))

%% poly plot( deprecated)
% %% Fit
% [xData, yData, zData] = prepareSurfaceData( airVecPos, E0, finalE );
% 
% % Set up fittype and options.
% ft = fittype( 'poly25' );
% 
% % Fit model to data.
% [fitresult, gof] = fit( [xData, yData], zData, ft );
% 
% fData = fitresult(airVecPos,E0'); 
% fData(fData<0)=0;

%% Estimate effect of kapton
kaptonThickness_um = 8;
kaptonThickness_cm = 1e-4*kaptonThickness_um;
Ekap_loss = nan(10,1);
Ekaps = 0.5:0.25:10;
for iKE = 1:numel(Ekaps)
    E0kap_vector = energyStoppingPowerKapton(Ekaps(iKE), [0 kaptonThickness_cm]);
    E0kap = E0kap_vector(2);
    Ekap_loss(iKE) = E0kap_vector(1) - E0kap_vector(2);
end
plot(Ekaps, Ekap_loss, 'o')
Fkap = fit(Ekaps', Ekap_loss, 'Power1', 'Robust', 'LAR')
aw = Fkap.a;
pw = Fkap.b;

EfunK = @(E0,z) real(((E0 - aw.*E0.^pw).^p - z./a).^(1/p))
%% Calculate fData
fData = Efun(E0', airVecPos); 