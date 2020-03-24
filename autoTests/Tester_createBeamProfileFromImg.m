%% Crear vacío
clear all
close all
dx = 0.1;
dy = 0.1;
sizeX = 20;
sizeY = 20;

%la imagen la estoy sacando del perfil gaussiano
emptyProfile = createEmptyCG2D(dx, dy, sizeX, sizeY);
centralPtX = emptyProfile.NX/2;
centralPtY = emptyProfile.NY/2;

emptyProfile.data(centralPtX:(centralPtX+1),centralPtY:(centralPtY+1)) = 0.25;

sigmaXY=7;
[beamProfile, beamProfile2] = createBeamProfileFromImg(dx, dy, sizeX, sizeY, sigmaXY, emptyProfile.data);
Xvalues = [beamProfile.minX]:[beamProfile.dx]:[beamProfile.maxX];
Yvalues = [beamProfile.minY]:[beamProfile.dy]:[beamProfile.maxY];
figure(1)
subplot(2,2,1);
emptyProfile.plotSlice
subplot(2,2,3)
beamProfile.plotSlice
title('Gaussian Profile');
xlabel('X(cm)');
ylabel('Y(cm)');
subplot(2,2,2);
Xsum = sum(beamProfile.data, 2);
plot(Xvalues, Xsum, 'o')
xlabel('X');
FX = fit(Xvalues', Xsum, 'gauss1')
sigmaX = FX.c1 / sqrt(2)
hold on
plot(Xvalues, FX(Xvalues));
subplot(2,2,4);
Ysum = sum(beamProfile.data, 1);
plot(Yvalues, Ysum, 'o')
xlabel('Y');
FY = fit(Yvalues', Ysum', 'gauss1')
sigmaY = FY.c1 / sqrt(2)
hold on
plot(Yvalues, FY(Yvalues));

figure(2)
subplot(2,2,1);
emptyProfile.plotSlice
subplot(2,2,3)
beamProfile2.plotSlice
title('Gaussian Profile 2');
xlabel('X(cm)');
ylabel('Y(cm)');
subplot(2,2,2);
Xsum = sum(beamProfile2.data, 2);
plot(Xvalues, Xsum, 'o')
xlabel('X');
FX = fit(Xvalues', Xsum, 'gauss1')
sigmaX = FX.c1 / sqrt(2)
hold on
plot(Xvalues, FX(Xvalues));
subplot(2,2,4);
Ysum = sum(beamProfile2.data, 1);
plot(Yvalues, Ysum, 'o')
xlabel('Y');
FY = fit(Yvalues', Ysum', 'gauss1')
sigmaY = FY.c1 / sqrt(2)
hold on
plot(Yvalues, FY(Yvalues));