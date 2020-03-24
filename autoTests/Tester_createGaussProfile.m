clear all
close all
dx=0.01;
dy=0.01;
sizeX=10;
sizeY=10;
sigmaX=1;
sigmaY=5;
beamProfile = createGaussProfile(dx, dy, sizeX, sizeY, sigmaX, sigmaY);
Xvalues = [beamProfile.minX]:[beamProfile.dx]:[beamProfile.maxX];
Yvalues = [beamProfile.minY]:[beamProfile.dy]:[beamProfile.maxY];
subplot(1,2,1);
beamProfile.plotSlice
title('Gaussian Profile');
xlabel('X(cm)');
ylabel('Y(cm)');
subplot(2,2,2);
Ysum = sum(beamProfile.data,1);
plot(Xvalues, Ysum', 'o')
FX = fit(Xvalues', Ysum', 'gauss1')
sigmaX = FX.c1 / sqrt(2)
hold on
plot(Xvalues, FX(Xvalues));
subplot(2,2,4);
Xsum = sum(beamProfile.data,2);
plot(Yvalues, Xsum, 'o')
FY = fit(Yvalues', Xsum, 'gauss1')
sigmaY = FY.c1 / sqrt(2)
hold on
plot(Yvalues, FY(Yvalues));