clear all
close all
%Perfil a emborronar
dx=0.01;
dy=0.01;
sizeX=10;
sizeY=10;
sigmaX=5;
sigmaY=5;
beamProfile = createGaussProfile(dx,dy, sizeX, sizeY, sigmaX, sigmaY);
imgPath=beamProfile.data;
%Sigma de emborronamiento
sigmaXY=7;
%Emborronamiento mediante el perfil gaussiano
beamProfileGauss = createGaussProfile(dx,dy, sizeX, sizeY, sigmaXY, sigmaXY);
%Emborronamiento mediante el flitrado de imagen
beamProfileIm = createBeamProfileFromImg(dx,dy, sizeX, sizeY, sigmaXY,imgPath);
%Valores de X e Y
Xvalues = [beamProfile.minX]:[beamProfile.dx]:[beamProfile.maxX];
Yvalues = [beamProfile.minY]:[beamProfile.dy]:[beamProfile.maxY];
%Figuras
figure
subplot(1,3,1)
imagesc(Xvalues,Yvalues,imgPath);
title('Original Profile');
xlabel('X(cm)');
ylabel('Y(cm)');
subplot(1,3,2);
imagesc(Xvalues,Yvalues,beamProfileGauss.data);
title('Gaussian Blurred Profile');
xlabel('X(cm)');
ylabel('Y(cm)');
subplot(1,3,3);
imagesc(Xvalues,Yvalues,beamProfileIm.data);
title('Imgaussfilt Blurred');
xlabel('X(cm)');
ylabel('Y(cm)');