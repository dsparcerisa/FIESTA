basePath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_07_08 - Estudio quenching CMAM';

fullPath = '/Users/dani/Google Drive/Trabajo/UCM/01-Proyectos/24-Experimentos CMAM/04-Radiocromicas/20_07_08 - Estudio quenching CMAM/scan-dani-6sept/EBT2/10MeV';

imgPath = [fullPath filesep 'Rad18.tif'];
I = imread(imgPath);

% El plan es:
% # X (cm), Y (cm), Z(cm), t (s)
% 5, 0, 0, 0.05
% 5, 0, -2, 0.05
% 5, 0, -4, 0.05
% 5, 0, -6, 0.05
% 5, 0, -8, 0.05
% 3, 2.5, 0, 0.05
% 2, 1.5, -1, 0.05
% 3, 0.5, -2, 0.05
% 2, -0.5, -3, 0.05
% 3, -1.5, -4, 0.05
% 2, -2.5, -5, 0.05
% 1, 2.5, -6, 0.05
% 0, 1.5, -7, 0.05
% 1, 0.5, -8, 0.05
% 0, -0.5, -9, 0.05
% 1, -1.5, -10, 0.05
% -1, 2.5, 0, 0.25
% -2, 1.5, -1, 0.25
% -1, 0.5, -2, 0.25
% -2, -0.5, -3, 0.25
% -1, -1.5, -4, 0.25
% -2, -2.5, -5, 0.25
% -3, 2.5, -6, 0.25
% -4, 1.5, -7, 0.25
% -3, 0.5, -8, 0.25
% -4, -0.5, -9, 0.25
% -3, -1.5, -10, 0.25

% Posici�n 1
[X1,Y1] = getMaxCoordinates(I);

% Posici�n 3
[X3,Y3] = getMaxCoordinates(I);

% Posici�n 27
[X27,Y27] = getMaxCoordinates(I);

% XY vectors
Yvec = [X3-X1, Y3-Y1] / 4;
Xvec = [X27-X1, Y27-Y1] / 6;

% Prueba vectores
xy2 = [X1, Y1] + 2*Yvec;
xy7 = [X1, Y1] + 2*Xvec;
xy21 = [X1, Y1] + 4*Xvec;

imshow(I); hold on
plot([X1 X3 X27 xy2(1) xy7(1) xy21(1)], [Y1 Y3 Y27 xy2(2) xy7(2) xy21(2)] , 'bo');

% Resto
Xpositions_cm = 3 - [3 2 3 2 3 2 1 0 1 0 1 -1 -2 -1 -2 -1 -2 -3 -4 -3 -4 -3]';
Ypositions_cm = 2.5 - [2.5 1.5 0.5 -0.5 -1.5 -2.5 2.5 1.5 0.5 -0.5 -1.5 2.5 1.5 0.5 -0.5 -1.5 -2.5 2.5 1.5 0.5 -0.5 -1.5]' ;

PxPositions =  Xpositions_cm * Xvec + Ypositions_cm * Yvec + [X1,Y1];

% Prueba
imshow(I); hold on
plot(PxPositions(:,1), PxPositions(:,2), 'bo');

%% Do the cropping
%Square size in pixels: 100 x 100;
SSiP = 100;

Npoints = size(PxPositions, 1);
figure(1); imshow(I)
figure(2);
allI = {};
for i = 1:Npoints
    allI{i} = imcrop(I, [PxPositions(i,1)-SSiP/2 PxPositions(i,2)-SSiP/2 SSiP SSiP]);
    subplot(4,6,i);
    imshow(allI{i});
    title(sprintf('%i',i));
end

%% 