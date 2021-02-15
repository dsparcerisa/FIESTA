%% Para comenzar desde aquí (específico de Rad3)
radName = 'Rad3';
saveFile = [radName '.mat']
load(saveFile);

distanciaBase = 5.3;
Zpos = (0:10) + distanciaBase;
% Corregir el error de las primeras irradiaciones
Zpos2 = Zpos;
Zpos(5) = Zpos(3);

%% Método 1: obtener las sigmas directamente de la imagen antes de hacer la
% conversión a dosis, usando meanAndCenterMass de mala manera

sigmasX_M1_R = nan(Npoints,1);
sigmasY_M1_R = nan(Npoints,1);
sigmasX_M1_G = nan(Npoints,1);
sigmasY_M1_G = nan(Npoints,1);
sigmasX_M1_B = nan(Npoints,1);
sigmasY_M1_B = nan(Npoints,1);

radius_pixels = radius_cm*pixCM;

for i=1:Npoints    
    [~, ~, ~, sigmaXR, sigmaYR, ~, ~, ~, deltasigmaXR, deltasigmaYR, ~] = meanAndCenterMass(double(allI{i}(:,:,1)),radius_pixels);
    [~, ~, ~, sigmaXG, sigmaYG, ~, ~, ~, deltasigmaXG, deltasigmaYG, ~] = meanAndCenterMass(double(allI{i}(:,:,2)),radius_pixels);
    [~, ~, ~, sigmaXB, sigmaYB, ~, ~, ~, deltasigmaXB, deltasigmaYB, ~] = meanAndCenterMass(double(allI{i}(:,:,3)),radius_pixels);
    sigmasX_M1_R(i) = 10*sigmaXR/pixCM;
    sigmasY_M1_R(i) = 10*sigmaYR/pixCM;
    sigmasX_M1_G(i) = 10*sigmaXG/pixCM;
    sigmasY_M1_G(i) = 10*sigmaYG/pixCM;
    sigmasX_M1_B(i) = 10*sigmaXB/pixCM;
    sigmasY_M1_B(i) = 10*sigmaYB/pixCM;
end
sigmasX_M1 = (sigmasX_M1_R + sigmasX_M1_G + sigmasX_M1_R)*(1/3);
sigmasY_M1 = (sigmasY_M1_R + sigmasY_M1_G + sigmasY_M1_R)*(1/3);

%% Plot 1
figure(1);
plot(Zpos, sigmasX_M1_R(1:11),'ro'); hold on
plot(Zpos, sigmasY_M1_R(1:11),'rx');
plot(Zpos, sigmasX_M1_G(1:11),'go');
plot(Zpos, sigmasY_M1_G(1:11),'gx');
plot(Zpos, sigmasX_M1_B(1:11),'bo');
plot(Zpos, sigmasY_M1_B(1:11),'bx');
plot(Zpos, sigmasX_M1(1:11),'ko');
plot(Zpos, sigmasY_M1(1:11),'kx');
grid on
% --> Conclusión: el azul no sale suficientemente bien, vamos al rojo y verde solamente:

%% Plot 2
sigmasX_M1b = 0.5 * (sigmasX_M1_G + sigmasX_M1_R)
sigmasY_M1b = 0.5 * (sigmasY_M1_G + sigmasY_M1_R)

figure(2);
plot(Zpos, sigmasX_M1_R(1:11),'ro'); hold on
plot(Zpos, sigmasY_M1_R(1:11),'rx');
plot(Zpos, sigmasX_M1_G(1:11),'go');
plot(Zpos, sigmasY_M1_G(1:11),'gx');
plot(Zpos, sigmasX_M1b(1:11),'ko-');
plot(Zpos, sigmasY_M1b(1:11),'kx-');
grid on
% --> Conclusión: este método es regulero

%% Plot 3: Comparar alta dosis y baja dosis
figure(3)
plot(Zpos, sigmasX_M1b(1:11), 'bo'); hold on
plot(Zpos, sigmasY_M1b(1:11), 'bx');
plot(Zpos2, sigmasX_M1b(12:22), 'ro');
plot(Zpos2, sigmasY_M1b(12:22), 'rx');
ylim([0 3]); grid on
% --> Conclusión: se mide más obviamente para la dosis alta

%% Método 2: obtener las sigmas de la dosis, usando meanAndCenterMass solamente

sigmasX_M2 = nan(Npoints,1);
sigmasY_M2 = nan(Npoints,1);

for i=1:Npoints  
    [dose, varMat] = getDoseMicke(double(allI{i}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
    [mask, xcentre, ycentre, sigmaX, sigmaY, meanvalue, ~, ~, deltasigmaX, deltasigmaY, stdDose] = meanAndCenterMass(-dose.data,radius_pixels);
    %[sigmaX, sigmaY, deltasigmaX, deltasigmaY, xCenter, yCenter] = getSigmas(10*dose.getAxisValues('X'),10*dose.getAxisValues('Y'),dose.data,0.05*dose.data+0.01);
    sigmasX_M2(i) = 10*sigmaX/pixCM;
    sigmasY_M2(i) = 10*sigmaY/pixCM;
end

%% Plot 4
figure(4)
plot(Zpos, sigmasX_M2(1:11), 'bo'); hold on
plot(Zpos, sigmasY_M2(1:11), 'bx');
plot(Zpos2, sigmasX_M2(12:22), 'ro');
plot(Zpos2, sigmasY_M2(12:22), 'rx');
ylim([0 3]); grid on
% --> Conclusión: el método de la imagen sobreestima el valor de sigma,
% mejor usar siempre la dosis. Parece mucho más estable para dosis altas
% (punto más pequeño) que para dosis bajas - puntos grandes.

%% Método 3: obtenerlo de la dosis y usando getSigmas con un error dado

sigmasX_M3 = nan(Npoints,1);
sigmasY_M3 = nan(Npoints,1);

for i=1:Npoints  
    [dose, varMat] = getDoseMicke(double(allI{i}), CoefR1, CoefG1, CoefB1, pixCM, deltas, maxBits);
    %[mask, xcentre, ycentre, sigmaX, sigmaY, meanvalue, ~, ~, deltasigmaX, deltasigmaY, stdDose] = meanAndCenterMass(-dose.data,radius_pixels);
    [sigmaX, sigmaY, deltasigmaX, deltasigmaY, xCenter, yCenter] = getSigmas(10*dose.getAxisValues('X'),10*dose.getAxisValues('Y'),dose.data,0.1*ones(size(dose.data)),4);
    sigmasX_M3(i) = sigmaX;
    sigmasY_M3(i) = sigmaY;
    %pause
end

%% Plot 5
figure(5)
plot(Zpos, sigmasX_M3(1:11), 'bo'); hold on
plot(Zpos, sigmasY_M3(1:11), 'bx');
plot(Zpos2, sigmasX_M3(12:22), 'ro');
plot(Zpos2, sigmasY_M3(12:22), 'rx');
ylim([0 3]); grid on
% --> Conclusión: el método de la imagen sobreestima el valor de sigma,
% mejor usar siempre la dosis. Parece mucho más estable para dosis altas
% (punto más pequeño) que para dosis bajas - puntos grandes.