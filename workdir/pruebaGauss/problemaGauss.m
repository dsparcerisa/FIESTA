% Problema: dado un haz con una determinada sigma (medida en X a partir de una integral),
% ¿cómo varía el número de partículas en un determinado radio R? ¿Es esto
% equivalente a la integral definida en 1D entre 0 y R?

%% 1. Hacer una simulación
% CG2D beamProfile = createGaussProfile(dx, dy, sizeX, sizeY, sigmaX, sigmaY)

dx = 0.01;
dy = 0.01;
sizeX = 3;
sizeY = 3;
sigmaX = 0.1; % cm
sigmaY = 0.1; % cm
beamProfile = createGaussProfile(dx, dy, sizeX, sizeY, sigmaX, sigmaY);
beamProfile.plotSlice

%% 2. Comprobar que la simulación es OK
radius = 0.25; % cm
radiusInVoxels = radius / dx;
[mask, xcentre_v, ycentre_v, sigmaX_vx, sigmaY_vx, meanValue, mX, mY] = meanAndCenterMass(-beamProfile.data, radiusInVoxels);

%% 3. Comprobar el valor numéricamente
sumInsideRadius = sum(beamProfile.data(mask)) / sum(beamProfile.data(:))

%% 4. Hacer la simulación completa para distintos Z
%sigFromZ = @(z) 0.125.*z - 0.026;
sigFromZ = @(z) 0.0127*z.^2 + 0.143.*z; 
Zval = 2:0.1:10;
sigVal = sigFromZ(Zval);
fracInsideRadius = nan(numel(sigVal), 1);

for i=1:numel(sigVal)
    beamProfile = createGaussProfile(dx, dy, sizeX, sizeY, sigVal(i)/10, sigVal(i)/10);
    mask = meanAndCenterMass(-beamProfile.data, radiusInVoxels);
    fracInsideRadius(i) = sum(beamProfile.data(mask)) / sum(beamProfile.data(:));
end

subplot(2,1,1);
plot(Zval, sigVal, 'r-')
grid on
xlabel('Air depth [cm]');
ylabel('Sigma [mm]');
subplot(2,1,2);
plot(Zval, fracInsideRadius, 'r-')
grid on
xlabel('Air depth [cm]');
ylabel('Beam fraction hitting water surface');

