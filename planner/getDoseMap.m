function doseMap = getDoseMap(doseCanvas, E0, z, dz, Nprot, targetTh, targetSPR, sigmaPolyX, sigmaPolyY)
% CG2D doseMap = getDoseMap(double E0, double z, double dz, double Nprot, 
%  double targetTh, double targetDens, double targetSPR, CG2D N0)

[~, ~, stoppingPowerWater] = energyStoppingPower(E0, 0:dz:z);

Sw_z = stoppingPowerWater(end);
sigmaX = polyval(sigmaPolyX, z) / 10;
sigmaY = polyval(sigmaPolyY, z) / 10;

NZ = createGaussProfile(doseCanvas.dx, doseCanvas.dy, doseCanvas.dx*doseCanvas.NX, doseCanvas.dy*doseCanvas.NY, sigmaX, sigmaY);

MeV2J = 1.602e-13;

E_Mev = Nprot*NZ.data*(Sw_z*targetSPR)*targetTh;
E_J = E_Mev * MeV2J;

m_1voxel_g = (doseCanvas.dx * doseCanvas.dy * targetTh);
m_1voxel_kg = m_1voxel_g / 1000;

D = E_J ./ m_1voxel_kg;

doseMap = CartesianGrid2D(doseCanvas);  
doseMap.data = D;

end

