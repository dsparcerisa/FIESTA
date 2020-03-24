function doseMap = getDoseMap(E0, z, dz, Nprot, targetTh, targetSPR, N0, sigmaPoly)
% CG2D doseMap = getDoseMap(double E0, double z, double dz, double Nprot, 
%  double targetTh, double targetDens, double targetSPR, CG2D N0)

[sigma, Sw_z] = getSandE(E0, z, dz);

if exist('sigmaPoly')
    sigma = polyval(sigmaPoly, z) / 10;
    NZ = createGaussProfile(N0.dx, N0.dy, N0.dx*N0.NX, N0.dy*N0.NY, sigma, sigma);
else
    NZ = convFluence(N0, sigma);
end

MeV2J = 1.602e-13;

E_Mev = Nprot*NZ.data*(Sw_z*targetSPR)*targetTh;
E_J = E_Mev * MeV2J;

m_1voxel_g = (N0.dx * N0.dy * targetTh);
m_1voxel_kg = m_1voxel_g / 1000;

D = E_J ./ m_1voxel_kg;

doseMap = CartesianGrid2D(N0);
doseMap.data = D;

end

