function doseMap = getDoseMap_ThickTarget(E0, Nprot, materialList, N0, sigmaPoly)
% CG2D doseMap = getDoseMap_ThickTarget(E0, Nprot, materialList, N0, sigmaPoly)
% Returns dose distribution in LAST MATERIAL

ETRACK_THRESHOLD = 0.1; % MeV
ELOSS_THRESHOLD = 0.03; % 1%

% Load tables
airTable = csvread('air10MeV.csv');
airTable = airTable(:, 1:2);
waterTable = csvread('water10MeV.csv');
airTable = airTable(:, 1:2);
copperTable = csvread('copper10MeV.csv');
copperTable = copperTable(:, 1:2);
goldTable = csvread('gold10MeV.csv');
goldTable = goldTable(:, 1:2);
aluminumTable = csvread('aluminum10MeV.csv');
aluminumTable = aluminumTable(:, 1:2);

% Aire
energy_Air = airTable(:,1); %Energía en MeV
stoppingPower_Air = airTable(:,2); %Poder de frenado en MeV cm2/g
rho_Air = 1.20479E-3; %densidad del aire
stoppingPower_Air = stoppingPower_Air*rho_Air; %obtenemos el poder de frenado en (MeV/cm)

% Water
energy_Water = waterTable(:,1);
stoppingPower_Water = waterTable(:,2);
rho_Water = 1; %densidad del aire
stoppingPower_Water = stoppingPower_Water*rho_Water;

% Aluminum
energy_Al = aluminumTable(:,1);
stoppingPower_Al = aluminumTable(:,2);
rho_Al = 2.6989; %densidad del aire
stoppingPower_Al = stoppingPower_Al*rho_Al;

% gold22qt
rho_Gold = 19.32;
rho_Copper = 8.96;
purityWt = 22/24;
rho_Au22 = 1/( purityWt/rho_Gold + (1-purityWt)/rho_Copper);
A_Au = 196.96657;
A_Cu = 63.546;
purityAt = 1 / (1 + (1/purityWt - 1)*A_Au/A_Cu);
energy_Au22 = goldTable(:,1);
stoppingPower_Au = goldTable(:,2);
stoppingPower_Cu = copperTable(:,2);
stoppingPower_Au22 = purityAt*stoppingPower_Au + (1-purityAt)*stoppingPower_Cu;
stoppingPower_Au22 = stoppingPower_Au22*rho_Au22;

% Copper
energy_Cu = copperTable(:,1);
stoppingPower_Cu = stoppingPower_Cu * rho_Copper;

% 1. Buscar E y Edep hasta la última fila de materialList

E = E0;
Nsteps = size(materialList,1);
exitE = nan(Nsteps, 1);
for i=1:Nsteps
    if E<=0 
        exitE(i) = 0;
        break;
    end
    
    Z = 0; % Z at current slab    
    if strcmp(materialList.Material{i},'Air')
        E_thisMat = energy_Air;
        S_thisMat = stoppingPower_Air;
    elseif strcmp(materialList.Material{i},'Water')
        E_thisMat = energy_Water;
        S_thisMat = stoppingPower_Water;
    elseif strcmp(materialList.Material{i},'Aluminum')
        E_thisMat = energy_Al;
        S_thisMat = stoppingPower_Al;
    elseif strcmp(materialList.Material{i},'Gold22qt')
        E_thisMat = energy_Au22;
        S_thisMat = stoppingPower_Au22;
    elseif strcmp(materialList.Material{i},'Copper')
        E_thisMat = energy_Cu;
        S_thisMat = stoppingPower_Cu;
    else
        error('Material %s not recognized', materialList.Material{i});
    end
    while (true)
        tentativeEdep = E*ELOSS_THRESHOLD;
        S = interp1(E_thisMat, S_thisMat, E);
        tentativeStep = tentativeEdep / S;
        maxStep = materialList.Thickness_cm(i) - Z;
        if tentativeStep <= maxStep
            Z = Z + tentativeStep;
            E = E - tentativeEdep;
            if (E <= ETRACK_THRESHOLD)
                E = 0;
                exitE(i) = 0; % Deposit final energy locally
                break;
            end
        else
            tentativeStep = maxStep;
            E = E - tentativeStep * S;
            exitE(i) = E;
            break;
        end
    end    
end

EdepInLastSlab = exitE(end-1) - exitE(end);

sigma = polyval(sigmaPoly, sum(materialList.Thickness_cm(1:(end-1)))) / 10;
NZ = createGaussProfile(N0.dx, N0.dy, N0.dx*N0.NX, N0.dy*N0.NY, sigma, sigma);

MeV2J = 1.602e-13;
 
E_Mev = Nprot*NZ.data*EdepInLastSlab;
E_J = E_Mev * MeV2J;
 
m_1voxel_g = (N0.dx * N0.dy * materialList.Thickness_cm(end));
m_1voxel_kg = m_1voxel_g / 1000;
 
D = E_J ./ m_1voxel_kg;

doseMap = CartesianGrid2D(N0);
doseMap.data = D;

end

