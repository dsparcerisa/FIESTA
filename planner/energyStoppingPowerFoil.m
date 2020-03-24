function [energyA, stoppingPowerAir, stoppingPowerWater] = energyStoppingPowerFoil(E0, Z, entrance2Foil, Au_um, Al_um)

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

%Creamos unos vectores iniciales para la energía y el poder de frenado que sean del tamaño de Z
energiesAir_Z = nan(size(Z));
stoppingPowerAir_Z = nan(size(Z));
energiesAir_Z(1) = E0;
stoppingPowerAir_Z(1) = interp1(energy_Air, stoppingPower_Air, E0); %interpolamos para que el valor del poder de frenado se corresponda con la energía
stoppingPowerWater_Z = nan(size(Z));
stoppingPowerWater_Z(1) = interp1(energy_Water, stoppingPower_Water, E0);

if (Au_um+Al_um)>0 && entrance2Foil > 0 && entrance2Foil < Z(end)
    
    % Fase 1: avanzar entrance2Foil hasta el Foil:
    NpreFoil = find(Z<entrance2Foil, 1, 'last');
    NFoil = NpreFoil + 1;
    for i=2:NpreFoil %hasta el número de elementos de Z
        stoppingPowerAir_Z(i) =  interp1(energy_Air,stoppingPower_Air, energiesAir_Z(i-1));
        energiesAir_Z(i) = energiesAir_Z(i-1) - stoppingPowerAir_Z(i-1)*(Z(i) - Z(i-1)); %la energía a una distancia determinada será igual a la energía inicial (3 Mev) menos la energía que pierde al desplazarse de Z1 a Z2 ((dE/dz)*z)
        stoppingPowerWater_Z(i) = interp1(energy_Air,stoppingPower_Water, energiesAir_Z(i-1)); %energía negativas
    end
    
    % Asumimos que el foil cae siempre al principio del voxel
    % Fase2a: Au
    EInFoil = energiesAir_Z(NpreFoil);
    stoppingPowerAu22_Z = interp1(energy_Au22, stoppingPower_Au22, EInFoil);
    EInFoil = EInFoil - stoppingPowerAu22_Z*Au_um/10000;
    % Fase2b: Al (en steps de 1 um)
    Al_um = round(Al_um);
    for i=1:Al_um
        stoppingPowerAl_Z = interp1(energy_Al, stoppingPower_Al, EInFoil);
        EInFoil = EInFoil - stoppingPowerAl_Z/10000;
    end
    % Fase 2c: remamente de aire
    remainingZ_Air = Z(NpreFoil + 1) - Z(NpreFoil) - (Al_um+Au_um)/10000;
    
    if remainingZ_Air>0
        
        stoppingPowerAir_Z(NFoil) =  interp1(energy_Air,stoppingPower_Air, EInFoil);
        energiesAir_Z(NFoil) = EInFoil;
        stoppingPowerWater_Z(NFoil) =  interp1(energy_Water,stoppingPower_Water, EInFoil);
        
        energiesAir_Z(NFoil+1) = EInFoil - remainingZ_Air*stoppingPowerAir_Z(NFoil);
        stoppingPowerAir_Z(NFoil+1) =  interp1(energy_Air,stoppingPower_Air, energiesAir_Z(NFoil+1));
        stoppingPowerWater_Z(NFoil+1) =  interp1(energy_Water,stoppingPower_Water, energiesAir_Z(NFoil+1));
        
    else remainingZ_Air<0
        error('Situation not yet coded');
    end
    
    % Fase 3: avanzar hasta el final
    for i=(NFoil+2):numel(Z)
        stoppingPowerAir_Z(i) =  interp1(energy_Air,stoppingPower_Air, energiesAir_Z(i-1));
        energiesAir_Z(i) = energiesAir_Z(i-1) - stoppingPowerAir_Z(i-1)*(Z(i) - Z(i-1)); %la energía a una distancia determinada será igual a la energía inicial (3 Mev) menos la energía que pierde al desplazarse de Z1 a Z2 ((dE/dz)*z)
        stoppingPowerWater_Z(i) = interp1(energy_Air,stoppingPower_Water, energiesAir_Z(i-1)); %energía negativas
    end
    
else
    
    for i=2:numel(Z)
        stoppingPowerAir_Z(i) =  interp1(energy_Air,stoppingPower_Air, energiesAir_Z(i-1));
        energiesAir_Z(i) = energiesAir_Z(i-1) - stoppingPowerAir_Z(i-1)*(Z(i) - Z(i-1)); %la energía a una distancia determinada será igual a la energía inicial (3 Mev) menos la energía que pierde al desplazarse de Z1 a Z2 ((dE/dz)*z)
        stoppingPowerWater_Z(i) = interp1(energy_Air,stoppingPower_Water, energiesAir_Z(i-1)); %energía negativas
    end
    
end

% figure(1)
% plot(Z,stoppingPowerAir_Z);
% xlabel('z (cm)');
% ylabel('Stopping Power in Air (MeV/cm)');
%
% figure(2)
% plot(Z,energiesAir_Z);
% xlabel('z (cm)');
% ylabel('Energy (MeV)');
% title('Energy (Air)');
%
% figure(3)
% plot(Z,stoppingPowerWater_Z);
% xlabel('z (cm)');
% ylabel('Stopping Power in Water (MeV/cm)');

stoppingPowerAir = stoppingPowerAir_Z;
stoppingPowerWater = stoppingPowerWater_Z;
energyA = energiesAir_Z;

end

