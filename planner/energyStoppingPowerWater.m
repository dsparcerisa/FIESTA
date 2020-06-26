function [energyW, stoppingPowerWater] = energyStoppingPowerWater(E0, Z)

waterTable = csvread('water10MeV.csv');
waterTable = waterTable(:, 1:2);

energy_Water = waterTable(:,1);
stoppingPower_Water = waterTable(:,2);

rho_Water = 1; %densidad del aire
stoppingPower_Water = stoppingPower_Water*rho_Water;

%Creamos unos vectores iniciales para la energía y el poder de frenado que sean del tamaño de Z
energiesWater_Z = nan(size(Z));
stoppingPowerWater_Z = nan(size(Z));
energiesWater_Z(1) = E0;
stoppingPowerWater_Z(1) = interp1(energy_Water, stoppingPower_Water, E0); %interpolamos para que el valor del poder de frenado se corresponda con la energía

for i=2:numel(Z) %hasta el número de elementos de Z    
    stoppingPowerWater_Z(i) =  interp1(energy_Water,stoppingPower_Water, energiesWater_Z(i-1));
    energiesWater_Z(i) = energiesWater_Z(i-1) - stoppingPowerWater_Z(i-1)*(Z(i) - Z(i-1)); %la energía a una distancia determinada será igual a la energía inicial (3 Mev) menos la energía que pierde al desplazarse de Z1 a Z2 ((dE/dz)*z)    
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

stoppingPowerWater = stoppingPowerWater_Z;
energyW = energiesWater_Z;   

end

