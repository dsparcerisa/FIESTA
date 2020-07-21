function [energyK, stoppingPowerKapton] = energyStoppingPowerKapton(E0, Z)

kaptonTable = csvread('kapton10MeV.csv');
kaptonTable = kaptonTable(:, 1:2);

energy_kapton = kaptonTable(:,1);
stoppingPower_kapton = kaptonTable(:,2);

rho_kapton = 1.42; %densidad del aire
stoppingPower_kapton = stoppingPower_kapton*rho_kapton;

%Creamos unos vectores iniciales para la energía y el poder de frenado que sean del tamaño de Z
energiesKapton_Z = nan(size(Z));
stoppingPowerKapton_Z = nan(size(Z));
energiesKapton_Z(1) = E0;
stoppingPowerKapton_Z(1) = interp1(energy_kapton, stoppingPower_kapton, E0); %interpolamos para que el valor del poder de frenado se corresponda con la energía

for i=2:numel(Z) %hasta el número de elementos de Z    
    stoppingPowerKapton_Z(i) =  interp1(energy_kapton,stoppingPower_kapton, energiesKapton_Z(i-1));
    energiesKapton_Z(i) = energiesKapton_Z(i-1) - stoppingPowerKapton_Z(i-1)*(Z(i) - Z(i-1)); %la energía a una distancia determinada será igual a la energía inicial (3 Mev) menos la energía que pierde al desplazarse de Z1 a Z2 ((dE/dz)*z)    
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

stoppingPowerKapton = stoppingPowerKapton_Z;
energyK = energiesKapton_Z;   

end

