E0 = 3;
Z = 0:0.005:20;

[energyA, stoppingPowerAir, stoppingPowerWater] = energyStoppingPower(E0, Z);

[energyA2, stoppingPowerAir2, stoppingPowerWater2] = energyStoppingPowerFoil(E0, Z, Z(2), 0.4, 16);

[energyA3, stoppingPowerAir3, stoppingPowerWater3] = energyStoppingPowerFoil(E0, Z, Z(2), 0, 16);

[energyA4, stoppingPowerAir4, stoppingPowerWater4] = energyStoppingPowerFoil(E0, Z, Z(2), 0, 32);

plot(Z, stoppingPowerAir); hold on; 
plot(Z, stoppingPowerAir2);
plot(Z, stoppingPowerAir3);
plot(Z, stoppingPowerAir4);