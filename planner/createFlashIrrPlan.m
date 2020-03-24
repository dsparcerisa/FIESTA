function [plan, totTime, doseRate] = createFlashIrrPlan( plateDose, doseSlice_1pC, I_muestra, deltaXY, miniSpotsPerSpot, shotTime_ms)
% creates a (FLASH rate) plan spreading the for each of the wanted Slots 
% from a given doseSlice calculated for 1 pC

if all(size(plateDose) == [8 12])
    well2wellDist_cm = 0.899;
    NX = 12; NY = 8;
elseif all(size(plateDose) == [4 2])
    well2wellDist_cm = 1.125;
    NX = 2; NY = 4;
else
    error('plateDose definido de un tamaño desconocido');
end

% Usamos el mismo well para la cubeta
wellDiam = 6.35; % mm
wellRadius_cm = 0.1 * wellDiam / 2;

% Shoot and create a 3x3 / 2x3x2 / 2x2 slice and compare with well to get the pC2Gy
% factor:
if miniSpotsPerSpot == 9
    deltaMatrix = ...
       [-deltaXY -deltaXY; ...
        -deltaXY 0; ...
        -deltaXY deltaXY; ...
        0 deltaXY; ...
        0 0; ...
        0 -deltaXY; ...
        deltaXY -deltaXY; ...
        deltaXY 0; ...
        deltaXY deltaXY ];
    relWeight = ones(9, 1);
    relWeight(5) = 0.9;

elseif miniSpotsPerSpot == 7
    deltaMatrix = ...
       [-deltaXY/2 -deltaXY*cosd(30); ...
        +deltaXY/2 -deltaXY*cosd(30); ...
        +deltaXY/2 +deltaXY*cosd(30); ...
        -deltaXY/2 +deltaXY*cosd(30); ...
        -deltaXY 0; ...
        0 0; ...
        deltaXY 0];        
    relWeight = ones(9, 1);
    relWeight(6) = 0.82;
    
elseif miniSpotsPerSpot == 4
    deltaMatrix = [-deltaXY/2 deltaXY/2; ...
                -deltaXY/2 -deltaXY/2; ...
                deltaXY/2 -deltaXY/2; ...
                deltaXY/2 deltaXY/2];
    relWeight = ones(4, 1);
    
elseif miniSpotsPerSpot == 1
    deltaMatrix = [0 0];
    relWeight = 1;
else
    error('Wrong number of miniSpots per spot');
end

doseSliceInWell = CartesianGrid2D(doseSlice_1pC);

for i=1:miniSpotsPerSpot    
    doseSliceMoveable = doseSlice_1pC.copy;
    doseSliceMoveable.shift(deltaMatrix(i, :));
    doseSliceMoveable.data = doseSliceMoveable.data * relWeight(i);
    doseSliceInWell = doseSliceInWell + doseSliceMoveable;
end

%% Enfocamos para un pocillo también en las microcubetas
well0 = getWell(CartesianGrid2D(doseSlice_1pC), wellRadius_cm, [0 0]);
wellDoses = getStats(well0, doseSliceInWell);
meanWellDose_1pC = mean(wellDoses);
doseRate = meanWellDose_1pC * I_muestra * 1000;

% Positions in reference with the center of the first spot
Xpos = well2wellDist_cm*(0:(-1):(-(NX-1)));
Ypos = well2wellDist_cm*(0:(-1):(-(NY-1)));
[x,y] = meshgrid(Xpos, Ypos);

% total number of spots
totalWells = sum(~isnan(plateDose(:)));
totalSpots = miniSpotsPerSpot * totalWells;
plan.X = nan(totalSpots, 1);
plan.Y = nan(totalSpots, 1);
plan.Q = nan(totalSpots, 1);
plan.numSpots = totalSpots;

j = 1;
nans = 0;
% Create plan 
for i=1:numel(plateDose)
    
    if isnan(plateDose(i))
        nans = nans + 1;
        continue
    end
    
    for k = 1:miniSpotsPerSpot
        plan.X(j) = x(i) + deltaMatrix(k, 1);
        plan.Y(j) = y(i) + deltaMatrix(k, 2);
        plan.Q(j) = relWeight(k) * plateDose(i) / meanWellDose_1pC;
        j = j+1;        
    end
    
end

% Calcular # shots
Qshot_pC = shotTime_ms * I_muestra;
plan.Nshots = round(plan.Q / Qshot_pC);
plan.Q = Qshot_pC * double(plan.Nshots);
meanWellDose_1shot = Qshot_pC * meanWellDose_1pC

%% Calculate irradiation time for Flash

% % Calculate TOTAL IRRADIATION TIME
% % NORMAL (3x3)
% % T_irr + T_between wells (Asumimos 1.2 + 0.66*x)
distanceX = abs(diff(plan.X));
distanceY = abs(diff(plan.Y));
if miniSpotsPerSpot==1
     distance = sum(distanceX) + sum(distanceY);
else
     distance = sum(distanceX(distanceX>deltaXY)) + sum(distanceY(distanceY>deltaXY));
end
irrTime = 0.32 * sum(plan.Nshots);

if miniSpotsPerSpot==9
    innerMovementTime = 10*totalWells;
elseif miniSpotsPerSpot==7
    innerMovementTime = 6*totalWells;
elseif miniSpotsPerSpot==4
    innerMovementTime = 3*totalWells; % Inventado
elseif miniSpotsPerSpot == 1
    innerMovementTime = 0;
else
    error('miniSpots per Spot undefined');
end

outerMovementTime = 1.2 * 0.66*distance;
totTime = innerMovementTime + irrTime + outerMovementTime;

end

