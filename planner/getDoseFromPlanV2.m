function dose = getDoseFromPlanV2(doseCanvas, plan, dz, targetTh, targetSPR, factorImuestra, airDepthAtPos0, sigmaPolyX, sigmaPolyY)
% CG2D dose = getDoseFromPlanV2(doseCanvas, plan, dz, targetTh, targetSPR, N0, factorImuestra, airDepthAtPos0, sigmaPolyX, sigmaPolyY)
thickTarget = false;

if strcmp(plan.mode, 'FLASH')
    plan.Q = plan.I * factorImuestra * plan.tRendija * double(plan.Nshots);
elseif strcmp(plan.mode, 'CONV')    
    if strcmp(plan.codFiltro, 'PP100')
        factorPP = 6.7e-4;
    elseif strcmp(plan.codFiltro, 'PP25')
        factorPP = 8.8e-5;      
    elseif strcmp(plan.codFiltro, '1')
        factorPP = 1;   
    elseif strcmp(plan.codFiltro, 'PP2capas')
        factorPP = 1/25; % Tentativo  
    end
    plan.Q = plan.I * factorImuestra * factorPP * plan.t_s * 1000;
else
    error('Plan mode %s not recognized', plan.mode);
end

% Calculate plan.Q

pC2pNumber = 1e7/1.602176634;
doseCanvas.data(:) = 0;
limits = [doseCanvas.minX doseCanvas.maxX doseCanvas.minY doseCanvas.maxY];

depths = airDepthAtPos0 - plan.Z;

doseKernel = getDoseMap(doseCanvas, plan.E, depths(1), dz, pC2pNumber, targetTh, targetSPR, sigmaPolyX, sigmaPolyY);
for i=1:plan.numSpots
    if i>1 && depths(i)~=depths(i-1)        
        doseKernel = getDoseMap(doseCanvas, plan.E, depths(1), dz, pC2pNumber, targetTh, targetSPR, sigmaPolyX, sigmaPolyY);
    end
    Di = doseKernel.copy;
    Di.data = plan.Q(i) * Di.data;
    Di.shift([plan.X(i) plan.Y(i)]);
    
    doseCanvas = doseCanvas + Di;
end
doseCanvas.crop(limits);
dose = doseCanvas;

end

