function dose = getDoseFromPlan_ThickTarget(doseCanvas, plan, N0, factorImuestra, materialTable, sigmaPoly)
% CG2D dose = getDoseFromPlan_ThickTarget(doseCanvas, plan, N0, factorImuestra, materialTable, sigmaPoly)

if strcmp(plan.mode, 'FLASH')
    plan.Q = plan.I * factorImuestra * plan.tRendija * double(plan.Nshots);
    % sigmaPoly = [0.0097 0.24 0]; % sin PP settings 5 feb
    % sigmaPoly = [0.0086 0.23 0]; % settings 6 feb (tras corregir la z)

elseif strcmp(plan.mode, 'CONV')    
    if strcmp(plan.codFiltro, 'PP100')
        factorPP = 6.7e-4;
    elseif strcmp(plan.codFiltro, 'PP25')
        factorPP = 8.8e-5;      
    elseif strcmp(plan.codFiltro, '1')
        factorPP = 1;
        % sigmaPoly = [0.0097 0.24 0]; sin PP settings 5 feb
        % sigmaPoly = [0.0086 0.23 0]; % settings 6 feb (tras corregir la z)    
    elseif strcmp(plan.codFiltro, 'PP2capas')
        factorPP = 1/25; % Tentativo
        %sigmaPoly = [0 0.26 0.63]; % con PP settings 5 feb
        % sigmaPoly = [0 0.24 0.76]; % con PP settings 6 feb        
    else
        error('Pepperpot %s not recognized', plan.codFiltro);
    end

    plan.Q = plan.I * factorImuestra * factorPP * plan.t_s * 1000;

else
    error('Plan mode %s not recognized', plan.mode);
end

% Calculate plan.Q

pC2pNumber = 1e7/1.602176634;
doseCanvas.data(:) = 0;
limits = [doseCanvas.minX doseCanvas.maxX doseCanvas.minY doseCanvas.maxY];

% SOLO VALE PARA UNA Z EN CADA PUNTO!

doseKernel = getDoseMap_ThickTarget(plan.E, pC2pNumber, materialTable, N0, sigmaPoly);

for i=1:plan.numSpots
    Di = doseKernel.copy;
    Di.data = plan.Q(i) * Di.data;
    Di.shift([plan.X(i) plan.Y(i)]);    
    doseCanvas = doseCanvas + Di;
end
doseCanvas.crop(limits);
dose = doseCanvas;

end

