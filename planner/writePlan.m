function writePlan(plan, outFile)
ID = fopen(outFile, 'w');

% Write name
fprintf(ID, '# %s\n', plan.name);

% Write type
fprintf(ID, '%s\n', plan.mode);

if strcmp(plan.mode,'FLASH')
    % Write header
    fprintf(ID, '# E (MeV), I (nA), t_rendija (ms), numSpots\n');
    fprintf(ID, '%3.2f, %3.3f, %3.3f, %u\n', plan.E, plan.I, plan.tRendija, plan.numSpots);
    
    % Write table
    fprintf(ID, '# X (cm), Y (cm), Z (cm), #shots\n');
    for i=1:plan.numSpots
        fprintf(ID, '%3.2f, %3.2f, %3.2f, %u\n', plan.X(i), plan.Y(i), plan.Z(i), plan.Nshots(i));
    end
    
elseif strcmp(plan.mode,'CONV')
    
    % Write header
    fprintf(ID, '# E (MeV), I (nA), codFiltro, numSpots\n');
    fprintf(ID, '%3.2f, %3.3f, %s, %u\n', plan.E, plan.I, plan.codFiltro, plan.numSpots);
    
    % Write table
    fprintf(ID, '# X (cm), Y (cm), Z (cm), t(s)\n');
    for i=1:plan.numSpots
        fprintf(ID, '%3.2f, %3.2f, %3.2f, %3.3f\n', plan.X(i), plan.Y(i), plan.Z(i), plan.t_s(i));
    end    
    
else
    fclose(ID);
    error('Unrecognized plan type (%s)', plan.mode);
end

fclose(ID);

end

