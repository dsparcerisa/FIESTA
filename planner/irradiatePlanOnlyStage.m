function irradiatePlanOnlyStage(Fcfg, COMstage, thePlan, vector2startPoint, stageLimits, numIterations)
% void irradiatePlanOnlyStage(Fcfg, COMstage, thePlan, vector2startPoint, stageLimits)
% Launches plan irradiation
if (exist('app'))
    app.logLine('Starting plan irradiation...');
else
    disp('Starting plan irradiation...');
end

%% Comprobar que el plan es válido antes de irradiar
planValid = thePlan;
planValid.X = planValid.X + vector2startPoint(1);
planValid.Y = planValid.Y + vector2startPoint(2);
planValid.Z = planValid.Z + vector2startPoint(3);

if planIsInvalid(planValid, stageLimits)
    error('Plan is invalid. Cannot irradiate');
end

%% Create log file
logFileName = datestr(now,'irrLog_yyyy_mm_dd_HH_MM_SS.log');
logFile = fopen(fullfile(Fcfg.logPath,logFileName),'w');

fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') 'Starting plan irradiation...' '\n']);
fprintf(logFile, '\tPlan name: %s\n', thePlan.name);
comment = input('Enter plan description text for log: ', 's');
fprintf(logFile, '\tPlan comment: %s\n', comment);
fprintf(logFile, '\tAbsolut shift from X to (0,0) well: [%3.3f %3.3f %3.3f]\n', vector2startPoint(1), vector2startPoint(2), vector2startPoint(3));

%% Irradiate plan

for j=1:numIterations
    msg = sprintf('Starting iteration %i', j);
    fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
    for i=1:numel(thePlan.X)
        absPos = [thePlan.X(i) thePlan.Y(i) thePlan.Z(i)] + vector2startPoint;
        msg = sprintf('Moving to: [%3.3f %3.3f %3.3f]', absPos(1), absPos(2), absPos(3));
        fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
        [finished, finalPos] = stageControl_moveToAbsPos(COMstage, absPos);
        msg = sprintf('Arrived at: [%3.3f %3.3f %3.3f]', finalPos(1), finalPos(2), finalPos(3));
        if (any(absPos ~= finalPos))
            msg = [msg ' WARNING!'];
        end
        fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
        
        msg = sprintf('Wait %3.3fs', thePlan.t_s(i));
        fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
        pause(thePlan.t_s(i));
        msg = sprintf('Continue movement.\n');
        fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
        
    end
    msg = sprintf('Finished iteration %i\n', j);
    fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
end

%% Close log file and end
msg = sprintf('Plan finished.');       
fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
fclose(logFile);

end

