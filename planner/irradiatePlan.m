function irradiatePlan(COMstage, COMshutter, thePlan, vector2startPoint, stageLimits, airDepthAtPos0, app)
% void irradiatePlan(COMstage, COMshutter, plan, vector2startPoint, stageLimits, airDepthAtPos0, app)
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
logFile = fopen(fullfile('logs',logFileName),'w');

fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') 'Starting plan irradiation...' '\n']);
fprintf(logFile, '\tPlan name: %s\n', thePlan.name);
comment = input('Enter plan description text for log: ');
fprintf(logFile, '\tPlan comment: %s\n', comment);
fprintf(logFile, '\tPlan type: %s\n', thePlan.mode);
if strcmp(thePlan.mode, 'FLASH')
    fprintf(logFile, '\tExpected slit time: %s\n', thePlan.tRendija);
else
    fprintf(logFile, '\tPepperPot code: %s\n', thePlan.codFiltro);
end

fprintf(logFile, '\tAir depth at stage position Z=0: %3.2f cm\n', airDepthAtPos0);
fprintf(logFile, '\tAbsolut shift from X to (0,0) well: [%3.3f %3.3f %3.3f]\n', vector2startPoint(1), vector2startPoint(2), vector2startPoint(3));

%% Input intensity at FC1
I_nA = input('Input currrent intensity at FC1 (nA): ');
fprintf(logFile, '\tCurrent intensity: %3.3f nA (planned for %3.3f nA)\n\n', I_nA, thePlan.I);

absPos = [thePlan.X(1) thePlan.Y(1) thePlan.Z(1)] + vector2startPoint;
msg = sprintf('Moving to: [%3.3f %3.3f %3.3f]', absPos(1), absPos(2), absPos(3));
    fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
    finished = stageControl_moveToAbsPos(COMstage, absPos);
    
%% Irradiate plan
for i=1:numel(thePlan.X)
    
    absPos = [thePlan.X(i) thePlan.Y(i) thePlan.Z(i)] + vector2startPoint;
    msg = sprintf('Moving to: [%3.3f %3.3f %3.3f]', absPos(1), absPos(2), absPos(3));
    fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
    finished = stageControl_moveToAbsPos(COMstage, absPos);
    msg = sprintf('Arrived at: [%3.3f %3.3f %3.3f]', absPos(1), absPos(2), absPos(3));
    fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);

    if strcmp(thePlan.mode, 'FLASH')
        msg = sprintf('Deliver: %i FLASH shots', thePlan.Nshots(i));    
        fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
    
        Shutter(COMshutter, 'f', thePlan.Nshots(i));
        msg = sprintf('Finished.\n');       
        fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
 
    elseif strcmp(thePlan.mode, 'CONV')
        Configure_shutter(COMshutter, 't', thePlan.t_s(i))
        msg = sprintf('Open shutter for %3.3fs', thePlan.t_s(i));  
        fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);

        Shutter(COMshutter,'n',1);
        msg = sprintf('Shutter closed.\n');       
        fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
        
    else
        error('Unrecognized plan type (%s)', thePlan.mode);
    end
       
    fprintf('Irradiated spot %i / %i\n', i, thePlan.numSpots);
end


%% Close log file and end
msg = sprintf('Plan irradiation finished.');       
fprintf(logFile, [datestr(now,'[HH:MM:SS.FFF] ') msg '\n']);
fclose(logFile);

if (exist('app'))
    app.logLine('Plan irradiation finished.');
else
    disp('Plan irradiation finished.');
end

end

